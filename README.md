# From_feather_to_fur
Code used in the manuscript: From feather to fur: gull and mink H5N1 clade 2.3.4.4b high pathogenicity avian influenza viruses in their original hosts and their spillover and spillback potential

## Variant Calling 

Use the variant_calling.sh to filter raw fastq files then map to the reference with minimap2 and finnaly obtain variants with clair3. 

## Variant Filtering and Annotation Workflow

This workflow identifies variants that are either **only present in the animal sample** or that **differ by at least 25% in frequency** between the inoculum and the animal sample.  
It then annotates them using **FluSurver** and merges the results for downstream interpretation.

---

### 1. Filter variants with `process_vcf`
Run `process_vcf` to compare the **animal VCF** against the **inoculum VCF**:

```bash
sbatch run_process_vcf.sh ./annotated_variants_mink2.vcf ./annotated_variants_inoc.vcf filtered_variants_mink2.txt filtered_variants_mink2.vcf
```

---

### 2. Compress and index the filtered VCF
```bash
bgzip filtered_variants_mink2.vcf
bcftools index filtered_variants_mink2.vcf.gz
```

---

### 3. Create a consensus FASTA
Generate a consensus FASTA containing all the variants (`-s -` flag):

```bash
bcftools consensus -f ../../reference/gulls_isidore.fasta -s - -o qsp_mink2.fasta filtered_variants_mink2.vcf.gz
```

---

### 4. Check variants with FluSurver
1. Upload the **consensus FASTA** to https://flusurver.bii.a-star.edu.sg/  
2. Download the results from FluSurver (replace the URL with your own result link):

```bash
wget https://flusurver.bii.a-star.edu.sg/tmp/flusurver_result42593.txt -O flusurver_mink2.txt
```

---

### 5. Merge FluSurver output with filtered variants
Use Python (via an Anaconda environment) to merge **FluSurver annotations** with the filtered variant list.

```python
import pandas as pd

# File paths
flusurver_mink2_path = "path_to/flusurver_mink2.txt"
filtered_variants_mink2_path = "path_to/filtered_variants_mink2.txt"

# Read input files
flusurver_mink2_data = pd.read_csv(flusurver_mink2_path, sep="\t", skiprows=1)
filtered_variants_mink2_data = pd.read_csv(filtered_variants_mink2_path, sep="\t")

# Merge datasets
merged_mink2_data = pd.merge(
    filtered_variants_mink2_data,
    flusurver_mink2_data,
    left_on=['CHROM', 'PROTEIN_CHANGE_animal'],
    right_on=['Query', 'Mutation'],
    how='left'
)

# Handle missing CHROM values
merged_mink2_data['CHROM'] = merged_mink2_data['CHROM'].replace('Missing_Chrom', pd.NA)
merged_mink2_data.loc[merged_mink2_data['CHROM'].isnull(), 'CHROM'] = 'NA'

# Save merged file
merged_mink2_data.to_csv("merged_mink2_data.txt", sep="\t", index=False)

# Preview results
print(merged_mink2_data.head())
```

---

### Notes
- Replace `"path_to/"` with your actual file paths.
- The merge keys (`CHROM` and `PROTEIN_CHANGE_animal` vs. `Query` and `Mutation`) must match the column names in your files.
- `skiprows=1` in `pd.read_csv` is used because the FluSurver file contains a header comment in the first row.

