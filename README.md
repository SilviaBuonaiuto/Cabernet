# Cabernet
## Required tools
* R
* Python

## Analysis
### Process IDAT files 
IDAT files are raw intensity data files generated by Illumina SNP genotyping arrays that contain fluorescence intensity measurements for each probe on the array, serving as the primary output format before conversion to genotype calls. Each sample typically has two IDAT files: Red and Green, representing different fluorescence channels.
This script processes raw IDAT files from Illumina SNP genotyping arrays to extract Log R Ratio (LRR) values. The script identifies sample IDs and their corresponding barcodes, locates the appropriate IDAT files, and processes them using an R script.

#### 1. Setup
1. Create file containing Sample_ID and barcode IDs and positions
2. Define paths to IDAT files, table with IDs and [R script](https://github.com/SilviaBuonaiuto/Cabernet/blob/main/scripts/processIDAT.R) inside the python script

#### 2. Run [script](https://github.com/SilviaBuonaiuto/Cabernet/blob/main/scripts/processingIDAT.py) to process IDAT files 
```
python3 processingIDAT.py 
```


### Copy Number Variation (CNV) analysis 

Run R script to perform Copy Number Variation (CNV) analysis on paired embryonic and extra-embryonic tissue samples. It processes the raw data, calculates Log R Ratio (LRR), performs segmentation, and compares the samples to assess their concordance.

**Required command-line arguments :**
> `<embryonic_tissue_file>` - Path to the embryonic tissue data file (results table from the company)



## Concordance Scores Significance

### BAF Correlation
#### Significance:

* ≥0.95: Samples from same individual
* 0.5-0.7: Parent-child relationship
* <0.3: Unrelated individuals


**Interpretation:** Measures similarity of B-allele frequency patterns, reflecting genotype concordance


### Mean Segment Concordance
#### Significance:

* ≥0.9: Samples from same individual
* 0.7-0.85: Related individuals
* <0.6: Unrelated individuals


**Interpretation:** Measures similarity of copy number variations/segmentation patterns

### Overall Concordance
#### Significance:

* ≥0.9: Samples from same individual
* 0.6-0.8: Related individuals
* <0.5: Unrelated individuals


**Interpretation:** Combined metric that balances genotype and structural concordance


### Concordance Quality Categories

* **Excellent (≥0.95):** Strong evidence of same individual with no mosaicism
* **Very Good (≥0.90):** Same individual with minimal differences
* **Good (≥0.85):** Same individual with possible minor mosaicism
* **Moderate (≥0.80):** Some discordance, possible significant mosaicism
* **Poor (≥0.70):** Major differences, likely significant mosaicism
* **Very Poor (<0.70):** Different individuals or major technical issues