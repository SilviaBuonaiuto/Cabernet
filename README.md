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
> `<embryonic_tissue_file>` - Path to the embryonic tissue data file (genotyping results from the SNP array)

> `<extra_embryonic_tissue_file>` - Path to the extra-embryonic tissue data file (genotyping results from the SNP array)

> `<multi_sample_table>` - Path to multi-sample reference table (reference info for SNP positions, quality and sample specific genotyping data)

> `<embryonic_idat_file>` - Path to embryonic tissue file containing LRR values created in previous step

> `<extra_embryonic_idat_file>` - Path to extra-embryonic tissue file containing LRR values created in previous step

> `<output_directory>` - Directory where output files should be saved

> `<prefix_for_output_files>` - Prefix to use for all output files


#### 1. Run Rscript
The [R script](https://github.com/SilviaBuonaiuto/Cabernet/blob/main/scripts/cnv.R) relies on custom functions defined in [CNV_functions.R](https://github.com/SilviaBuonaiuto/Cabernet/blob/main/scripts/CNV_functions.R)

```
Rscript cnv.R examples/sample1Embryo.csv examples/sample1Extra.csv examples/multi_sample.tsv examples/sample1Embryo_lrr.tsv examples/sample1Extra_lrr.tsv out_directory out_prefix
```
#### Output Files

The script generates the following output files in the specified `<output_directory>`:

##### Segment Plots 
* `<prefix>_embryonic_segments.png` - Segmentation plot for embryonic tissue	
* `<prefix>_extra_embryonic_segments.png` - Segmentation plot for extra-embryonic tissue

##### Manhattan Plots
* `<prefix>_manhattan_panel.png` - Combined Manhattan plot comparing both samples

##### Concordance Analysis
* `<prefix>_baf_correlation.png` - B-allele frequency correlation between samples
* `<prefix>_segment_difference.png` - Visual representation of segment differences
* `<prefix>_combined_report.png` - Combined visualization of concordance metrics

##### Data Files
* `<prefix>_concordance_summary.csv` - Summary statistics of concordance analysis
* `<prefix>_segment_comparison.csv` - Detailed comparison of segments between samples

## Results

Results can be found here :
[Concordance analysis table](https://github.com/SilviaBuonaiuto/Cabernet/blob/main/data/cabernet_results.tsv)
[Segment plots](https://github.com/SilviaBuonaiuto/Cabernet/tree/main/plots/segment_plots)
[Manhattan plots](https://github.com/SilviaBuonaiuto/Cabernet/tree/main/plots/manhattan_plots)

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