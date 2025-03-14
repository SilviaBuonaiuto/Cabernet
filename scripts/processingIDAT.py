import os
import subprocess
import pandas as pd

# Define the base directory for IDAT files
BASE_DIR = "/path/to/IDAT_files"

# Path to your sample information file (table containing sample names and corresponding IDs in sequencer)
# Change this to your actual file path
SAMPLE_INFO_FILE = "/path/to/sample_info.tsv"  

# Path to your R script
R_SCRIPT_PATH = "/path/to/processIDAT.R"

# Path to the Rscript executable in your conda environment
# Replace this with the actual path to Rscript in your conda environment
RSCRIPT_PATH = "/path/to/your/conda/env/bin/Rscript"

# Directory for output files
# Change this to your desired output directory
OUTPUT_DIR = "output_data"

# Create output directory if it doesn't exist
import os
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Read the sample information table
# Adjust sep parameter if your file is CSV or has different separator
sample_data = pd.read_csv(SAMPLE_INFO_FILE, sep='\t')

# Function to process a sample
def process_sample(sample_id, barcode, position):
    # Define file paths
    idat_dir = os.path.join(BASE_DIR, str(barcode))
    red_idat = os.path.join(idat_dir, f"{barcode}_{position}_Red.idat")
    green_idat = os.path.join(idat_dir, f"{barcode}_{position}_Grn.idat")
    output_file = f"{sample_id}_lrr.tsv"
    
    # Print information for verification
    print(f"Processing sample: {sample_id}")
    print(f"  Red IDAT: {red_idat}")
    print(f"  Green IDAT: {green_idat}")
    print(f"  Output: {output_file}")
    
    # Check if files exist
    if not os.path.exists(red_idat):
        print(f"  Error: Red IDAT file not found: {red_idat}")
        return False
    
    if not os.path.exists(green_idat):
        print(f"  Error: Green IDAT file not found: {green_idat}")
        return False
    
    # Run the R script with full path to Rscript from conda environment
    try:
        cmd = [RSCRIPT_PATH, R_SCRIPT_PATH, red_idat, green_idat, output_file]
        subprocess.run(cmd, check=True)
        print(f"  Successfully processed {sample_id}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  Error running R script: {e}")
        return False

# Verify R script exists
if not os.path.exists(R_SCRIPT_PATH):
    print(f"Error: R script not found at {R_SCRIPT_PATH}")
    exit(1)

# Process each sample
successful = 0
failed = 0

for _, row in sample_data.iterrows():
    if process_sample(row['Sample_ID'], row['SentrixBarcode_A'], row['SentrixPosition_A']):
        successful += 1
    else:
        failed += 1

print(f"\nProcessing complete. Successfully processed {successful} samples, {failed} failed.")