library(tidyr)
library(dplyr)
library(DNAcopy)
library(bracer)
library(magrittr)
library(gridExtra)
library(ggplot2)

# Source the functions file
source("/path/to/CNV_functions.R")

# Parse command-line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop('specify seven arguments: 
       <embryonic_tissue_file> 
       <extra_embryonic_tissue_file> 
       <multi_sample_table> 
       <embryonic_idat_file> 
       <extra_embryonic_idat_file>
       <output_directory>
       <prefix_for_output_files>', call.=FALSE)
}

# Assign arguments to named variables for clarity
embryonic_file <- args[1]
extra_embryonic_file <- args[2]
multi_sample_path <- args[3]
embryonic_idat_path <- args[4]
extra_embryonic_idat_path <- args[5]
output_dir <- args[6]
output_prefix <- args[7]

# Ensure output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read shared table once
multi_sample_table <- read.table(multi_sample_path, sep = "\t", header = TRUE, dec = ",", check.names = FALSE)

# Read idat tables
embryonic_idat <- read.table(embryonic_idat_path, header = TRUE)
extra_embryonic_idat <- read.table(extra_embryonic_idat_path, header = TRUE)

# ===== PROCESS EMBRYONIC SAMPLE =====
cat("Processing embryonic sample...\n")

# Step 1: Process sample data
embryonic_data <- process_sample(embryonic_file, multi_sample_table)
embryonic_id <- attr(embryonic_data, "sample_id")

# Step 2: Calculate LRR
embryonic_lrr <- calculate_lrr(embryonic_data, embryonic_idat, "embryonic_lrr")

# Step 3: Perform segmentation
embryonic_segmentation <- perform_segmentation(embryonic_lrr, "embryonic_lrr")

# Step 4: Create plots
embryonic_segment_plot <- create_segment_plot(
  embryonic_segmentation$filtered_segments, 
  embryonic_id
)

embryonic_manhattan <- create_manhattan_plot(
  embryonic_segmentation$filtered_segments, 
  embryonic_id
)


# Save plots
ggsave(
  file.path(output_dir, paste0(output_prefix, "_embryonic_segments.png")), 
  plot = embryonic_segment_plot, 
  width = 15, 
  height = 12
)


#ggsave(
  #file.path(output_dir, paste0(output_prefix, "_embryonic_manhattan.png")), 
  #plot = embryonic_manhattan$plot, 
  #width = 18, 
  #height = 4
#)

# ===== PROCESS EXTRA-EMBRYONIC SAMPLE =====
cat("Processing extra-embryonic sample...\n")

# Step 1: Process sample data
extra_embryonic_data <- process_sample(extra_embryonic_file, multi_sample_table)
extra_embryonic_id <- attr(extra_embryonic_data, "sample_id")

# Step 2: Calculate LRR
extra_embryonic_lrr <- calculate_lrr(extra_embryonic_data, extra_embryonic_idat, "extra_embryonic_lrr")

# Step 3: Perform segmentation
extra_embryonic_segmentation <- perform_segmentation(extra_embryonic_lrr, "extra_embryonic_lrr")

# Step 4: Create plots
extra_embryonic_segment_plot <- create_segment_plot(
  extra_embryonic_segmentation$filtered_segments, 
  extra_embryonic_id
)

extra_embryonic_manhattan <- create_manhattan_plot(
  extra_embryonic_segmentation$filtered_segments, 
  extra_embryonic_id
)


# Save plots
ggsave(
  file.path(output_dir, paste0(output_prefix, "_extra_embryonic_segments.png")), 
  plot = extra_embryonic_segment_plot, 
  width = 15, 
  height = 12
)


#ggsave(
  #file.path(output_dir, paste0(output_prefix, "_extra_embryonic_manhattan.png")), 
  #plot = extra_embryonic_manhattan$plot, 
  #width = 18, 
  #height = 4
#)

# ===== COMPARE SAMPLES =====
cat("Comparing samples...\n")

# Perform comparison
comparison_results <- compare_samples(
  embryonic_data,
  extra_embryonic_data,
  embryonic_id,
  extra_embryonic_id,
  embryonic_segmentation$filtered_segments,
  extra_embryonic_segmentation$filtered_segments
)

# Create comparison plots
concordance_plots <- create_concordance_plots(
  comparison_results,
  embryonic_id,
  extra_embryonic_id
)

# Create combined Manhattan plot
manhattan_panel <- create_manhattan_panel(
  embryonic_manhattan,
  extra_embryonic_manhattan,
  embryonic_id,
  extra_embryonic_id
)
# Save combined Manhattan plot
ggsave(
  file.path(output_dir, paste0(output_prefix, "_manhattan_panel.png")), 
  plot = manhattan_panel, 
  width = 18, 
  height = 8
)

# Save comparison results
write.csv(
  comparison_results$summary_report,
  file.path(output_dir, paste0(output_prefix, "_concordance_summary.csv")),
  row.names = FALSE
)

write.csv(
  comparison_results$segment_comparison,
  file.path(output_dir, paste0(output_prefix, "_segment_comparison.csv")),
  row.names = FALSE
)

# Save comparison plots
ggsave(
  file.path(output_dir, paste0(output_prefix, "_baf_correlation.png")), 
  plot = concordance_plots$baf_plot, 
  width = 8, 
  height = 6
)

ggsave(
  file.path(output_dir, paste0(output_prefix, "_segment_difference.png")), 
  plot = concordance_plots$segment_plot, 
  width = 10, 
  height = 6
)

# Create a combined report figure
combined_plot <- grid.arrange(
  concordance_plots$baf_plot,
  concordance_plots$segment_plot,
  ncol = 2
)

ggsave(
  file.path(output_dir, paste0(output_prefix, "_combined_report.png")), 
  plot = combined_plot, 
  width = 12, 
  height = 6
)

# Print summary to console
cat("\n===== CONCORDANCE ANALYSIS SUMMARY =====\n")
cat(paste("Samples:", embryonic_id, "vs", extra_embryonic_id, "\n"))
cat(paste("BAF Correlation:", round(comparison_results$baf_correlation, 4), "\n"))
cat(paste("Mean Segment Concordance:", 
          round(comparison_results$summary_report$mean_segment_concordance, 4), "\n"))
cat(paste("Overall Concordance:", 
          round(comparison_results$summary_report$overall_concordance, 4), "\n"))
cat("=====================================\n")
cat("Analysis completed successfully!\n")