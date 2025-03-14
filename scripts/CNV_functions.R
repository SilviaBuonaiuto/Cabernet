#' Create a panel with multiple Manhattan plots aligned
#' 
#' @param sample1_manhattan Manhattan plot from first sample
#' @param sample2_manhattan Manhattan plot from second sample
#' @param sample1_id ID of first sample
#' @param sample2_id ID of second sample
#' @return A grid.arrange object with aligned Manhattan plots
create_manhattan_panel <- function(sample1_manhattan, sample2_manhattan, sample1_id, sample2_id) {
  # Create the panel using grid.arrange
  panel <- grid.arrange(
    sample1_manhattan$plot + ggtitle(paste("Manhattan Plot -", sample1_id)),
    sample2_manhattan$plot + ggtitle(paste("Manhattan Plot -", sample2_id)),
    ncol = 1  # Stack plots vertically
  )
  
  return(panel)
}
library(tidyr)
library(dplyr)
library(DNAcopy)
library(bracer)
library(magrittr)
library(gridExtra)
library(ggplot2)

#' Process a sample from SNP array data
#' 
#' @param sample_file Path to the sample result table
#' @param multi_sample_table Data frame containing multiple samples
#' @param separator Separator character in the sample file (default: ";")
#' @return A processed data frame with BAF values

process_sample <- function(sample_file, multi_sample_table, separator = ";") {
  # Extract sample ID from filename
  sample_id <- strsplit(basename(sample_file), ".", fixed = TRUE) %>% 
               sapply(extract2, 1)
  
  # Read the sample data
  sample_df <- read.table(sample_file, sep = separator, header = FALSE, skip = 11)
  colnames(sample_df)<-c("SNP.Name", "Sample.ID", "GC.Score", "Allele1", "Allele2")
  
  # Extract sample from multi-sample table
  gtype_col <- paste0(sample_id, ".GType")
  sample_subset <- multi_sample_table %>% 
                   select(Address, Chr, Position, Name, all_of(gtype_col))
  
  # Merge sample data
  merged_data <- merge(sample_subset, sample_df, by.x = "Name", by.y = "SNP.Name")
  
  # Filter by GC score and calculate BAF
  processed_data <- merged_data %>% 
                    filter(GC.Score >= 0.6) %>% 
                    mutate(BAF = ifelse(.data[[gtype_col]] == "AA", 0, 
                                        ifelse(.data[[gtype_col]] == "BB", 1, 0.5)))
  
  # Add the sample ID as an attribute for reference
  attr(processed_data, "sample_id") <- sample_id
  attr(processed_data, "gtype_col") <- gtype_col
  
  return(processed_data)
}

#' Calculate LRR for a sample
#' 
#' @param processed_data Processed sample data from process_sample()
#' @param idat_table Table containing idat data
#' @param lrr_col_name Name for the LRR column (default: "lrr")
#' @return Data frame with LRR values and filtered chromosomes
calculate_lrr <- function(processed_data, idat_table, lrr_col_name = "lrr") {
  # Extract sample ID
  sample_id <- attr(processed_data, "sample_id")
  
  # Normalize mean intensities and calculate LRR
  df <- idat_table %>% 
        mutate(NormalizedMean = Mean/median(Mean)) %>% 
        mutate(!!lrr_col_name := log2(NormalizedMean))
  
  # Merge processed data with LRR and filter out unwanted chromosomes
  final_data <- merge(processed_data, df, by.x = "Address", by.y = "ProbeID") %>% 
                filter(Chr != "0" & Chr != "MT" & Chr != "XY")
  
  return(final_data)
}

#' Perform segmentation on sample data
#' 
#' @param data Processed data with LRR values
#' @param lrr_col_name Name of the LRR column
#' @return List containing segmentation results
perform_segmentation <- function(data, lrr_col_name = "lrr") {
  # Create CNA object
  cna <- CNA(genomdat = data[[lrr_col_name]], 
             chrom = data$Chr, 
             maploc = data$Position, 
             data.type = "logratio")
  
  # Perform segmentation
  segment_cna <- segment(cna)
  segments <- segment_cna$output
  
  # Filter segments with more than 10 markers
  filtered_segments <- segments %>% 
                       filter(num.mark > 10) %>% 
                       mutate(abs_seg_mean = abs(seg.mean))
  
  # Format chromosome names
  filtered_segments$chrom <- paste0("chr", filtered_segments$chrom)
  filtered_segments$chrom <- factor(
    filtered_segments$chrom, 
    levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
               "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
               "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
               "chrX", "chrY")
  )
  
  return(list(
    cna = cna,
    segment_cna = segment_cna,
    segments = segments,
    filtered_segments = filtered_segments
  ))
}

#' Create segment plot
#' 
#' @param filtered_segments Filtered segments from perform_segmentation()
#' @param sample_id Sample identifier for the plot title
#' @return ggplot object
create_segment_plot <- function(filtered_segments, sample_id) {
  p <- ggplot(filtered_segments, 
              aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean)) + 
       geom_segment() + 
       geom_hline(yintercept = 0, color = "red") + 
       facet_wrap(~chrom, scales = "free_x") + 
       theme_bw() + 
       labs(title = "CNV Profile", x = "Genomic Position", y = "Log R Ratio") + 
       ggtitle(paste(sample_id))
  
  return(p)
}

#' Create Manhattan plot
#' 
#' @param filtered_segments Filtered segments from perform_segmentation()
#' @param sample_id Sample identifier for the plot title
#' @return List containing the plot and chromosome length data
create_manhattan_plot <- function(filtered_segments, sample_id) {
  # Calculate chromosome length
  chrom_lengths <- filtered_segments %>% 
                  mutate(length = loc.end - loc.start) %>% 
                  select(chrom, length) %>% 
                  group_by(chrom) %>% 
                  summarise(chr_len = sum(length))
  
  # Convert to numeric
  chrom_lengths$chr_len <- as.numeric(chrom_lengths$chr_len)
  
  # Calculate positions
  chrom_lengths <- chrom_lengths %>% 
                  mutate(cumulative_length = cumsum(chr_len)) %>% 
                  mutate(start_position = lag(cumulative_length, default = 0)) %>% 
                  mutate(center_position = (start_position + chr_len / 2))
  
  # Join with segments
  myf <- left_join(filtered_segments, chrom_lengths, by = "chrom")
  myf <- myf %>% mutate(adjusted_position = start_position + (loc.start + loc.end) / 2)
  
  # Create plot
  p <- myf %>% 
       ggplot(aes(adjusted_position, seg.mean)) + 
       geom_point(aes(color = as.factor(chrom)), alpha = 0.8, size = 2) + 
       scale_x_continuous(label = myf$chrom, breaks = myf$center_position) + 
       scale_y_continuous(expand = c(0, 0), limits = c(-2.5, 2)) + 
       theme_bw() + 
       theme(legend.position = "none", 
             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 12), 
             axis.text.y = element_text(size = 12), 
             axis.title = element_text(size = 14), 
             plot.title = element_text(size = 14)) + 
       scale_color_manual(values = rep(c("#050C9C", "#FFA823"), 12)) + 
       labs(x = "Genomic Position", y = "Log R Ratio") + 
       ggtitle(paste(sample_id))
  
  return(list(
    plot = p,
    chrom_lengths = chrom_lengths,
    manhattan_data = myf
  ))
}

#' Compare two samples and calculate concordance metrics
#' 
#' @param sample1 Processed data for sample 1
#' @param sample2 Processed data for sample 2
#' @param sample1_id ID of sample 1
#' @param sample2_id ID of sample 2
#' @param sample1_segments Segments from sample 1
#' @param sample2_segments Segments from sample 2
#' @return List containing concordance metrics and comparison data
compare_samples <- function(sample1, sample2, sample1_id, sample2_id, 
                           sample1_segments, sample2_segments) {
  # Calculate genotype concordance
  merged_data <- merge(
    sample1 %>% select(Name, BAF), 
    sample2 %>% select(Name, BAF),
    by = "Name", 
    suffixes = c("_sample1", "_sample2")
  )
  
  # Calculate BAF correlation
  baf_correlation <- cor(merged_data$BAF_sample1, 
                         merged_data$BAF_sample2, 
                         use = "pairwise.complete.obs")
  
  # Compare segments by chromosome
  chromosomes <- unique(as.character(sample1_segments$chrom))
  segment_comparison <- data.frame()
  
  for (chr in chromosomes) {
    # Get segments for this chromosome
    sample1_chr_segs <- sample1_segments %>% filter(chrom == chr)
    sample2_chr_segs <- sample2_segments %>% filter(chrom == chr)
    
    # Calculate difference metrics (handling case where segment counts differ)
    # This is a simplified approach - in reality you might need more complex matching
    if (nrow(sample1_chr_segs) > 0 && nrow(sample2_chr_segs) > 0) {
      # Calculate mean segment values for each sample on this chromosome
      sample1_mean <- mean(sample1_chr_segs$seg.mean)
      sample2_mean <- mean(sample2_chr_segs$seg.mean)
      diff <- abs(sample1_mean - sample2_mean)
    } else {
      diff <- NA
    }
    
    # Add to comparison table
    segment_comparison <- rbind(segment_comparison, data.frame(
      chromosome = chr,
      sample1_segments = nrow(sample1_chr_segs),
      sample2_segments = nrow(sample2_chr_segs),
      mean_segment_difference = diff
    ))
  }
  
  # Calculate overall segment concordance
  segment_comparison$concordance_score <- 1 - (segment_comparison$mean_segment_difference / 2)
  
  # Create summary report
  overall_concordance <- (baf_correlation + mean(segment_comparison$concordance_score, na.rm = TRUE)) / 2
  
  # Add concordance quality assessment
  concordance_quality <- case_when(
    overall_concordance >= 0.95 ~ "Excellent",
    overall_concordance >= 0.90 ~ "Very Good",
    overall_concordance >= 0.85 ~ "Good",
    overall_concordance >= 0.80 ~ "Moderate",
    overall_concordance >= 0.70 ~ "Poor",
    TRUE ~ "Very Poor"
  )

  summary_report <- data.frame(
    sample1 = sample1_id,
    sample2 = sample2_id,
    baf_correlation = baf_correlation,
    mean_segment_concordance = mean(segment_comparison$concordance_score, na.rm = TRUE),
    overall_concordance = (baf_correlation + mean(segment_comparison$concordance_score, na.rm = TRUE)) / 2,
    concordance_quality = concordance_quality
  )
  
  # Return results
  return(list(
    merged_data = merged_data,
    baf_correlation = baf_correlation,
    segment_comparison = segment_comparison,
    summary_report = summary_report
  ))
}

#' Create concordance plots for visual comparison
#' 
#' @param comparison_results Results from compare_samples()
#' @param sample1_id ID of sample 1
#' @param sample2_id ID of sample 2
#' @return List of plot objects
create_concordance_plots <- function(comparison_results, sample1_id, sample2_id) {
  # BAF Correlation plot
  baf_plot <- ggplot(comparison_results$merged_data, 
                    aes(x = BAF_sample1, y = BAF_sample2)) +
              geom_point(alpha = 0.3) +
              geom_smooth(method = "lm", color = "red") +
              theme_bw() +
              labs(title = paste("BAF Correlation (r =", 
                                round(comparison_results$baf_correlation, 3), ")"),
                   x = paste("BAF -", sample1_id),
                   y = paste("BAF -", sample2_id))
  
  # Segment difference plot
  segment_plot <- ggplot(comparison_results$segment_comparison, 
                        aes(x = chromosome, y = mean_segment_difference)) +
                  geom_bar(stat = "identity", fill = "steelblue") +
                  theme_bw() +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                  labs(title = "Segment Difference by Chromosome",
                       x = "Chromosome",
                       y = "Mean |Segment Difference|")
  
  return(list(
    baf_plot = baf_plot,
    segment_plot = segment_plot
  ))
}