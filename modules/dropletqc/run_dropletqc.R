library(DropletQC)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description='Run DropletQC nuclear fraction analysis')
parser$add_argument('--bam_file', type='character', required=TRUE,
                    help='Path to BAM file')
parser$add_argument('--barcodes_file', type='character', required=TRUE,
                    help='Path to barcodes file (TSV or CSV)')
parser$add_argument('--sample_name', type='character', required=TRUE,
                    help='Sample name for output files')
parser$add_argument('--cores', type='integer', default=1,
                    help='Number of cores to use (default: 1)')

args <- parser$parse_args()

# Set up paths
bam_file <- args$bam_file
barcodes_file <- args$barcodes_file

cat("BAM file:", bam_file, "\n")
cat("Barcodes file:", barcodes_file, "\n")

# Check if files exist
if (!file.exists(bam_file)) {
  stop("BAM file not found: ", bam_file)
}
if (!file.exists(barcodes_file)) {
  stop("Barcodes file not found: ", barcodes_file)
}

# Run nuclear fraction analysis
cat("Running nuclear_fraction_tags with", args$cores, "cores...\n")
nf_result <- nuclear_fraction_tags(
  bam = bam_file,
  barcodes = barcodes_file,
  cores = args$cores
)

cat("Nuclear fraction analysis completed\n")
cat("Number of cells analyzed:", nrow(nf_result), "\n")
cat("Mean nuclear fraction:", round(mean(nf_result$nuclear_fraction, na.rm = TRUE), 3), "\n")
cat("Median nuclear fraction:", round(median(nf_result$nuclear_fraction, na.rm = TRUE), 3), "\n")

# Add cell barcodes as a column
nf_result$cell_barcode <- rownames(nf_result)

# Save results
write.csv(nf_result, paste0(args$sample_name, "_dropletqc_metrics.csv"), row.names = FALSE)

# Create summary file
summary_stats <- data.frame(
  metric = c("total_cells", "mean_nuclear_fraction", "median_nuclear_fraction", 
             "min_nuclear_fraction", "max_nuclear_fraction", "sd_nuclear_fraction"),
  value = c(
    nrow(nf_result),
    round(mean(nf_result$nuclear_fraction, na.rm = TRUE), 4),
    round(median(nf_result$nuclear_fraction, na.rm = TRUE), 4),
    round(min(nf_result$nuclear_fraction, na.rm = TRUE), 4),
    round(max(nf_result$nuclear_fraction, na.rm = TRUE), 4),
    round(sd(nf_result$nuclear_fraction, na.rm = TRUE), 4)
  )
)

writeLines(c(
  paste("DropletQC Summary for", args$sample_name),
  paste("====================================="),
  paste("Total cells analyzed:", nrow(nf_result)),
  paste("Mean nuclear fraction:", round(mean(nf_result$nuclear_fraction, na.rm = TRUE), 3)),
  paste("Median nuclear fraction:", round(median(nf_result$nuclear_fraction, na.rm = TRUE), 3)),
  paste("Standard deviation:", round(sd(nf_result$nuclear_fraction, na.rm = TRUE), 3)),
  paste("Min nuclear fraction:", round(min(nf_result$nuclear_fraction, na.rm = TRUE), 3)),
  paste("Max nuclear fraction:", round(max(nf_result$nuclear_fraction, na.rm = TRUE), 3)),
  paste("Cores used:", args$cores)
), paste0(args$sample_name, "_dropletqc_summary.txt"))

cat("DropletQC analysis completed successfully\n")
