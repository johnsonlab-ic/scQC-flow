library(SingleCellExperiment)
library(scDblFinder)
library(Matrix)
library(argparse)
library(Seurat)

# Parse command line arguments
parser <- ArgumentParser(description='Run scDblFinder doublet detection')
parser$add_argument('--h5_file', type='character', required=TRUE,
                    help='Path to H5 counts matrix file')
parser$add_argument('--sample_name', type='character', required=TRUE,
                    help='Sample name for output files')

args <- parser$parse_args()

# Set up paths
h5_file <- args$h5_file

cat("Loading data from H5 file:", h5_file, "\n")

# Check if H5 file exists
if (!file.exists(h5_file)) {
  stop("H5 file not found: ", h5_file)
}

# Load H5 data
counts <- Read10X_h5(h5_file)
cat("Loaded", ncol(counts), "cells and", nrow(counts), "features\n")

# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = counts))

# Add basic cell and gene info
sce$cell_barcode <- colnames(sce)
rownames(sce) <- rownames(counts)

cat("Created SingleCellExperiment object\n")
cat("Running scDblFinder doublet detection...\n")

# Run scDblFinder
sce <- scDblFinder(sce)

cat("scDblFinder analysis completed\n")

# Extract metadata
metadata <- colData(sce)
doublet_results <- data.frame(
  cell_barcode = metadata$cell_barcode,
  doublet_score = metadata$scDblFinder.score,
  doublet_class = metadata$scDblFinder.class,
  stringsAsFactors = FALSE
)

cat("Number of cells analyzed:", nrow(doublet_results), "\n")
cat("Number of doublets detected:", sum(doublet_results$doublet_class == "doublet"), "\n")
cat("Doublet rate:", round(mean(doublet_results$doublet_class == "doublet") * 100, 2), "%\n")

# Save results
write.csv(doublet_results, paste0(args$sample_name, "_scdbl_metrics.csv"), row.names = FALSE)

# Create summary file
doublet_count <- sum(doublet_results$doublet_class == "doublet")
singlet_count <- sum(doublet_results$doublet_class == "singlet")
doublet_rate <- round((doublet_count / nrow(doublet_results)) * 100, 2)

writeLines(c(
  paste("scDblFinder Summary for", args$sample_name),
  paste("====================================="),
  paste("Total cells analyzed:", nrow(doublet_results)),
  paste("Singlets detected:", singlet_count),
  paste("Doublets detected:", doublet_count),
  paste("Doublet rate:", paste0(doublet_rate, "%")),
  paste("Mean doublet score:", round(mean(doublet_results$doublet_score), 4)),
  paste("Median doublet score:", round(median(doublet_results$doublet_score), 4))
), paste0(args$sample_name, "_scdbl_summary.txt"))

cat("scDblFinder analysis completed successfully\n")
