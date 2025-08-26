library(SingleCellExperiment)
library(scDblFinder)
library(Matrix)
library(argparse)
library(Seurat)

# Parse command line arguments
parser <- ArgumentParser(description='Run scDblFinder doublet detection')
parser$add_argument('--mapping_dir', type='character', required=TRUE,
                    help='Path to Cell Ranger mapping directory')
parser$add_argument('--sample_name', type='character', required=TRUE,
                    help='Sample name for output files')

args <- parser$parse_args()

# Set up paths
data_path <- file.path(args$mapping_dir, "outs/filtered_feature_bc_matrix")

cat("Loading data from:", data_path, "\n")

# Check if data directory exists
if (!dir.exists(data_path)) {
  stop("Data directory not found: ", data_path)
}

# Load Cell Ranger data
counts <- Read10X(data.dir = data_path)
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
