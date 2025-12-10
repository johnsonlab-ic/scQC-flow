#!/usr/bin/env Rscript
# extract_modalities.R
# Extract Gene Expression and ATAC modalities from multiome H5 and write as separate H5 files
# GEX H5 is used for CellBender, ATAC H5 is used later for Seurat multiome object creation

suppressPackageStartupMessages({
    library(argparse)
    library(Seurat)
    library(DropletUtils)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Extract modalities from multiome H5")
parser$add_argument("--input_h5", type = "character", required = TRUE,
                    help = "Input multiome raw H5 file")
parser$add_argument("--output_gex_h5", type = "character", required = TRUE,
                    help = "Output Gene Expression H5 file")
parser$add_argument("--output_atac_h5", type = "character", required = TRUE,
                    help = "Output ATAC/Peaks H5 file")
parser$add_argument("--sample_name", type = "character", required = TRUE,
                    help = "Sample name for logging")

args <- parser$parse_args()

cat("=== Extracting Modalities from Multiome H5 ===\n")
cat("Sample:", args$sample_name, "\n")
cat("Input H5:", args$input_h5, "\n")
cat("Output GEX H5:", args$output_gex_h5, "\n")
cat("Output ATAC H5:", args$output_atac_h5, "\n\n")

# Read the multiome H5 - returns a list with Gene Expression and Peaks
cat("Reading multiome H5 file...\n")
counts <- Read10X_h5(args$input_h5)

# Check if it's a list (multiome) or single matrix
if (!is.list(counts)) {
    stop("Expected multiome data with multiple modalities, but got single matrix")
}

cat("Detected multiome data with modalities:", paste(names(counts), collapse = ", "), "\n")

# Extract Gene Expression modality
if ("Gene Expression" %in% names(counts)) {
    gex_counts <- counts$`Gene Expression`
    cat("Extracted Gene Expression modality:", nrow(gex_counts), "genes x", ncol(gex_counts), "barcodes\n")
} else {
    stop("No 'Gene Expression' modality found in H5 file")
}

# Extract ATAC/Peaks modality
if ("Peaks" %in% names(counts)) {
    atac_counts <- counts$Peaks
    cat("Extracted Peaks modality:", nrow(atac_counts), "peaks x", ncol(atac_counts), "barcodes\n")
} else {
    stop("No 'Peaks' modality found in H5 file")
}

# Write Gene Expression H5 using DropletUtils
cat("\nWriting Gene Expression H5 file...\n")
write10xCounts(
    path = args$output_gex_h5,
    x = gex_counts,
    type = "HDF5",
    genome = "GRCh38",
    version = "3",
    overwrite = TRUE
)
cat("Successfully wrote Gene Expression H5 to:", args$output_gex_h5, "\n")

# Write ATAC/Peaks H5 using DropletUtils
cat("\nWriting ATAC/Peaks H5 file...\n")
write10xCounts(
    path = args$output_atac_h5,
    x = atac_counts,
    type = "HDF5",
    genome = "GRCh38",
    version = "3",
    overwrite = TRUE
)
cat("Successfully wrote ATAC/Peaks H5 to:", args$output_atac_h5, "\n")

cat("\n=== Modality extraction complete ===\n")
