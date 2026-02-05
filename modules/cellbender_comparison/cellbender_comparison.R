#!/usr/bin/env Rscript

# CellBender vs Cell Ranger droplet comparison analysis script
# Compares droplet calling between Cell Ranger and CellBender
# Generates metrics CSV and knee-plot visualization

suppressPackageStartupMessages({
    library(argparse)
    library(Matrix)
    library(ggplot2)
    library(dplyr)
})

# Parse command-line arguments
parser <- ArgumentParser(description = "Compare droplet calling between Cell Ranger and CellBender")
parser$add_argument("--raw_h5", type = "character", required = TRUE,
                    help = "Path to Cell Ranger raw_feature_bc_matrix.h5")
parser$add_argument("--filtered_h5", type = "character", required = TRUE,
                    help = "Path to Cell Ranger filtered_feature_bc_matrix.h5")
parser$add_argument("--cellbender_h5", type = "character", required = FALSE, default = NULL,
                    help = "Path to CellBender output H5 (optional)")
parser$add_argument("--output_metrics", type = "character", required = TRUE,
                    help = "Output CSV file for droplet metrics")
parser$add_argument("--output_plot", type = "character", required = TRUE,
                    help = "Output PNG file for knee-plot")
parser$add_argument("--sample_name", type = "character", required = TRUE,
                    help = "Sample name for labeling")

args <- parser$parse_args()

cat("Starting droplet calling comparison analysis...\n")
cat("Sample:", args$sample_name, "\n")
cat("Raw H5:", args$raw_h5, "\n")
cat("Filtered H5:", args$filtered_h5, "\n")
if (!is.null(args$cellbender_h5)) {
    cat("CellBender H5:", args$cellbender_h5, "\n")
}

# Load H5 files using rhdf5
library(rhdf5)

# Function to read H5 and count barcodes
read_barcodes_from_h5 <- function(h5_path) {
    tryCatch({
        h5file <- h5open(h5_path)
        # Look for barcodes under different possible paths
        if (exists("barcodes", where = h5file)) {
            barcodes <- h5read(h5file, "barcodes")
        } else if (exists("matrix/barcodes", where = h5file)) {
            barcodes <- h5read(h5file, "matrix/barcodes")
        } else {
            # Try to find any barcodes dataset
            barcodes_list <- h5ls(h5file)
            barcode_indices <- grep("barcode", tolower(barcodes_list$name), ignore.case = TRUE)
            if (length(barcode_indices) > 0) {
                barcode_path <- barcodes_list$name[barcode_indices[1]]
                barcodes <- h5read(h5file, barcode_path)
            } else {
                stop("Could not find barcodes in H5 file")
            }
        }
        h5close(h5file)
        return(barcodes)
    }, error = function(e) {
        cat("Error reading H5:", conditionMessage(e), "\n")
        return(NULL)
    })
}

# Read barcodes from all sources
raw_barcodes <- read_barcodes_from_h5(args$raw_h5)
filtered_barcodes <- read_barcodes_from_h5(args$filtered_h5)
cellbender_barcodes <- if (!is.null(args$cellbender_h5)) {
    read_barcodes_from_h5(args$cellbender_h5)
} else {
    NULL
}

# Count droplets
n_raw_droplets <- length(raw_barcodes)
n_cellranger_cells <- length(filtered_barcodes)
n_cellbender_cells <- if (!is.null(cellbender_barcodes)) length(cellbender_barcodes) else NA

cat("\n=== Droplet Counts ===\n")
cat("Raw droplets (all):", n_raw_droplets, "\n")
cat("Cell Ranger called cells:", n_cellranger_cells, "\n")
if (!is.na(n_cellbender_cells)) {
    cat("CellBender called cells:", n_cellbender_cells, "\n")
}

# Create metrics data frame
metrics_df <- data.frame(
    metric = c("Total Raw Droplets", "CellRanger Called Cells", "CellBender Called Cells"),
    count = c(n_raw_droplets, n_cellranger_cells, n_cellbender_cells),
    percentage_of_raw = c(100.0, 
                          (n_cellranger_cells / n_raw_droplets) * 100,
                          if (!is.na(n_cellbender_cells)) (n_cellbender_cells / n_raw_droplets) * 100 else NA)
)

cat("\n=== Metrics Summary ===\n")
print(metrics_df)

# Write metrics to CSV
write.csv(metrics_df, file = args$output_metrics, row.names = FALSE)
cat("\nMetrics saved to:", args$output_metrics, "\n")

# ============================================================================
# Create Knee-Plot Comparison
# ============================================================================

# Read counts for UMI/rank analysis
read_counts_from_h5 <- function(h5_path) {
    tryCatch({
        h5file <- h5open(h5_path)
        
        # Read data matrix
        if (exists("data", where = h5file)) {
            data <- h5read(h5file, "data")
        } else if (exists("matrix/data", where = h5file)) {
            data <- h5read(h5file, "matrix/data")
        } else {
            data_list <- h5ls(h5file)
            data_indices <- grep("data", tolower(data_list$name), ignore.case = TRUE)
            if (length(data_indices) > 0) {
                data_path <- data_list$name[data_indices[1]]
                data <- h5read(h5file, data_path)
            } else {
                stop("Could not find data in H5 file")
            }
        }
        
        h5close(h5file)
        return(data)
    }, error = function(e) {
        cat("Error reading counts:", conditionMessage(e), "\n")
        return(NULL)
    })
}

# Try to build UMI rank data for knee plot
tryCatch({
    cat("\nGenerating knee-plot data...\n")
    
    # Read raw counts matrix to get UMI per barcode
    raw_counts <- read_counts_from_h5(args$raw_h5)
    
    if (!is.null(raw_counts)) {
        # Sum UMIs per barcode (column sum)
        if (is.list(raw_counts)) {
            # Sparse matrix format
            umi_per_barcode <- colSums(as.matrix(raw_counts))
        } else if (nrow(raw_counts) != length(raw_barcodes)) {
            # Transpose if needed
            umi_per_barcode <- rowSums(raw_counts)
        } else {
            umi_per_barcode <- colSums(raw_counts)
        }
        
        # Sort in decreasing order for knee plot
        umi_sorted <- sort(umi_per_barcode, decreasing = TRUE)
        rank <- seq_along(umi_sorted)
        
        knee_data <- data.frame(
            rank = rank,
            umi_count = umi_sorted,
            source = "All Droplets"
        )
        
        # Mark Cell Ranger cells
        cr_mask <- names(umi_sorted) %in% filtered_barcodes
        knee_data$source[cr_mask] <- "CellRanger Called"
        
        # Mark CellBender cells if available
        if (!is.null(cellbender_barcodes)) {
            cb_mask <- names(umi_sorted) %in% cellbender_barcodes
            knee_data$source[cb_mask] <- "CellBender Called"
        }
        
        # Create knee plot
        p <- ggplot(knee_data, aes(x = rank, y = umi_count, color = source)) +
            geom_point(size = 1, alpha = 0.6) +
            scale_x_log10(name = "Barcode Rank (log10)") +
            scale_y_log10(name = "UMI Count (log10)") +
            scale_color_manual(
                values = c("All Droplets" = "#CCCCCC", 
                          "CellRanger Called" = "#2E86AB",
                          "CellBender Called" = "#A23B72"),
                na.value = "#CCCCCC"
            ) +
            ggtitle(paste0("Droplet Calling Comparison - ", args$sample_name)) +
            theme_minimal() +
            theme(
                plot.title = element_text(face = "bold", size = 14),
                legend.position = "top",
                legend.title = element_blank(),
                panel.grid.major = element_line(color = "#E8E8E8"),
                panel.grid.minor = element_blank()
            )
        
        ggsave(args$output_plot, plot = p, width = 8, height = 6, dpi = 300)
        cat("Knee-plot saved to:", args$output_plot, "\n")
        
    } else {
        cat("Warning: Could not generate knee-plot (counts matrix not accessible)\n")
    }
    
}, error = function(e) {
    cat("Warning: Could not generate knee-plot:", conditionMessage(e), "\n")
})

cat("\nDroplet calling comparison analysis completed!\n")
