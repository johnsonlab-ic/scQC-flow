#!/usr/bin/env Rscript

# CellBender vs Cell Ranger droplet comparison analysis script
# Compares droplet calling between Cell Ranger and CellBender
# Generates metrics CSV and knee-plot visualization

suppressPackageStartupMessages({
    library(argparse)
    library(Matrix)
    library(ggplot2)
    library(dplyr)
    library(Seurat)
    library(patchwork)
    library(scales)
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
parser$add_argument("--output_plot_cr", type = "character", required = FALSE, default = NULL,
                    help = "Output PNG file for Cell Ranger knee-plot")
parser$add_argument("--output_plot_cb", type = "character", required = FALSE, default = NULL,
                    help = "Output PNG file for CellBender knee-plot")
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

# Read H5 via Seurat (uses hdf5r under the hood)
read_matrix_from_h5 <- function(h5_path) {
    tryCatch({
        obj <- Seurat::Read10X_h5(h5_path)
        if (is.list(obj)) {
            if ("Gene Expression" %in% names(obj)) {
                return(obj[["Gene Expression"]])
            }
            return(obj[[1]])
        }
        return(obj)
    }, error = function(e) {
        cat("Error reading H5:", conditionMessage(e), "\n")
        return(NULL)
    })
}

read_barcodes_from_h5 <- function(h5_path) {
    mat <- read_matrix_from_h5(h5_path)
    if (is.null(mat)) {
        return(NULL)
    }
    return(colnames(mat))
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
percentage_of_raw <- function(count, total) {
    if (is.na(total) || total == 0) {
        return(NA)
    }
    return((count / total) * 100)
}

metrics_df <- data.frame(
    metric = c("Total Raw Droplets", "CellRanger Called Cells", "CellBender Called Cells"),
    count = c(n_raw_droplets, n_cellranger_cells, n_cellbender_cells),
    percentage_of_raw = c(
        if (n_raw_droplets > 0) 100.0 else NA,
        percentage_of_raw(n_cellranger_cells, n_raw_droplets),
        if (!is.na(n_cellbender_cells)) percentage_of_raw(n_cellbender_cells, n_raw_droplets) else NA
    )
)

cat("\n=== Metrics Summary ===\n")
print(metrics_df)

# Write metrics to CSV
write.csv(metrics_df, file = args$output_metrics, row.names = FALSE)
cat("\nMetrics saved to:", args$output_metrics, "\n")

# ============================================================================
# Create Knee-Plot Comparison
# ============================================================================

# Try to build UMI rank data for knee plot
tryCatch({
    cat("\nGenerating knee-plot data...\n")
    
    # Read raw counts matrix to get UMI per barcode
    raw_counts <- read_matrix_from_h5(args$raw_h5)
    
    if (!is.null(raw_counts)) {
        # Sum UMIs per barcode (column sum)
        umi_per_barcode <- Matrix::colSums(raw_counts)
        names(umi_per_barcode) <- colnames(raw_counts)
        
        # Sort in decreasing order for knee plot
        umi_sorted <- sort(umi_per_barcode, decreasing = TRUE)
        rank <- seq_along(umi_sorted)
        
        # Base data for all droplets
        knee_data <- data.frame(
            rank = rank,
            umi_count = umi_sorted,
            barcode = names(umi_sorted),
            stringsAsFactors = FALSE
        )
        
        # Create Cell Ranger knee plot
        knee_data_cr <- knee_data
        knee_data_cr$called <- ifelse(knee_data_cr$barcode %in% filtered_barcodes, 
                                       "Cell Ranger Called", "Background")
        
        p_cr <- ggplot(knee_data_cr, aes(x = rank, y = umi_count, color = called)) +
            geom_point(size = 0.8, alpha = 0.6) +
            scale_x_log10(name = "Barcode Rank (log10)", 
                         labels = scales::comma) +
            scale_y_log10(name = "UMI Count (log10)",
                         labels = scales::comma) +
            scale_color_manual(
                values = c("Background" = "#E0E0E0", 
                          "Cell Ranger Called" = "#2E86AB"),
                breaks = c("Cell Ranger Called", "Background")
            ) +
            ggtitle(paste0("Cell Ranger Droplet Calling - ", args$sample_name)) +
            theme_minimal(base_size = 14) +
            theme(
                plot.title = element_text(face = "bold", size = 16),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size = 12),
                panel.grid.major = element_line(color = "#F0F0F0"),
                panel.grid.minor = element_blank()
            )
        
        # Save Cell Ranger plot
        if (!is.null(args$output_plot_cr)) {
            ggsave(args$output_plot_cr, plot = p_cr, width = 8, height = 6, dpi = 300)
            cat("Cell Ranger knee-plot saved to:", args$output_plot_cr, "\n")
        }
        
        # Create CellBender knee plot if CellBender data available
        if (!is.null(cellbender_barcodes)) {
            knee_data_cb <- knee_data
            knee_data_cb$called <- ifelse(knee_data_cb$barcode %in% cellbender_barcodes, 
                                          "CellBender Called", "Background")
            
            p_cb <- ggplot(knee_data_cb, aes(x = rank, y = umi_count, color = called)) +
                geom_point(size = 0.8, alpha = 0.6) +
                scale_x_log10(name = "Barcode Rank (log10)",
                             labels = scales::comma) +
                scale_y_log10(name = "UMI Count (log10)",
                             labels = scales::comma) +
                scale_color_manual(
                    values = c("Background" = "#E0E0E0", 
                              "CellBender Called" = "#A23B72"),
                    breaks = c("CellBender Called", "Background")
                ) +
                ggtitle(paste0("CellBender Droplet Calling - ", args$sample_name)) +
                theme_minimal(base_size = 14) +
                theme(
                    plot.title = element_text(face = "bold", size = 16),
                    axis.title = element_text(size = 14),
                    axis.text = element_text(size = 12),
                    legend.position = "top",
                    legend.title = element_blank(),
                    legend.text = element_text(size = 12),
                    panel.grid.major = element_line(color = "#F0F0F0"),
                    panel.grid.minor = element_blank()
                )
            
            # Save CellBender plot
            if (!is.null(args$output_plot_cb)) {
                ggsave(args$output_plot_cb, plot = p_cb, width = 8, height = 6, dpi = 300)
                cat("CellBender knee-plot saved to:", args$output_plot_cb, "\n")
            }
            
            # Create combined plot (side-by-side)
            p_combined <- p_cr + p_cb + plot_layout(ncol = 2)
            ggsave(args$output_plot, plot = p_combined, width = 18, height = 10, dpi = 300)
            cat("Combined knee-plot saved to:", args$output_plot, "\n")
        } else {
            # If no CellBender, just save the Cell Ranger plot
            ggsave(args$output_plot, plot = p_cr, width = 10, height = 10, dpi = 300)
            cat("Knee-plot saved to:", args$output_plot, "\n")
        }
        
    } else {
        cat("Warning: Could not generate knee-plot (counts matrix not accessible)\n")
    }
    
}, error = function(e) {
    cat("Warning: Could not generate knee-plot:", conditionMessage(e), "\n")
})

cat("\nDroplet calling comparison analysis completed!\n")
