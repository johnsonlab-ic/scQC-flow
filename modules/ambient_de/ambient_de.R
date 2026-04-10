#!/usr/bin/env Rscript
# ambient_de.R
#
# Differential expression analysis on empty droplets vs. cells to identify
# ambient RNA genes, following the scprocess mlm_hvgs workflow.
#
# Takes:
#   - Per-sample raw H5 files (barcode_matrix_*.h5) with all droplets
#   - Per-sample knee CSV files (knee_plot_data_*.csv) with empty/cell boundaries
#   - Per-sample decontX-filtered H5 files (filt_counts_*.h5) for cell barcodes
#   - GTF file for gene annotation (SYMBOL_ENSEMBL mapping)
#
# Builds:
#   - Pseudobulk matrix for empty droplets (genes x samples)
#   - Pseudobulk matrix for cells (genes x samples)
#
# Outputs:
#   - edger_dt.csv.gz: Differential expression results (gene_id, log2FC, p-value, etc.)
#   - pb_empties.rds: SummarizedExperiment with pseudobulk empty droplet counts
#
# Requirements:
#   - library(argparse)
#   - library(data.table)
#   - library(Matrix)
#   - library(rhdf5)
#   - library(SummarizedExperiment)
#   - library(edgeR)

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(Matrix)
    library(rhdf5)
    library(SummarizedExperiment)
    library(edgeR)
})

# =============================================================================
# H5 I/O Helpers
# =============================================================================

.read_h5_csc <- function(h5_path) {
    # Read a stacked S+U+A H5 into a CSC sparse matrix
    tryCatch({
        indptr   <- rhdf5::h5read(h5_path, "matrix/indptr")
        indices  <- rhdf5::h5read(h5_path, "matrix/indices")
        data     <- rhdf5::h5read(h5_path, "matrix/data")
        features <- rhdf5::h5read(h5_path, "matrix/features/name")
        barcodes <- rhdf5::h5read(h5_path, "matrix/barcodes")
        shape    <- rhdf5::h5read(h5_path, "matrix/shape")
        rhdf5::h5closeAll()
        
        mat <- Matrix::sparseMatrix(
            i = indices + 1L,  # Convert 0-indexed to 1-indexed
            p = indptr,
            x = as.numeric(data),
            dims = as.integer(shape),
            giveCsparse = TRUE
        )
        return(list(mat = mat, features = features, barcodes = barcodes))
    }, error = function(e) {
        stop(sprintf("Error reading H5 file %s: %s", h5_path, e$message))
    })
}

.sum_sua <- function(mat_csc, features) {
    # Sum spliced (_S), unspliced (_U), and ambiguous (_A) rows.
    # Returns list(summed_csc, gene_ids) where gene_ids have the _S/_U/_A suffix removed.
    types <- c('_S', '_U', '_A')
    mats <- list()
    
    for (t in types) {
        idx <- which(grepl(sprintf("%s$", t), features))
        if (length(idx) > 0) {
            mats[[t]] <- mat_csc[idx, ]
        }
    }
    
    # Get gene names from _S rows
    s_idx <- which(grepl("_S$", features))
    genes <- sub(r"(_[SUA]$)", "", features[s_idx])
    
    summed <- mats[['_S']]
    if (!is.null(mats[['_U']])) {
        summed <- summed + mats[['_U']]
    }
    if (!is.null(mats[['_A']])) {
        summed <- summed + mats[['_A']]
    }
    
    return(list(summed = summed, genes = genes))
}

# =============================================================================
# GTF Parsing
# =============================================================================

.parse_gtf <- function(gtf_path) {
    # Parse GTF file and return data.table with columns:
    # ensembl_id, gene_name, gene_type
    cat("Parsing GTF:", gtf_path, "\n")
    
    # Decompress if needed
    if (endsWith(gtf_path, ".gz")) {
        gtf_path_tmp <- tempfile(fileext = ".gtf")
        system(sprintf("gzip -cd %s > %s", shQuote(gtf_path), shQuote(gtf_path_tmp)))
        on.exit(unlink(gtf_path_tmp), add = TRUE)
        gtf_path <- gtf_path_tmp
    }
    
    # Read GTF — filter for 'gene' features only
    gtf_dt <- fread(
        cmd = sprintf("grep '^[^#]' %s", shQuote(gtf_path)),
        sep = "\t",
        header = FALSE,
        select = c(1, 3, 9),
        col.names = c("seqname", "feature", "attributes")
    )
    gtf_dt <- gtf_dt[feature == "gene", ]
    
    # Parse attributes: extract ensembl_id, gene_name, gene_type
    parse_attr <- function(attrs, key) {
        match_str <- sprintf('%s "([^"]+)"', key)
        result <- sub(sprintf('^.*%s "([^"]+)".*$', key), "\\1", attrs)
        result[!grepl(match_str, attrs)] <- NA_character_
        return(result)
    }
    
    gtf_dt[, ensembl_id := parse_attr(attributes, "gene_id")]
    gtf_dt[, gene_name := parse_attr(attributes, "gene_name")]
    gtf_dt[, gene_type := parse_attr(attributes, "gene_type")]
    
    # Handle missing gene_type by checking for gene_biotype
    missing_type <- is.na(gtf_dt$gene_type) | gtf_dt$gene_type == ""
    if (any(missing_type)) {
        gtf_dt[missing_type, gene_type := parse_attr(attributes[missing_type], "gene_biotype")]
    }
    
    # Keep unique genes by ensembl_id (take first occurrence)
    gtf_dt <- gtf_dt[, .SD[1], by = ensembl_id]
    
    cat(sprintf("Parsed %d genes from GTF\n", nrow(gtf_dt)))
    if (anyNA(gtf_dt$gene_name)) {
        cat(sprintf("Warning: %d genes missing gene_name\n", sum(is.na(gtf_dt$gene_name))))
    }
    
    return(gtf_dt[, .(ensembl_id, gene_name, gene_type)])
}

# =============================================================================
# Main Workflow
# =============================================================================

# Parse command-line arguments
parser <- ArgumentParser(description = "Run differential expression on empty vs. cell droplets")
parser$add_argument("--raw_h5_pattern", type = "character", default = "barcode_matrix_*.h5",
                    help = "Pattern to match raw H5 files (default: barcode_matrix_*.h5)")
parser$add_argument("--knee_pattern", type = "character", default = "knee_plot_data_*.csv",
                    help = "Pattern to match knee CSV files (default: knee_plot_data_*.csv)")
parser$add_argument("--filt_h5_pattern", type = "character", default = "filt_counts_*.h5",
                    help = "Pattern to match filtered H5 files (default: filt_counts_*.h5)")
parser$add_argument("--gtf", type = "character", required = TRUE,
                    help = "Path to GTF file (required)")
parser$add_argument("--out_de", type = "character", default = "edger_dt.csv.gz",
                    help = "Output DE results (default: edger_dt.csv.gz)")
parser$add_argument("--out_pb_empties", type = "character", default = "pb_empties.rds",
                    help = "Output pseudobulk empties SE (default: pb_empties.rds)")

args <- parser$parse_args()

cat("=== AMBIENT_DE ===\n")
cat("raw_h5_pattern:", args$raw_h5_pattern, "\n")
cat("knee_pattern:  ", args$knee_pattern, "\n")
cat("filt_h5_pattern:", args$filt_h5_pattern, "\n")
cat("gtf:           ", args$gtf, "\n")
cat("out_de:        ", args$out_de, "\n")
cat("out_pb_empties:", args$out_pb_empties, "\n\n")

# ------------------------------------------------------------------
# Locate input files by pattern
# ------------------------------------------------------------------
raw_h5_files <- Sys.glob(args$raw_h5_pattern)
knee_files <- Sys.glob(args$knee_pattern)
filt_h5_files <- Sys.glob(args$filt_h5_pattern)

cat(sprintf("Found %d raw H5 files\n", length(raw_h5_files)))
cat(sprintf("Found %d knee CSV files\n", length(knee_files)))
cat(sprintf("Found %d filtered H5 files\n", length(filt_h5_files)))

if (length(raw_h5_files) == 0 || length(knee_files) == 0 || length(filt_h5_files) == 0) {
    stop("Not all required input files found. Check patterns in current directory.")
}

# Extract sample IDs from filenames
extract_sample_id <- function(fname) {
    # Try common patterns: *_<sampleId>.* or *_<sampleId>.csv
    basename_only <- basename(fname)
    
    # Pattern 1: barcode_matrix_<sampleId>.h5
    m <- regexpr("barcode_matrix_(.+?)\\.h5$", basename_only)
    if (m > 0) return(regmatches(basename_only, m, invert = FALSE))
    
    # Pattern 2: knee_plot_data_<sampleId>.csv
    m <- regexpr("knee_plot_data_(.+?)\\.csv$", basename_only)
    if (m > 0) return(regmatches(basename_only, m, invert = FALSE))
    
    # Pattern 3: filt_counts_<sampleId>.h5
    m <- regexpr("filt_counts_(.+?)\\.h5$", basename_only)
    if (m > 0) return(regmatches(basename_only, m, invert = FALSE))
    
    return(NA_character_)
}

# Build a mapping from sampleId to files
get_sample_id <- function(fname) {
    basename_only <- basename(fname)
    if (grepl("barcode_matrix_", basename_only)) {
        return(sub("^barcode_matrix_(.+?)\\.h5$", "\\1", basename_only))
    }
    if (grepl("knee_plot_data_", basename_only)) {
        return(sub("^knee_plot_data_(.+?)\\.csv$", "\\1", basename_only))
    }
    if (grepl("filt_counts_", basename_only)) {
        return(sub("^filt_counts_(.+?)\\.h5$", "\\1", basename_only))
    }
    stop("Cannot extract sample ID from filename:", fname)
}

raw_h5_map <- setNames(raw_h5_files, sapply(raw_h5_files, get_sample_id))
knee_map <- setNames(knee_files, sapply(knee_files, get_sample_id))
filt_h5_map <- setNames(filt_h5_files, sapply(filt_h5_files, get_sample_id))

all_sample_ids <- unique(c(names(raw_h5_map), names(knee_map), names(filt_h5_map)))
cat(sprintf("Found %d unique sample(s): %s\n\n", length(all_sample_ids), paste(all_sample_ids, collapse = ", ")))

# Parse GTF for gene annotation
gtf_dt <- .parse_gtf(args$gtf)

# ------------------------------------------------------------------
# Build symbol-to-ensembl mapping
# ------------------------------------------------------------------
sym_map <- gtf_dt[!is.na(gene_name), .(sym = gene_name, ens = ensembl_id)]
setkey(sym_map, sym)

# ------------------------------------------------------------------
# Per-sample pseudobulk loading
# ------------------------------------------------------------------
pb_empties_list <- list()
pb_cells_list <- list()
gene_union <- NULL

for (sample_id in all_sample_ids) {
    cat(sprintf("\n=== Processing sample: %s ===\n", sample_id))
    
    if (!(sample_id %in% names(raw_h5_map))) {
        warning(sprintf("No raw H5 found for %s, skipping", sample_id))
        next
    }
    if (!(sample_id %in% names(knee_map))) {
        warning(sprintf("No knee CSV found for %s, skipping", sample_id))
        next
    }
    if (!(sample_id %in% names(filt_h5_map))) {
        warning(sprintf("No filtered H5 found for %s, skipping", sample_id))
        next
    }
    
    # Load raw H5 and identify empty barcodes
    cat(sprintf("  Loading raw H5: %s\n", basename(raw_h5_map[sample_id])))
    raw_data <- .read_h5_csc(raw_h5_map[sample_id])
    raw_mat <- raw_data$mat
    raw_features <- raw_data$features
    raw_barcodes <- raw_data$barcodes
    
    # Load knee CSV to identify empty droplets
    cat(sprintf("  Loading knee CSV: %s\n", basename(knee_map[sample_id])))
    knee_dt <- fread(knee_map[sample_id])
    
    # Knee CSV should have columns like "in_empty_plateau" or "is_cell" or similar
    # Assuming structure from DropletUtils::barcodeRanks: low-rank barcodes are empty
    # We'll use barcodes marked as in_plateau=TRUE or is_cell=FALSE
    if ("in_empty_plateau" %in% names(knee_dt)) {
        empty_barcodes <- knee_dt[in_empty_plateau == TRUE, barcode]
    } else if ("is_cell" %in% names(knee_dt)) {
        empty_barcodes <- knee_dt[is_cell == FALSE, barcode]
    } else {
        # Fallback: assume the CSV has a 'barcode' column and some threshold
        # Most likely structure: use a low rank cutoff
        cat("    Warning: Cannot determine empty droplet column; using lowest-ranked barcodes.\n")
        knee_dt <- knee_dt[order(rank)]
        n_empty <- max(1, nrow(knee_dt) %/% 10)  # Bottom 10%
        empty_barcodes <- knee_dt[1:n_empty, barcode]
    }
    
    cat(sprintf("    Empty barcodes: %d\n", length(empty_barcodes)))
    
    # Extract empty barcode indices from raw H5
    empty_idx <- which(raw_barcodes %in% empty_barcodes)
    if (length(empty_idx) > 0) {
        raw_empty_mat <- raw_mat[, empty_idx]
        cat(sprintf("    Raw empty matrix: %d genes x %d barcodes\n", nrow(raw_empty_mat), ncol(raw_empty_mat)))
        
        # Sum S+U+A per gene
        sua_result <- .sum_sua(raw_empty_mat, raw_features)
        pb_empty_summed <- as.vector(Matrix::rowSums(sua_result$summed))
        names(pb_empty_summed) <- sua_result$genes
        
        pb_empties_list[[sample_id]] <- pb_empty_summed
    } else {
        cat("    WARNING: No empty barcodes found in raw H5\n")
    }
    
    # Load decontX-filtered H5 and sum across cells
    cat(sprintf("  Loading filtered H5: %s\n", basename(filt_h5_map[sample_id])))
    filt_data <- .read_h5_csc(filt_h5_map[sample_id])
    filt_mat <- filt_data$mat
    filt_features <- filt_data$features
    filt_barcodes <- filt_data$barcodes
    
    cat(sprintf("    Filtered matrix: %d genes x %d barcodes\n", nrow(filt_mat), ncol(filt_mat)))
    
    # Sum S+U+A per gene
    sua_result <- .sum_sua(filt_mat, filt_features)
    pb_cell_summed <- as.vector(Matrix::rowSums(sua_result$summed))
    names(pb_cell_summed) <- sua_result$genes
    
    pb_cells_list[[sample_id]] <- pb_cell_summed
    
    # Accumulate union of gene IDs
    if (is.null(gene_union)) {
        gene_union <- unique(c(names(pb_empty_summed), names(pb_cell_summed)))
    } else {
        gene_union <- unique(c(gene_union, names(pb_empty_summed), names(pb_cell_summed)))
    }
    
    cat(sprintf("    Gene union so far: %d genes\n", length(gene_union)))
}

if (length(pb_empties_list) == 0) {
    stop("No samples processed successfully.")
}

# ------------------------------------------------------------------
# Build aligned pseudobulk matrices (genes x samples)
# ------------------------------------------------------------------
cat(sprintf("\n=== Building pseudobulk matrices ===\n"))
cat(sprintf("Total unique genes: %d\n", length(gene_union)))

pb_empty_mat <- matrix(0, nrow = length(gene_union), ncol = length(pb_empties_list))
rownames(pb_empty_mat) <- gene_union
colnames(pb_empty_mat) <- names(pb_empties_list)

pb_cell_mat <- matrix(0, nrow = length(gene_union), ncol = length(pb_cells_list))
rownames(pb_cell_mat) <- gene_union
colnames(pb_cell_mat) <- names(pb_cells_list)

for (i in seq_along(names(pb_empties_list))) {
    sample_id <- names(pb_empties_list)[i]
    pb_genes <- names(pb_empties_list[[sample_id]])
    pb_empty_mat[pb_genes, i] <- pb_empties_list[[sample_id]]
}

for (i in seq_along(names(pb_cells_list))) {
    sample_id <- names(pb_cells_list)[i]
    pb_genes <- names(pb_cells_list[[sample_id]])
    pb_cell_mat[pb_genes, i] <- pb_cells_list[[sample_id]]
}

cat("Empty pseudobulk matrix:", paste(dim(pb_empty_mat), collapse = " x "), "\n")
cat("Cell pseudobulk matrix: ", paste(dim(pb_cell_mat), collapse = " x "), "\n")

# ------------------------------------------------------------------
# Save pseudobulk empty SummarizedExperiment
# ------------------------------------------------------------------
cat("\n=== Saving pseudobulk objects ===\n")
pb_empties_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = pb_empty_mat),
    colData = data.frame(sample_id = colnames(pb_empty_mat), row.names = colnames(pb_empty_mat))
)
saveRDS(pb_empties_se, args$out_pb_empties)
cat("Saved pseudobulk empties SE:", args$out_pb_empties, "\n")

# ------------------------------------------------------------------
# EdgeR differential expression: empty vs. cells
# ------------------------------------------------------------------
cat("\n=== Running EdgeR analysis ===\n")

# Combine counts: rbind(empties, cells) as groups
combined_counts <- cbind(pb_empty_mat, pb_cell_mat)
group <- factor(c(rep("empty", ncol(pb_empty_mat)), rep("cell", ncol(pb_cell_mat))))

# Create DGEList
dge <- edgeR::DGEList(counts = combined_counts, group = group)
cat("DGEList created:", nrow(dge), "genes x", ncol(dge), "samples\n")

# Filter lowly-expressed genes (require at least 1 CPM in at least 2 samples)
keep <- edgeR::filterByExpr(dge, min.count = 5)
dge <- dge[keep, ]
cat("After filtering:", nrow(dge), "genes\n")

# Calculate normalization factors
dge <- edgeR::calcNormFactors(dge)

# Estimate dispersions
dge <- edgeR::estimateDisp(dge)

# Fit GLM
design <- model.matrix(~group)
fit <- edgeR::glmQLFit(dge, design)

# Test: compare empty vs. cell (lfc = 0)
res <- edgeR::glmQLFTest(fit, coef = 2)

# Extract results
res_dt <- as.data.table(res$table, keep.rownames = TRUE)
setnames(res_dt, "rn", "gene_id")

# Add FDR
res_dt[, FDR := p.adjust(PValue, method = "BH")]

# Map gene IDs to SYMBOL_ENSEMBL format
res_dt[, ensembl_id := gene_id]
res_dt[, symbol := NA_character_]

# Lookup symbols from GTF
for (i in 1:nrow(res_dt)) {
    ens_id <- res_dt[i, ensembl_id]
    sym_match <- gtf_dt[ensembl_id == ens_id, gene_name]
    if (length(sym_match) > 0 && !is.na(sym_match)) {
        res_dt[i, symbol := sym_match]
    }
}

# Create SYMBOL_ENSEMBL format
res_dt[is.na(symbol), symbol_ensembl := ensembl_id]
res_dt[!is.na(symbol), symbol_ensembl := sprintf("%s_%s", symbol, ensembl_id)]

# Reorder columns for output
res_out <- res_dt[, .(gene_id, symbol_ensembl, logFC, logCPM, LR, PValue, FDR)]

# Sort by FDR, then by logFC (descending)
res_out <- res_out[order(FDR, -abs(logFC))]

# Write output
cat(sprintf("Writing %d DE genes to %s\n", nrow(res_out), args$out_de))
fwrite(res_out, file = args$out_de, sep = ",")
system(sprintf("gzip -f %s", args$out_de))
cat("Completed: .gz version written\n")

# Summary statistics
ambient_genes <- res_out[FDR < 0.05 & logFC > 0, ]
cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("Ambient genes (FDR<0.05, logFC>0): %d\n", nrow(ambient_genes)))
if (nrow(ambient_genes) > 0) {
    cat(sprintf("Top ambient gene: %s (logFC=%.2f)\n", ambient_genes[1, symbol_ensembl], ambient_genes[1, logFC]))
}

cat("AMBIENT_DE complete.\n")
