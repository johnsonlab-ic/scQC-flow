#!/usr/bin/env Rscript
# run_integration.R
#
# Seurat + Harmony integration mirroring scprocess's integration pipeline.
#
# Takes the HVG count matrix (hvg_counts.h5) from HVG_SELECTION and an
# optional user-provided metadata CSV, normalises (LogNormalize 10k), scales,
# and runs:
#   PCA -> Harmony (if batch vars supplied) -> FindNeighbors ->
#   FindClusters/Leiden (multiple resolutions) -> UMAP
#
# Outputs:
#   integration_dt.csv.gz  — embedding + UMAP1/2 + PCA/Harmony dims + cluster
#                            assignments per cell

suppressPackageStartupMessages({
  library(argparse)
  library(rhdf5)
  library(Matrix)
  library(Seurat)
  library(harmony)
  library(data.table)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

parser <- ArgumentParser(description = "Run Seurat + Harmony integration")
parser$add_argument("--hvg_h5", type = "character", required = TRUE,
  help = "HVG count matrix H5 from hvg_selection.py")
parser$add_argument("--metadata_csv", type = "character", default = NULL,
  help = "Path to metadata CSV (optional)")
parser$add_argument("--metadata_id_col", type = "character", default = "sample_id",
  help = "Column in metadata CSV that maps to sample IDs")
parser$add_argument("--metadata_vars", type = "character", default = "",
  help = "Space-separated metadata columns for Harmony (e.g. 'brainregion condition')")
parser$add_argument("--n_dims", type = "integer", default = 30,
  help = "Number of PCA dimensions")
parser$add_argument("--theta", type = "double", default = 2.0,
  help = "Harmony theta")
parser$add_argument("--leiden_res", type = "character", default = "0.3 0.5 1.0",
  help = "Space-separated Leiden resolutions")
parser$add_argument("--out_csv", type = "character", default = "integration_dt.csv.gz",
  help = "Output integration CSV.gz path")

args <- parser$parse_args()

hvg_h5          <- args$hvg_h5
metadata_csv    <- args$metadata_csv
metadata_id_col <- args$metadata_id_col
metadata_vars   <- if (nzchar(args$metadata_vars)) strsplit(trimws(args$metadata_vars), "\\s+")[[1]] else character(0)
n_dims          <- args$n_dims
theta           <- args$theta
res_ls          <- strsplit(trimws(args$leiden_res), "\\s+")[[1]]
out_csv         <- args$out_csv

cat("=== INTEGRATION ===\n")

# ---------------------------------------------------------------------------
# 1. Load HVG count matrix (CSC genes x cells)
# ---------------------------------------------------------------------------
cat(sprintf("Loading HVG matrix: %s\n", hvg_h5))

indptr     <- h5read(hvg_h5, "matrix/indptr")
indices    <- h5read(hvg_h5, "matrix/indices")
vals       <- h5read(hvg_h5, "matrix/data")
shape      <- h5read(hvg_h5, "matrix/shape")
features   <- h5read(hvg_h5, "matrix/features/name")
barcodes   <- h5read(hvg_h5, "matrix/barcodes")
sample_ids <- h5read(hvg_h5, "matrix/sample_ids")

# build CSC sparse matrix (genes x cells)
mat_csc <- sparseMatrix(
  i    = as.integer(indices) + 1L,
  p    = as.integer(indptr),
  x    = as.numeric(vals),
  dims = as.integer(shape)
)
colnames(mat_csc) <- barcodes
rownames(mat_csc) <- features
cat(sprintf("  %d genes x %d cells\n", nrow(mat_csc), ncol(mat_csc)))

# ---------------------------------------------------------------------------
# 2. Create Seurat object
# ---------------------------------------------------------------------------
seu <- CreateSeuratObject(counts = mat_csc, project = "integration")
seu$cell_id   <- barcodes
seu$sample_id <- sample_ids

# ---------------------------------------------------------------------------
# 3. Join metadata CSV on sample_id
# ---------------------------------------------------------------------------
if (!is.null(metadata_csv) && length(metadata_vars) > 0) {
  cat(sprintf("Loading metadata: %s\n", metadata_csv))
  meta <- fread(metadata_csv)

  if (!metadata_id_col %in% names(meta)) {
    stop(sprintf("metadata_id_col '%s' not found in %s", metadata_id_col, metadata_csv))
  }
  setnames(meta, metadata_id_col, "sample_id")
  meta[, sample_id := as.character(sample_id)]

  missing_vars <- setdiff(metadata_vars, names(meta))
  if (length(missing_vars) > 0) {
    stop(sprintf("metadata_vars not found in CSV: %s", paste(missing_vars, collapse = ", ")))
  }

  keep_cols <- c("sample_id", metadata_vars)
  meta      <- unique(meta[, ..keep_cols])

  # join to Seurat metadata
  seu_meta  <- data.table(cell_id = seu$cell_id, sample_id = seu$sample_id)
  seu_meta  <- merge(seu_meta, meta, by = "sample_id", all.x = TRUE)
  setkey(seu_meta, cell_id)
  for (v in metadata_vars) {
    seu@meta.data[[v]] <- seu_meta[colnames(seu)][[v]]
  }
  cat(sprintf("  Metadata vars joined: %s\n", paste(metadata_vars, collapse = ", ")))
} else {
  cat("  No metadata vars — running without Harmony batch correction\n")
  metadata_vars <- character(0)
}

# ---------------------------------------------------------------------------
# 4. Normalise (LogNormalize = CPM10k + log1p)
# ---------------------------------------------------------------------------
cat("Normalising (LogNormalize, scale.factor=10000) ...\n")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4,
  verbose = FALSE)

# ---------------------------------------------------------------------------
# 5. Scale (all HVG features already selected)
# ---------------------------------------------------------------------------
cat("Scaling ...\n")
seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)

# ---------------------------------------------------------------------------
# 6. PCA
# ---------------------------------------------------------------------------
cat(sprintf("PCA (n_dims=%d) ...\n", n_dims))
seu <- RunPCA(seu, npcs = n_dims, features = rownames(seu), verbose = FALSE)

# ---------------------------------------------------------------------------
# 7. Harmony (if metadata_vars provided and >1 batch)
# ---------------------------------------------------------------------------
n_batches  <- length(unique(seu$sample_id))
sel_embed  <- "pca"

if (length(metadata_vars) > 0 && n_batches > 1) {
  harmony_key <- if (length(metadata_vars) == 1) metadata_vars else metadata_vars
  cat(sprintf("Harmony integration (key=%s, theta=%.1f) ...\n",
    paste(harmony_key, collapse = "+"), theta))
  seu       <- RunHarmony(seu,
    group.by.vars   = harmony_key,
    theta           = rep(theta, length(harmony_key)),
    reduction       = "pca",
    reduction.save  = "harmony",
    verbose         = FALSE)
  sel_embed <- "harmony"
} else {
  cat("Skipping Harmony (single batch or no metadata vars)\n")
}

# ---------------------------------------------------------------------------
# 8. Neighbors -> Leiden -> UMAP
# ---------------------------------------------------------------------------
cat("Finding neighbors ...\n")
seu <- FindNeighbors(seu, reduction = sel_embed, dims = seq_len(n_dims), verbose = FALSE)

cat(sprintf("Leiden clustering (resolutions: %s) ...\n", paste(res_ls, collapse = " ")))
for (res in res_ls) {
  seu <- FindClusters(seu,
    resolution  = as.numeric(res),
    algorithm   = 4,        # Leiden
    random.seed = 42,
    verbose     = FALSE)
  # FindClusters stores active.ident and RNA_snn_res.<res> in meta.data
}

cat("UMAP ...\n")
seu <- RunUMAP(seu, reduction = sel_embed, dims = seq_len(n_dims), verbose = FALSE)

# ---------------------------------------------------------------------------
# 9. Build output data.frame
# ---------------------------------------------------------------------------
cat("Building output table ...\n")

out_dt <- data.table(
  embedding = if (sel_embed == "harmony") "harmony" else "pca",
  cell_id   = seu$cell_id,
  sample_id = seu$sample_id
)

# metadata vars
for (v in metadata_vars) {
  out_dt[[v]] <- seu@meta.data[[v]]
}

# cluster columns
for (res in res_ls) {
  col          <- paste0("RNA_snn_res.", res)
  leiden_col   <- paste0("leiden_", res)
  cl_vals      <- as.character(seu@meta.data[[col]])
  out_dt[[col]]        <- cl_vals
  out_dt[[leiden_col]] <- cl_vals
}

# UMAP
umap_mat       <- Embeddings(seu, "umap")
out_dt$UMAP1   <- umap_mat[, 1]
out_dt$UMAP2   <- umap_mat[, 2]

# PCA / Harmony dims
emb_mat <- Embeddings(seu, sel_embed)
prefix  <- if (sel_embed == "harmony") "hmny_pca" else "pca"
for (i in seq_len(ncol(emb_mat))) {
  out_dt[[sprintf("%s_%02d", prefix, i)]] <- emb_mat[, i]
}

# ---------------------------------------------------------------------------
# 10. Save
# ---------------------------------------------------------------------------
fwrite(out_dt, file = out_csv, compress = "gzip")
cat(sprintf("Written integration results: %s  (%d cells)\n", out_csv, nrow(out_dt)))
cat("=== INTEGRATION done ===\n")
