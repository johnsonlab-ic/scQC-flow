#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(scuttle)
  library(data.table)
})

parser <- ArgumentParser(description = "Build a per-sample query SCE for annotation workflows")
parser$add_argument("--sample_id", type = "character", required = TRUE,
  help = "Sample identifier")
parser$add_argument("--integration_csv", type = "character", required = TRUE,
  help = "Path to integration_dt.csv.gz")
parser$add_argument("--utils_r", type = "character", required = TRUE,
  help = "Path to export_utils.R")
parser$add_argument("--h5_file", type = "character", required = TRUE,
  help = "Filtered H5 file for the sample")
parser$add_argument("--output_rds", type = "character", required = TRUE,
  help = "Output query SCE RDS path")

opts <- parser$parse_args()

source(opts$utils_r)

drop_export_bloat_cols <- function(dt) {
  bloat_cols <- grep("^(hmny_pca_|harmony_pca_)", names(dt), value = TRUE)
  if (length(bloat_cols) > 0) {
    dt[, (bloat_cols) := NULL]
  }
  dt
}

filter_annotation_cells <- function(dt) {
  keep_idx <- rep(TRUE, nrow(dt))

  if ("embedding" %in% names(dt)) {
    embedding_vals <- as.character(dt$embedding)
    keep_idx <- !is.na(embedding_vals) & nzchar(embedding_vals)
  } else if (all(c("UMAP1", "UMAP2") %in% names(dt))) {
    keep_idx <- !is.na(dt$UMAP1) & !is.na(dt$UMAP2)
  } else {
    cluster_cols <- grep("^RNA_snn_res\\.", names(dt), value = TRUE)
    if (length(cluster_cols) > 0) {
      keep_idx <- Reduce(`|`, lapply(cluster_cols, function(cluster_col) {
        cluster_vals <- as.character(dt[[cluster_col]])
        !is.na(cluster_vals) & nzchar(cluster_vals)
      }))
    } else if (all(c("is_dbl", "in_dbl_cl") %in% names(dt))) {
      keep_idx <- (!as.logical(dt$is_dbl)) & (!as.logical(dt$in_dbl_cl))
      keep_idx[is.na(keep_idx)] <- FALSE
    }
  }

  n_keep <- sum(keep_idx)
  if (n_keep == 0) {
    stop("No clean pass-2 integration cells available for annotation", call. = FALSE)
  }

  n_drop <- sum(!keep_idx)
  if (n_drop > 0) {
    message(sprintf(
      "Retaining %d clean integrated cells and dropping %d pass-1-only / filtered rows before annotation",
      n_keep,
      n_drop
    ))
  }

  dt[keep_idx]
}

message("Loading integration metadata: ", opts$integration_csv)
obs_dt <- fread(opts$integration_csv)
obs_dt[, sample_id := as.character(sample_id)]
obs_dt[, cell_id := as.character(cell_id)]
obs_dt <- obs_dt[sample_id == opts$sample_id]
obs_dt <- unique(obs_dt, by = "cell_id")
obs_dt <- drop_export_bloat_cols(obs_dt)
obs_dt <- filter_annotation_cells(obs_dt)

if (nrow(obs_dt) == 0) {
  stop(sprintf("No integration rows found for sample_id '%s'", opts$sample_id), call. = FALSE)
}

message("Building count matrix for sample: ", opts$sample_id)
counts_mat <- read_sparse_h5_matrix_export(opts$h5_file)
counts_mat <- sum_sua_counts_export(counts_mat)

global_ids <- paste(opts$sample_id, colnames(counts_mat), sep = "_")
index_map <- setNames(seq_along(global_ids), global_ids)
wanted_ids <- obs_dt$cell_id
missing_ids <- setdiff(wanted_ids, names(index_map))
if (length(missing_ids) > 0) {
  stop(sprintf(
    "Filtered H5 for sample %s is missing cells: %s",
    opts$sample_id,
    paste(head(missing_ids, 5), collapse = ", ")
  ))
}

sel_idx <- unname(index_map[wanted_ids])
counts_mat <- counts_mat[, sel_idx, drop = FALSE]
colnames(counts_mat) <- wanted_ids
gene_keys <- rownames(counts_mat)

obs_dt <- obs_dt[match(wanted_ids, cell_id)]
if (anyNA(obs_dt$cell_id)) {
  stop("Failed to align integration metadata to exported counts", call. = FALSE)
}

sce <- SingleCellExperiment(
  assays = list(counts = counts_mat),
  colData = S4Vectors::DataFrame(as.data.frame(obs_dt, stringsAsFactors = FALSE))
)

rowData(sce)$gene_key <- gene_keys

if (all(c("UMAP1", "UMAP2") %in% names(obs_dt))) {
  umap_mat <- as.matrix(obs_dt[, .(UMAP1, UMAP2)])
  storage.mode(umap_mat) <- "double"
  reducedDim(sce, "X_umap") <- umap_mat
}

sce <- scuttle::logNormCounts(sce)

dir.create(dirname(opts$output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, opts$output_rds)

message(sprintf("Wrote annotation query SCE for %s: %s (%d genes x %d cells)", opts$sample_id, opts$output_rds, nrow(sce), ncol(sce)))