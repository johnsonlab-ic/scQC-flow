#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(scuttle)
  library(data.table)
})

parser <- ArgumentParser(description = "Build a query SCE for annotation workflows")
parser$add_argument("--integration_csv", type = "character", required = TRUE,
  help = "Path to integration_dt.csv.gz")
parser$add_argument("--genome_gtf", type = "character", required = TRUE,
  help = "Path to genome GTF")
parser$add_argument("--utils_r", type = "character", required = TRUE,
  help = "Path to export_utils.R")
parser$add_argument("--h5_pattern", type = "character", default = "filt_counts_*.h5",
  help = "Glob for filtered H5 files")
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

message("Loading integration metadata: ", opts$integration_csv)
obs_dt <- fread(opts$integration_csv)
obs_dt[, cell_id := as.character(cell_id)]
obs_dt <- unique(obs_dt, by = "cell_id")
obs_dt <- drop_export_bloat_cols(obs_dt)

message("Building count matrix from filtered H5 files")
counts_obj <- build_export_counts(opts$h5_pattern, obs_dt)
gtf_dt <- parse_gtf_annotations_export(opts$genome_gtf)
var_dt <- build_gene_metadata_export(counts_obj$gene_keys, gtf_dt)

obs_dt <- obs_dt[match(counts_obj$ordered_ids, cell_id)]
if (anyNA(obs_dt$cell_id)) {
  stop("Failed to align integration metadata to exported counts", call. = FALSE)
}

counts_mat <- counts_obj$counts
rownames(counts_mat) <- counts_obj$gene_keys
colnames(counts_mat) <- counts_obj$ordered_ids

sce <- SingleCellExperiment(
  assays = list(counts = counts_mat),
  colData = S4Vectors::DataFrame(as.data.frame(obs_dt, stringsAsFactors = FALSE))
)

rowData(sce)$gene_key <- var_dt$gene_key
rowData(sce)$gene_id <- var_dt$gene_id
rowData(sce)$symbol <- var_dt$symbol
rowData(sce)$gene_type <- var_dt$gene_type

if (all(c("UMAP1", "UMAP2") %in% names(obs_dt))) {
  umap_mat <- as.matrix(obs_dt[, .(UMAP1, UMAP2)])
  storage.mode(umap_mat) <- "double"
  reducedDim(sce, "X_umap") <- umap_mat
}

sce <- scuttle::logNormCounts(sce)

dir.create(dirname(opts$output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, opts$output_rds)

message(sprintf("Wrote annotation query SCE: %s (%d genes x %d cells)", opts$output_rds, nrow(sce), ncol(sce)))