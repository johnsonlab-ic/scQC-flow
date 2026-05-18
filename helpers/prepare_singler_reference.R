#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(data.table)
})

parser <- ArgumentParser(description = "Prepare a SingleR reference SCE")
parser$add_argument("--input_rds", type = "character", required = TRUE,
  help = "Input SingleCellExperiment reference RDS")
parser$add_argument("--output_rds", type = "character", required = TRUE,
  help = "Output cleaned SingleCellExperiment reference RDS")
parser$add_argument("--required_label_cols", type = "character", default = "",
  help = "Space-separated label columns that must be non-missing")
parser$add_argument("--symbol_col", type = "character", default = "symbol",
  help = "rowData column to use for symbols")
parser$add_argument("--gene_id_col", type = "character", default = "gene_id",
  help = "rowData column to use for gene ids")

opts <- parser$parse_args()

message("Loading input reference: ", opts$input_rds)
sce <- readRDS(opts$input_rds)
if (!inherits(sce, "SingleCellExperiment")) {
  stop("Input reference must be a SingleCellExperiment", call. = FALSE)
}
if (!"counts" %in% assayNames(sce)) {
  stop("Input reference must contain a 'counts' assay", call. = FALSE)
}

required_label_cols <- trimws(strsplit(opts$required_label_cols, "\\s+")[[1]])
required_label_cols <- required_label_cols[nzchar(required_label_cols)]
if (length(required_label_cols) > 0) {
  missing_label_cols <- setdiff(required_label_cols, colnames(colData(sce)))
  if (length(missing_label_cols) > 0) {
    stop(sprintf("Missing required label columns: %s", paste(missing_label_cols, collapse = ", ")), call. = FALSE)
  }
  keep_idx <- rep(TRUE, ncol(sce))
  for (label_col in required_label_cols) {
    vals <- as.character(colData(sce)[[label_col]])
    keep_idx <- keep_idx & !is.na(vals) & nzchar(vals)
  }
  sce <- sce[, keep_idx]
}

if (!"logcounts" %in% assayNames(sce)) {
  sce <- scuttle::logNormCounts(sce)
}

if (opts$gene_id_col %in% colnames(rowData(sce))) {
  gene_ids <- as.character(rowData(sce)[[opts$gene_id_col]])
} else {
  gene_ids <- rownames(sce)
}
if (opts$symbol_col %in% colnames(rowData(sce))) {
  symbols <- as.character(rowData(sce)[[opts$symbol_col]])
} else {
  symbols <- gene_ids
}

keep_genes <- !is.na(gene_ids) & nzchar(gene_ids)
gene_ids <- gene_ids[keep_genes]
symbols <- symbols[keep_genes]
sce <- sce[keep_genes, ]

rownames(sce) <- gene_ids
rowData(sce)$gene_id <- gene_ids
rowData(sce)$symbol <- ifelse(is.na(symbols) | !nzchar(symbols), gene_ids, symbols)

dir.create(dirname(opts$output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, opts$output_rds)

message(sprintf("Wrote cleaned SingleR reference: %s (%d genes x %d cells)", opts$output_rds, nrow(sce), ncol(sce)))
if (length(required_label_cols) > 0) {
  label_summary <- rbindlist(lapply(required_label_cols, function(label_col) {
    vals <- as.character(colData(sce)[[label_col]])
    data.table(label_col = label_col, label = vals)[, .(n_cells = .N), by = .(label_col, label)][order(label_col, -n_cells)]
  }))
  print(label_summary)
}