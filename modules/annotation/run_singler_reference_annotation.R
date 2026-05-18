#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(SingleR)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(data.table)
})

option_list <- list(
  make_option("--query_rds", type = "character", help = "Path to prepared query SCE RDS"),
  make_option("--reference_rds", type = "character", help = "Path to SingleR reference RDS"),
  make_option("--reference_label_col", type = "character", help = "Reference label column"),
  make_option("--method_id", type = "character", help = "Method identifier"),
  make_option("--reference_name", type = "character", help = "Reference display name"),
  make_option("--output_cells_csv", type = "character", help = "Per-cell output CSV.GZ"),
  make_option("--output_export_csv", type = "character", help = "Export metadata CSV.GZ"),
  make_option("--output_cluster_csv", type = "character", help = "Cluster summary CSV.GZ"),
  make_option("--ncores", type = "integer", default = 1L, help = "BiocParallel workers [default %default]"),
  make_option("--bp_type", type = "character", default = "multicore", help = "BiocParallel backend (multicore|snow) [default %default]"),
  make_option("--fine_tune", type = "character", default = "false", help = "Enable SingleR fine-tuning [default %default]"),
  make_option("--prune", type = "character", default = "true", help = "Enable SingleR pruning [default %default]")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

required <- c(
  "query_rds", "reference_rds", "reference_label_col", "method_id",
  "reference_name", "output_cells_csv", "output_export_csv", "output_cluster_csv"
)
missing <- required[vapply(required, function(key) is.null(opts[[key]]) || !nzchar(opts[[key]]), logical(1))]
if (length(missing) > 0) {
  print_help(parser)
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

choose_bp <- function(ncores, bp_type = "multicore") {
  if (ncores > 1) {
    if (bp_type == "multicore" && .Platform$OS.type != "windows") {
      return(BiocParallel::MulticoreParam(workers = ncores))
    }
    return(BiocParallel::SnowParam(workers = ncores))
  }
  BiocParallel::SerialParam()
}

align_genes <- function(query_sce, ref_sce) {
  common <- intersect(rownames(query_sce), rownames(ref_sce))
  if (length(common) >= 200) {
    return(list(query = query_sce[common, ], ref = ref_sce[common, ]))
  }

  query_uc <- toupper(rownames(query_sce))
  ref_uc <- toupper(rownames(ref_sce))
  common_uc <- intersect(query_uc, ref_uc)
  common_uc <- common_uc[!duplicated(common_uc)]
  if (length(common_uc) < 200) {
    stop("Too few overlapping genes between query and SingleR reference", call. = FALSE)
  }

  query_idx <- match(common_uc, query_uc)
  ref_idx <- match(common_uc, ref_uc)
  query_aligned <- query_sce[query_idx, ]
  ref_aligned <- ref_sce[ref_idx, ]
  rownames(query_aligned) <- rownames(ref_aligned)
  list(query = query_aligned, ref = ref_aligned)
}

sanitize_prefix <- function(x) {
  gsub("[^A-Za-z0-9_]", "_", x)
}

make_cluster_prediction_summary <- function(cells_dt) {
  cluster_cols <- grep("^(RNA_snn_res\\.|leiden_)", names(cells_dt), value = TRUE)
  if (length(cluster_cols) == 0) {
    return(data.table(
      cluster_col = character(),
      cluster_resolution = character(),
      cluster = character(),
      n_cells = integer(),
      top_label = character(),
      cs = numeric(),
      mean_score = numeric()
    ))
  }

  rbindlist(lapply(cluster_cols, function(cluster_col) {
    dt <- cells_dt[!is.na(get(cluster_col)) & nzchar(as.character(get(cluster_col))), .(
      cluster = as.character(get(cluster_col)),
      label = as.character(label),
      score = as.numeric(score)
    )]
    if (nrow(dt) == 0) {
      return(NULL)
    }

    counts_dt <- dt[, .N, by = .(cluster, label)]
    totals_dt <- dt[, .(
      n_cells = .N,
      mean_score = mean(score, na.rm = TRUE)
    ), by = cluster]

    summary_dt <- merge(counts_dt, totals_dt, by = "cluster", all.x = TRUE, sort = FALSE)
    summary_dt[, cs := N / n_cells]
    setorder(summary_dt, cluster, -cs, -N, label)
    summary_dt <- summary_dt[, .SD[1], by = cluster]
    summary_dt[, `:=`(
      cluster_col = cluster_col,
      cluster_resolution = sub("^(RNA_snn_res\\.|leiden_)", "", cluster_col),
      top_label = label
    )]
    summary_dt[, .(cluster_col, cluster_resolution, cluster, n_cells, top_label, cs, mean_score)]
  }), use.names = TRUE, fill = TRUE)
}

message("Loading query SCE: ", opts$query_rds)
query_sce <- readRDS(opts$query_rds)
if (!"logcounts" %in% assayNames(query_sce)) {
  query_sce <- scuttle::logNormCounts(query_sce)
}

message("Loading SingleR reference: ", opts$reference_rds)
ref_sce <- readRDS(opts$reference_rds)
if (!opts$reference_label_col %in% colnames(SummarizedExperiment::colData(ref_sce))) {
  stop(sprintf("SingleR reference is missing label column '%s'", opts$reference_label_col), call. = FALSE)
}
if (!"logcounts" %in% assayNames(ref_sce)) {
  ref_sce <- scuttle::logNormCounts(ref_sce)
}

aligned <- align_genes(query_sce, ref_sce)
query_sce <- aligned$query
ref_sce <- aligned$ref

labels <- SummarizedExperiment::colData(ref_sce)[[opts$reference_label_col]]
bp <- choose_bp(opts$ncores, opts$bp_type)

pred <- SingleR::SingleR(
  test = query_sce,
  ref = ref_sce,
  labels = labels,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  fine.tune = as_flag(opts$fine_tune),
  prune = as_flag(opts$prune),
  BPPARAM = bp
)

per_cell_labels <- pred$labels
per_cell_scores <- apply(pred$scores, 1, max, na.rm = TRUE)

cell_ids <- colnames(query_sce)
query_meta <- as.data.table(as.data.frame(SummarizedExperiment::colData(query_sce), stringsAsFactors = FALSE))
query_meta[, cell_id := cell_ids]

cells_dt <- copy(query_meta)
cells_dt[, `:=`(
  sample_id = as.character(sample_id),
  cell_id = cell_ids,
  label = as.character(per_cell_labels),
  score = as.numeric(per_cell_scores),
  method_id = opts$method_id,
  engine = "singler",
  reference_name = opts$reference_name,
  reference_label_col = opts$reference_label_col
)]

if ("X_umap" %in% reducedDimNames(query_sce)) {
  umap_mat <- reducedDim(query_sce, "X_umap")
  if (!is.null(umap_mat) && ncol(umap_mat) >= 2) {
    cells_dt[, `:=`(
      UMAP1 = as.numeric(umap_mat[, 1]),
      UMAP2 = as.numeric(umap_mat[, 2])
    )]
  }
}

cluster_summary_dt <- make_cluster_prediction_summary(cells_dt)

front_cols <- c(
  "sample_id", "cell_id", "label", "score", "method_id",
  "engine", "reference_name", "reference_label_col"
)
setcolorder(cells_dt, c(front_cols, setdiff(names(cells_dt), front_cols)))

prefix <- sanitize_prefix(opts$method_id)
export_dt <- data.table(cell_id = cells_dt$cell_id)
export_dt[[paste0(prefix, "_annotation")]] <- cells_dt$label
export_dt[[paste0(prefix, "_annotation_score")]] <- cells_dt$score

fwrite(cells_dt, opts$output_cells_csv)
fwrite(cluster_summary_dt, opts$output_cluster_csv)
fwrite(export_dt, opts$output_export_csv)

message(sprintf("Wrote SingleR annotation outputs for %s", opts$method_id))