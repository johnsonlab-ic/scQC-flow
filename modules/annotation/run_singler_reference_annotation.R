#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(SingleR)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(data.table)
})

parser <- ArgumentParser(description = "Run per-cell SingleR annotation for one reference")
parser$add_argument("--query_rds", type = "character", required = TRUE,
  help = "Path to prepared query SCE RDS")
parser$add_argument("--reference_rds", type = "character", required = TRUE,
  help = "Path to SingleR reference RDS")
parser$add_argument("--reference_label_col", type = "character", required = TRUE,
  help = "Reference label column")
parser$add_argument("--method_id", type = "character", required = TRUE,
  help = "Method identifier")
parser$add_argument("--reference_name", type = "character", required = TRUE,
  help = "Reference display name")
parser$add_argument("--output_cells_csv", type = "character", required = TRUE,
  help = "Per-cell output CSV.GZ")
parser$add_argument("--output_export_csv", type = "character", required = TRUE,
  help = "Export metadata CSV.GZ")
parser$add_argument("--output_cluster_csv", type = "character", required = TRUE,
  help = "Cluster summary CSV.GZ")
parser$add_argument("--ncores", type = "integer", default = 1L,
  help = "BiocParallel workers")
parser$add_argument("--bp_type", type = "character", default = "multicore",
  help = "BiocParallel backend (multicore|snow)")
parser$add_argument("--fine_tune", type = "character", default = "false",
  help = "Enable SingleR fine-tuning")
parser$add_argument("--prune", type = "character", default = "true",
  help = "Enable SingleR pruning")

opts <- parser$parse_args()

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

extract_terminal_ensg <- function(ids) {
  ids <- as.character(ids)
  has_ensg <- grepl("(^|[_-])(ENS(?:MUS)?G[^._-]+(?:\\.[0-9]+)?)$", ids)
  out <- rep(NA_character_, length(ids))
  out[has_ensg] <- sub("^.*[_-](ENS(?:MUS)?G[^._-]+(?:\\.[0-9]+)?)$", "\\1", ids[has_ensg])
  out
}

extract_symbol_prefix <- function(ids) {
  ids <- as.character(ids)
  has_ensg <- grepl("^.+[_-](ENS(?:MUS)?G[^._-]+(?:\\.[0-9]+)?)$", ids)
  out <- rep(NA_character_, length(ids))
  out[has_ensg] <- sub("^(.+)[_-](ENS(?:MUS)?G[^._-]+(?:\\.[0-9]+)?)$", "\\1", ids[has_ensg])
  out
}

build_gene_candidates <- function(sce) {
  raw_ids <- as.character(rownames(sce))
  raw_novers <- sub("\\..*$", "", raw_ids)

  rd <- as.data.frame(SummarizedExperiment::rowData(sce), stringsAsFactors = FALSE)
  gene_id <- if ("gene_id" %in% names(rd)) as.character(rd$gene_id) else raw_ids
  symbol <- if ("symbol" %in% names(rd)) as.character(rd$symbol) else raw_ids

  gene_id_novers <- sub("\\..*$", "", gene_id)
  raw_tail_ensg <- extract_terminal_ensg(raw_ids)
  raw_tail_ensg_novers <- sub("\\..*$", "", raw_tail_ensg)
  raw_symbol_prefix <- extract_symbol_prefix(raw_ids)

  list(
    raw = raw_ids,
    raw_uc = toupper(raw_ids),
    raw_novers = raw_novers,
    gene_id = gene_id,
    gene_id_novers = gene_id_novers,
    symbol = symbol,
    symbol_uc = toupper(symbol),
    raw_tail_ensg = raw_tail_ensg,
    raw_tail_ensg_novers = raw_tail_ensg_novers,
    raw_symbol_prefix = raw_symbol_prefix,
    raw_symbol_prefix_uc = toupper(raw_symbol_prefix),
    symbol_ensg = ifelse(
      !is.na(symbol) & nzchar(symbol) & !is.na(gene_id_novers) & nzchar(gene_id_novers),
      paste0(symbol, "-", gene_id_novers),
      NA_character_
    )
  )
}

pair_overlap <- function(query_ids, ref_ids) {
  valid_query <- !is.na(query_ids) & nzchar(query_ids)
  valid_ref <- !is.na(ref_ids) & nzchar(ref_ids)
  if (!any(valid_query) || !any(valid_ref)) {
    return(integer())
  }

  query_counts <- table(query_ids[valid_query])
  ref_counts <- table(ref_ids[valid_ref])
  intersect(names(query_counts)[query_counts == 1L], names(ref_counts)[ref_counts == 1L])
}

align_genes <- function(query_sce, ref_sce) {
  query_candidates <- build_gene_candidates(query_sce)
  ref_candidates <- build_gene_candidates(ref_sce)

  candidate_pairs <- list(
    c("raw", "raw"),
    c("raw_novers", "raw_novers"),
    c("raw_uc", "raw_uc"),
    c("raw_tail_ensg", "raw"),
    c("raw_tail_ensg_novers", "raw_novers"),
    c("raw_tail_ensg", "gene_id"),
    c("raw_tail_ensg_novers", "gene_id_novers"),
    c("raw_symbol_prefix", "symbol"),
    c("raw_symbol_prefix_uc", "symbol_uc"),
    c("symbol", "symbol"),
    c("symbol_uc", "symbol_uc"),
    c("symbol_ensg", "raw"),
    c("symbol_ensg", "raw_novers")
  )

  overlap_dt <- rbindlist(lapply(candidate_pairs, function(pair) {
    common <- pair_overlap(query_candidates[[pair[[1]]]], ref_candidates[[pair[[2]]]])
    data.table(query_scheme = pair[[1]], ref_scheme = pair[[2]], overlap = length(common))
  }))
  best <- overlap_dt[which.max(overlap)]

  if (nrow(best) == 0 || best$overlap[[1]] < 200) {
    message("SingleR gene alignment overlap summary:")
    print(overlap_dt[order(-overlap)])
    stop("Too few overlapping genes between query and SingleR reference", call. = FALSE)
  }

  common <- pair_overlap(
    query_candidates[[best$query_scheme[[1]]]],
    ref_candidates[[best$ref_scheme[[1]]]]
  )
  query_idx <- match(common, query_candidates[[best$query_scheme[[1]]]])
  ref_idx <- match(common, ref_candidates[[best$ref_scheme[[1]]]])
  query_aligned <- query_sce[query_idx, ]
  ref_aligned <- ref_sce[ref_idx, ]
  rownames(query_aligned) <- common
  rownames(ref_aligned) <- common
  message(sprintf(
    "Aligned SingleR genes using query scheme '%s' and reference scheme '%s' (%d genes)",
    best$query_scheme[[1]],
    best$ref_scheme[[1]],
    length(common)
  ))
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