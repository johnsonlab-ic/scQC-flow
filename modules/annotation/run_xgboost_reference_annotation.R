#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(xgboost)
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
})

option_list <- list(
  make_option("--query_rds", type = "character", help = "Path to prepared query SCE RDS"),
  make_option("--model_rds", type = "character", help = "Path to pretrained xgboost RDS"),
  make_option("--class_csv", type = "character", help = "Path to class lookup CSV"),
  make_option("--cluster_col", type = "character", default = "", help = "Optional query cluster column [default disabled]"),
  make_option("--method_id", type = "character", help = "Method identifier"),
  make_option("--reference_name", type = "character", help = "Reference display name"),
  make_option("--output_cells_csv", type = "character", help = "Per-cell output CSV.GZ"),
  make_option("--output_export_csv", type = "character", help = "Export metadata CSV.GZ"),
  make_option("--output_cluster_csv", type = "character", help = "Cluster summary CSV.GZ"),
  make_option("--chunk_size", type = "integer", default = 10000L, help = "Prediction chunk size [default %default]"),
  make_option("--scale_factor", type = "double", default = 1e4, help = "Library-size scaling factor [default %default]")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

required <- c(
  "query_rds", "model_rds", "class_csv", "method_id", "reference_name",
  "output_cells_csv", "output_export_csv", "output_cluster_csv"
)
missing <- required[vapply(required, function(key) is.null(opts[[key]]) || !nzchar(opts[[key]]), logical(1))]
if (length(missing) > 0) {
  print_help(parser)
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

sanitize_prefix <- function(x) {
  gsub("[^A-Za-z0-9_]", "_", x)
}

choose_feature_ids <- function(ids, hvgs) {
  ids <- as.character(ids)
  ids_novers <- sub("\\..*$", "", ids)

  candidates <- list(
    raw = ids,
    no_version = ids_novers,
    ensg_dash = sub("_ENSG", "-ENSG", ids, fixed = TRUE)
  )

  ensg_frac <- mean(grepl("^ENSG", ids_novers))
  if (!is.na(ensg_frac) && ensg_frac > 0.5) {
    if (!requireNamespace("ensembldb", quietly = TRUE) || !requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      stop("Packages 'ensembldb' and 'EnsDb.Hsapiens.v86' are required for ENSG feature conversion", call. = FALSE)
    }

    map_dt <- data.table::as.data.table(
      ensembldb::select(
        EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
        keys = unique(ids_novers[grepl("^ENSG", ids_novers)]),
        keytype = "GENEID",
        columns = c("GENEID", "SYMBOL")
      )
    )
    map_dt <- map_dt[!is.na(SYMBOL) & nzchar(SYMBOL)]
    map_dt <- unique(map_dt, by = "GENEID")
    map_vec <- setNames(map_dt$SYMBOL, map_dt$GENEID)
    symbols <- unname(map_vec[ids_novers])
    candidates$symbol <- symbols
    candidates$symbol_ensg <- ifelse(!is.na(symbols) & nzchar(symbols), paste0(symbols, "-", ids_novers), NA_character_)
  }

  overlap_dt <- rbindlist(lapply(names(candidates), function(nm) {
    feat <- candidates[[nm]]
    keep <- !is.na(feat) & nzchar(feat)
    data.table(scheme = nm, overlap = length(intersect(unique(feat[keep]), hvgs)))
  }))

  best <- overlap_dt[which.max(overlap)]
  if (nrow(best) == 0 || best$overlap[[1]] == 0) {
    stop("No overlap between query feature IDs and XGBoost model HVGs", call. = FALSE)
  }
  candidates[[best$scheme[[1]]]]
}

predict_chunks <- function(xgb_obj, cls_dt, counts_mat, hvgs, chunk_size, scale_factor) {
  if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- as(counts_mat, "dgCMatrix")
  }

  cell_ids <- colnames(counts_mat)
  n_cells <- ncol(counts_mat)
  hvg_idx <- match(hvgs, rownames(counts_mat))
  have_hvg <- !is.na(hvg_idx)
  chunk_start <- seq.int(1L, n_cells, by = chunk_size)
  out <- vector("list", length(chunk_start))

  for (i in seq_along(chunk_start)) {
    start_idx <- chunk_start[[i]]
    end_idx <- min(start_idx + chunk_size - 1L, n_cells)
    idx <- start_idx:end_idx
    lib_sizes <- Matrix::colSums(counts_mat[, idx, drop = FALSE])
    lib_sizes[lib_sizes <= 0] <- 1

    dense_chunk <- matrix(0, nrow = length(idx), ncol = length(hvgs))
    if (any(have_hvg)) {
      sub_counts <- counts_mat[hvg_idx[have_hvg], idx, drop = FALSE]
      dense_chunk[, have_hvg] <- as.matrix(Matrix::t(sub_counts))
    }

    dense_chunk <- dense_chunk / as.numeric(lib_sizes)
    dense_chunk <- log1p(dense_chunk * scale_factor)

    probs <- predict(xgb_obj, dense_chunk, reshape = TRUE)
    if (is.null(dim(probs))) {
      probs <- matrix(probs, nrow = length(idx), byrow = TRUE)
    }
    colnames(probs) <- cls_dt$cluster
    pred_idx <- max.col(probs, ties.method = "first")

    out[[i]] <- data.table(
      cell_id = cell_ids[idx],
      label = cls_dt$cluster[pred_idx],
      score = probs[cbind(seq_len(length(idx)), pred_idx)]
    )
  }

  rbindlist(out)
}

message("Loading query SCE: ", opts$query_rds)
query_sce <- readRDS(opts$query_rds)
if (!"counts" %in% assayNames(query_sce)) {
  stop("Query SCE is missing the 'counts' assay", call. = FALSE)
}

message("Loading pretrained XGBoost model: ", opts$model_rds)
xgb_obj <- readRDS(opts$model_rds)
cls_dt <- fread(opts$class_csv)
if (!"cluster" %in% names(cls_dt)) {
  stop("XGBoost class CSV must contain a 'cluster' column", call. = FALSE)
}
hvgs <- variable.names(xgb_obj)

counts_mat <- assay(query_sce, "counts")
chosen_ids <- choose_feature_ids(rownames(counts_mat), hvgs)
keep_rows <- !is.na(chosen_ids) & nzchar(chosen_ids)
counts_mat <- counts_mat[keep_rows, , drop = FALSE]
rownames(counts_mat) <- chosen_ids[keep_rows]
counts_mat <- counts_mat[!duplicated(rownames(counts_mat)), , drop = FALSE]

preds_dt <- predict_chunks(
  xgb_obj = xgb_obj,
  cls_dt = cls_dt,
  counts_mat = counts_mat,
  hvgs = hvgs,
  chunk_size = opts$chunk_size,
  scale_factor = opts$scale_factor
)

query_meta <- as.data.table(as.data.frame(SummarizedExperiment::colData(query_sce), stringsAsFactors = FALSE))
query_meta[, cell_id := colnames(query_sce)]
cells_dt <- merge(query_meta, preds_dt, by = "cell_id", all.x = TRUE, sort = FALSE)
cells_dt[, `:=`(
  method_id = opts$method_id,
  engine = "xgboost",
  reference_name = opts$reference_name
)]

cluster_col <- trimws(opts$cluster_col)
if (nzchar(cluster_col) && !cluster_col %in% names(cells_dt)) {
  stop(sprintf("cluster_col '%s' not found in query metadata", cluster_col), call. = FALSE)
}

if ("X_umap" %in% reducedDimNames(query_sce)) {
  umap_mat <- reducedDim(query_sce, "X_umap")
  if (!is.null(umap_mat) && ncol(umap_mat) >= 2) {
    cells_dt[, `:=`(
      UMAP1 = as.numeric(umap_mat[, 1]),
      UMAP2 = as.numeric(umap_mat[, 2])
    )]
  }
}

cluster_summary_dt <- if (nzchar(cluster_col)) {
  cells_dt[, .(
    n_cells = .N,
    top_label = names(sort(table(label), decreasing = TRUE))[1],
    top_score = mean(score, na.rm = TRUE)
  ), by = .(cluster = as.character(get(cluster_col)))]
} else {
  cells_dt[, .(
    n_cells = .N,
    top_label = names(sort(table(label), decreasing = TRUE))[1],
    top_score = mean(score, na.rm = TRUE)
  ), by = .(cluster = label)]
}

prefix <- sanitize_prefix(opts$method_id)
export_dt <- data.table(cell_id = cells_dt$cell_id)
export_dt[[paste0(prefix, "_annotation")]] <- as.character(cells_dt$label)
export_dt[[paste0(prefix, "_annotation_score")]] <- as.numeric(cells_dt$score)

keep_front <- c("sample_id", "cell_id")
if (nzchar(cluster_col)) {
  setnames(cells_dt, cluster_col, "cluster")
  keep_front <- c(keep_front, "cluster")
} else if (!"cluster" %in% names(cells_dt)) {
  cells_dt[, cluster := NA_character_]
  keep_front <- c(keep_front, "cluster")
}
setcolorder(cells_dt, c(unique(c(keep_front, "label", "score", "method_id", "engine", "reference_name")), setdiff(names(cells_dt), unique(c(keep_front, "label", "score", "method_id", "engine", "reference_name")))))

fwrite(cells_dt, opts$output_cells_csv)
fwrite(cluster_summary_dt, opts$output_cluster_csv)
fwrite(export_dt, opts$output_export_csv)

message(sprintf("Wrote XGBoost annotation outputs for %s", opts$method_id))