#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
  library(xgboost)
})

option_list <- list(
  make_option("--reference_rds", type = "character", help = "Labelled SingleCellExperiment reference"),
  make_option("--label_col", type = "character", help = "Reference label column to predict"),
  make_option("--out_model", type = "character", help = "Output XGBoost model RDS"),
  make_option("--out_classes", type = "character", help = "Output class lookup CSV"),
  make_option("--n_hvgs", type = "integer", default = 4000L, help = "Number of HVGs to select [default %default]"),
  make_option("--seed", type = "integer", default = 1L, help = "Random seed [default %default]")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

required <- c("reference_rds", "label_col", "out_model", "out_classes")
missing <- required[vapply(required, function(key) is.null(opts[[key]]) || !nzchar(opts[[key]]), logical(1))]
if (length(missing) > 0) {
  print_help(parser)
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

set.seed(opts$seed)

message("Loading reference: ", opts$reference_rds)
sce <- readRDS(opts$reference_rds)
if (!opts$label_col %in% colnames(colData(sce))) {
  stop(sprintf("Label column not found in reference: %s", opts$label_col), call. = FALSE)
}
if (!"counts" %in% assayNames(sce)) {
  stop("Reference SCE must contain a 'counts' assay", call. = FALSE)
}

labels_chr <- as.character(colData(sce)[[opts$label_col]])
keep_cells <- !is.na(labels_chr) & nzchar(labels_chr)
sce <- sce[, keep_cells]
labels_chr <- labels_chr[keep_cells]

counts <- assay(sce, "counts")
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

gene_keep <- !is.na(rownames(counts)) & nzchar(rownames(counts))
counts <- counts[gene_keep, , drop = FALSE]
lib_sizes <- Matrix::colSums(counts)
lib_sizes[lib_sizes <= 0] <- 1

n_cells_total <- ncol(counts)
gene_means <- Matrix::rowSums(counts) / n_cells_total
gene_means_sq <- Matrix::rowSums(counts^2) / n_cells_total
vars <- gene_means_sq - gene_means^2
ord <- order(vars, decreasing = TRUE)
hvgs <- rownames(counts)[head(ord, min(opts$n_hvgs, length(ord)))]

logm_hvg <- counts[hvgs, , drop = FALSE]
logm_hvg <- Matrix::t(Matrix::t(logm_hvg) / lib_sizes) * 1e4
logm_hvg@x <- log1p(logm_hvg@x)

classes <- sort(unique(labels_chr))
cls_dt <- data.table(cluster = classes)
dir.create(dirname(opts$out_classes), recursive = TRUE, showWarnings = FALSE)
fwrite(cls_dt, opts$out_classes)

y <- match(labels_chr, classes) - 1L
num_class <- length(classes)

n <- ncol(logm_hvg)
idx <- sample.int(n)
n_train <- floor(0.9 * n)
tr_idx <- idx[seq_len(n_train)]
va_idx <- idx[(n_train + 1):n]

x_train <- Matrix::t(logm_hvg[, tr_idx, drop = FALSE])
x_valid <- Matrix::t(logm_hvg[, va_idx, drop = FALSE])
colnames(x_train) <- hvgs
colnames(x_valid) <- hvgs

dtrain <- xgb.DMatrix(data = x_train, label = y[tr_idx])
dvalid <- xgb.DMatrix(data = x_valid, label = y[va_idx])
setinfo(dtrain, "feature_name", hvgs)
setinfo(dvalid, "feature_name", hvgs)

params <- list(
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = num_class,
  eta = 0.1,
  max_depth = 8,
  min_child_weight = 2,
  subsample = 0.8,
  colsample_bytree = 0.8,
  tree_method = "hist"
)

message(sprintf("Training XGBoost model with %d cells, %d HVGs, %d classes", ncol(logm_hvg), nrow(logm_hvg), num_class))
fit <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 600,
  watchlist = list(train = dtrain, valid = dvalid),
  early_stopping_rounds = 30,
  verbose = 1
)

dir.create(dirname(opts$out_model), recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, opts$out_model)

message("Wrote model: ", opts$out_model)
message("Wrote classes: ", opts$out_classes)