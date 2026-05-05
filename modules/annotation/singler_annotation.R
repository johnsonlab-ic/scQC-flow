#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(ExperimentHub)
  library(AnnotationHub)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--seurat", type = "character",
              default = "/data/harmony_integration/harmony_bonn_ucl_9jan/integrated_harmony_seurat_bonn_ucl.rds",
              help = "Path to integrated Seurat .rds [default %default]"),
  make_option("--output", type = "character",
              default = "/data/singleR_annotation/singleR_bonn_ucl_9jan",
              help = "Output directory for SingleR annotations [default %default]"),
  make_option("--default-level", type = "character", default = "mid.cell.class",
              help = "Annotation level stored in `cell_type` [default %default]"),
  make_option("--individual-col", type = "character", default = "Individual_ID",
              help = "Meta-data column copied into `individual` [default %default]"),
  make_option("--ncores", type = "integer", default = 1L,
              help = "Cores for SingleR / BiocParallel"),
  make_option("--assay", type = "character", default = "RNA",
              help = "Seurat assay to operate on [default %default]"),
  make_option("--cluster-col", type = "character", default = "harmony_leiden",
              help = "Meta-data column for cluster-level SingleR; set to '' to disable [default %default]"),
  make_option("--levels", type = "character",
              default = "cell.type2,superfine.cell.class,fine.cell.class,mid.cell.class,broad.cell.class",
              help = "Comma-separated annotation levels to run [default %default]"),
  make_option("--fine-tune", action = "store_true", default = FALSE,
              help = "Enable SingleR fine-tuning (slower) [default %default]"),
  make_option("--prune", action = "store_true", default = TRUE,
              help = "Enable SingleR pruning [default %default]"),
  make_option("--bp-type", type = "character", default = "multicore",
              help = "BiocParallel backend (multicore|snow) [default %default]")
)

parser <- OptionParser(
  option_list = option_list,
  description = "Annotate an integrated Seurat object with SingleR using the humanHippocampus2024 reference."
)
opts <- parse_args(parser)

if (is.null(opts$seurat) || is.null(opts$output)) {
  print_help(parser)
  stop("Missing required arguments: --seurat and --output", call. = FALSE)
}

dir.create(opts$output, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# Load integrated data and convert to SCE
# ----------------------------------------------------------------------
message("Loading Seurat object: ", opts$seurat)
seu <- readRDS(opts$seurat)
DefaultAssay(seu) <- opts$assay

# Join layered assays (Seurat v5 multi-sample) into a single counts layer for downstream conversion
layer_names <- tryCatch(Layers(seu[[opts$assay]]), error = function(e) character(0))
if (length(layer_names) > 0 && any(grepl("^counts[.]", layer_names))) {
  message("Joining layered assay '", opts$assay, "' into a single counts matrix for SingleR ...")
  seu <- SeuratObject::JoinLayers(seu, assay = opts$assay)
  DefaultAssay(seu) <- opts$assay
}

to_sce <- function(seu, assay = "RNA") {
  sce <- tryCatch({
    as.SingleCellExperiment(seu)
  }, error = function(e) NULL)

  if (is.null(sce)) {
    counts <- tryCatch({
      Seurat::GetAssayData(seu[[assay]], layer = "counts")
    }, error = function(e) {
      Seurat::GetAssayData(seu[[assay]], slot = "counts")
    })
    sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)))
    colData(sce) <- DataFrame(seu@meta.data)
  }

  if (!"logcounts" %in% assayNames(sce)) {
    sce <- scuttle::logNormCounts(sce)
  }
  sce
}

sce_test <- to_sce(seu, assay = opts$assay)
message("Test SCE dims: ", paste(dim(sce_test), collapse = " x "))

# ----------------------------------------------------------------------
# Load reference (humanHippocampus2024 EH9606)
# ----------------------------------------------------------------------
message("Fetching humanHippocampus2024 reference (EH9606) ...")
eh <- ExperimentHub::ExperimentHub()
files <- AnnotationHub::query(eh, "humanHippocampus2024")
sce_ref <- files[["EH9606"]]
message("Reference SCE dims: ", paste(dim(sce_ref), collapse = " x "))

if (!"logcounts" %in% assayNames(sce_ref)) {
  sce_ref <- scuttle::logNormCounts(sce_ref)
}

# ----------------------------------------------------------------------
# Align genes between query and reference
# ----------------------------------------------------------------------
align_genes <- function(test_sce, ref_sce) {
  common <- intersect(rownames(test_sce), rownames(ref_sce))
  if (length(common) >= 200) {
    return(list(test = test_sce[common, ], ref = ref_sce[common, ]))
  }

  test_uc <- toupper(rownames(test_sce))
  ref_uc  <- toupper(rownames(ref_sce))
  common_uc <- intersect(test_uc, ref_uc)
  common_uc <- common_uc[!duplicated(common_uc)]

  if (length(common_uc) < 200) {
    stop("Too few overlapping genes between the dataset and reference. Ensure gene identifiers are compatible.")
  }

  test_idx <- match(common_uc, test_uc)
  ref_idx  <- match(common_uc, ref_uc)
  test_sce <- test_sce[test_idx, ]
  ref_sce  <- ref_sce[ref_idx, ]
  rownames(test_sce) <- rownames(ref_sce)
  list(test = test_sce, ref = ref_sce)
}

aligned <- align_genes(sce_test, sce_ref)
sce_test_al <- aligned$test
sce_ref_al  <- aligned$ref

# ----------------------------------------------------------------------
# Run SingleR at requested annotation resolutions
# ----------------------------------------------------------------------
choose_bp <- function(ncores, bp_type = "multicore") {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("Package 'BiocParallel' is required.")
  }
  if (ncores > 1) {
    if (bp_type == "multicore" && .Platform$OS.type != "windows") {
      return(BiocParallel::MulticoreParam(workers = ncores))
    } else {
      return(BiocParallel::SnowParam(workers = ncores))
    }
  }
  BiocParallel::SerialParam()
}

run_singler <- function(test_sce, ref_sce, labels, clusters = NULL,
                        ncores = 1L, fine_tune = FALSE, prune = TRUE, bp_type = "multicore") {
  if (!requireNamespace("SingleR", quietly = TRUE))
    stop("Package 'SingleR' is required.")
  bp <- choose_bp(ncores, bp_type)
  SingleR::SingleR(
    test = test_sce,
    ref = ref_sce,
    labels = labels,
    clusters = clusters,
    assay.type.test = "logcounts",
    assay.type.ref = "logcounts",
    fine.tune = fine_tune,
    prune = prune,
    BPPARAM = bp
  )
}

cluster_vec <- NULL
if (nzchar(opts$`cluster-col`) &&
    opts$`cluster-col` %in% colnames(SummarizedExperiment::colData(sce_test_al))) {
  cluster_vec <- as.factor(SummarizedExperiment::colData(sce_test_al)[[opts$`cluster-col`]])
  message("Using cluster-level SingleR on column: ", opts$`cluster-col`)
} else if (nzchar(opts$`cluster-col`)) {
  warning("Requested cluster-col '", opts$`cluster-col`, "' not found; running per-cell SingleR.")
}

label_cols <- trimws(unlist(strsplit(opts$levels, ",", fixed = TRUE)))
valid_cols <- c("cell.type2","superfine.cell.class","fine.cell.class","mid.cell.class","broad.cell.class")
if (length(setdiff(label_cols, valid_cols))) {
  stop("Unknown levels in --levels; valid levels: ", paste(valid_cols, collapse = ", "))
}
missing_cols <- setdiff(label_cols, colnames(SummarizedExperiment::colData(sce_ref_al)))
if (length(missing_cols)) {
  stop("Reference is missing requested label columns: ", paste(missing_cols, collapse = ", "))
}

expand_cluster_result <- function(pred, cluster_vec) {
  cl_map <- pred$labels
  names(cl_map) <- rownames(pred)
  per_cell_labels <- cl_map[as.character(cluster_vec)]
  max_scores <- apply(pred$scores, 1, max, na.rm = TRUE)
  names(max_scores) <- rownames(pred)
  per_cell_scores <- max_scores[as.character(cluster_vec)]
  list(labels = per_cell_labels, scores = per_cell_scores)
}

preds <- list()
for (col in label_cols) {
  message("Running SingleR for level: ", col)
  ref_labels <- SummarizedExperiment::colData(sce_ref_al)[[col]]
  pred <- run_singler(
    test_sce = sce_test_al,
    ref_sce  = sce_ref_al,
    labels   = ref_labels,
    clusters = cluster_vec,
    ncores   = opts$ncores,
    fine_tune = opts$`fine-tune`,
    prune     = opts$prune,
    bp_type   = opts$`bp-type`
  )
  if (!is.null(cluster_vec)) {
    expanded <- expand_cluster_result(pred, cluster_vec)
    preds[[col]] <- list(
      labels = expanded$labels,
      score_vec = expanded$scores,
      cluster_result = pred
    )
  } else {
    score_vec <- apply(pred$scores, 1, max, na.rm = TRUE)
    preds[[col]] <- list(
      labels = pred$labels,
      score_vec = score_vec,
      cluster_result = pred
    )
  }
}

add_pred_to_meta <- function(seu, pred_list, prefix) {
  md <- seu@meta.data
  md[[paste0("singler_", prefix)]] <- pred_list$labels
  md[[paste0("singler_", prefix, "_score")]] <- pred_list$score_vec
  seu@meta.data <- md
  seu
}

for (col in names(preds)) {
  key <- gsub("[.]", "_", col)
  seu <- add_pred_to_meta(seu, preds[[col]], key)
}


# ----------------------------------------------------------------------
# Ensure required meta columns
# ----------------------------------------------------------------------
default_key <- gsub("[.]", "_", opts$`default-level`)
default_col <- paste0("singler_", default_key)
if (!default_col %in% colnames(seu@meta.data)) {
  stop("Default level column not found: ", default_col)
}
seu$cell_type <- seu[[default_col]][, 1, drop = TRUE]

ind_source <- opts$`individual-col`
if (!ind_source %in% colnames(seu@meta.data)) {
  warning("Requested individual column '", ind_source, "' not found. Using Sample_ID if available.")
  ind_source <- if ("Sample_ID" %in% colnames(seu@meta.data)) "Sample_ID" else NULL
}
seu$individual <- if (!is.null(ind_source)) seu[[ind_source]][, 1, drop = TRUE] else "unknown"

# ----------------------------------------------------------------------
# Save annotated Seurat object
# ----------------------------------------------------------------------
annot_path <- file.path(opts$output, "integrated_with_singleR.rds")
saveRDS(seu, annot_path)
message("Saved annotated Seurat object: ", annot_path)

# ----------------------------------------------------------------------
# Export per-level predictions and counts
# ----------------------------------------------------------------------
for (nm in names(preds)) {
  key <- gsub("[.]", "_", nm)
  pred <- preds[[nm]]
  score_vec <- pred$score_vec
  per_cell <- data.frame(
    cell_barcode = rownames(seu@meta.data),
    label = pred$labels,
    score = score_vec,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    per_cell,
    file = file.path(opts$output, paste0("singleR_cells_", key, ".csv")),
    row.names = FALSE
  )

  counts <- as.data.frame(sort(table(pred$labels), decreasing = TRUE))
  colnames(counts) <- c("label", "n_cells")
  utils::write.csv(
    counts,
    file = file.path(opts$output, paste0("singleR_counts_", key, ".csv")),
    row.names = FALSE
  )

  if (!is.null(cluster_vec) && !is.null(pred$cluster_result)) {
    cl_map <- data.frame(
      cluster_id = rownames(pred$cluster_result),
      label = pred$cluster_result$labels,
      max_score = apply(pred$cluster_result$scores, 1, max, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      cl_map,
      file = file.path(opts$output, paste0("singleR_cluster_labels_", key, ".csv")),
      row.names = FALSE
    )

    cl_counts <- as.data.frame(sort(table(pred$cluster_result$labels), decreasing = TRUE))
    colnames(cl_counts) <- c("cluster_label", "n_clusters")
    utils::write.csv(
      cl_counts,
      file = file.path(opts$output, paste0("singleR_cluster_counts_", key, ".csv")),
      row.names = FALSE
    )
  }
}

# ----------------------------------------------------------------------
# UMAP plots by annotation level
# ----------------------------------------------------------------------
if ("umap" %in% Reductions(seu)) {
  for (nm in names(preds)) {
    key <- gsub("[.]", "_", nm)
    coln <- paste0("singler_", key)
    plot_path <- file.path(opts$output, paste0("umap_singler_", key, ".png"))
    p <- DimPlot(
      seu,
      reduction = "umap",
      group.by = coln,
      label = TRUE,
      label.size = 3,
      repel = TRUE,
      raster = TRUE
    ) +
      ggtitle(paste0("UMAP: SingleR ", nm)) +
      theme_minimal(base_size = 12)
    ggsave(plot_path, plot = p, width = 10, height = 8, dpi = 300)
    message("Saved UMAP: ", plot_path)
  }
} else {
  message("UMAP reduction not found; skipping UMAP plots.")
}


# ----------------------------------------------------------------------
# Summary log
# ----------------------------------------------------------------------
summary_file <- file.path(opts$output, "singleR_summary.txt")
summary_lines <- c(
  paste0("SingleR annotation completed: ", Sys.time()),
  paste0("Seurat input: ", opts$seurat),
  paste0("Output dir:  ", normalizePath(opts$output)),
  paste0("Default level assigned to `cell_type`: ", opts$`default-level`),
  paste0("`individual` column sourced from: ", opts$`individual-col`),
  paste0("Cluster column: ", ifelse(nzchar(opts$`cluster-col`), opts$`cluster-col`, "<cell-level>")),
  paste0("Levels evaluated: ", opts$levels),
  paste0("Fine-tune: ", opts$`fine-tune`, "; prune: ", opts$prune),
  paste0("BiocParallel: ", opts$`bp-type`, " (ncores=", opts$ncores, ")")
)
writeLines(summary_lines, summary_file)
message("Summary written to: ", summary_file)
