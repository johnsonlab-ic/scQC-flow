#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

parse_args <- function(args) {
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs")
  }
  out <- list()
  idx <- seq(1, length(args), by = 2)
  for (i in idx) {
    key <- sub("^--", "", args[[i]])
    out[[key]] <- args[[i + 1]]
  }
  out
}

add_reduction_if_present <- function(seu, meta_dt, sel_cols, reduction_name, key, dim_names = NULL) {
  if (length(sel_cols) == 0) {
    return(seu)
  }
  emb <- as.matrix(meta_dt[, ..sel_cols])
  storage.mode(emb) <- "double"
  rownames(emb) <- rownames(meta_dt)
  if (!is.null(dim_names)) {
    colnames(emb) <- dim_names
  }
  seu[[reduction_name]] <- CreateDimReducObject(
    embeddings = emb,
    key = key,
    assay = DefaultAssay(seu)
  )
  seu
}

build_seurat_object <- function(counts_mat, obs_dt, var_dt, is_all_cells) {
  meta_dt <- as.data.frame(obs_dt)
  rownames(meta_dt) <- meta_dt$cell_id
  meta_dt$cell_id <- NULL

  seu <- CreateSeuratObject(
    counts = counts_mat,
    project = "scqcflow",
    meta.data = meta_dt
  )

  feature_meta <- as.data.frame(var_dt[, .(gene_id, symbol, gene_type)])
  rownames(feature_meta) <- var_dt$gene_key
  seu[[DefaultAssay(seu)]]@meta.features <- feature_meta
  seu <- NormalizeData(seu, verbose = FALSE)

  if (is_all_cells && all(c("dbl_UMAP1", "dbl_UMAP2") %in% names(meta_dt))) {
    emb <- as.matrix(meta_dt[, c("dbl_UMAP1", "dbl_UMAP2"), drop = FALSE])
    storage.mode(emb) <- "double"
    rownames(emb) <- rownames(meta_dt)
    colnames(emb) <- c("UMAP_1", "UMAP_2")
    seu[["umap"]] <- CreateDimReducObject(
      embeddings = emb,
      key = "UMAP_",
      assay = DefaultAssay(seu)
    )
  } else if (all(c("UMAP1", "UMAP2") %in% names(meta_dt))) {
    emb <- as.matrix(meta_dt[, c("UMAP1", "UMAP2"), drop = FALSE])
    storage.mode(emb) <- "double"
    rownames(emb) <- rownames(meta_dt)
    colnames(emb) <- c("UMAP_1", "UMAP_2")
    seu[["umap"]] <- CreateDimReducObject(
      embeddings = emb,
      key = "UMAP_",
      assay = DefaultAssay(seu)
    )
  }

  hmny_cols <- grep("^hmny_pca_[0-9]+$", names(meta_dt), value = TRUE)
  if (length(hmny_cols) > 0 && !is_all_cells) {
    ord <- order(as.integer(sub("^hmny_pca_", "", hmny_cols)))
    hmny_cols <- hmny_cols[ord]
    dim_names <- paste0("HARMONY_", seq_along(hmny_cols))
    seu <- add_reduction_if_present(seu, meta_dt, hmny_cols, "harmony", "HARMONY_", dim_names)
  }

  pca_cols <- grep("^pca_[0-9]+$", names(meta_dt), value = TRUE)
  if (length(pca_cols) > 0 && !is_all_cells) {
    ord <- order(as.integer(sub("^pca_", "", pca_cols)))
    pca_cols <- pca_cols[ord]
    dim_names <- paste0("PC_", seq_along(pca_cols))
    seu <- add_reduction_if_present(seu, meta_dt, pca_cols, "pca", "PC_", dim_names)
  }

  seu
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  required <- c("h5_pattern", "qc_pattern", "integration_csv", "annotation_csv", "genome_gtf", "utils_r", "out_all", "out_clean")
  missing <- setdiff(required, names(args))
  if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
  }

  source(args$utils_r)

  obs_dt <- load_export_metadata(args$integration_csv, args$qc_pattern, args$annotation_csv)
  gtf_dt <- parse_gtf_annotations_export(args$genome_gtf)

  all_counts_obj <- build_export_counts(args$h5_pattern, obs_dt)
  all_obs <- obs_dt[match(all_counts_obj$ordered_ids, cell_id)]
  var_dt <- build_gene_metadata_export(all_counts_obj$gene_keys, gtf_dt)
  all_counts <- all_counts_obj$counts
  rownames(all_counts) <- all_counts_obj$gene_keys
  colnames(all_counts) <- all_counts_obj$ordered_ids
  all_seu <- build_seurat_object(all_counts, all_obs, var_dt, is_all_cells = TRUE)
  saveRDS(all_seu, args$out_all)

  clean_obs <- all_obs[is_dbl == FALSE & in_dbl_cl == FALSE]
  clean_counts_obj <- build_export_counts(args$h5_pattern, clean_obs)
  clean_obs <- clean_obs[match(clean_counts_obj$ordered_ids, cell_id)]
  clean_counts <- clean_counts_obj$counts
  rownames(clean_counts) <- clean_counts_obj$gene_keys
  colnames(clean_counts) <- clean_counts_obj$ordered_ids
  clean_seu <- build_seurat_object(clean_counts, clean_obs, var_dt, is_all_cells = FALSE)
  saveRDS(clean_seu, args$out_clean)
}

main()