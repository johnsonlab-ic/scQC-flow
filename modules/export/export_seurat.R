#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(assertthat)
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
  emb <- as.matrix(meta_dt[, sel_cols, drop = FALSE])
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
  seu[[DefaultAssay(seu)]] <- AddMetaData(
    object = seu[[DefaultAssay(seu)]],
    metadata = feature_meta
  )

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

merge_saved_seurat_objects <- function(rds_paths) {
  if (length(rds_paths) == 0) {
    return(NULL)
  }

  message(sprintf("  Incrementally merging %d saved Seurat objects", length(rds_paths)))
  merged <- readRDS(rds_paths[[1]])

  if (length(rds_paths) > 1) {
    for (path in rds_paths[-1]) {
      next_obj <- readRDS(path)
      merged <- merge(merged, next_obj, add.cell.ids = NULL, merge.data = FALSE)
      rm(next_obj)
      gc(verbose = FALSE)
    }
  }

  merged <- merge_with_bpcells(list(merged))
  merged
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  required <- c("h5_pattern", "qc_pattern", "integration_csv", "annotation_csv", "genome_gtf", "utils_r", "out_dir", "write_combined")
  missing <- setdiff(required, names(args))
  if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
  }

  source(args$utils_r)

  write_combined <- as.logical(args$write_combined)
  out_dir <- args$out_dir
  dir.create(file.path(out_dir, "seurat"), showWarnings = FALSE, recursive = TRUE)

  message("Loading metadata and annotations...")
  obs_dt <- load_export_metadata(args$integration_csv, args$qc_pattern, args$annotation_csv)
  gtf_dt <- parse_gtf_annotations_export(args$genome_gtf)

  # Get unique samples
  unique_samples <- unique(obs_dt$sample_id)
  message(sprintf("Processing %d samples", length(unique_samples)))

  # Per-sample processing
  all_rds_paths <- character()
  clean_rds_paths <- character()

  for (sid in unique_samples) {
    message(sprintf("Processing sample: %s", sid))

    # All cells for this sample
    sample_obs_all <- obs_dt[sample_id == sid]
    if (nrow(sample_obs_all) == 0) {
      message(sprintf("  Skipping %s: no cells found", sid))
      next
    }

    counts_obj_all <- build_export_counts(args$h5_pattern, sample_obs_all)
    sample_obs_all <- sample_obs_all[match(counts_obj_all$ordered_ids, cell_id)]
    var_dt <- build_gene_metadata_export(counts_obj_all$gene_keys, gtf_dt)

    counts_all <- counts_obj_all$counts
    rownames(counts_all) <- counts_obj_all$gene_keys
    colnames(counts_all) <- counts_obj_all$ordered_ids

    seu_all <- build_seurat_object(counts_all, sample_obs_all, var_dt, is_all_cells = TRUE)
    all_rds_path <- file.path(out_dir, "seurat", sprintf("sample_%s_all.rds", sid))
    saveRDS(seu_all, all_rds_path)
    message(sprintf("  Saved all cells: %d cells", ncol(seu_all)))

    if (write_combined) {
      all_rds_paths <- c(all_rds_paths, all_rds_path)
    }

    # Clean cells for this sample
    # Only create clean subset if doublet columns are present
    has_doublet_cols <- all(c("is_dbl", "in_dbl_cl") %in% colnames(sample_obs_all))
    
    if (has_doublet_cols) {
      sample_obs_clean <- sample_obs_all[is_dbl == FALSE & in_dbl_cl == FALSE]
      if (nrow(sample_obs_clean) > 0) {
        keep_clean <- match(sample_obs_clean$cell_id, colnames(counts_all))
        counts_clean <- counts_all[, keep_clean, drop = FALSE]
        colnames(counts_clean) <- sample_obs_clean$cell_id

        seu_clean <- build_seurat_object(counts_clean, sample_obs_clean, var_dt, is_all_cells = FALSE)
        clean_rds_path <- file.path(out_dir, "seurat", sprintf("sample_%s_clean.rds", sid))
        saveRDS(seu_clean, clean_rds_path)
        message(sprintf("  Saved clean cells: %d cells", ncol(seu_clean)))

        if (write_combined) {
          clean_rds_paths <- c(clean_rds_paths, clean_rds_path)
        }

        rm(seu_clean, counts_clean, sample_obs_clean)
        gc(verbose = FALSE)
      }
    } else {
      message(sprintf("  Skipping clean subset: doublet columns not found in metadata"))
    }

    rm(seu_all, counts_all, counts_obj_all, sample_obs_all, var_dt)
    gc(verbose = FALSE)
  }

  # Write combined objects if requested
  if (write_combined && length(all_rds_paths) > 0) {
    message("Building combined objects with BPCells backing...")

    # Combine all cells with BPCells
    message("  Combining all cells...")
    combined_all <- merge_saved_seurat_objects(all_rds_paths)
    saveRDS(combined_all, file.path(out_dir, "seurat", "combined_all.rds"))
    message(sprintf("  Saved combined all: %d cells", ncol(combined_all)))
    rm(combined_all)
    gc(verbose = FALSE)

    # Combine clean cells with BPCells
    if (length(clean_rds_paths) > 0) {
      message("  Combining clean cells...")
      combined_clean <- merge_saved_seurat_objects(clean_rds_paths)
      saveRDS(combined_clean, file.path(out_dir, "seurat", "combined_clean.rds"))
      message(sprintf("  Saved combined clean: %d cells", ncol(combined_clean)))
      rm(combined_clean)
      gc(verbose = FALSE)
    }
  }

  message("Seurat export complete!")
}

main()