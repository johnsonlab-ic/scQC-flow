#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript make_seurat.R <sample_name> <mapping_dir> <dropletqc_csv> <scdbl_csv>")
}

sample_name <- args[1]
mapping_dir_arg <- args[2]
dropletqc_file <- args[3]
scdbl_file <- args[4]

cat("make_seurat.R: sample:", sample_name, "mapping_dir_arg:", mapping_dir_arg, "\n")

# Try to find data path: prefer a local symlink in the workdir that contains outs/filtered_feature_bc_matrix
workdir_candidates <- list.files(path = getwd(), full.names = TRUE)
data_path <- NULL
for (cand in workdir_candidates) {
  if (dir.exists(file.path(cand, "outs/filtered_feature_bc_matrix"))) {
    data_path <- file.path(cand, "outs/filtered_feature_bc_matrix")
    break
  }
}
# fallback to provided mapping_dir_arg
if (is.null(data_path)) {
  possible <- file.path(mapping_dir_arg, "outs/filtered_feature_bc_matrix")
  if (dir.exists(possible)) data_path <- possible
}

if (is.null(data_path) || !dir.exists(data_path)) {
  stop(paste("Data directory not found for sample", sample_name, "tried:", paste(c(workdir_candidates, mapping_dir_arg), collapse = ", ")))
}

cat("Loading counts from:", data_path, "\n")
counts <- Read10X(data.dir = data_path)

seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Add DropletQC nuclear fraction if available
if (file.exists(dropletqc_file)) {
  dq <- tryCatch(read.csv(dropletqc_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(dq)) {
    seurat_bcs <- colnames(seurat_obj)
    clean_seurat_bcs <- sub("-.*$", "", seurat_bcs)
    # Try common barcode columns
    bc_col <- NULL
    if ("cell_barcode" %in% colnames(dq)) bc_col <- "cell_barcode"
    else if ("barcode" %in% colnames(dq)) bc_col <- "barcode"
    if (!is.null(bc_col) && ("nuclear_fraction" %in% colnames(dq))) {
      clean_dq_bcs <- sub("-.*$", "", dq[[bc_col]])
      nf <- rep(NA_real_, length(seurat_bcs))
      match_idx <- match(clean_seurat_bcs, clean_dq_bcs)
      matched <- which(!is.na(match_idx))
      nf[matched] <- dq$nuclear_fraction[match_idx[matched]]
      seurat_obj[["nuclear_fraction"]] <- nf
      cat("Added nuclear_fraction for", length(matched), "cells\n")
    }
  }
}

# Add scDbl metadata if available
if (file.exists(scdbl_file)) {
  sd <- tryCatch(read.csv(scdbl_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(sd)) {
    seurat_bcs <- colnames(seurat_obj)
    clean_seurat_bcs <- sub("-.*$", "", seurat_bcs)
    bc_col <- NULL
    if ("cell_barcode" %in% colnames(sd)) bc_col <- "cell_barcode"
    else if ("barcode" %in% colnames(sd)) bc_col <- "barcode"
    if (!is.null(bc_col)) {
      clean_sd_bcs <- sub("-.*$", "", sd[[bc_col]])
      match_idx <- match(clean_seurat_bcs, clean_sd_bcs)
      if ("doublet_score" %in% colnames(sd)) {
        ds <- rep(NA_real_, length(seurat_bcs))
        ds[which(!is.na(match_idx))] <- sd$doublet_score[match_idx[which(!is.na(match_idx))]]
        seurat_obj[["doublet_score"]] <- ds
      }
      if ("doublet_class" %in% colnames(sd)) {
        dc <- rep(NA_character_, length(seurat_bcs))
        dc[which(!is.na(match_idx))] <- as.character(sd$doublet_class[match_idx[which(!is.na(match_idx))]])
        seurat_obj[["doublet_class"]] <- dc
      }
      cat("Added scDbl metadata for", sum(!is.na(match_idx)), "cells\n")
    }
  }
}

# Write pre-QC object
pre_path <- file.path(getwd(), paste0(sample_name, "_seurat_object.rds"))
saveRDS(seurat_obj, pre_path)
cat("Saved pre-QC Seurat object to:", pre_path, "\n")

# Simple post-QC filtering
post_obj <- seurat_obj
min_features <- 200
min_counts <- 500
max_mito <- 20
keep <- rep(TRUE, ncol(post_obj))
if ("nFeature_RNA" %in% colnames(post_obj@meta.data)) keep <- keep & (post_obj@meta.data$nFeature_RNA >= min_features)
if ("nCount_RNA" %in% colnames(post_obj@meta.data)) keep <- keep & (post_obj@meta.data$nCount_RNA >= min_counts)
if ("percent.mt" %in% colnames(post_obj@meta.data)) keep <- keep & (post_obj@meta.data$percent.mt <= max_mito)
if ("doublet_class" %in% colnames(post_obj@meta.data)) {
  dc <- tolower(as.character(post_obj@meta.data$doublet_class))
  keep <- keep & (dc != "doublet")
}
post_obj <- subset(post_obj, cells = colnames(post_obj)[which(keep)])
post_path <- file.path(getwd(), paste0(sample_name, "_seurat_object_postqc.rds"))
saveRDS(post_obj, post_path)
cat("Saved post-QC Seurat object to:", post_path, "\n")

invisible(list(pre = pre_path, post = post_path))
