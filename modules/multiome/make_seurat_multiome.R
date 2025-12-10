#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(dplyr)
  library(argparse)
})

# Set up argument parser
parser <- ArgumentParser(description = 'Create Seurat objects with QC filtering for Multiome data')
parser$add_argument('sample_name', help = 'Sample name')
parser$add_argument('mapping_dir', help = 'Path to mapping directory')
parser$add_argument('--h5_path', type = 'character', default = NULL, help = 'Path to H5 counts file (CellBender or 10X) for Gene Expression')
parser$add_argument('--atac_h5_path', type = 'character', default = NULL, help = 'Path to ATAC/Peaks H5 file for ChromatinAssay')
parser$add_argument('dropletqc_csv', help = 'Path to DropletQC metrics CSV')
parser$add_argument('scdbl_csv', help = 'Path to scDbl metrics CSV')
parser$add_argument('--max_mito', type = 'double', default = 20.0, 
                   help = 'Maximum mitochondrial percentage threshold (default: 20)')
parser$add_argument('--min_nuclear', type = 'double', default = 0.0,
                   help = 'Minimum nuclear fraction threshold (default: 0.0)')
parser$add_argument('--metadata', type = 'character', default = NULL,
                   help = 'Optional metadata CSV file (sampleid column must match sample_name)')

args <- parser$parse_args()


sample_name <- args$sample_name
mapping_dir_arg <- args$mapping_dir
h5_path <- args$h5_path
atac_h5_path <- args$atac_h5_path
dropletqc_file <- args$dropletqc_csv
scdbl_file <- args$scdbl_csv
max_mito <- args$max_mito
min_nuclear <- args$min_nuclear
metadata_file <- args$metadata


cat("make_seurat_multiome.R: sample:", sample_name, "mapping_dir_arg:", mapping_dir_arg, "\n")
cat("QC thresholds: max_mito =", max_mito, ", min_nuclear =", min_nuclear, "\n")
cat("GEX H5 path:", h5_path, "\n")
cat("ATAC H5 path:", atac_h5_path, "\n")
if (!is.null(metadata_file)) cat("Using metadata file:", metadata_file, "\n")


# Load counts from H5 if provided, else fallback to directory
if (!is.null(h5_path) && file.exists(h5_path)) {
  cat("Loading counts from H5 file:", h5_path, "\n")
  counts_raw <- Read10X_h5(h5_path)
  
  # Extract Gene Expression modality for multiome data
  if (is.list(counts_raw)) {
    cat("Multiome H5 detected - extracting Gene Expression modality\n")
    counts <- counts_raw$`Gene Expression`
  } else {
    counts <- counts_raw
  }
} else {
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
  counts_raw <- Read10X(data.dir = data_path)
  
  # Extract Gene Expression modality for multiome data
  if (is.list(counts_raw)) {
    cat("Multiome data detected - extracting Gene Expression modality\n")
    counts <- counts_raw$`Gene Expression`
  } else {
    counts <- counts_raw
  }
}

seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200, assay = "RNA")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# =============================================================================
# Add ATAC ChromatinAssay if ATAC H5 is provided
# =============================================================================
if (!is.null(atac_h5_path) && file.exists(atac_h5_path)) {
  cat("\n=== Adding ATAC ChromatinAssay ===\n")
  cat("Loading ATAC counts from:", atac_h5_path, "\n")
  
  # Read ATAC counts
  atac_counts <- Read10X_h5(atac_h5_path)
  
  # If it's a list (shouldn't be since we extracted Peaks only), get Peaks
  if (is.list(atac_counts)) {
    if ("Peaks" %in% names(atac_counts)) {
      atac_counts <- atac_counts$Peaks
    } else {
      stop("No 'Peaks' modality found in ATAC H5 file")
    }
  }
  
  cat("ATAC matrix dimensions (before filtering):", nrow(atac_counts), "peaks x", ncol(atac_counts), "barcodes\n")
  
  # Get barcodes from the GEX Seurat object (these are the CellBender-filtered cells)
  gex_barcodes <- colnames(seurat_obj)
  cat("GEX cells (CellBender-filtered):", length(gex_barcodes), "\n")
  
  # Match barcodes between ATAC and GEX
  # Clean barcodes (remove suffix like -1) for matching
  clean_gex_bcs <- sub("-.*$", "", gex_barcodes)
  clean_atac_bcs <- sub("-.*$", "", colnames(atac_counts))
  
  # Find common barcodes
  common_bcs_clean <- intersect(clean_gex_bcs, clean_atac_bcs)
  cat("Common barcodes after cleaning:", length(common_bcs_clean), "\n")
  
  # Get the original barcode names that correspond to common barcodes
  gex_match_idx <- match(common_bcs_clean, clean_gex_bcs)
  atac_match_idx <- match(common_bcs_clean, clean_atac_bcs)
  
  gex_common_bcs <- gex_barcodes[gex_match_idx]
  atac_common_bcs <- colnames(atac_counts)[atac_match_idx]
  
  # Subset ATAC counts to only include cells present in GEX
  atac_counts_subset <- atac_counts[, atac_match_idx, drop = FALSE]
  
  # Rename ATAC barcodes to match GEX barcodes exactly
  colnames(atac_counts_subset) <- gex_common_bcs
  
  cat("ATAC matrix dimensions (after filtering to GEX cells):", nrow(atac_counts_subset), "peaks x", ncol(atac_counts_subset), "barcodes\n")
  
  # Subset the Seurat object to only include cells with ATAC data
  if (length(gex_common_bcs) < length(gex_barcodes)) {
    cat("Subsetting Seurat object to cells with ATAC data\n")
    seurat_obj <- subset(seurat_obj, cells = gex_common_bcs)
    cat("Seurat object now has", ncol(seurat_obj), "cells\n")
  }
  
  # Try to find fragments file for ChromatinAssay
  fragpath <- file.path(mapping_dir_arg, "outs/atac_fragments.tsv.gz")
  if (!file.exists(fragpath)) {
    # Try alternative locations
    fragpath_alt <- list.files(path = getwd(), pattern = "atac_fragments.tsv.gz", 
                                full.names = TRUE, recursive = TRUE)
    if (length(fragpath_alt) > 0) {
      fragpath <- fragpath_alt[1]
    } else {
      fragpath <- NULL
    }
  }
  
  # Create ChromatinAssay
  # Note: Without annotation and fragments, we create a basic assay
  if (!is.null(fragpath) && file.exists(fragpath)) {
    cat("Creating ChromatinAssay with fragments file:", fragpath, "\n")
    atac_assay <- CreateChromatinAssay(
      counts = atac_counts_subset,
      sep = c(":", "-"),
      fragments = fragpath,
      min.cells = 1,
      min.features = 0
    )
  } else {
    cat("Creating ChromatinAssay without fragments file\n")
    atac_assay <- CreateChromatinAssay(
      counts = atac_counts_subset,
      sep = c(":", "-"),
      min.cells = 1,
      min.features = 0
    )
  }
  
  # Add ATAC assay to Seurat object
  seurat_obj[["ATAC"]] <- atac_assay
  cat("Added ATAC assay to Seurat object\n")
  cat("ATAC assay features:", nrow(seurat_obj[["ATAC"]]), "\n")
  cat("=== ATAC ChromatinAssay added successfully ===\n\n")
} else {
  cat("No ATAC H5 path provided or file not found, skipping ATAC assay creation\n")
}

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


# Add user metadata if provided
if (!is.null(metadata_file) && file.exists(metadata_file)) {
  meta <- tryCatch(read.csv(metadata_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(meta) && "sampleid" %in% colnames(meta)) {
    meta_row <- meta[meta$sampleid == sample_name, , drop = FALSE]
    if (nrow(meta_row) == 1) {
      # Add all columns except sampleid as metadata (repeat for all cells)
      meta_fields <- setdiff(colnames(meta_row), "sampleid")
      for (field in meta_fields) {
        value <- meta_row[[field]][1]
        seurat_obj[[field]] <- rep(value, ncol(seurat_obj))
      }
      cat("Added user metadata fields:", paste(meta_fields, collapse=", "), "\n")
    } else {
      cat("No matching row in metadata for sampleid:", sample_name, "\n")
    }
  } else {
    cat("Metadata file missing sampleid column or could not be read\n")
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
keep <- rep(TRUE, ncol(post_obj))
if ("nFeature_RNA" %in% colnames(post_obj@meta.data)) keep <- keep & (post_obj@meta.data$nFeature_RNA >= min_features)
if ("nCount_RNA" %in% colnames(post_obj@meta.data)) keep <- keep & (post_obj@meta.data$nCount_RNA >= min_counts)
if ("percent.mt" %in% colnames(post_obj@meta.data)) keep <- keep & (post_obj@meta.data$percent.mt <= max_mito)
if ("nuclear_fraction" %in% colnames(post_obj@meta.data)) keep <- keep & (post_obj@meta.data$nuclear_fraction >= min_nuclear)
if ("doublet_class" %in% colnames(post_obj@meta.data)) {
  dc <- tolower(as.character(post_obj@meta.data$doublet_class))
  keep <- keep & (dc != "doublet")
}
cat("Applied QC filters: max_mito =", max_mito, ", min_nuclear =", min_nuclear, "\n")
cat("Cells before filtering:", ncol(seurat_obj), ", after filtering:", sum(keep), "\n")
post_obj <- subset(post_obj, cells = colnames(post_obj)[which(keep)])

# Add user metadata to post-QC object as well
if (!is.null(metadata_file) && file.exists(metadata_file)) {
  meta <- tryCatch(read.csv(metadata_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(meta) && "sampleid" %in% colnames(meta)) {
    meta_row <- meta[meta$sampleid == sample_name, , drop = FALSE]
    if (nrow(meta_row) == 1) {
      meta_fields <- setdiff(colnames(meta_row), "sampleid")
      for (field in meta_fields) {
        value <- meta_row[[field]][1]
        post_obj[[field]] <- rep(value, ncol(post_obj))
      }
      cat("Added user metadata fields to post-QC object:\n")
    }
  }
}

post_path <- file.path(getwd(), paste0(sample_name, "_seurat_object_postqc.rds"))
saveRDS(post_obj, post_path)
cat("Saved post-QC Seurat object to:", post_path, "\n")

invisible(list(pre = pre_path, post = post_path))
