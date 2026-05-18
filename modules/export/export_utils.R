suppressPackageStartupMessages({
  library(assertthat)
  library(data.table)
  library(Matrix)
  library(rhdf5)
  library(stringr)
  library(Seurat)
  tryCatch(library(BPCells), error = function(e) {
    warning("BPCells not available; using standard matrices for combined objects")
  })
})

parse_gtf_annotations_export <- function(gtf_f) {
  decomp_cmd <- if (grepl("\\.gz$", gtf_f)) {
    paste0("zcat ", shQuote(gtf_f))
  } else {
    paste0("cat ", shQuote(gtf_f))
  }

  gtf_lines <- fread(
    cmd = paste0(decomp_cmd, " | grep -v '^#'"),
    sep = "\t",
    header = FALSE,
    select = c(3, 9),
    col.names = c("feature", "attributes")
  )
  gtf_lines <- gtf_lines[feature == "gene"]

  gtf_lines[, gene_id := str_match(attributes, 'gene_id "([^"]+)"')[, 2]]
  gtf_lines[, symbol := str_match(attributes, 'gene_name "([^"]+)"')[, 2]]
  gtf_lines[, gene_type := str_match(attributes, 'gene_type "([^"]+)"')[, 2]]
  na_type <- is.na(gtf_lines$gene_type)
  if (any(na_type)) {
    gtf_lines[na_type, gene_type := str_match(attributes, 'gene_biotype "([^"]+)"')[, 2]]
  }

  out_dt <- unique(gtf_lines[, .(gene_id, symbol, gene_type)], by = "gene_id")
  out_dt[is.na(symbol) | symbol == "", symbol := gene_id]
  out_dt[is.na(gene_type) | gene_type == "", gene_type := "unknown"]
  out_dt[, gene_key := paste(symbol, gene_id, sep = "_")]
  out_dt
}

read_sparse_h5_matrix_export <- function(h5_f) {
  h5 <- H5Fopen(h5_f, flags = "H5F_ACC_RDONLY")
  on.exit(H5Fclose(h5))

  mat <- sparseMatrix(
    i = as.vector(h5$matrix$indices + 1L),
    p = as.vector(h5$matrix$indptr),
    x = as.vector(h5$matrix$data),
    repr = "C",
    dims = h5$matrix$shape
  )
  colnames(mat) <- as.character(h5$matrix$barcodes)
  rownames(mat) <- as.character(h5$matrix$features$name)
  mat
}

sum_sua_counts_export <- function(stacked_mat) {
  suffixes <- c("S", "U", "A")
  mats <- lapply(suffixes, function(suffix) {
    idx <- grep(paste0("_", suffix, "$"), rownames(stacked_mat))
    stacked_mat[idx, , drop = FALSE]
  })
  names(mats) <- suffixes

  if (any(vapply(mats, nrow, integer(1)) == 0L)) {
    stop("Expected stacked S/U/A rows in filtered H5 matrix")
  }

  counts_mat <- Reduce(`+`, mats)
  rownames(counts_mat) <- sub("_[SUA]$", "", rownames(mats[[1]]))
  counts_mat
}

infer_sample_id_from_h5_export <- function(h5_f) {
  sample_id <- basename(h5_f)
  sample_id <- sub("^filt_counts_", "", sample_id)
  sub("\\.h5$", "", sample_id)
}

combine_duplicate_columns_export <- function(dt, suffix) {
  dup_cols <- grep(paste0(suffix, "$"), names(dt), value = TRUE)
  for (dup_col in dup_cols) {
    base_col <- sub(paste0(suffix, "$"), "", dup_col)
    if (base_col %in% names(dt)) {
      dt[is.na(get(base_col)), (base_col) := get(dup_col)]
      dt[, (dup_col) := NULL]
    } else {
      setnames(dt, dup_col, base_col)
    }
  }
  dt
}

load_qc_metadata_export <- function(qc_pattern) {
  qc_files <- Sys.glob(qc_pattern)
  assert_that(length(qc_files) > 0, msg = sprintf("No QC files matched '%s'", qc_pattern))
  qc_dt <- rbindlist(lapply(qc_files, fread), use.names = TRUE, fill = TRUE)
  qc_dt[, sample_id := as.character(sample_id)]
  qc_dt[, cell_id := paste(sample_id, as.character(cell_id), sep = "_")]
  qc_dt[, `:=`(
    n_counts = sum,
    n_features = detected,
    splice_fraction = (total_spliced + 1) / (total_spliced + total_unspliced + 2),
    mito_fraction = mito_pct
  )]
  qc_dt
}

load_annotation_metadata_export <- function(annotation_pattern) {
  ann_files <- Sys.glob(annotation_pattern)
  ann_files <- ann_files[basename(ann_files) != "NO_FILE"]
  if (length(ann_files) == 0) {
    return(NULL)
  }

  ann_list <- lapply(ann_files, function(path) {
    ann_dt <- fread(path)
    if (!"cell_id" %in% names(ann_dt)) {
      stop(sprintf("Annotation metadata file is missing cell_id: %s", path))
    }
    ann_dt[, cell_id := as.character(cell_id)]
    unique(ann_dt, by = "cell_id")
  })

  merged_dt <- Reduce(function(x, y) {
    out_dt <- merge(x, y, by = "cell_id", all = TRUE, suffixes = c("", "__ann"), sort = FALSE)
    combine_duplicate_columns_export(out_dt, "__ann")
  }, ann_list)

  merged_dt
}

load_export_metadata <- function(integration_csv, qc_pattern, annotation_pattern) {
  int_dt <- fread(integration_csv)
  int_dt[, cell_id := as.character(cell_id)]
  int_dt <- unique(int_dt, by = "cell_id")

  qc_dt <- load_qc_metadata_export(qc_pattern)
  merged_dt <- merge(int_dt, qc_dt, by = "cell_id", all.x = TRUE, suffixes = c("", "__qc"), sort = FALSE)
  merged_dt <- combine_duplicate_columns_export(merged_dt, "__qc")

  ann_dt <- load_annotation_metadata_export(annotation_pattern)
  if (!is.null(ann_dt)) {
    merged_dt <- merge(merged_dt, ann_dt, by = "cell_id", all.x = TRUE, suffixes = c("", "__ann"), sort = FALSE)
    merged_dt <- combine_duplicate_columns_export(merged_dt, "__ann")
  }

  if (!"sample_id" %in% names(merged_dt)) {
    stop("Merged export metadata is missing sample_id")
  }
  if (any(is.na(merged_dt$sample_id))) {
    missing_ids <- merged_dt[is.na(sample_id), head(cell_id, 5)]
    stop(sprintf("sample_id missing for exported cells: %s", paste(missing_ids, collapse = ", ")))
  }

  merged_dt[, sample_id := as.character(sample_id)]
  merged_dt
}

build_export_counts <- function(h5_pattern, obs_dt) {
  h5_files <- Sys.glob(h5_pattern)
  assert_that(length(h5_files) > 0, msg = sprintf("No filtered H5 files matched '%s'", h5_pattern))

  sample_ids <- unique(as.character(obs_dt$sample_id))
  sample_ids <- trimws(sample_ids)
  sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]
  if (length(sample_ids) == 0L) {
    stop("No valid sample_id values found for export after filtering NA/blank IDs")
  }

  h5_map <- setNames(h5_files, vapply(h5_files, infer_sample_id_from_h5_export, character(1)))
  missing_samples <- setdiff(sample_ids, names(h5_map))
  if (length(missing_samples) > 0L) {
    stop(sprintf(
      "Missing filtered H5 files for export samples: %s",
      paste(head(missing_samples, 10), collapse = ", ")
    ))
  }

  gene_keys <- NULL
  mats <- vector("list", length(sample_ids))
  ordered_ids <- character()

  for (idx in seq_along(sample_ids)) {
    sample_id_val <- sample_ids[[idx]]
    counts_mat <- read_sparse_h5_matrix_export(h5_map[[sample_id_val]])
    counts_mat <- sum_sua_counts_export(counts_mat)

    if (is.null(gene_keys)) {
      gene_keys <- rownames(counts_mat)
    } else {
      assert_that(identical(gene_keys, rownames(counts_mat)), msg = "Gene keys differ across filtered H5 files")
    }

    sample_obs <- obs_dt[sample_id == sample_id_val]
    global_ids <- paste(sample_id_val, colnames(counts_mat), sep = "_")
    index_map <- setNames(seq_along(global_ids), global_ids)
    wanted_ids <- sample_obs$cell_id
    missing_ids <- setdiff(wanted_ids, names(index_map))
    if (length(missing_ids) > 0) {
      stop(sprintf("Filtered H5 for sample %s is missing cells: %s", sample_id_val, paste(head(missing_ids, 5), collapse = ", ")))
    }

    sel_idx <- unname(index_map[wanted_ids])
    mats[[idx]] <- counts_mat[, sel_idx, drop = FALSE]
    ordered_ids <- c(ordered_ids, wanted_ids)
  }

  combined <- Reduce(cbind2, mats)
  colnames(combined) <- ordered_ids
  list(counts = combined, gene_keys = gene_keys, ordered_ids = ordered_ids)
}

build_gene_metadata_export <- function(gene_keys, gtf_dt) {
  var_dt <- data.table(gene_key = gene_keys)
  var_dt <- merge(var_dt, gtf_dt[, .(gene_key, gene_id, symbol, gene_type)], by = "gene_key", all.x = TRUE, sort = FALSE)
  missing_idx <- is.na(var_dt$gene_id)
  if (any(missing_idx)) {
    fallback <- str_match(var_dt$gene_key[missing_idx], "^(.*)_(ENSG[^_]+|ENSMUSG[^_]+)$")
    var_dt[missing_idx, symbol := ifelse(is.na(fallback[, 2]), gene_key, fallback[, 2])]
    var_dt[missing_idx, gene_id := ifelse(is.na(fallback[, 3]), gene_key, fallback[, 3])]
    var_dt[missing_idx, gene_type := "unknown"]
  }
  var_dt
}

convert_seurat_to_bpcells <- function(seu) {
  merged <- seu

  # Seurat v5: collapse per-object layers into unified layers when available.
  tryCatch({
    if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
      merged <- SeuratObject::JoinLayers(merged)
    }
  }, error = function(e) {
    message(sprintf("  Note: JoinLayers skipped (%s)", e$message))
  })

  # Attempt to use BPCells for memory efficiency if available
  tryCatch({
    if (requireNamespace("BPCells", quietly = TRUE)) {
      message("  Converting counts to BPCells format for memory efficiency...")
      counts_mat <- GetAssayData(merged, assay = DefaultAssay(merged), slot = "counts")

      # Create temporary BPCells matrix directory
      bpcells_dir <- tempfile(pattern = "bpcells_")
      dir.create(bpcells_dir)

      # Convert sparse matrix to BPCells and read back
      bpcells_mat <- BPCells::write_matrix_dir(
        mat = counts_mat,
        dir = file.path(bpcells_dir, "counts"),
        overwrite = TRUE
      )
      bpcells_mat_read <- BPCells::read_matrix_dir(file.path(bpcells_dir, "counts"))

      # Replace counts with BPCells matrix
      merged[[DefaultAssay(merged)]]@counts <- bpcells_mat_read

      message("  BPCells conversion complete")
    }
  }, error = function(e) {
    message(sprintf("  Note: BPCells not available or conversion failed; using standard matrices (%s)", e$message))
  })

  merged
}

merge_with_bpcells <- function(seu_list) {
  assert_that(length(seu_list) > 0, msg = "No Seurat objects provided for merging")

  if (length(seu_list) == 1) {
    return(convert_seurat_to_bpcells(seu_list[[1]]))
  }

  message(sprintf("  Merging %d Seurat objects", length(seu_list)))

  # Standard merge to combine metadata and reductions
  merged <- seu_list[[1]]
  for (i in seq(2, length(seu_list))) {
    merged <- merge(merged, seu_list[[i]], add.cell.ids = NULL, merge.data = TRUE)
  }

  convert_seurat_to_bpcells(merged)
}
