#!/usr/bin/env Rscript
# doublet_detection.R
#
# Per-sample doublet detection with scDblFinder.
# Inputs:  decontX-filtered H5 (stacked S+U+A) + genome GTF
# Outputs: scdblfinder_<sampleId>.csv.gz
#
# GTF is parsed here solely to exclude rRNA / Mt_rRNA genes before
# running scDblFinder — their high expression would confound doublet scoring.

suppressPackageStartupMessages({
  library(magrittr)
  library(rhdf5)
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(strex)
  library(BiocParallel)
})

# Container isn't CPU-capped (docker --cpu-shares, not --cpus), so
# scDblFinder's default parallel backend can see the host's full core
# count and over-subscribe when several samples run concurrently.
# Cap workers to what this process is actually allocated.
register(MulticoreParam(workers = 4))

# ---------------------------------------------------------------------------
# Read H5 (stacked S+U+A)
# ---------------------------------------------------------------------------
.get_h5_mx <- function(af_mat_f) {
  h5   <- H5Fopen(af_mat_f, flags = "H5F_ACC_RDONLY")
  mat  <- sparseMatrix(
    i    = as.vector(h5$matrix$indices + 1L),
    p    = as.vector(h5$matrix$indptr),
    x    = as.vector(h5$matrix$data),
    repr = "C",
    dims = h5$matrix$shape
  )
  colnames(mat) <- as.character(h5$matrix$barcodes)
  rownames(mat) <- as.character(h5$matrix$features$name)
  H5Fclose(h5)
  mat
}

# ---------------------------------------------------------------------------
# Parse GTF — returns ensembl_id / symbol / gene_type table
# ---------------------------------------------------------------------------
.parse_gtf_annotations <- function(gtf_f) {
  decomp_cmd <- if (grepl("\\.gz$", gtf_f)) {
    paste0("zcat ", shQuote(gtf_f))
  } else {
    paste0("cat ", shQuote(gtf_f))
  }
  gtf_lines <- fread(
    cmd       = paste0(decomp_cmd, " | grep -v '^#'"),
    sep       = "\t",
    header    = FALSE,
    select    = c(3, 9),
    col.names = c("feature", "attributes")
  )
  gtf_lines <- gtf_lines[feature == "gene"]

  gtf_lines[, ensembl_id := str_match(attributes, 'gene_id "([^"]+)"')[, 2]]
  gtf_lines[, symbol     := str_match(attributes, 'gene_name "([^"]+)"')[, 2]]
  gtf_lines[, gene_type  := str_match(attributes, 'gene_type "([^"]+)"')[, 2]]

  na_type <- is.na(gtf_lines$gene_type)
  if (any(na_type)) {
    gtf_lines[na_type, gene_type := str_match(attributes, 'gene_biotype "([^"]+)"')[, 2]]
  }

  gtf_lines[, attributes := NULL]
  gtf_lines[, feature    := NULL]
  unique(gtf_lines, by = "ensembl_id")
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
run_doublet_detection <- function() {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3) {
    stop("Usage: doublet_detection.R <sampleId> <h5_file> <genome_gtf>")
  }

  sample_id <- args[1]
  h5_file   <- args[2]
  gtf_file  <- args[3]

  message("=== DOUBLET_DETECTION: ", sample_id, " ===")
  message("H5:  ", h5_file)
  message("GTF: ", gtf_file)

  # -------------------------------------------------------------------------
  # 1. Load filtered H5
  # -------------------------------------------------------------------------
  mat <- .get_h5_mx(h5_file)
  message("Loaded H5: ", nrow(mat), " features x ", ncol(mat), " barcodes")

  # -------------------------------------------------------------------------
  # 2. Sum S + U + A into a single counts matrix
  # -------------------------------------------------------------------------
  splice_ns  <- c("S", "U", "A")
  usa_ls     <- lapply(splice_ns, function(s) {
    idx <- grep(paste0("_", s, "$"), rownames(mat))
    mat[idx, , drop = FALSE]
  }) %>% setNames(splice_ns)

  counts_mat <- Reduce("+", usa_ls)
  rownames(counts_mat) <- sub("_[SUA]$", "", rownames(usa_ls[["S"]]))

  # -------------------------------------------------------------------------
  # 3. Exclude rRNA / Mt_rRNA genes (high expression confounds doublet scoring)
  # -------------------------------------------------------------------------
  gene_annots <- .parse_gtf_annotations(gtf_file)
  rrna_ids    <- gene_annots[gene_type == "rRNA",    ensembl_id]
  mt_rrna_ids <- gene_annots[gene_type == "Mt_rRNA", ensembl_id]
  exclude_gs  <- rownames(counts_mat) %in% c(rrna_ids, mt_rrna_ids)
  message("Excluding ", sum(exclude_gs), " rRNA / Mt_rRNA genes")
  counts_mat  <- counts_mat[!exclude_gs, ]

  # -------------------------------------------------------------------------
  # 4. Drop zero-count cells
  # -------------------------------------------------------------------------
  cell_counts <- Matrix::colSums(counts_mat)
  cell_feats  <- Matrix::colSums(counts_mat > 0)
  keep_cells  <- cell_counts > 0
  counts_mat  <- counts_mat[, keep_cells, drop = FALSE]
  cell_feats  <- cell_feats[keep_cells]
  message("Non-zero cells: ", ncol(counts_mat))

  # -------------------------------------------------------------------------
  # 5. scDblFinder (require >= 100 detected genes per cell)
  # -------------------------------------------------------------------------
  keep_for_dbl <- cell_feats >= 100
  n_dbl_input  <- sum(keep_for_dbl)
  message("Running scDblFinder on ", n_dbl_input, " cells ...")

  if (n_dbl_input >= 100) {
    dbl_cells <- colnames(counts_mat)[keep_for_dbl]
    dbl_mat   <- counts_mat[, dbl_cells, drop = FALSE]
    dbl_sce   <- SingleCellExperiment(assays = list(counts = dbl_mat))
    dbl_res   <- scDblFinder(dbl_sce, returnType = "table",
                             multiSampleMode = "singleModel", verbose = FALSE,
                             BPPARAM = MulticoreParam(workers = 4))
    dbl_dt    <- as.data.table(dbl_res, keep.rownames = "cell_id")
    dbl_dt    <- dbl_dt[type == "real"]
    dbl_dt[, sample_id := sample_id]
    setnames(dbl_dt, "class", "scdbl_class")
    message("scDblFinder: ", sum(dbl_dt$scdbl_class == "doublet"),
            " doublets / ", nrow(dbl_dt), " cells")
  } else {
    message("Too few cells (", n_dbl_input, ") for scDblFinder — skipping")
    dbl_dt <- data.table(
      cell_id     = colnames(counts_mat),
      sample_id   = sample_id,
      scdbl_class = "singlet",
      scdbl_score = NA_real_
    )
  }

  out_f <- sprintf("scdblfinder_%s.csv.gz", sample_id)
  fwrite(dbl_dt, out_f)
  message("Written: ", out_f)
  message("=== DOUBLET_DETECTION done ===")
}

run_doublet_detection()
