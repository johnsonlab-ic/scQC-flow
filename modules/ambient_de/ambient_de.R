#!/usr/bin/env Rscript
# ambient_de.R
#
# Pseudobulk empty-droplet profiles + EdgeR DE (empty vs cells).
# Mirrors scprocess mlm_empties analysis exactly.
#
# Inputs (staged by Nextflow into the work directory):
#   barcode_matrix_<sampleId>.h5    — raw stacked S+U+A H5 (BARCODE_ESTIMATION)
#   knee_plot_data_<sampleId>.csv   — empty plateau flags (in_empty_plateau column)
#   filt_counts_<sampleId>.h5       — decontX-filtered S+U+A H5 (DECONTX)
#   <genome.gtf[.gz]>               — genome annotation (for SYMBOL_ENSEMBL mapping)
#
# Outputs:
#   edger_dt.csv.gz   — EdgeR DE: gene_id (SYMBOL_ENSEMBL), logFC, unshrunk.logFC,
#                        logCPM, PValue, FDR, mean_logcpm.cells, mean_logcpm.empties,
#                        is_ambient
#   pb_empties.rds    — SummarizedExperiment: pseudobulk empty profiles
#                        (rows = SYMBOL_ENSEMBL genes, cols = samples)

suppressPackageStartupMessages({
  library(magrittr)
  library(rhdf5)
  library(Matrix)
  library(data.table)
  library(SummarizedExperiment)
  library(edgeR)
  library(strex)
})

# ---------------------------------------------------------------------------
# Read H5 (stacked S+U+A format — mirrors .get_h5_mx in decontx.R)
# ---------------------------------------------------------------------------
.get_h5_mx <- function(h5_path) {
  h5  <- H5Fopen(h5_path, flags = "H5F_ACC_RDONLY")
  mat <- sparseMatrix(
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
# Sum _S + _U + _A rows into per-gene counts.
# Strips the suffix from each modality, aligns by gene, and sums.
# Returns a sparse matrix with plain Ensembl IDs as rownames.
# ---------------------------------------------------------------------------
.sum_sua <- function(mat) {
  result <- NULL
  for (sfx in c("_S", "_U", "_A")) {
    idx <- which(endsWith(rownames(mat), sfx))
    if (length(idx) == 0L) next
    m   <- mat[idx, , drop = FALSE]
    rownames(m) <- sub(paste0("\\", sfx, "$"), "", rownames(m))
    result <- if (is.null(result)) m else result + m[rownames(result), , drop = FALSE]
  }
  result
}

# ---------------------------------------------------------------------------
# Parse GTF — identical to doublet_detection.R
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
    gtf_lines[na_type, gene_type :=
      str_match(attributes[na_type], 'gene_biotype "([^"]+)"')[, 2]]
  }
  gtf_lines[, c("attributes", "feature") := NULL]
  unique(gtf_lines, by = "ensembl_id")
}

# ===========================================================================
# Main
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(1L, 2L))) stop("Usage: ambient_de.R <genome_gtf> [subset_qc_csv_glob]")
gtf_f <- args[1]
subset_qc_pattern <- if (length(args) == 2L) args[2] else NULL

# --- Locate input files ---
raw_h5_files   <- sort(Sys.glob("barcode_matrix_*.h5"))
knee_files     <- sort(Sys.glob("knee_plot_data_*.csv"))
filt_h5_files  <- sort(Sys.glob("filt_counts_*.h5"))
empty_bc_files <- sort(Sys.glob("empty_barcodes_*.csv"))   # is_empty from CELL_CALLING

if (length(raw_h5_files)   == 0L) stop("No barcode_matrix_*.h5 files found")
if (length(knee_files)     == 0L) stop("No knee_plot_data_*.csv files found")
if (length(filt_h5_files)  == 0L) stop("No filt_counts_*.h5 files found")
if (length(empty_bc_files) == 0L) stop("No empty_barcodes_*.csv files found")

subset_qc_dt <- NULL
if (!is.null(subset_qc_pattern)) {
  subset_qc_files <- sort(Sys.glob(subset_qc_pattern))
  if (length(subset_qc_files) == 0L) {
    stop("No subset QC files matched pattern: ", subset_qc_pattern)
  }
  subset_qc_dt <- rbindlist(lapply(subset_qc_files, fread), use.names = TRUE, fill = TRUE)
  subset_qc_dt[, sample_id := as.character(sample_id)]
  subset_qc_dt[, cell_id := as.character(cell_id)]
}

# Extract sample IDs from raw H5 filenames (canonical order)
sample_ids <- sub("\\.h5$", "", sub("^barcode_matrix_", "", basename(raw_h5_files)))
message("Samples: ", paste(sample_ids, collapse = ", "))

# --- GTF → SYMBOL_ENSEMBL lookup table ---
message("Parsing GTF: ", gtf_f)
gtf_dt  <- .parse_gtf_annotations(gtf_f)
sym_map <- gtf_dt[, .(
  ensembl_id,
  symbol_ensembl = paste0(
    ifelse(is.na(symbol) | symbol == "", ensembl_id, symbol),
    "_", ensembl_id
  )
)]
setkey(sym_map, ensembl_id)
message("GTF: ", nrow(sym_map), " genes parsed")

# ---------------------------------------------------------------------------
# Per-sample pseudobulk construction
# ---------------------------------------------------------------------------
empty_cols <- list()
cell_cols  <- list()

for (i in seq_along(sample_ids)) {
  sid    <- sample_ids[i]
  raw_f  <- raw_h5_files[i]
  knee_f <- grep(sid, knee_files,    value = TRUE, fixed = TRUE)[1]
  filt_f <- grep(sid, filt_h5_files, value = TRUE, fixed = TRUE)[1]

  if (is.na(knee_f)) stop("No knee CSV for sample: ", sid)
  if (is.na(filt_f)) stop("No filtered H5 for sample: ", sid)

  message("=== ", sid, " ===")

  # --- Empties pseudobulk (caller-independent: knee empty plateau) ---
  # scprocess (.get_one_empty_pb) defines empties as in_empty_plateau==TRUE from the
  # knee data — a low-UMI plateau that is stable across cell-callers and sample
  # quality. This is deliberately NOT the cell-caller's is_empty set (which extends
  # to the GMM cell_floor and drifts on bad samples).
  raw_mat   <- .get_h5_mx(raw_f)
  knee_dt   <- fread(knee_f)
  if (!"in_empty_plateau" %in% names(knee_dt))
    stop("knee CSV lacks in_empty_plateau column: ", knee_f)
  empty_bcs <- intersect(knee_dt[in_empty_plateau == TRUE, barcode], colnames(raw_mat))
  message("  Empty-plateau barcodes: ", length(empty_bcs))
  if (length(empty_bcs) == 0L) stop("No empty-plateau barcodes for sample: ", sid)

  raw_sua   <- .sum_sua(raw_mat[, empty_bcs, drop = FALSE])
  empty_cols[[sid]] <- Matrix::rowSums(raw_sua)
  rm(raw_mat, raw_sua); gc()

  # --- Cells pseudobulk ---
  filt_mat  <- .get_h5_mx(filt_f)
  if (!is.null(subset_qc_dt)) {
    sel_bcs <- intersect(subset_qc_dt[sample_id == sid, unique(cell_id)], colnames(filt_mat))
    message("  Selected zoom cells: ", length(sel_bcs))
    if (length(sel_bcs) == 0L) {
      # A zoom subset can legitimately contain no cells from some samples. Skip the
      # whole sample rather than aborting the run, and drop the empties column added
      # above so the empties/cells pseudobulk matrices stay column-aligned.
      message("  No zoom cells for sample ", sid, " — skipping from ambient DE")
      empty_cols[[sid]] <- NULL
      rm(filt_mat); gc()
      next
    }
    filt_mat <- filt_mat[, sel_bcs, drop = FALSE]
  }
  filt_sua  <- .sum_sua(filt_mat)
  cell_cols[[sid]] <- Matrix::rowSums(filt_sua)
  rm(filt_mat, filt_sua); gc()
}

# ---------------------------------------------------------------------------
# Build aligned pseudobulk matrices (genes × samples)
# ---------------------------------------------------------------------------
if (length(cell_cols) == 0L)
  stop("No samples had cells in this subset — cannot build cells pseudobulk for ambient DE")

all_genes <- Reduce(union, lapply(empty_cols, names))
message("Total genes (union): ", length(all_genes))

make_mat <- function(col_list, genes) {
  m <- matrix(0L, nrow = length(genes), ncol = length(col_list),
               dimnames = list(genes, names(col_list)))
  for (sid in names(col_list)) {
    v <- col_list[[sid]]
    shared <- intersect(genes, names(v))
    m[shared, sid] <- as.integer(round(v[shared]))
  }
  m
}

empty_mat <- make_mat(empty_cols, all_genes)
cell_mat  <- make_mat(cell_cols,  all_genes)

# Apply SYMBOL_ENSEMBL rownames (fall back to ensembl_id for unmapped genes)
sym_vec <- sym_map[J(all_genes)]$symbol_ensembl
sym_vec[is.na(sym_vec)] <- all_genes[is.na(sym_vec)]
rownames(empty_mat) <- sym_vec
rownames(cell_mat)  <- sym_vec

# ---------------------------------------------------------------------------
# Save pb_empties SummarizedExperiment
# ---------------------------------------------------------------------------
pb_se <- SummarizedExperiment(
  assays  = list(counts = empty_mat),
  colData = DataFrame(sample_id = sample_ids, row.names = sample_ids)
)
saveRDS(pb_se, "pb_empties.rds")
message("Written pb_empties.rds  (", nrow(pb_se), " genes x ", ncol(pb_se), " samples)")

# ---------------------------------------------------------------------------
# EdgeR DE: empty vs cells
#
# Combined matrix columns named empty_<sid> / cell_<sid> to avoid duplicates.
# group is releveled so "cell" is the reference → groupempty coef = empty − cell.
# Positive logFC → higher in empty droplets → ambient gene.
# ---------------------------------------------------------------------------
n_s          <- length(sample_ids)
combined_mat <- cbind(empty_mat, cell_mat)
colnames(combined_mat) <- c(paste0("empty_", sample_ids), paste0("cell_", sample_ids))

group  <- factor(c(rep("empty", n_s), rep("cell", n_s)))
group  <- relevel(group, ref = "cell")

dge    <- DGEList(counts = combined_mat, group = group)
keep   <- filterByExpr(dge, group = group)
dge    <- dge[keep, , keep.lib.sizes = FALSE]
message("Genes after filterByExpr: ", nrow(dge))

dge    <- calcNormFactors(dge)
design <- model.matrix(~ group)
dge    <- estimateDisp(dge, design)   # pass design so trended dispersion is used

fit    <- glmQLFit(dge, design)
qlf    <- glmQLFTest(fit, coef = "groupempty")

# topTags returns logFC, [unshrunk.logFC], logCPM, F, PValue, FDR
tt_dt  <- as.data.table(
  topTags(qlf, n = Inf, sort.by = "none")$table,
  keep.rownames = "gene_id"
)

# Older EdgeR versions omit unshrunk.logFC — fall back to logFC
if (!"unshrunk.logFC" %in% names(tt_dt)) {
  tt_dt[, unshrunk.logFC := logFC]
}

# Per-group mean log-CPM
logcpm    <- cpm(dge, log = TRUE)
empty_idx <- grep("^empty_", colnames(logcpm))
cell_idx  <- grep("^cell_",  colnames(logcpm))

m_idx <- match(tt_dt$gene_id, rownames(logcpm))
tt_dt[, mean_logcpm.cells   := rowMeans(logcpm[m_idx, cell_idx,  drop = FALSE])]
tt_dt[, mean_logcpm.empties := rowMeans(logcpm[m_idx, empty_idx, drop = FALSE])]

# is_ambient: FDR < 0.01 AND logFC > 0 (higher in empty droplets)
tt_dt[, is_ambient := FDR < 0.01 & logFC > 0]

# Final column order matching scprocess edger_empty_genes output
setcolorder(tt_dt, c("gene_id", "logFC", "unshrunk.logFC", "logCPM",
                      "PValue", "FDR", "mean_logcpm.cells", "mean_logcpm.empties",
                      "is_ambient"))

message("Ambient genes (FDR<0.01, logFC>0): ", sum(tt_dt$is_ambient))
fwrite(tt_dt, "edger_dt.csv.gz")
message("Written edger_dt.csv.gz  (", nrow(tt_dt), " genes)")
message("=== AMBIENT_DE done ===")
