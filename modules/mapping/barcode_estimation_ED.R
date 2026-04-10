#!/usr/bin/env Rscript
# barcode_estimation_ED.R
#
# Like barcode_estimation.R but uses EmptyDrops (DropletUtils) for cell calling
# instead of a simple barcode-rank knee cut.
#
# Outputs:
#   - af_counts_mat.h5          (S+U+A stacked matrix, 10x v3 — same format for decontX)
#   - ed_plot_data_<id>.csv     (per-barcode knee data + EmptyDrops stats)
#   - ambient_params_<id>.env   (KEY=VALUE; CB_EXPECTED_CELLS = EmptyDrops cell count)
#
# EmptyDrops setup (mirrors scprocess ambient.R):
#   - known.empty = barcodes in the empty plateau (from 4-point knee estimation)
#   - retain      = knee1 UMI value (barcodes above this are auto-called as cells)
#   - niters      = 10000
#   - FDR threshold for testing: 0.001
#
# Usage: Rscript barcode_estimation_ED.R <sampleId> <quantDir> <h5Out> <edOut> <ambientEnvOut>

suppressPackageStartupMessages({
  library(magrittr)
  library(fishpond)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(BiocParallel)
  library(data.table)
  library(parallel)
  library(strex)
})

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

run_barcode_estimation_ED <- function() {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5) {
    stop("Usage: barcode_estimation_ED.R <sampleId> <quantDir> <h5Out> <edOut> <ambientEnvOut>")
  }

  sample_id       <- args[1]
  quant_dir       <- args[2]
  h5_out          <- args[3]
  ed_out          <- args[4]
  ambient_env_out <- args[5]

  message("Sample:   ", sample_id)
  message("Quant:    ", quant_dir)
  message("H5 out:   ", h5_out)
  message("ED out:   ", ed_out)
  message("Params:   ", ambient_env_out)

  # -------------------------------------------------------------------------
  # Load alevin-fry output (S, U, A modes)
  # -------------------------------------------------------------------------

  sce <- loadFry(
    quant_dir,
    outputFormat = list(S = c("S"), U = c("U"), A = c("A"))
  )

  # Gene-level S+U+A sum — input for emptyDrops
  total_counts <- assay(sce, "S") + assay(sce, "U") + assay(sce, "A")

  # Stacked S/U/A matrix with prefixed row names — saved as H5 for decontX
  mat <- assayNames(sce) %>% lapply(function(n) {
    m            <- assay(sce, n)
    rownames(m)  <- paste0(rownames(m), "_", n)
    m
  }) %>% do.call(rbind, .)

  # Drop all-zero barcodes (consistent with barcode_estimation.R)
  keep_bcs     <- colSums(mat) > 0
  mat          <- mat[, keep_bcs]
  total_counts <- total_counts[, keep_bcs]
  message("Barcodes retained (total > 0): ", ncol(mat))

  # Save H5 (stacked S/U/A for downstream decontX)
  write10xCounts(h5_out, mat, version = "3", overwrite = TRUE)
  message("Written H5: ", h5_out)

  # Per-barcode spliced / unspliced counts
  splice_dt <- data.table(
    barcode   = colnames(sce)[keep_bcs],
    spliced   = colSums(assay(sce, "S")[, keep_bcs]),
    unspliced = colSums(assay(sce, "U")[, keep_bcs])
  )

  # -------------------------------------------------------------------------
  # 4-point knee estimation — provides empty plateau + knee1 for emptyDrops
  # colSums(stacked mat) == colSums(total_counts) so barcodeRanks is equivalent
  # -------------------------------------------------------------------------

  bender_ps <- calc_ambient_params(
    split_mat = mat,
    run       = sample_id,
    run_var   = "sample_id"
  )

  knee1_val  <- as.integer(unique(bender_ps$knee1))
  empty_bcs  <- bender_ps[in_empty_plateau == TRUE, barcode]
  empty_idx  <- which(colnames(total_counts) %in% empty_bcs)

  message("Knee1 (retain threshold for emptyDrops): ", knee1_val)
  message("Empty plateau barcodes (known.empty):     ", length(empty_idx))

  # -------------------------------------------------------------------------
  # EmptyDrops cell calling
  # -------------------------------------------------------------------------

  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
  bpparam <- MulticoreParam(workers = max(1L, n_cores), progressbar = FALSE)

  message("Running emptyDrops (retain=", knee1_val, ", niters=10000, ncores=", n_cores, ") ...")
  ed_res <- emptyDrops(
    m           = total_counts,
    niters      = 10000L,
    BPPARAM     = bpparam,
    known.empty = empty_idx,
    retain      = knee1_val
  )

  # Convert to data.table; prefix column names to avoid collision with knee data
  ed_dt <- as.data.table(as.data.frame(ed_res), keep.rownames = TRUE)
  setnames(ed_dt, "rn", "barcode")
  setnames(ed_dt,
           c("Total",    "LogProb",    "PValue",    "Limited",    "FDR"),
           c("ED_Total", "ED_LogProb", "ED_PValue", "ED_Limited", "ED_FDR"))

  # is_cell:
  #   - barcodes with Total >= knee1 passed to retain → FDR = NA, always cells
  #   - barcodes below retain that pass the FDR test
  ed_dt[, is_cell := (!is.na(ED_FDR) & ED_FDR <= 0.001) |
                     (!is.na(ED_Total) & ED_Total >= knee1_val)]

  n_cells <- sum(ed_dt$is_cell, na.rm = TRUE)
  message("Cells called by EmptyDrops (FDR <= 0.001 or Total >= knee1): ", n_cells)

  # -------------------------------------------------------------------------
  # Merge knee data + EmptyDrops results + splice stats; override expected_cells
  # -------------------------------------------------------------------------

  out_dt <- merge(
    bender_ps,
    ed_dt[, .(barcode, ED_Total, ED_LogProb, ED_PValue, ED_Limited, ED_FDR, is_cell)],
    by = "barcode", all.x = TRUE
  )
  out_dt[, expected_cells := n_cells]   # override knee-based estimate with ED count

  out_dt <- merge(out_dt, splice_dt, by = "barcode")
  setorder(out_dt, rank)

  fwrite(out_dt, file = ed_out)
  message("Written ED plot data: ", ed_out)

  # -------------------------------------------------------------------------
  # Ambient params .env  (CB_EXPECTED_CELLS = EmptyDrops cell count)
  # -------------------------------------------------------------------------

  amb <- list(
    run                        = sample_id,
    cb_total_droplets_included = as.integer(unique(bender_ps$total_droplets_included)),
    cb_expected_cells          = n_cells,
    cb_low_count_threshold     = as.integer(unique(bender_ps$low_count_threshold)),
    knee1                      = as.integer(unique(bender_ps$knee1)),
    shin1                      = as.integer(unique(bender_ps$shin1)),
    knee2                      = as.integer(unique(bender_ps$knee2)),
    shin2                      = as.integer(unique(bender_ps$shin2))
  )

  writeLines(c(
    sprintf("RUN=%s",                        amb$run),
    sprintf("CB_TOTAL_DROPLETS_INCLUDED=%s", amb$cb_total_droplets_included),
    sprintf("CB_EXPECTED_CELLS=%s",          amb$cb_expected_cells),
    sprintf("CB_LOW_COUNT_THRESHOLD=%s",     amb$cb_low_count_threshold),
    sprintf("KNEE1=%s",                      amb$knee1),
    sprintf("SHIN1=%s",                      amb$shin1),
    sprintf("KNEE2=%s",                      amb$knee2),
    sprintf("SHIN2=%s",                      amb$shin2)
  ), ambient_env_out)
  message("Written ambient params: ", ambient_env_out)
}


# ===========================================================================
# Functions — verbatim from barcode_estimation.R
# ===========================================================================

calc_ambient_params <- function(split_mat, run,
                                min_umis_empty      = 5,
                                min_umis_cells      = NULL,
                                rank_empty_plateau  = NULL,
                                low_count_threshold = "shin2",
                                expected_cells      = NA,
                                total_included      = NA,
                                run_var             = "sample_id",
                                knee1 = NA, shin1 = NA,
                                knee2 = NA, shin2 = NA) {
  if (is.character(low_count_threshold)) {
    if (!(low_count_threshold %in% c("knee2", "shin2"))) {
      stop('low_count_threshold must be "knee2", "shin2", or an integer')
    }
  }

  knee1_ls <- .get_knee_and_shin_1(split_mat, min_umis_cells, knee1, shin1, knee2)
  knee2_ls <- .get_knee_and_shin_2(split_mat, knee1_ls$ranks_dt,
                                   rank_empty_plateau, min_umis_empty,
                                   knee1_ls$shin1_x, knee2, shin2)
  params_ls <- .get_params_ls(knee1_ls, knee2_ls, low_count_threshold,
                               expected_cells, total_included)

  bender_ps <- knee1_ls$ranks_dt %>%
    .[, (run_var) := run] %>%
    .[, `:=`(
      knee1                   = knee1_ls$sel_knee["knee"],
      shin1                   = knee1_ls$sel_knee["shin"],
      knee2                   = knee2_ls$sel_knee["knee"],
      shin2                   = knee2_ls$sel_knee["shin"],
      total_droplets_included = params_ls$total_included,
      low_count_threshold     = params_ls$lc,
      expected_cells          = params_ls$expected_cells
    )]

  bender_ps <- .get_empty_plateau(
    knee_df        = bender_ps,
    shin1          = knee1_ls$sel_knee["shin"],
    total_included = params_ls$total_included,
    knee2          = knee2_ls$sel_knee["knee"]
  )

  return(bender_ps)
}


.get_empty_plateau <- function(knee_df, shin1, total_included, knee2) {
  shin1_idx  <- which.min(abs(knee_df$total - shin1))[1]
  shin1_x    <- knee_df[shin1_idx, rank]

  empty_start <- copy(knee_df)[, n := .I] %>%
    .[rank %between% c(shin1_x, total_included), n] %>%
    log10() %>% mean() %>% (function(x) 10^x)()

  empty_end <- copy(knee_df)[total == knee2, unique(rank)]

  knee_df[, in_empty_plateau :=
    fifelse(rank %between% c(empty_start, empty_end), TRUE, FALSE)]
  return(knee_df)
}


.get_knee_and_shin_1 <- function(split_mat, min_umis_cells,
                                  knee1 = NA, shin1 = NA, knee2 = NA) {
  if (all(sapply(c(knee1, shin1, knee2), function(p) !is.na(p)))) {
    low       <- median(c(knee2, shin1))
    ranks_obj <- barcodeRanks(split_mat, lower = low)
    sel_knee  <- c(shin = shin1, knee = knee1)

  } else if (!is.null(min_umis_cells)) {
    ranks_obj <- barcodeRanks(split_mat, lower = min_umis_cells)
    sel_knee  <- c(
      shin = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )

  } else {
    ranks_ls <- lapply(seq(1000, 100, by = -100),
                       function(x) barcodeRanks(split_mat, lower = x))
    ks_strs  <- lapply(ranks_ls, function(x)
      paste0(
        as.integer(round(as.numeric(as.character(metadata(x)$inflection)))),
        "_",
        as.integer(round(as.numeric(as.character(metadata(x)$knee))))
      )) %>% unlist()

    votes_tbl <- table(ks_strs)
    sel_cut   <- names(votes_tbl)[which.max(votes_tbl)]
    sel_knee  <- strsplit(sel_cut, "_")[[1]] %>% as.integer() %>%
      setNames(c("shin", "knee"))
    ranks_obj <- ranks_ls[[which(ks_strs == sel_cut)[1]]]
  }

  ranks_dt  <- ranks_obj %>% as.data.frame() %>%
    as.data.table(keep.rownames = TRUE) %>%
    setnames("rn", "barcode") %>%
    .[order(rank)]

  shin1_idx <- which.min(abs(ranks_dt$total - sel_knee[1]))[1]
  shin1_x   <- ranks_dt[shin1_idx, rank]

  list(ranks_dt = ranks_dt, sel_knee = sel_knee, shin1_x = shin1_x)
}


.get_knee_and_shin_2 <- function(split_mat, ranks_dt, rank_empty_plateau,
                                  min_umis_empty, shin1_x, knee2, shin2) {
  if (all(sapply(c(knee2, shin2), function(p) !is.na(p)))) {
    shin2_corr <- ranks_dt[which.min(abs(ranks_dt$total - shin2))[1], total]
    knee2_corr <- ranks_dt[which.min(abs(ranks_dt$total - knee2))[1], total]
    sel_knee   <- c(shin = shin2_corr, knee = knee2_corr)

  } else if (!is.null(rank_empty_plateau)) {
    ranks_smol <- ranks_dt[rank > rank_empty_plateau, barcode]
    ranks_obj  <- barcodeRanks(split_mat[, ranks_smol], lower = min_umis_empty)
    sel_knee   <- c(
      shin = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )

  } else {
    cuts     <- .calc_small_knee_cuts_ls(ranks_dt, min_umis_empty, shin1_x)
    ranks_ls <- lapply(cuts, function(this_cut) {
      smol      <- ranks_dt[rank > this_cut, barcode]
      obj       <- barcodeRanks(split_mat[, smol], lower = min_umis_empty)
      sh2       <- as.integer(round(as.numeric(as.character(metadata(obj)$inflection))))
      kn2       <- as.integer(round(as.numeric(as.character(metadata(obj)$knee))))
      paste0(sh2, "_", kn2)
    }) %>% unlist()

    shin_tbl  <- str_before_first(ranks_ls, "_") %>% table()
    sel_i     <- names(shin_tbl)[which.max(shin_tbl)]

    match_idx <- grepl(paste0("^", sel_i, "_"), ranks_ls)
    match_ks  <- str_after_last(ranks_ls[match_idx], "_") %>% as.numeric()
    med_val   <- median(match_ks)
    sel_k     <- match_ks[which.min(abs(match_ks - med_val))[1]]

    sel_knee  <- c(shin = as.integer(sel_i), knee = as.integer(sel_k))
  }

  knee2_x <- ranks_dt[total == sel_knee["knee"], rank] %>% unique()
  list(sel_knee = sel_knee, knee2_x = knee2_x)
}


.calc_small_knee_cuts_ls <- function(ranks_dt, min_umis_empty, shin_x) {
  last   <- tail(unique(ranks_dt[total > min_umis_empty, total]), n = 3)[1]
  last_x <- ranks_dt[total == last, rank][1]

  middle <- copy(ranks_dt) %>%
    .[, n := .I] %>%
    .[rank %between% c(shin_x, last_x), n] %>%
    log10() %>% mean() %>% `^`(10, .)

  10^seq(log10(shin_x), log10(middle), length.out = 10)
}


.get_params_ls <- function(knee1_ls, knee2_ls, low_count_threshold,
                            expected_cells = NA, total_included = NA) {
  ranks_dt <- knee1_ls$ranks_dt
  shin1_x  <- knee1_ls$shin1_x
  knee2_x  <- knee2_ls$knee2_x

  if (is.na(expected_cells)) {
    expected_cells <- shin1_x
  }

  if (is.na(total_included)) {
    total_included <- copy(ranks_dt) %>%
      .[, n := .I] %>%
      .[rank %between% c(shin1_x, knee2_x), n] %>%
      log10() %>% mean() %>% `^`(10, .) %>% round()
  }

  if (is.character(low_count_threshold)) {
    lc <- if (low_count_threshold == "knee2") knee2_ls$sel_knee["knee"]
          else                                knee2_ls$sel_knee["shin"]
  } else {
    expected_x <- ranks_dt[which.min(abs(rank - expected_cells)), total]
    total_x    <- ranks_dt[which.min(abs(rank - total_included)), total]
    if (!((low_count_threshold < expected_x) & (low_count_threshold < total_x))) {
      stop("low_count_threshold must be less than expected_cells and total_droplets_included")
    }
    lc <- low_count_threshold
  }

  list(expected_cells = expected_cells, total_included = total_included, lc = lc)
}


run_barcode_estimation_ED()
