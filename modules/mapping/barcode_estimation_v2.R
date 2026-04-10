#!/usr/bin/env Rscript
# barcode_estimation_v2.R
#
# Like barcode_estimation.R but replaces every call to the installed
# barcodeRanks() with a standalone reimplementation of DropletUtils 1.26.0.
# This exactly mirrors the scprocess conda environment behaviour.
#
# Produces (no H5 — already written by barcode_estimation.R):
#   - knee_plot_data_v2_<id>.csv
#   - ambient_params_v2_<id>.env

suppressPackageStartupMessages({
  library(magrittr)
  library(fishpond)
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
  library(strex)
})

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------

run_barcode_estimation_v2 <- function() {

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: barcode_estimation_v2.R <sampleId> <quantDir> <kneeOut> <ambientEnvOut>")
}

sample_id       <- args[1]
quant_dir       <- args[2]
knee_out        <- args[3]
ambient_env_out <- args[4]

message("Sample:   ", sample_id)
message("Quant:    ", quant_dir)
message("Knee out: ", knee_out)
message("Params:   ", ambient_env_out)

# ---------------------------------------------------------------------------
# Load alevin-fry output (S, U, A modes)
# ---------------------------------------------------------------------------

sce <- loadFry(
  quant_dir,
  outputFormat = list(S = c("S"), U = c("U"), A = c("A"))
)

mat <- assayNames(sce) %>% lapply(function(n) {
  m            <- assay(sce, n)
  rownames(m)  <- paste0(rownames(m), "_", n)
  m
}) %>% do.call(rbind, .)

mat <- mat[, colSums(mat) > 0]
message("Barcodes retained (total > 0): ", ncol(mat))

splice_dt <- data.table(
  barcode   = colnames(sce),
  spliced   = colSums(assay(sce, "S")),
  unspliced = colSums(assay(sce, "U"))
)

# ---------------------------------------------------------------------------
# Knee / ambient parameter estimation using DropletUtils 1.26 logic
# ---------------------------------------------------------------------------

bender_ps <- calc_ambient_params_v2(
  split_mat = mat,
  run       = sample_id,
  run_var   = "sample_id"
)

bender_ps <- merge(bender_ps, splice_dt, by = "barcode") %>%
  .[order(rank)]

fwrite(bender_ps, file = knee_out)
message("Written knee data (v2): ", knee_out)

amb <- list(
  run                        = sample_id,
  cb_total_droplets_included = as.integer(unique(bender_ps$total_droplets_included)),
  cb_expected_cells          = as.integer(unique(bender_ps$expected_cells)),
  cb_low_count_threshold     = as.integer(unique(bender_ps$low_count_threshold)),
  knee1                      = as.integer(unique(bender_ps$knee1)),
  shin1                      = as.integer(unique(bender_ps$shin1)),
  knee2                      = as.integer(unique(bender_ps$knee2)),
  shin2                      = as.integer(unique(bender_ps$shin2))
)
writeLines(c(
  sprintf("RUN=%s",                           amb$run),
  sprintf("CB_TOTAL_DROPLETS_INCLUDED=%s",    amb$cb_total_droplets_included),
  sprintf("CB_EXPECTED_CELLS=%s",             amb$cb_expected_cells),
  sprintf("CB_LOW_COUNT_THRESHOLD=%s",        amb$cb_low_count_threshold),
  sprintf("KNEE1=%s",                         amb$knee1),
  sprintf("SHIN1=%s",                         amb$shin1),
  sprintf("KNEE2=%s",                         amb$knee2),
  sprintf("SHIN2=%s",                         amb$shin2)
), ambient_env_out)
message("Written ambient params (v2): ", ambient_env_out)
}


# ===========================================================================
# DropletUtils 1.26.0 barcodeRanks — standalone reimplementation
#
# Key differences from installed 1.30:
#   - inflection = right endpoint of the steepest-derivative zone (not midpoint)
#   - knee       = minimum signed curvature of smooth.spline on that zone
#   - returns a plain list (not a BarcodeRanks DataFrame)
# ===========================================================================

.find_curve_bounds_v126 <- function(x, y, exclude.from) {
  d1n  <- diff(y) / diff(x)
  skip <- min(length(d1n) - 1L, sum(x <= log10(exclude.from)))
  d1n  <- tail(d1n, length(d1n) - skip)
  right.edge <- which.min(d1n)
  left.edge  <- which.max(d1n[seq_len(right.edge)])
  c(left = left.edge, right = right.edge) + skip
}

barcodeRanks_v126 <- function(m, lower = 100, exclude.from = 50, df = 20) {
  totals_raw <- Matrix::colSums(m)
  bc_names   <- names(totals_raw)            # preserve before as.numeric() strips them
  totals     <- as.numeric(totals_raw)       # as.numeric() required for rle()
  o          <- order(totals, decreasing = TRUE)

  stuff      <- rle(totals[o])
  run.rank   <- cumsum(stuff$lengths) - (stuff$lengths - 1) / 2
  run.totals <- stuff$values

  keep <- run.totals > lower
  keep[run.rank <= exclude.from] <- FALSE
  if (sum(keep) < 2L) stop("insufficient unique points for v126 barcodeRanks")

  y <- log10(run.totals[keep])
  x <- log10(run.rank[keep])

  edge.out   <- .find_curve_bounds_v126(x = x, y = y, exclude.from = exclude.from)
  inflection <- 10^(y[edge.out["right"]])

  new.keep <- edge.out["left"]:edge.out["right"]
  if (length(new.keep) >= 4L) {
    fit       <- smooth.spline(x[new.keep], y[new.keep], df = df)
    curvature <- predict(fit, deriv = 2)$y / (1 + predict(fit, deriv = 1)$y^2)^1.5
    knee      <- 10^(y[new.keep][which.min(curvature)])
  } else {
    knee <- 10^(y[new.keep[1]])
  }

  # Reconstruct full-length rank/total in original barcode order
  all.ranks  <- rep(run.rank,   stuff$lengths)
  all.totals <- rep(run.totals, stuff$lengths)
  result_ranks  <- numeric(length(totals))
  result_totals <- numeric(length(totals))
  result_ranks[o]  <- all.ranks
  result_totals[o] <- all.totals

  list(
    knee       = knee,
    inflection = inflection,
    rank       = setNames(result_ranks,  bc_names),
    total      = setNames(result_totals, bc_names)
  )
}


# ---------------------------------------------------------------------------
# Helper: convert barcodeRanks_v126 list output to ranked data.table
# ---------------------------------------------------------------------------

.v126_to_ranks_dt <- function(res) {
  data.table(
    barcode = names(res$rank),
    rank    = as.numeric(res$rank),
    total   = as.numeric(res$total)
  )[order(rank)]
}


# ===========================================================================
# Functions mirroring scprocess mapping.R (all barcodeRanks calls → v126)
# ===========================================================================

calc_ambient_params_v2 <- function(split_mat, run,
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

  knee1_ls  <- .get_knee_and_shin_1_v2(split_mat, min_umis_cells, knee1, shin1, knee2)
  knee2_ls  <- .get_knee_and_shin_2_v2(split_mat, knee1_ls$ranks_dt,
                                        rank_empty_plateau, min_umis_empty,
                                        knee1_ls$shin1_x, knee2, shin2)
  params_ls <- .get_params_ls_v2(knee1_ls, knee2_ls, low_count_threshold,
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

  bender_ps <- .get_empty_plateau_v2(
    knee_df        = bender_ps,
    shin1          = knee1_ls$sel_knee["shin"],
    total_included = params_ls$total_included,
    knee2          = knee2_ls$sel_knee["knee"]
  )

  return(bender_ps)
}


.get_empty_plateau_v2 <- function(knee_df, shin1, total_included, knee2) {
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


.get_knee_and_shin_1_v2 <- function(split_mat, min_umis_cells,
                                     knee1 = NA, shin1 = NA, knee2 = NA) {
  if (all(sapply(c(knee1, shin1, knee2), function(p) !is.na(p)))) {
    low <- median(c(knee2, shin1))
    res <- barcodeRanks_v126(split_mat, lower = low)
    sel_knee <- c(shin = shin1, knee = knee1)

  } else if (!is.null(min_umis_cells)) {
    res <- barcodeRanks_v126(split_mat, lower = min_umis_cells)
    sel_knee <- c(
      shin = as.integer(round(as.numeric(res$inflection))),
      knee = as.integer(round(as.numeric(res$knee)))
    )

  } else {
    results <- lapply(seq(1000, 100, by = -100),
                      function(x) barcodeRanks_v126(split_mat, lower = x))
    ks_strs <- sapply(results, function(r)
      paste0(
        as.integer(round(as.numeric(r$inflection))), "_",
        as.integer(round(as.numeric(r$knee)))
      ))

    votes_tbl <- table(ks_strs)
    sel_cut   <- names(votes_tbl)[which.max(votes_tbl)]
    sel_knee  <- strsplit(sel_cut, "_")[[1]] %>% as.integer() %>%
      setNames(c("shin", "knee"))
    res <- results[[which(ks_strs == sel_cut)[1]]]
  }

  ranks_dt  <- .v126_to_ranks_dt(res)
  shin1_idx <- which.min(abs(ranks_dt$total - sel_knee[1]))[1]
  shin1_x   <- ranks_dt[shin1_idx, rank]

  list(ranks_dt = ranks_dt, sel_knee = sel_knee, shin1_x = shin1_x)
}


.get_knee_and_shin_2_v2 <- function(split_mat, ranks_dt, rank_empty_plateau,
                                     min_umis_empty, shin1_x, knee2, shin2) {
  if (all(sapply(c(knee2, shin2), function(p) !is.na(p)))) {
    shin2_corr <- ranks_dt[which.min(abs(ranks_dt$total - shin2))[1], total]
    knee2_corr <- ranks_dt[which.min(abs(ranks_dt$total - knee2))[1], total]
    sel_knee   <- c(shin = shin2_corr, knee = knee2_corr)

  } else if (!is.null(rank_empty_plateau)) {
    smol <- ranks_dt[rank > rank_empty_plateau, barcode]
    res  <- barcodeRanks_v126(split_mat[, smol], lower = min_umis_empty)
    sel_knee <- c(
      shin = as.integer(round(as.numeric(res$inflection))),
      knee = as.integer(round(as.numeric(res$knee)))
    )

  } else {
    cuts     <- .calc_small_knee_cuts_ls_v2(ranks_dt, min_umis_empty, shin1_x)
    ranks_ls <- lapply(cuts, function(this_cut) {
      smol <- ranks_dt[rank > this_cut, barcode]
      res  <- barcodeRanks_v126(split_mat[, smol], lower = min_umis_empty)
      sh2  <- as.integer(round(as.numeric(res$inflection)))
      kn2  <- as.integer(round(as.numeric(res$knee)))
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


.calc_small_knee_cuts_ls_v2 <- function(ranks_dt, min_umis_empty, shin_x) {
  last   <- tail(unique(ranks_dt[total > min_umis_empty, total]), n = 3)[1]
  last_x <- ranks_dt[total == last, rank][1]

  middle <- copy(ranks_dt) %>%
    .[, n := .I] %>%
    .[rank %between% c(shin_x, last_x), n] %>%
    log10() %>% mean() %>% `^`(10, .)

  10^seq(log10(shin_x), log10(middle), length.out = 10)
}


.get_params_ls_v2 <- function(knee1_ls, knee2_ls, low_count_threshold,
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
          else                                 knee2_ls$sel_knee["shin"]
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


run_barcode_estimation_v2()
