#!/usr/bin/env Rscript
# run_barcode_estimation.R
#
# Loads a simpleaf/alevin-fry quant directory and produces four outputs:
#
#   1. <sample_id>_counts.h5           — unfiltered total (S+U+A) count matrix in
#                                        10x v3 format; input for CellBender or decontX.
#   2. <sample_id>_nuclear_fraction.csv — per-barcode nuclear_fraction (U / (S+U)),
#                                        spliced and unspliced totals kept for reference.
#   3. <sample_id>_knee_data.csv       — per-barcode knee estimation data (rank, total,
#                                        knee1, shin1, knee2, shin2, expected_cells,
#                                        total_droplets_included, low_count_threshold,
#                                        in_empty_plateau, spliced, unspliced); used
#                                        directly by the Quarto barcode report.
#   4. <sample_id>_knee_params.txt     — single CSV line: expected_cells,total_droplets,
#                                        low_count_threshold; read by the Nextflow workflow
#                                        to pass CellBender / decontX parameters as
#                                        channel values (no YAML hand-off needed).
#
# Knee estimation mirrors scprocess exactly:
#   - First knee (knee1/shin1): voting across barcodeRanks lower=[100..1000] runs;
#     the inflection/knee pair with the most votes wins.
#   - Second knee (knee2/shin2): voting across filtered subsets of barcodes below
#     the first inflection; same majority-vote logic.
#   - expected_cells         = rank of shin1 (first inflection)
#   - total_droplets_included = geometric mean of ranks between shin1 and knee2
#   - low_count_threshold    = shin2 UMI value (default)
#
# Usage:
#   Rscript run_barcode_estimation.R <sample_id> <quant_dir> <h5_out>
#                                    <nf_csv_out> <knee_data_out> <params_out>

suppressPackageStartupMessages({
    library(fishpond)
    library(SingleCellExperiment)
    library(DropletUtils)
    library(Matrix)
    library(data.table)
})

# ===========================================================================
# Helper functions — ported from scprocess/scripts/mapping.R
# Must be defined before the main execution block.
# ===========================================================================

# ---------------------------------------------------------------------------
# .calc_ambient_params
#   Orchestrates the four-point knee estimation and returns a per-barcode
#   data.table with all knee/shin values and CellBender/decontX parameters.
#
#   mat              : sparse matrix (genes x barcodes) — only colSums are used
#   run              : sample identifier (written into the sample_id column)
#   min_umis_empty   : minimum UMI count to include in empty-droplet analysis
#   low_count_threshold : how to derive the low-count threshold; 'shin2' or 'knee2'
# ---------------------------------------------------------------------------
.calc_ambient_params <- function(mat, run, min_umis_empty = 5,
                                  low_count_threshold = "shin2") {
    # First knee/inflection — voting across multiple lower parameters
    knee1_ls <- .get_knee_and_shin_1(mat)

    # Second knee/inflection — voting across filtered barcode subsets
    knee2_ls <- .get_knee_and_shin_2(mat, knee1_ls$ranks_dt,
                                      min_umis_empty, knee1_ls$shin1_x)

    # Derive CellBender / decontX parameters
    params_ls <- .get_params_ls(knee1_ls, knee2_ls,
                                 low_count_threshold = low_count_threshold)

    # Build per-barcode result data.table
    bender_ps <- copy(knee1_ls$ranks_dt)[
        , sample_id                := run
    ][, `:=`(
        knee1                   = knee1_ls$sel_knee["knee"],
        shin1                   = knee1_ls$sel_knee["shin"],
        knee2                   = knee2_ls$sel_knee["knee"],
        shin2                   = knee2_ls$sel_knee["shin"],
        total_droplets_included = params_ls$total_included,
        low_count_threshold     = params_ls$lc,
        expected_cells          = params_ls$expected_cells
    )]

    # Label barcodes within the empty droplet plateau
    bender_ps <- .get_empty_plateau(
        knee_df        = bender_ps,
        shin1          = knee1_ls$sel_knee["shin"],
        total_included = params_ls$total_included,
        knee2          = knee2_ls$sel_knee["knee"]
    )

    return(bender_ps)
}


# ---------------------------------------------------------------------------
# .get_knee_and_shin_1
#   Runs barcodeRanks with lower = seq(1000, 100, by = -100) and picks the
#   knee+inflection pair that receives the most votes.
# ---------------------------------------------------------------------------
.get_knee_and_shin_1 <- function(mat, min_umis_cells = NULL, knee1 = NA,
                                  shin1 = NA, knee2 = NA) {
    if (all(sapply(c(knee1, shin1, knee2), function(p) !is.na(p)))) {
        # Custom values provided
        low      <- median(c(knee2, shin1))
        ranks_obj <- barcodeRanks(mat, lower = low)
        sel_knee  <- c(shin = shin1, knee = knee1)

    } else if (!is.null(min_umis_cells)) {
        ranks_obj <- barcodeRanks(mat, lower = min_umis_cells)
        sel_knee  <- c(
            shin = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
            knee = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
        )

    } else {
        # Vote across lower = 1000, 900, ..., 100
        ranks_ls <- lapply(seq(1000, 100, by = -100), function(x) {
            barcodeRanks(mat, lower = x)
        })

        ranks_ls_ks <- sapply(ranks_ls, function(x) {
            paste0(
                as.integer(round(as.numeric(as.character(metadata(x)$inflection)))),
                "_",
                as.integer(round(as.numeric(as.character(metadata(x)$knee))))
            )
        })

        votes_tbl <- table(ranks_ls_ks)
        sel_cut   <- names(votes_tbl)[which.max(votes_tbl)]

        parts    <- strsplit(sel_cut, "_")[[1]]
        sel_knee <- c(shin = as.integer(parts[1]), knee = as.integer(parts[2]))

        ranks_obj <- ranks_ls[[which(ranks_ls_ks == sel_cut)[1]]]
    }

    # Convert barcodeRanks to data.table
    ranks_dt <- as.data.table(as.data.frame(ranks_obj), keep.rownames = TRUE)
    setnames(ranks_dt, "rn", "barcode")
    setorder(ranks_dt, rank)

    # X-coordinate (rank) of shin1 (first inflection)
    shin1_idx <- which.min(abs(ranks_dt$total - sel_knee["shin"]))[1]
    shin1_x   <- ranks_dt[shin1_idx, rank]

    list(ranks_dt = ranks_dt, sel_knee = sel_knee, shin1_x = shin1_x)
}


# ---------------------------------------------------------------------------
# .get_knee_and_shin_2
#   Estimates the second knee/inflection (lower boundary of real cells /
#   upper boundary of empty droplets) using barcodes below the first inflection.
# ---------------------------------------------------------------------------
.get_knee_and_shin_2 <- function(mat, ranks_dt, min_umis_empty, shin1_x,
                                  knee2 = NA, shin2 = NA) {
    if (all(sapply(c(knee2, shin2), function(p) !is.na(p)))) {
        # Custom values — snap to nearest actual UMI counts in the data
        shin2_idx  <- which.min(abs(ranks_dt$total - shin2))[1]
        shin2_corr <- ranks_dt[shin2_idx, total]
        knee2_idx  <- which.min(abs(ranks_dt$total - knee2))[1]
        knee2_corr <- ranks_dt[knee2_idx, total]
        sel_knee   <- c(shin = shin2_corr, knee = knee2_corr)

    } else {
        # Vote across multiple cut-points below shin1
        cuts <- .calc_small_knee_cuts_ls(ranks_dt, min_umis_empty, shin1_x)

        ranks_ls <- sapply(cuts, function(this_cut) {
            ranks_smol <- ranks_dt[rank > this_cut, barcode]
            if (length(ranks_smol) < 3) return(NA_character_)
            tryCatch({
                ro     <- barcodeRanks(mat[, ranks_smol, drop = FALSE],
                                       lower = min_umis_empty)
                shin2v <- as.integer(round(as.numeric(as.character(metadata(ro)$inflection))))
                knee2v <- as.integer(round(as.numeric(as.character(metadata(ro)$knee))))
                paste0(shin2v, "_", knee2v)
            }, error = function(e) NA_character_)
        })

        ranks_ls <- ranks_ls[!is.na(ranks_ls)]
        if (length(ranks_ls) == 0) {
            # Fall back: use the lower bound of the barcode rank range
            last_bc    <- tail(ranks_dt[total > min_umis_empty, total], 1)
            sel_knee   <- c(shin = as.integer(last_bc), knee = as.integer(last_bc))
        } else {
            # Pick inflection with most votes; among tied inflections pick knee
            # closest to median
            shin_vals <- sub("_.*", "", ranks_ls)
            shin_tbl  <- table(shin_vals)
            sel_i     <- names(shin_tbl)[which.max(shin_tbl)]

            match_idx  <- grepl(paste0("^", sel_i, "_"), ranks_ls)
            match_ks   <- as.numeric(sub(".*_", "", ranks_ls[match_idx]))
            med_val    <- median(match_ks)
            sel_k      <- match_ks[which.min(abs(match_ks - med_val))[1]]

            sel_knee   <- c(shin = as.integer(sel_i), knee = as.integer(sel_k))
        }
    }

    # Rank of the second knee
    knee2_x <- unique(ranks_dt[total == sel_knee["knee"], rank])

    list(sel_knee = sel_knee, knee2_x = knee2_x)
}


# ---------------------------------------------------------------------------
# .calc_small_knee_cuts_ls
#   Generates 10 cut-point ranks spanning the log-scale range between shin1
#   and the lower tail of the barcode rank distribution.
# ---------------------------------------------------------------------------
.calc_small_knee_cuts_ls <- function(ranks_dt, min_umis_empty, shin_x) {
    # Last rank that still has enough unique totals for barcodeRanks (needs >= 3)
    last   <- tail(unique(ranks_dt[total > min_umis_empty, total]), n = 3)[1]
    last_x <- ranks_dt[total == last, rank][1]

    # Geometric midpoint between shin1 and the lower tail
    middle <- copy(ranks_dt)[, n := .I][
        rank %between% c(shin_x, last_x), n
    ] |> log10() |> mean() |> (\(x) 10^x)()

    cuts <- 10^seq(log10(shin_x), log10(middle), length.out = 10)
    cuts
}


# ---------------------------------------------------------------------------
# .get_params_ls
#   Derives CellBender / decontX parameters from the two knee estimates.
#
#   expected_cells        = rank at shin1 (first inflection)
#   total_droplets_included = geometric mean rank between shin1 and knee2
#   low_count_threshold   = shin2 UMI value (by default)
# ---------------------------------------------------------------------------
.get_params_ls <- function(knee1_ls, knee2_ls, low_count_threshold = "shin2",
                            expected_cells = NA, total_included = NA) {
    ranks_dt <- knee1_ls$ranks_dt
    shin1_x  <- knee1_ls$shin1_x
    knee2_x  <- knee2_ls$knee2_x

    if (is.na(expected_cells)) {
        expected_cells <- shin1_x
    }

    if (is.na(total_included)) {
        total_included <- copy(ranks_dt)[, n := .I][
            rank %between% c(shin1_x, knee2_x), n
        ] |> log10() |> mean() |> (\(x) round(10^x))()
    }

    if (is.character(low_count_threshold)) {
        if (low_count_threshold == "knee2") {
            lc <- knee2_ls$sel_knee["knee"]
        } else {
            lc <- knee2_ls$sel_knee["shin"]   # default: shin2
        }
    } else {
        lc <- low_count_threshold
    }

    list(expected_cells = expected_cells,
         total_included = total_included,
         lc             = lc)
}


# ---------------------------------------------------------------------------
# .get_empty_plateau
#   Labels barcodes that fall within the "empty droplet plateau" — the flat
#   region between shin1 and total_droplets_included on the rank plot.
# ---------------------------------------------------------------------------
.get_empty_plateau <- function(knee_df, shin1, total_included, knee2) {
    shin1_idx <- which.min(abs(knee_df$total - shin1))[1]
    shin1_x   <- knee_df[shin1_idx, rank]

    # Geometric mean row index between shin1_x and total_included (on log scale)
    empty_start <- copy(knee_df)[, n := .I][
        rank %between% c(shin1_x, total_included), n
    ] |> log10() |> mean() |> (\(x) 10^x)()

    empty_end <- knee_df[total == knee2, unique(rank)]

    knee_df[, in_empty_plateau := fifelse(
        rank %between% c(empty_start, empty_end), TRUE, FALSE
    )]
    knee_df
}

# ===========================================================================
# Main execution
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
    stop("Usage: run_barcode_estimation.R <sample_id> <quant_dir> <h5_out> ",
         "<nf_csv_out> <knee_data_out> <params_out>")
}

sample_id     <- args[1]
quant_dir     <- args[2]
h5_out        <- args[3]
nf_csv_out    <- args[4]
knee_data_out <- args[5]
params_out    <- args[6]

message("=== BARCODE ESTIMATION ===")
message("Sample:    ", sample_id)
message("Quant dir: ", quant_dir)

fry_dir <- file.path(quant_dir, "af_quant")
if (!dir.exists(fry_dir)) {
    stop("af_quant directory not found inside quant_dir: ", fry_dir)
}

# ---------------------------------------------------------------------------
# 1. Load with individual spliced (S), unspliced (U), ambiguous (A) assays
# ---------------------------------------------------------------------------
message("Loading alevin-fry quantification ...")
sce <- loadFry(fry_dir, outputFormat = list(S = c("S"), U = c("U"), A = c("A")))

# Total counts = S + U + A (sum all three assays)
total_counts <- assay(sce, "S") + assay(sce, "U") + assay(sce, "A")

# Drop barcodes with zero total counts
keep <- colSums(total_counts) > 0
sce  <- sce[, keep]
total_counts <- total_counts[, keep]

message("Barcodes (non-zero): ", ncol(sce))
message("Genes:               ", nrow(sce))

# ---------------------------------------------------------------------------
# 2. Write unfiltered total count H5 (input for CellBender / decontX)
# ---------------------------------------------------------------------------
message("Writing unfiltered count H5: ", h5_out)
write10xCounts(h5_out, total_counts, version = "3", overwrite = TRUE)

# ---------------------------------------------------------------------------
# 3. Compute nuclear_fraction per barcode
# ---------------------------------------------------------------------------
# nuclear_fraction = U / (S + U), pseudocount of 1 in numerator + denominator
# to handle barcodes where S+U == 0 (all counts are ambiguous).
message("Computing nuclear_fraction ...")

S <- colSums(assay(sce, "S"))
U <- colSums(assay(sce, "U"))
nuclear_fraction <- (U + 1) / (S + U + 2)

nf_df <- data.frame(
    barcode          = colnames(sce),
    nuclear_fraction = round(nuclear_fraction, 6),
    spliced          = as.integer(S),
    unspliced        = as.integer(U)
)
write.csv(nf_df, nf_csv_out, row.names = FALSE, quote = FALSE)
message("Nuclear fraction written: ", nf_csv_out)
message("  Median nuclear_fraction: ", round(median(nuclear_fraction), 3))

# ---------------------------------------------------------------------------
# 4. Four-point knee estimation (knee1, shin1, knee2, shin2)
#    Ported from scprocess/scripts/mapping.R — voting approach for robustness.
# ---------------------------------------------------------------------------
message("Running 4-point knee estimation ...")

bender_ps <- .calc_ambient_params(total_counts, run = sample_id)

# Merge spliced / unspliced counts into knee_data
splice_dt <- data.table(
    barcode   = colnames(sce),
    spliced   = as.integer(S),
    unspliced = as.integer(U)
)
bender_ps <- merge(bender_ps, splice_dt, by = "barcode")
setorder(bender_ps, rank)

fwrite(bender_ps, file = knee_data_out)
message("Knee data written: ", knee_data_out)

expected_cells <- unique(bender_ps$expected_cells)
total_droplets <- unique(bender_ps$total_droplets_included)
low_count_thr  <- unique(bender_ps$low_count_threshold)

message("  knee1 (upper cell boundary):  ", unique(bender_ps$knee1))
message("  shin1 (first inflection):     ", unique(bender_ps$shin1))
message("  knee2 (lower cell boundary):  ", unique(bender_ps$knee2))
message("  shin2 (empty droplet bottom): ", unique(bender_ps$shin2))
message("  expected_cells:               ", expected_cells)
message("  total_droplets_included:      ", total_droplets)
message("  low_count_threshold:          ", low_count_thr)

writeLines(paste(
    as.integer(round(expected_cells)),
    as.integer(round(total_droplets)),
    as.integer(round(low_count_thr)),
    sep = ","
), params_out)
message("Knee params written: ", params_out)
message("=== DONE ===")
