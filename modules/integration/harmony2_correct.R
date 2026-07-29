#!/usr/bin/env Rscript
# Harmony2 (R harmony >=2.0) batch correction for one embedding, called as a subprocess from
# run_integration.py (--harmony_impl harmony2). Pure binary/CSV round-trip, no rpy2:
#   in : embedding.bin (float64, row-major, n_cells x n_dims), dims.txt ("n_cells,n_dims"),
#        meta.csv.gz (batch columns), args: vars (comma-sep), theta, ncores
#   out: corrected.bin (float64, row-major, n_cells x n_dims)
suppressWarnings(suppressMessages({
  .libPaths(c(Sys.getenv("HARMONY2_LIB", "/mnt/data/harmony2_Rlib"), .libPaths()))
  library(harmony)
  library(data.table)
}))

a <- commandArgs(trailingOnly = TRUE)
emb_bin <- a[1]; dims_txt <- a[2]; meta_csv <- a[3]; vars_str <- a[4]
theta   <- suppressWarnings(as.numeric(a[5])); ncores <- as.integer(a[6]); out_bin <- a[7]

dims   <- as.integer(strsplit(readLines(dims_txt), ",")[[1]])
n_cells <- dims[1]; n_dims <- dims[2]
x <- readBin(emb_bin, "double", n = n_cells * n_dims, size = 8)
emb <- matrix(x, nrow = n_cells, ncol = n_dims, byrow = TRUE)   # cells x dims

meta <- as.data.frame(fread(meta_csv))
vars <- strsplit(vars_str, ",")[[1]]
vars <- vars[nzchar(vars)]
stopifnot(nrow(meta) == n_cells, all(vars %in% names(meta)))
for (v in vars) meta[[v]] <- as.factor(as.character(meta[[v]]))

cat(sprintf("Harmony2 %s: %d cells x %d dims, vars=[%s], theta=%s, ncores=%d\n",
            as.character(packageVersion("harmony")), n_cells, n_dims,
            paste(vars, collapse = ","), ifelse(is.na(theta), "default", theta), ncores))

harmony_args <- list(data_mat = emb, meta_data = meta, vars_use = vars,
                     ncores = ncores, return_object = FALSE, verbose = TRUE)
if (!is.na(theta)) harmony_args$theta <- rep(theta, length(vars))   # Harmony2 needs theta per variable

corrected <- do.call(harmony::RunHarmony, harmony_args)
corrected <- as.matrix(corrected)
if (nrow(corrected) == n_dims && ncol(corrected) == n_cells) corrected <- t(corrected)  # dims x cells -> cells x dims
stopifnot(nrow(corrected) == n_cells, ncol(corrected) == n_dims)

writeBin(as.double(as.vector(t(corrected))), out_bin, size = 8)   # row-major to match python
cat("HARMONY2_DONE\n")
