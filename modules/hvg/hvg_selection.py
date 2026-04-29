#!/usr/bin/env python3
# hvg_selection.py
#
# Seurat VST HVG selection — matches scprocess hvgs.py logic.
#
# Two-pass design to avoid loading all samples simultaneously:
#
#   Pass 1 (per sample, one at a time):
#     Load H5 → sum S+U+A → filter to QC singlets → compute VST stats → free matrix.
#     Matches scprocess calculate_std_var_stats_for_sample pattern.
#
#   Aggregate: concat per-sample stats → rank HVGs (multi-batch Seurat VST).
#
#   Pass 2 (per sample, one at a time):
#     Load H5 → sum S+U+A → extract HVG rows only → append to output lists → free.
#     hstack of HVG-only slices is safe: ~4k genes x sparse cells << full gene matrix.
#     Matches scprocess create_hvg_matrix / create_doublets_matrix pattern.
#
# Peak memory: O(n_genes x max_cells_per_sample), independent of number of samples.
# Output: hvg_stats.csv.gz, hvg_counts.h5 (singlets), dbl_hvg_counts.h5 (doublets)

import argparse
import glob
import gzip
import re
import sys
import warnings

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import csc_matrix, csr_matrix, hstack
from skmisc.loess import loess


# ---------------------------------------------------------------------------
# GTF parser — produces SYMBOL_ENSEMBL lookup
# ---------------------------------------------------------------------------

def _parse_gtf(gtf_path):
    """Parse gene entries from GTF. Returns dict: ensembl_id -> symbol_ensembl."""
    opener = gzip.open if gtf_path.endswith('.gz') else open
    sym_map = {}
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attrs = fields[8]
            m_id  = re.search(r'gene_id "([^"]+)"', attrs)
            m_sym = re.search(r'gene_name "([^"]+)"', attrs)
            if m_id is None:
                continue
            ens_id = m_id.group(1)
            symbol = m_sym.group(1) if m_sym else ''
            sym_ensembl = f"{symbol}_{ens_id}" if symbol else f"{ens_id}_{ens_id}"
            sym_map[ens_id] = sym_ensembl
    print(f"GTF parsed: {len(sym_map)} genes")
    return sym_map


# ---------------------------------------------------------------------------
# H5 I/O helpers
# ---------------------------------------------------------------------------

def _read_h5_csc(h5_path):
    """Read a stacked S+U+A H5 into a CSC sparse matrix (genes x cells)."""
    with h5py.File(h5_path, 'r') as f:
        indptr   = f['matrix/indptr'][:]
        indices  = f['matrix/indices'][:]
        data     = f['matrix/data'][:]
        features = f['matrix/features/name'][:].astype(str)
        barcodes = f['matrix/barcodes'][:].astype(str)
        shape    = tuple(f['matrix/shape'][:])
    mat = csc_matrix((data, indices, indptr), shape=shape)
    return mat, features, barcodes


def _write_h5_csc(h5_path, mat_csc, features, barcodes):
    """Write a CSC sparse matrix to H5 (genes x cells)."""
    mat_csc = mat_csc.tocsc()
    with h5py.File(h5_path, 'w') as f:
        f.create_dataset('matrix/indptr',        data=mat_csc.indptr,  compression='gzip')
        f.create_dataset('matrix/indices',       data=mat_csc.indices, compression='gzip')
        f.create_dataset('matrix/data',          data=mat_csc.data,    compression='gzip')
        f.create_dataset('matrix/shape',         data=np.array(mat_csc.shape))
        f.create_dataset('matrix/features/name', data=np.array(features, dtype='S'), compression='gzip')
        f.create_dataset('matrix/barcodes',      data=np.array(barcodes, dtype='S'), compression='gzip')


# ---------------------------------------------------------------------------
# Sum S+U+A layers (identical to scprocess sum_SUA)
# ---------------------------------------------------------------------------

def _sum_sua(mat_csc, features):
    """Sum spliced (_S), unspliced (_U), and ambiguous (_A) rows."""
    row_names = np.array(features, dtype=str)
    types = ['_S$', '_U$', '_A$']
    mats = []
    for t in types:
        indices = [i for i, name in enumerate(row_names) if re.search(t, name)]
        mats.append(mat_csc[indices, :])

    # gene names from stripping suffix
    genes = [re.sub(r'_[SUA]$', '', name) for name in row_names]
    uniq_genes = list(dict.fromkeys(genes))
    uniq_genes = np.array(uniq_genes)

    mats_sum = sum(mats)
    return mats_sum, uniq_genes


# ---------------------------------------------------------------------------
# Loess fitting (identical to scprocess safe_loess)
# ---------------------------------------------------------------------------

def _safe_loess(x, y, span, initial_amount=1e-16, max_attempts=5, seed=1234):
    """Fit loess with increasing jitter on failure (matches scprocess)."""
    attempts = 0
    amount = initial_amount

    while attempts < max_attempts:
        try:
            if attempts == 0:
                x_jitter = x
            else:
                np.random.seed(seed)
                jitter = np.random.uniform(-amount, amount, len(x))
                x_jitter = x + jitter

            model = loess(x_jitter, y, span=span, degree=2)
            model.fit()
            return model

        except Exception as e:
            warnings.warn(f"Attempt {attempts + 1} failed: {e}. Increasing jitter to {amount * 10}.")
            amount *= 10

        attempts += 1

    raise RuntimeError(f"Failed to fit loess model after {max_attempts} attempts.")


# ---------------------------------------------------------------------------
# Seurat VST variance calculation (identical to scprocess)
# ---------------------------------------------------------------------------

def _calculate_feature_stats(sparse_csr, features):
    """Calculate mean and variance per gene (row) from genes×cells CSR matrix."""
    n_cells = sparse_csr.shape[1]
    # Use scanpy's efficient implementation if available, else manual
    mean = np.array(sparse_csr.mean(axis=1)).flatten()
    # Var = E[X^2] - E[X]^2
    sq_mean = np.array(sparse_csr.multiply(sparse_csr).mean(axis=1)).flatten()
    var = sq_mean - mean**2
    # Bessel correction
    var = var * n_cells / (n_cells - 1) if n_cells > 1 else var

    return mean, var


def _calculate_estimated_vars(mean, var, n_cells, span=0.3):
    """
    Fit loess to log10(mean) vs log10(var) to get estimated variance.
    Returns (reg_std, clip_val, est_var) — matching scprocess logic.
    """
    not_const = var > 0
    y = np.log10(var[not_const])
    x = np.log10(mean[not_const])

    # Fit loess (matching scprocess: unique x,y pairs)
    xy_df = pd.DataFrame({'x': x, 'y': y}).drop_duplicates()
    model = _safe_loess(xy_df['x'].values, xy_df['y'].values, span=span)

    # Predict estimated variance for all genes
    est_var = np.zeros(len(mean), dtype=np.float64)
    est_var[not_const] = model.predict(x).values

    # Regularized std and clip value (Seurat approach)
    reg_std = np.sqrt(10**est_var)
    vmax = np.sqrt(n_cells)
    clip_val = reg_std * vmax + mean

    return reg_std, clip_val, est_var


def _calculate_regularized_variance(sparse_csr, reg_std, clip_val, mean, n_cells):
    """
    Calculate standardized variance with clipping (identical to scprocess).
    Uses np.bincount for vectorized accumulation — matches scprocess numba.prange logic.
    """
    n_genes = sparse_csr.shape[0]
    csc = sparse_csr.tocsc()
    indices = csc.indices
    data = csc.data.astype(np.float64)

    clipped = np.minimum(data, clip_val[indices])
    sq_counts_sum = np.bincount(indices, weights=clipped ** 2, minlength=n_genes)
    counts_sum    = np.bincount(indices, weights=clipped,      minlength=n_genes)

    variances_norm = (
        (n_cells * mean**2 + sq_counts_sum - 2 * counts_sum * mean) /
        ((n_cells - 1) * reg_std**2)
    )
    return variances_norm


# ---------------------------------------------------------------------------
# HVG ranking with ambient exclusion (identical to scprocess)
# ---------------------------------------------------------------------------

def _rank_hvgs(gene_ids, variances_norm, ambient_genes, n_hvgs):
    """
    Rank genes by standardized variance, excluding ambient genes.
    Returns DataFrame with gene_id, variances_norm, highly_variable, rank.
    """
    df = pd.DataFrame({
        'gene_id': gene_ids,
        'variances_norm': variances_norm,
    })

    # Ambient genes get NaN rank (excluded from selection)
    df['variance_tmp'] = df['variances_norm']
    df.loc[df['gene_id'].isin(ambient_genes), 'variance_tmp'] = np.nan

    # Rank by variance (descending, so rank 1 = highest variance)
    df['highly_variable_rank'] = df['variance_tmp'].rank(ascending=False, method='average')

    # Verify ambient genes are excluded
    amb_ranked = df[df['gene_id'].isin(ambient_genes) & df['highly_variable_rank'].notna()]
    if len(amb_ranked) > 0:
        raise ValueError("Some ambient genes received a rank — logic error")

    # Sort by variance descending
    df = df.sort_values('variances_norm', ascending=False)

    # Top N (by rank) are HVGs
    df['highly_variable'] = df['highly_variable_rank'] <= n_hvgs

    # Add batch count (1 for single-batch; multi-batch handled separately)
    df['highly_variable_nbatches'] = 1

    return df


def _rank_hvgs_multi_batch(per_sample_stats, ambient_genes, n_hvgs):
    """
    Multi-sample HVG ranking (identical to scprocess _process_multiple_groups).
    Ranks per sample, then aggregates by median rank and nbatches.
    """
    # Exclude ambient genes from ranking
    per_sample_stats['variance_tmp'] = per_sample_stats['variances_norm']
    per_sample_stats.loc[
        per_sample_stats['gene_id'].isin(ambient_genes), 'variance_tmp'
    ] = np.nan

    # Rank within each sample
    per_sample_stats['rank_within'] = (
        per_sample_stats.groupby('sample_id')['variance_tmp']
        .rank(ascending=False, method='average')
    )

    # Aggregate across samples
    agg = per_sample_stats.groupby('gene_id').agg(
        highly_variable_nbatches=('rank_within', lambda x: (x <= n_hvgs).sum()),
        highly_variable_rank=('rank_within', 'median'),
        variances_norm=('variances_norm', 'mean'),
    ).reset_index()

    # Sort: most batches first, then by median rank
    agg = agg.sort_values(
        ['highly_variable_nbatches', 'highly_variable_rank'],
        ascending=[False, True]
    )

    # Top N are HVGs
    agg['highly_variable'] = np.arange(len(agg)) < n_hvgs
    if agg['highly_variable'].sum() != n_hvgs:
        raise ValueError(f"Expected {n_hvgs} HVGs but got {agg['highly_variable'].sum()}")

    return agg


# ---------------------------------------------------------------------------
# Main HVG selection function
# ---------------------------------------------------------------------------

def run_hvg_selection(h5_files, qc_csv_files, n_top_genes, out_stats, out_h5,
                      out_dbl_h5, gtf_path, edger_csv=None):
    print("=== HVG_SELECTION (Seurat VST) ===")
    print(f"H5 files:    {h5_files}")
    print(f"QC CSVs:     {qc_csv_files}")
    print(f"n_top_genes: {n_top_genes}")

    # ------------------------------------------------------------------
    # 0. Load auxiliary data
    # ------------------------------------------------------------------
    sym_map = _parse_gtf(gtf_path) if gtf_path else None

    ambient_genes = set()
    if edger_csv:
        print(f"Loading ambient genes from: {edger_csv}")
        edger_df = pd.read_csv(edger_csv)
        ambient_genes = set(edger_df[edger_df['is_ambient'] == True]['gene_id'].tolist())
        print(f"  Ambient genes to exclude: {len(ambient_genes)}")

    # ------------------------------------------------------------------
    # 1. Load QC metrics — build per-sample singlet/doublet cell sets
    # ------------------------------------------------------------------
    qc_frames = [pd.read_csv(f) for f in qc_csv_files]
    qc_all = pd.concat(qc_frames, ignore_index=True)
    qc_all['cell_id'] = qc_all['cell_id'].astype(str)
    qc_all['sample_id'] = qc_all['sample_id'].astype(str)

    qc_by_sample = {str(sid): grp for sid, grp in qc_all.groupby('sample_id')}

    print(f"QC-passing singlets: {int(((qc_all['keep'] == True) & (qc_all['scdbl_class'] != 'doublet')).sum())}")
    print(f"Doublets:            {int((qc_all['scdbl_class'] == 'doublet').sum())}")

    # ------------------------------------------------------------------
    # Pass 1: per-sample VST stats — one sample at a time, no combined matrix.
    # Matches scprocess calculate_std_var_stats_for_sample pattern.
    # ------------------------------------------------------------------
    print("\n--- Pass 1: per-sample VST stats ---")
    h5_files_sorted = sorted(h5_files)
    per_sample_stats = []
    gene_list = None
    sym_ensembl_ids = None

    for h5f in h5_files_sorted:
        sample_id = re.sub(r'^(?:.*/)?filt_counts_', '', h5f).replace('.h5', '')

        mat, features, barcodes = _read_h5_csc(h5f)
        summed, genes = _sum_sua(mat, features)
        del mat

        if gene_list is None:
            gene_list = genes
            sym_ensembl_ids = (
                np.array([sym_map.get(g, f"{g}_{g}") for g in gene_list])
                if sym_map is not None else gene_list
            )
        elif not np.array_equal(gene_list, genes):
            raise ValueError(f"Gene lists differ between samples: {h5f}")

        qc_s = qc_by_sample.get(sample_id)
        if qc_s is None:
            raise ValueError(f"No QC CSV for sample '{sample_id}'. "
                             f"Available: {list(qc_by_sample.keys())}")

        singlet_cells_s = set(qc_s.loc[
            (qc_s['keep'] == True) & (qc_s['scdbl_class'] != 'doublet'), 'cell_id'
        ])
        sing_idx = np.array([i for i, bc in enumerate(barcodes) if bc in singlet_cells_s])

        if len(sing_idx) == 0:
            print(f"  {sample_id}: no singlets, skipping")
            del summed
            continue

        # Column-slice to singlets, free the full sample matrix immediately
        sing_mat = summed[:, sing_idx]
        del summed

        csr = sing_mat.tocsr()
        del sing_mat
        n_cells = csr.shape[1]

        mean, var = _calculate_feature_stats(csr, sym_ensembl_ids)
        reg_std, clip_val, _ = _calculate_estimated_vars(mean, var, n_cells, span=0.3)
        variances_norm = _calculate_regularized_variance(csr, reg_std, clip_val, mean, n_cells)
        del csr

        per_sample_stats.append(pd.DataFrame({
            'gene_id':        sym_ensembl_ids,
            'sample_id':      sample_id,
            'mean':           mean,
            'variance':       var,
            'variances_norm': variances_norm,
            'n_cells':        n_cells,
        }))
        print(f"  {sample_id}: {n_cells} singlets — VST stats done")

    if not per_sample_stats:
        raise ValueError("No singlet cells remain after QC filtering.")

    # ------------------------------------------------------------------
    # Aggregate stats and rank HVGs
    # ------------------------------------------------------------------
    print("\n--- Ranking HVGs ---")
    all_stats = pd.concat(per_sample_stats, ignore_index=True)
    del per_sample_stats
    n_samples = all_stats['sample_id'].nunique()
    print(f"Aggregating stats across {n_samples} samples...")

    if n_samples == 1:
        hvg_df = _rank_hvgs(
            sym_ensembl_ids, all_stats['variances_norm'].values, ambient_genes, n_top_genes
        )
    else:
        hvg_df = _rank_hvgs_multi_batch(all_stats, ambient_genes, n_top_genes)
    del all_stats

    # ------------------------------------------------------------------
    # Save HVG stats CSV
    # ------------------------------------------------------------------
    hvg_out = hvg_df[['gene_id', 'variances_norm', 'highly_variable',
                       'highly_variable_nbatches', 'highly_variable_rank']].copy()
    hvg_out = hvg_out.sort_values(
        ['highly_variable_nbatches', 'highly_variable_rank'],
        ascending=[False, True], na_position='last'
    )
    with gzip.open(out_stats, 'wt') as fh:
        hvg_out.to_csv(fh, index=False)
    print(f"Written HVG stats: {out_stats}  ({int(hvg_out['highly_variable'].sum())} HVGs)")

    # ------------------------------------------------------------------
    # Pass 2: build HVG-only count matrices — one sample at a time.
    # HVG rows only (4k << 60k genes), so the final hstack is small.
    # Matches scprocess create_hvg_matrix / create_doublets_matrix pattern.
    # ------------------------------------------------------------------
    print("\n--- Pass 2: building HVG count matrices ---")
    hvg_gene_ids = set(hvg_df.loc[hvg_df['highly_variable'] == True, 'gene_id'].tolist())
    hvg_mask = np.array([g in hvg_gene_ids for g in sym_ensembl_ids])
    hvg_genes = sym_ensembl_ids[hvg_mask]

    singlet_hvg_mats = []
    all_bcs_singlet = []
    all_sample_ids_singlet = []
    doublet_hvg_mats = []
    all_bcs_doublet = []

    for h5f in h5_files_sorted:
        sample_id = re.sub(r'^(?:.*/)?filt_counts_', '', h5f).replace('.h5', '')

        mat, features, barcodes = _read_h5_csc(h5f)
        summed, _ = _sum_sua(mat, features)
        del mat

        qc_s = qc_by_sample[sample_id]
        singlet_cells_s = set(qc_s.loc[
            (qc_s['keep'] == True) & (qc_s['scdbl_class'] != 'doublet'), 'cell_id'
        ])
        doublet_cells_s = set(qc_s.loc[qc_s['scdbl_class'] == 'doublet', 'cell_id'])

        sing_idx = np.array([i for i, bc in enumerate(barcodes) if bc in singlet_cells_s])
        dbl_idx  = np.array([i for i, bc in enumerate(barcodes) if bc in doublet_cells_s])

        # Convert to CSR once for efficient row-slicing to HVGs, then column-slice per group
        summed_csr = summed.tocsr()
        del summed
        hvg_summed = summed_csr[hvg_mask, :]   # (n_hvgs x all_cells) — small
        del summed_csr

        if len(sing_idx) > 0:
            singlet_hvg_mats.append(hvg_summed[:, sing_idx])
            all_bcs_singlet.extend([f"{sample_id}_{bc}" for bc in barcodes[sing_idx]])
            all_sample_ids_singlet.extend([sample_id] * len(sing_idx))
            print(f"  {sample_id}: {len(sing_idx)} singlets")

        if len(dbl_idx) > 0:
            doublet_hvg_mats.append(hvg_summed[:, dbl_idx])
            all_bcs_doublet.extend([f"{sample_id}_{bc}" for bc in barcodes[dbl_idx]])
            print(f"  {sample_id}: {len(dbl_idx)} doublets")

        del hvg_summed

    if not singlet_hvg_mats:
        raise ValueError("No singlet cells after pass 2.")

    # hstack HVG-only slices — safe: n_hvgs x total_cells is sparse and small
    combined_singlet = hstack(singlet_hvg_mats, format='csc')
    all_bcs_singlet_arr = np.array(all_bcs_singlet)
    all_sample_ids_arr  = np.array(all_sample_ids_singlet)

    _write_h5_csc(out_h5, combined_singlet, hvg_genes, all_bcs_singlet_arr)
    with h5py.File(out_h5, 'a') as f:
        f.create_dataset('matrix/sample_ids',
                         data=np.array(all_sample_ids_arr, dtype='S'),
                         compression='gzip')
    print(f"Written singlet HVG matrix: {out_h5}  "
          f"({combined_singlet.shape[0]} genes x {combined_singlet.shape[1]} cells)")

    if doublet_hvg_mats:
        combined_doublet = hstack(doublet_hvg_mats, format='csc')
        _write_h5_csc(out_dbl_h5, combined_doublet, hvg_genes, np.array(all_bcs_doublet))
        print(f"Written doublet HVG matrix: {out_dbl_h5}  "
              f"({combined_doublet.shape[0]} genes x {combined_doublet.shape[1]} cells)")
    else:
        empty_mat = csc_matrix((len(hvg_genes), 0), dtype=np.float64)
        _write_h5_csc(out_dbl_h5, empty_mat, hvg_genes, np.array([], dtype='U'))
        print(f"Written empty doublet HVG matrix: {out_dbl_h5} (no doublets)")

    print("=== HVG_SELECTION done ===")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Seurat VST HVG selection (scprocess-identical)')
    parser.add_argument('--h5_pattern',   type=str, default='filt_counts_*.h5')
    parser.add_argument('--qc_pattern',   type=str, default='qc_metrics_*.csv.gz')
    parser.add_argument('--n_top_genes',  type=int, default=4000)
    parser.add_argument('--out_stats',    type=str, default='hvg_stats.csv.gz')
    parser.add_argument('--out_h5',       type=str, default='hvg_counts.h5')
    parser.add_argument('--out_dbl_h5',   type=str, default='dbl_hvg_counts.h5')
    parser.add_argument('--gtf',          type=str, default=None)
    parser.add_argument('--edger_csv',    type=str, default=None,
                        help='edger_dt.csv.gz from AMBIENT_DE (for ambient gene exclusion)')
    args = parser.parse_args()

    h5_files  = sorted(glob.glob(args.h5_pattern))
    qc_files  = sorted(glob.glob(args.qc_pattern))

    if not h5_files:
        sys.exit(f"ERROR: no H5 files matched pattern '{args.h5_pattern}'")
    if not qc_files:
        sys.exit(f"ERROR: no QC CSV files matched pattern '{args.qc_pattern}'")

    run_hvg_selection(
        h5_files     = h5_files,
        qc_csv_files = qc_files,
        n_top_genes  = args.n_top_genes,
        out_stats    = args.out_stats,
        out_h5       = args.out_h5,
        out_dbl_h5   = args.out_dbl_h5,
        gtf_path     = args.gtf,
        edger_csv    = args.edger_csv,
    )
