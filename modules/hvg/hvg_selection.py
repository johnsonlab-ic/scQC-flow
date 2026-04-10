#!/usr/bin/env python3
# hvg_selection.py
#
# Per-dataset HVG selection mirroring scprocess's hvgs pipeline.
#
# Takes all per-sample filtered H5 files (stacked S+U+A, from decontX) and
# per-sample QC metrics CSVs, sums S+U+A counts per gene, filters to
# QC-passing singlets, and runs scanpy HVG selection with batch correction.
#
# Outputs:
#   hvg_stats.csv.gz   — per-gene stats (mean, variance, HVG flag, rank)
#   hvg_counts.h5      — raw SUA-summed counts for all cells (HVG genes only)
#                        in 10x-style CSC H5 format

import argparse
import glob
import gzip
import re
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc


# ---------------------------------------------------------------------------
# H5 I/O helpers (mirrors scprocess .get_h5_mx)
# ---------------------------------------------------------------------------

def _read_h5_csc(h5_path):
    """Read a stacked S+U+A H5 into a CSC sparse matrix."""
    with h5py.File(h5_path, 'r') as f:
        indptr   = f['matrix/indptr'][:]
        indices  = f['matrix/indices'][:]
        data     = f['matrix/data'][:]
        features = f['matrix/features/name'][:].astype(str)
        barcodes = f['matrix/barcodes'][:].astype(str)
        shape    = tuple(f['matrix/shape'][:])
    mat = sp.csc_matrix((data, indices, indptr), shape=shape)
    return mat, features, barcodes


def _write_h5_csc(h5_path, mat_csc, features, barcodes):
    """Write a CSC sparse matrix to an H5 file in 10x-style layout."""
    mat_csc = mat_csc.tocsc()
    with h5py.File(h5_path, 'w') as f:
        f.create_dataset('matrix/indptr',          data=mat_csc.indptr,  compression='gzip')
        f.create_dataset('matrix/indices',         data=mat_csc.indices, compression='gzip')
        f.create_dataset('matrix/data',            data=mat_csc.data,    compression='gzip')
        f.create_dataset('matrix/shape',           data=np.array(mat_csc.shape))
        f.create_dataset('matrix/features/name',   data=np.array(features, dtype='S'),  compression='gzip')
        f.create_dataset('matrix/barcodes',        data=np.array(barcodes, dtype='S'),  compression='gzip')


# ---------------------------------------------------------------------------
# Sum S+U+A layers (mirrors scprocess sum_SUA)
# ---------------------------------------------------------------------------

def _sum_sua(mat_csc, features):
    """
    Sum spliced (_S), unspliced (_U), and ambiguous (_A) rows.
    Returns (summed_csc, gene_ids) where gene_ids have the _S/_U/_A suffix removed.
    """
    types = ['_S', '_U', '_A']
    mats  = []
    for t in types:
        idx = np.array([i for i, n in enumerate(features) if n.endswith(t)])
        if len(idx) > 0:
            mats.append(mat_csc[idx, :])

    # gene names from _S rows (same ordering as _U and _A)
    s_idx  = np.array([i for i, n in enumerate(features) if n.endswith('_S')])
    genes  = np.array([re.sub(r'_[SUA]$', '', features[i]) for i in s_idx])

    summed = mats[0].copy()
    for m in mats[1:]:
        summed = summed + m

    return summed, genes


# ---------------------------------------------------------------------------
# Main HVG selection function
# ---------------------------------------------------------------------------

def run_hvg_selection(h5_files, qc_csv_files, n_top_genes, out_stats, out_h5):
    print(f"=== HVG_SELECTION ===")
    print(f"H5 files:    {h5_files}")
    print(f"QC CSVs:     {qc_csv_files}")
    print(f"n_top_genes: {n_top_genes}")

    # ------------------------------------------------------------------
    # 1. Load QC metrics — keep only keep=TRUE & scdbl_class='singlet'
    # ------------------------------------------------------------------
    qc_frames = [pd.read_csv(f) for f in qc_csv_files]
    qc_all    = pd.concat(qc_frames, ignore_index=True)
    qc_pass   = qc_all[(qc_all['keep'] == True) & (qc_all['scdbl_class'] == 'singlet')]
    keep_cells = set(qc_pass['cell_id'].astype(str))
    print(f"QC-passing singlets: {len(keep_cells)}")

    # ------------------------------------------------------------------
    # 2. Load each H5, sum SUA, filter to QC-passing cells
    # ------------------------------------------------------------------
    # Infer sample_id from filename: filt_counts_<sampleId>.h5
    sample_mats  = []
    sample_bcs   = []
    sample_ids   = []
    gene_list    = None

    for h5f in sorted(h5_files):
        sample_id = re.sub(r'^filt_counts_', '', h5f.replace('.h5', ''))
        # strip path if present
        sample_id = re.sub(r'^.*/filt_counts_', '', h5f).replace('.h5', '')

        mat, features, barcodes = _read_h5_csc(h5f)
        summed, genes = _sum_sua(mat, features)

        if gene_list is None:
            gene_list = genes
        else:
            if not np.array_equal(gene_list, genes):
                raise ValueError(f"Gene lists differ between samples: {h5f}")

        # filter barcodes to QC-passing cells
        keep_idx = np.array([i for i, bc in enumerate(barcodes) if bc in keep_cells])
        if len(keep_idx) == 0:
            print(f"  WARNING: no QC-passing cells found in {h5f}, skipping")
            continue

        summed_filt = summed[:, keep_idx]
        bcs_filt    = barcodes[keep_idx]

        sample_mats.append(summed_filt)
        sample_bcs.extend(bcs_filt.tolist())
        sample_ids.extend([sample_id] * len(bcs_filt))

        print(f"  {sample_id}: {len(keep_idx)} cells loaded")

    if len(sample_mats) == 0:
        raise ValueError("No cells remain after QC filtering.")

    # ------------------------------------------------------------------
    # 3. Concatenate samples (genes x cells CSC)
    # ------------------------------------------------------------------
    combined_csc = sp.hstack(sample_mats, format='csc')
    all_barcodes = np.array(sample_bcs)
    all_sample_ids = np.array(sample_ids)
    print(f"Combined matrix: {combined_csc.shape[0]} genes x {combined_csc.shape[1]} cells")

    # ------------------------------------------------------------------
    # 4. Build AnnData (cells x genes CSR) for scanpy
    # ------------------------------------------------------------------
    adata = sc.AnnData(
        X   = combined_csc.T.tocsr(),
        obs = pd.DataFrame({'sample_id': all_sample_ids}, index=all_barcodes),
        var = pd.DataFrame(index=gene_list)
    )
    adata.obs_names_make_unique()

    # ------------------------------------------------------------------
    # 5. Normalise (CPM10k + log1p) — required for seurat_v3 HVG
    # ------------------------------------------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ------------------------------------------------------------------
    # 6. HVG selection (seurat_v3 with batch_key)
    # ------------------------------------------------------------------
    n_batches = adata.obs['sample_id'].nunique()
    batch_key = 'sample_id' if n_batches > 1 else None
    print(f"Running HVG selection: n_top_genes={n_top_genes}, batch_key={batch_key}")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes  = n_top_genes,
        flavor       = 'seurat_v3',
        batch_key    = batch_key,
        span         = 0.3,
    )

    # ------------------------------------------------------------------
    # 7. Save HVG stats CSV
    # ------------------------------------------------------------------
    hvg_df = adata.var.copy().reset_index()
    hvg_df.rename(columns={'index': 'gene_id'}, inplace=True)
    # add rank (1 = most variable)
    hvg_df['hvg_rank'] = hvg_df['highly_variable_rank'].rank(method='first').astype('Int64')
    with gzip.open(out_stats, 'wt') as fh:
        hvg_df.to_csv(fh, index=False)
    print(f"Written HVG stats: {out_stats}  ({hvg_df['highly_variable'].sum()} HVGs)")

    # ------------------------------------------------------------------
    # 8. Save raw (un-normalised) HVG count matrix
    #    Use the original combined_csc, subset to HVG genes, in CSC format
    # ------------------------------------------------------------------
    hvg_mask  = adata.var['highly_variable'].values
    hvg_genes = gene_list[hvg_mask]
    hvg_csc   = combined_csc[hvg_mask, :]   # HVG genes x all cells
    _write_h5_csc(out_h5, hvg_csc, hvg_genes, all_barcodes)
    # also store sample_id per barcode as a separate dataset for integration
    with h5py.File(out_h5, 'a') as f:
        f.create_dataset('matrix/sample_ids',
                         data=np.array(all_sample_ids, dtype='S'),
                         compression='gzip')
    print(f"Written HVG count matrix: {out_h5}  ({hvg_csc.shape[0]} genes x {hvg_csc.shape[1]} cells)")
    print("=== HVG_SELECTION done ===")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='HVG selection from decontX H5 files')
    parser.add_argument('--h5_pattern',   type=str, default='filt_counts_*.h5',
                        help='Glob pattern for input H5 files (default: filt_counts_*.h5)')
    parser.add_argument('--qc_pattern',   type=str, default='qc_metrics_*.csv.gz',
                        help='Glob pattern for QC metrics CSVs (default: qc_metrics_*.csv.gz)')
    parser.add_argument('--n_top_genes',  type=int, default=4000)
    parser.add_argument('--out_stats',    type=str, default='hvg_stats.csv.gz')
    parser.add_argument('--out_h5',       type=str, default='hvg_counts.h5')
    args = parser.parse_args()

    h5_files  = sorted(glob.glob(args.h5_pattern))
    qc_files  = sorted(glob.glob(args.qc_pattern))

    if not h5_files:
        sys.exit(f"ERROR: no H5 files matched pattern '{args.h5_pattern}'")
    if not qc_files:
        sys.exit(f"ERROR: no QC CSV files matched pattern '{args.qc_pattern}'")

    run_hvg_selection(
        h5_files    = h5_files,
        qc_csv_files = qc_files,
        n_top_genes  = args.n_top_genes,
        out_stats    = args.out_stats,
        out_h5       = args.out_h5,
    )
