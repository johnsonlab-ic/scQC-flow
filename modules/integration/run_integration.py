#!/usr/bin/env python3
# run_integration.py
#
# Harmony integration mirroring scprocess's integration pipeline.
#
# Takes the HVG count matrix (hvg_counts.h5) from HVG_SELECTION and a
# user-provided metadata CSV, normalises (CPM10k + log1p), scales, and runs:
#   PCA -> Harmony -> neighbors -> Leiden (multiple resolutions) -> UMAP
#
# Outputs:
#   integration_dt.csv.gz  — UMAP1/2 + PCA dims + cluster assignments per cell

import argparse
import gzip
import re
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import harmonypy as hm
import scanpy as sc
import scanpy.external as sce


# ---------------------------------------------------------------------------
# H5 I/O helper
# ---------------------------------------------------------------------------

def _read_hvg_h5(h5_path):
    """Read the HVG count matrix written by hvg_selection.py."""
    with h5py.File(h5_path, 'r') as f:
        indptr     = f['matrix/indptr'][:]
        indices    = f['matrix/indices'][:]
        data       = f['matrix/data'][:]
        features   = f['matrix/features/name'][:].astype(str)
        barcodes   = f['matrix/barcodes'][:].astype(str)
        shape      = tuple(f['matrix/shape'][:])
        sample_ids = f['matrix/sample_ids'][:].astype(str)
    # stored as genes x cells CSC
    mat_csc = sp.csc_matrix((data, indices, indptr), shape=shape)
    return mat_csc, features, barcodes, sample_ids


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------

def run_integration(hvg_h5, metadata_csv, metadata_id_col, metadata_vars_str,
                    n_dims, theta, leiden_res_str, out_csv):

    print("=== INTEGRATION ===")
    res_ls = leiden_res_str.split()
    metadata_vars = metadata_vars_str.split() if metadata_vars_str else []

    # ------------------------------------------------------------------
    # 1. Load HVG count matrix
    # ------------------------------------------------------------------
    print(f"Loading HVG matrix: {hvg_h5}")
    mat_csc, features, barcodes, sample_ids = _read_hvg_h5(hvg_h5)
    print(f"  {mat_csc.shape[0]} genes x {mat_csc.shape[1]} cells")

    # ------------------------------------------------------------------
    # 2. Build AnnData (cells x genes)
    # ------------------------------------------------------------------
    obs_df = pd.DataFrame({'cell_id': barcodes, 'sample_id': sample_ids})
    obs_df.index = barcodes
    adata = sc.AnnData(
        X   = mat_csc.T.tocsr(),
        obs = obs_df,
        var = pd.DataFrame(index=features)
    )
    adata.obs_names_make_unique()

    # ------------------------------------------------------------------
    # 3. Join metadata CSV on sample_id
    # ------------------------------------------------------------------
    if metadata_csv and metadata_vars:
        print(f"Loading metadata: {metadata_csv}")
        meta = pd.read_csv(metadata_csv)
        if metadata_id_col not in meta.columns:
            raise KeyError(f"metadata_id_col '{metadata_id_col}' not found in {metadata_csv}")
        meta = meta.rename(columns={metadata_id_col: 'sample_id'})
        meta['sample_id'] = meta['sample_id'].astype(str)

        # check all requested vars exist
        missing = [v for v in metadata_vars if v not in meta.columns]
        if missing:
            raise KeyError(f"metadata_vars not found in CSV: {missing}")

        keep_cols = ['sample_id'] + metadata_vars
        meta = meta[keep_cols].drop_duplicates(subset='sample_id')

        adata.obs = adata.obs.merge(meta, on='sample_id', how='left')
        adata.obs.index = barcodes
        print(f"  Metadata vars joined: {metadata_vars}")
    else:
        print("  No metadata vars — running without Harmony batch correction")
        metadata_vars = []

    # ------------------------------------------------------------------
    # 4. Normalise (CPM10k + log1p)
    # ------------------------------------------------------------------
    print("Normalising (CPM10k + log1p) ...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ------------------------------------------------------------------
    # 5. Scale
    # ------------------------------------------------------------------
    print("Scaling ...")
    sc.pp.scale(adata, max_value=10)

    # ------------------------------------------------------------------
    # 6. PCA
    # ------------------------------------------------------------------
    print(f"PCA (n_dims={n_dims}) ...")
    sc.tl.pca(adata, n_comps=n_dims)

    # ------------------------------------------------------------------
    # 7. Harmony (if metadata_vars provided and >1 batch)
    # ------------------------------------------------------------------
    sel_embed = 'X_pca'
    n_batches = adata.obs['sample_id'].nunique()

    if metadata_vars and n_batches > 1:
        # harmonypy calls .describe().loc['unique'] which only exists for
        # object/string columns — cast all key columns to str to avoid KeyError
        for v in metadata_vars:
            adata.obs[v] = adata.obs[v].astype(str)
        harmony_key = metadata_vars if len(metadata_vars) > 1 else metadata_vars[0]
        print(f"Harmony integration (key={harmony_key}, theta={theta}) ...")
        ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, harmony_key, theta=theta)
        z = np.asarray(ho.Z_corr)
        # old harmonypy: Z_corr is (n_pcs, n_cells); new PyTorch harmonypy: (n_cells, n_pcs)
        adata.obsm['X_pca_harmony'] = z.T if z.shape[0] != adata.n_obs else z
        sel_embed = 'X_pca_harmony'
    else:
        print("Skipping Harmony (single batch or no metadata vars)")

    # ------------------------------------------------------------------
    # 8. Neighbors -> Leiden -> UMAP
    # ------------------------------------------------------------------
    print("Finding neighbors ...")
    if np.isnan(adata.obsm[sel_embed]).any():
        raise ValueError("NaN values detected in embedding — check input data")
    sc.pp.neighbors(adata, n_pcs=n_dims, use_rep=sel_embed)

    print(f"Leiden clustering (resolutions: {res_ls}) ...")
    for res in res_ls:
        key = f"RNA_snn_res.{res}"
        sc.tl.leiden(adata, key_added=key, resolution=float(res))

    print("UMAP ...")
    sc.tl.umap(adata, maxiter=750)

    # ------------------------------------------------------------------
    # 9. Build output DataFrame
    # ------------------------------------------------------------------
    print("Building output table ...")
    out_df = adata.obs[['cell_id', 'sample_id']].copy()

    # add cluster columns
    cl_cols = [f"RNA_snn_res.{r}" for r in res_ls]
    for col in cl_cols:
        out_df[col] = adata.obs[col].values
        out_df[col.replace('RNA_snn_res.', 'leiden_')] = adata.obs[col].values

    # add metadata vars
    for v in metadata_vars:
        if v in adata.obs.columns:
            out_df[v] = adata.obs[v].values

    # add UMAP
    umap = adata.obsm['X_umap']
    out_df['UMAP1'] = umap[:, 0]
    out_df['UMAP2'] = umap[:, 1]

    # add PCA / Harmony dims (prefix accordingly)
    pca_arr = adata.obsm[sel_embed]
    prefix  = 'hmny_pca' if sel_embed == 'X_pca_harmony' else 'pca'
    for i in range(pca_arr.shape[1]):
        out_df[f"{prefix}_{i+1:02d}"] = pca_arr[:, i]

    # ------------------------------------------------------------------
    # 10. Save
    # ------------------------------------------------------------------
    with gzip.open(out_csv, 'wt') as fh:
        out_df.to_csv(fh, index=False)
    print(f"Written integration results: {out_csv}  ({len(out_df)} cells)")
    print("=== INTEGRATION done ===")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Harmony integration of HVG counts')
    parser.add_argument('--hvg_h5',           type=str, required=True,
                        help='HVG count matrix H5 from hvg_selection.py')
    parser.add_argument('--metadata_csv',     type=str, default=None,
                        help='Path to metadata CSV')
    parser.add_argument('--metadata_id_col',  type=str, default='sample_id',
                        help='Column in metadata CSV that maps to sample IDs')
    parser.add_argument('--metadata_vars',    type=str, default='',
                        help='Space-separated metadata columns for Harmony (e.g. "brainregion condition")')
    parser.add_argument('--n_dims',           type=int, default=30)
    parser.add_argument('--theta',            type=float, default=2.0)
    parser.add_argument('--leiden_res',       type=str, default='0.3 0.5 1.0',
                        help='Space-separated Leiden resolutions')
    parser.add_argument('--out_csv',          type=str, default='integration_dt.csv.gz')
    args = parser.parse_args()

    run_integration(
        hvg_h5           = args.hvg_h5,
        metadata_csv     = args.metadata_csv,
        metadata_id_col  = args.metadata_id_col,
        metadata_vars_str = args.metadata_vars,
        n_dims           = args.n_dims,
        theta            = args.theta,
        leiden_res_str   = args.leiden_res,
        out_csv          = args.out_csv,
    )
