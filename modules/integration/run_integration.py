#!/usr/bin/env python3
# run_integration.py
#
# Two-pass Harmony integration — identical to scprocess's integration.py.
#
# Pass 1 (doublet detection):
#   - Load singlet + doublet HVG matrices, hstack
#   - Normalise using QC library sizes (sum column), CPM10k + log1p on sparse
#   - scale(max_value=10) -> PCA -> Harmony(theta=0, key=sample_id)
#   - Leiden at high resolution (dbl_res) -> UMAP
#   - Compute doublet proportion per cluster, flag enriched clusters
#
# Pass 2 (clean integration):
#   - Remove doublets AND cells in doublet-enriched clusters
#   - Re-use normalised matrix (slice to clean cells)
#   - scale -> PCA -> Harmony(theta=user_theta, key=batch_var) -> Leiden(res_ls) -> UMAP
#   - Rename clusters to cl01, cl02 by descending size
#
# Output:
#   integration_dt.csv.gz — full join of pass-2 results + pass-1 doublet data

import argparse
from contextlib import contextmanager
from datetime import datetime
import gc
import glob
import gzip
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import csc_matrix, hstack
import anndata as ad
import scanpy as sc
import harmonypy


def _log(message):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}", flush=True)


def _format_seconds(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, rem = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m {rem:.1f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h {int(minutes)}m {rem:.1f}s"


@contextmanager
def _timed_step(message):
    _log(message)
    start = time.perf_counter()
    try:
        yield
    except Exception:
        elapsed = time.perf_counter() - start
        _log(f"{message} failed after {_format_seconds(elapsed)}")
        raise
    else:
        elapsed = time.perf_counter() - start
        _log(f"{message} completed in {_format_seconds(elapsed)}")


def _describe_sparse_matrix(name, matrix):
    total = matrix.shape[0] * matrix.shape[1]
    density = (matrix.nnz / total) if total else 0.0
    _log(
        f"{name}: {matrix.shape[0]} genes x {matrix.shape[1]} cells; "
        f"nnz={matrix.nnz:,}; density={density:.4%}"
    )


def _run_leiden(adata, resolution):
    key_added = f"RNA_snn_res.{resolution}"
    try:
        sc.tl.leiden(
            adata,
            key_added=key_added,
            resolution=float(resolution),
            flavor='igraph',
            directed=False,
            n_iterations=2,
        )
        _log(f"  Leiden resolution={resolution} used flavor=igraph")
    except Exception as err:
        _log(
            "  WARNING: leiden flavor=igraph failed "
            f"({type(err).__name__}: {err}); falling back to the default backend"
        )
        sc.tl.leiden(
            adata,
            key_added=key_added,
            resolution=float(resolution),
        )
        _log(f"  Leiden resolution={resolution} used fallback backend")


# ---------------------------------------------------------------------------
# Harmony compatibility wrapper
# ---------------------------------------------------------------------------

def _run_harmony_integrate_compat(adata, key, basis='X_pca', adjusted_basis='X_pca_harmony', **kwargs):
    """Run harmonypy directly and store corrected PCs with the right orientation.

    scprocess' integration flow was written against older harmonypy releases where
    ``Z_corr`` was returned as PCs x cells. In the landmark image, harmonypy 0.2.x
    returns cells x PCs instead. Scanpy 1.12.1 still unconditionally transposes,
    which breaks pass 2 of the two-pass integration. Handle both layouts here so
    the pipeline remains compatible with the scprocess logic and with newer images.
    """
    if isinstance(key, str):
        key_cols = [key]
    else:
        key_cols = list(key)

    missing_cols = [col for col in key_cols if col not in adata.obs.columns]
    if missing_cols:
        raise KeyError(f"Harmony batch columns missing from AnnData: {missing_cols}")

    for col in key_cols:
        if adata.obs[col].isna().any():
            raise ValueError(f"Harmony batch column '{col}' contains missing values")

    x = np.asarray(adata.obsm[basis], dtype=np.float64)
    harmony_out = harmonypy.run_harmony(x, adata.obs, key, **kwargs)
    z_corr = np.asarray(harmony_out.Z_corr)

    expected_shape = x.shape
    if z_corr.shape == expected_shape:
        corrected = z_corr
    elif z_corr.shape == (expected_shape[1], expected_shape[0]):
        corrected = z_corr.T
    else:
        raise ValueError(
            "Harmony returned corrected PCs with unexpected shape "
            f"{z_corr.shape}; expected {expected_shape} or {(expected_shape[1], expected_shape[0])}"
        )

    adata.obsm[adjusted_basis] = corrected


# ---------------------------------------------------------------------------
# H5 I/O
# ---------------------------------------------------------------------------

def _read_h5(h5_path):
    """Read HVG count matrix (genes x cells CSC) from H5."""
    with h5py.File(h5_path, 'r') as f:
        indptr   = f['matrix/indptr'][:]
        indices  = f['matrix/indices'][:]
        data     = f['matrix/data'][:]
        features = f['matrix/features/name'][:].astype(str)
        barcodes = f['matrix/barcodes'][:].astype(str)
        shape    = tuple(f['matrix/shape'][:])
    mat_csc = csc_matrix((data, indices, indptr), shape=shape)
    return mat_csc, features, barcodes


# ---------------------------------------------------------------------------
# Load and combine matrices
# ---------------------------------------------------------------------------

def _get_hvg_mat(hvg_mat_f, dbl_hvg_mat_f):
    """Load singlet + doublet HVG matrices, hstack into one genes x cells CSC."""
    mat_singlet, features, bcs_passed = _read_h5(hvg_mat_f)

    bcs_dbl = []
    if dbl_hvg_mat_f:
        mat_dbl, _, bcs_dbl_arr = _read_h5(dbl_hvg_mat_f)
        bcs_dbl = bcs_dbl_arr.tolist()
        if mat_dbl.shape[1] > 0:
            all_hvg_mat = hstack([mat_singlet, mat_dbl])
        else:
            all_hvg_mat = mat_singlet
    else:
        all_hvg_mat = mat_singlet

    return all_hvg_mat, bcs_passed.tolist(), bcs_dbl


# ---------------------------------------------------------------------------
# Build cells_df from QC metrics
# ---------------------------------------------------------------------------

def _get_cells_df(qc_csv_files, bcs_passed, bcs_dbl, metadata_vars=None):
    """Build cells DataFrame matching matrix column order.

    Marks doublet cells with is_dbl_int column.
    Joins metadata if provided.
    """
    # Load and concat QC CSVs
    qc_frames = [pd.read_csv(f) for f in qc_csv_files]
    all_coldata = pd.concat(qc_frames, ignore_index=True)
    all_coldata['sample_id'] = all_coldata['sample_id'].astype(str)
    all_coldata['cell_id'] = all_coldata['cell_id'].astype(str)
    # Prefix cell_id with sample_id to match barcodes written by hvg_selection.py
    all_coldata['cell_id'] = all_coldata['sample_id'] + '_' + all_coldata['cell_id']

    # Identify passed (keep=True) and doublet cells
    passed_idx = all_coldata['keep'].astype(bool)

    # Doublet identification: scdbl_class == 'doublet'
    dbl_idx = (all_coldata['scdbl_class'] == 'doublet')

    # Validate barcodes
    passed_cells = set(all_coldata.loc[passed_idx, 'cell_id'])
    if not set(bcs_passed).issubset(passed_cells):
        raise ValueError("QC-passed barcodes from HVG matrix don't match QC cell_ids")

    if bcs_dbl:
        dbl_cells = set(all_coldata.loc[dbl_idx, 'cell_id'])
        if not set(bcs_dbl).issubset(dbl_cells):
            raise ValueError("Doublet barcodes from HVG matrix don't match QC cell_ids")

    # Add doublet label
    all_coldata['is_dbl_int'] = dbl_idx

    # Order cells_df to match matrix column order (barcodes are now globally unique)
    all_bcs = bcs_passed + bcs_dbl
    order_df = pd.DataFrame({'cell_id': all_bcs, '_col_idx': range(len(all_bcs))})
    cells_df = (order_df
                .merge(all_coldata, on='cell_id', how='left')
                .sort_values('_col_idx')
                .drop(columns='_col_idx')
                .reset_index(drop=True))

    # Validate no missing QC entries
    n_missing = cells_df['sum'].isna().sum()
    if n_missing > 0:
        missing_bcs = cells_df.loc[cells_df['sum'].isna(), 'cell_id'].tolist()[:5]
        raise ValueError(
            f"QC data missing for {n_missing} barcodes from HVG matrix "
            f"(first few: {missing_bcs})"
        )

    if metadata_vars:
        missing = [v for v in metadata_vars if v not in cells_df.columns]
        if missing:
            raise KeyError(f"metadata_vars not found in QC metadata columns: {missing}")

        missing_meta = cells_df[metadata_vars].isna().any(axis=1)
        if missing_meta.any():
            missing_samples = sorted(cells_df.loc[missing_meta, 'sample_id'].unique().tolist())
            raise ValueError(
                "Metadata values are missing for one or more run samples in QC metrics. "
                f"Missing sample_id values: {missing_samples}"
            )

    return cells_df


# ---------------------------------------------------------------------------
# Normalisation (identical to scprocess)
# ---------------------------------------------------------------------------

def _normalize_hvg_mat(hvg_mat, cells_df, exclude_mito, scale_f=10000):
    """Normalise genes x cells matrix using QC library sizes.

    Divides each cell's counts by its total UMI count (from QC 'sum' column),
    multiplies by scale_f, then log1p — all in-place on sparse data.
    """
    if exclude_mito:
        lib_sizes = (cells_df['sum'] - cells_df['mito_sum']).values.astype(np.float64)
    else:
        lib_sizes = cells_df['sum'].values.astype(np.float64)

    # Convert to CSR (rows = genes, cols = cells) for efficient per-cell ops
    hvg_mat = hvg_mat.tocsr()
    # In CSR with shape (genes, cells): indices are column (cell) indices
    hvg_mat.data = hvg_mat.data.astype(np.float64)
    hvg_mat.data /= lib_sizes[hvg_mat.indices]
    hvg_mat.data *= scale_f
    np.log1p(hvg_mat.data, out=hvg_mat.data)

    return hvg_mat


# ---------------------------------------------------------------------------
# Single integration pass
# ---------------------------------------------------------------------------

def _do_one_integration(adata, batch_var, n_dims, res_ls, theta):
    """Run one integration pass: scale -> PCA -> Harmony -> Leiden -> UMAP.

    Returns a DataFrame with cell_id, embedding coords, UMAP, clusters.
    """
    # Check whether we have >1 batch
    n_batches = adata.obs[batch_var].nunique()
    if n_batches == 1:
        this_embedding = 'pca'
    else:
        this_embedding = 'harmony'
    _log(
        f"Integration pass setup: batch_var='{batch_var}', levels={n_batches}, "
        f"n_dims={n_dims}, theta={theta}, resolutions={res_ls}"
    )

    # Scale
    with _timed_step('  Scaling expression matrix'):
        sc.pp.scale(adata, max_value=10)

    # PCA
    with _timed_step('  Running PCA'):
        sc.tl.pca(adata, n_comps=n_dims)

    sel_embed = 'X_pca'
    if this_embedding == 'harmony' and float(theta) != 0.0:
        with _timed_step(f"  Running Harmony on '{batch_var}'"):
            _run_harmony_integrate_compat(adata, key=batch_var, theta=theta)
        sel_embed = 'X_pca_harmony'
    elif this_embedding == 'harmony' and float(theta) == 0.0:
        # theta=0 means no diversity penalty → no correction; PCA is equivalent
        _log('  Harmony theta=0: skipping correction, using PCA')
        this_embedding = 'pca'
    else:
        _log(f"  Single-level batch variable '{batch_var}': using PCA only")

    # Neighbors
    _log(f"  Using embedding '{sel_embed}' for neighbor graph")
    if np.isnan(adata.obsm[sel_embed]).any():
        raise ValueError("NaN values in embedding — check input data")
    with _timed_step('  Building neighbor graph'):
        sc.pp.neighbors(adata, n_pcs=n_dims, use_rep=sel_embed)

    # Leiden clustering
    _log('  Finding clusters')
    if not isinstance(res_ls, list):
        res_ls = [res_ls]
    for res in res_ls:
        with _timed_step(f'    Leiden clustering at resolution={res}'):
            _run_leiden(adata, res)

    # UMAP
    with _timed_step('  Running UMAP'):
        sc.tl.umap(adata, maxiter=750)

    # Extract results
    with _timed_step('  Recording cluster assignments'):
        clusts_df = _get_clusts_from_adata(adata, this_embedding, batch_var)

    with _timed_step('  Extracting embeddings'):
        embeds_df = _get_embeddings_from_adata(adata, this_embedding, sel_embed)

    int_df = clusts_df.merge(embeds_df, on='cell_id', how='inner')
    _log(f"  Integration pass yielded {len(int_df)} rows")

    return int_df


# ---------------------------------------------------------------------------
# Extract clusters from AnnData (with cl01/cl02 renaming)
# ---------------------------------------------------------------------------

def _get_clusts_from_adata(adata, embedding, batch_var):
    """Extract cluster assignments, rename to cl01/cl02 by descending size."""
    clusts_df = adata.obs.copy()

    # Get clustering columns
    cl_vs = [col for col in clusts_df.columns if re.match(r'RNA_snn_res\..*', col)]
    sample_vs = list(set([batch_var, 'sample_id']) & set(clusts_df.columns))
    all_cols = ['cell_id'] + cl_vs + ['embedding'] if 'embedding' in clusts_df.columns else ['cell_id'] + cl_vs
    all_cols = list(set(['cell_id'] + cl_vs + sample_vs))

    # Add embedding column
    clusts_df['embedding'] = embedding
    all_cols.append('embedding')

    clusts_df = clusts_df[all_cols].copy()

    # Rename clusters: cl01, cl02, ... ordered by descending count
    for cl_v in cl_vs:
        counts = clusts_df[cl_v].value_counts().sort_values(ascending=False)
        rank_map = {}
        for rank, cluster_val in enumerate(counts.index, start=1):
            rank_map[cluster_val] = f"cl{rank:02d}"
        clusts_df[cl_v] = clusts_df[cl_v].map(rank_map)

    return clusts_df


# ---------------------------------------------------------------------------
# Extract embeddings from AnnData
# ---------------------------------------------------------------------------

def _get_embeddings_from_adata(adata, embedding, sel_embed):
    """Extract PCA/Harmony dims and UMAP coordinates."""
    pca_array = adata.obsm[sel_embed]
    n_dims = pca_array.shape[1]

    prefix = 'hmny_pca' if embedding == 'harmony' else 'pca'
    pca_col_names = [
        re.sub(r'_(\d)$', r'_0\1', f'{prefix}_{i}')
        for i in range(1, n_dims + 1)
    ]

    embeds_df = pd.DataFrame(pca_array, columns=pca_col_names)
    embeds_df['cell_id'] = adata.obs['cell_id'].values

    umap_array = adata.obsm['X_umap']
    embeds_df['UMAP1'] = umap_array[:, 0]
    embeds_df['UMAP2'] = umap_array[:, 1]

    return embeds_df


# ---------------------------------------------------------------------------
# Doublet data calculation
# ---------------------------------------------------------------------------

def _calc_dbl_data(int_dbl, cells_df, dbl_res, dbl_cl_prop):
    """Compute doublet proportion per cluster, flag enriched clusters."""
    dbl_ids = set(cells_df.loc[cells_df['is_dbl_int'] == True, 'cell_id'].tolist())

    dbl_clust_col = f"RNA_snn_res.{dbl_res}"

    dbl_data = int_dbl[['cell_id', 'UMAP1', 'UMAP2', dbl_clust_col]].copy()
    dbl_data = dbl_data.rename(columns={
        'UMAP1': 'dbl_UMAP1',
        'UMAP2': 'dbl_UMAP2',
        dbl_clust_col: 'dbl_cluster'
    })
    dbl_data['is_dbl'] = dbl_data['cell_id'].isin(dbl_ids)
    dbl_data['dbl_cluster'] = dbl_data['dbl_cluster'].astype(str)

    # Proportion of doublets per cluster
    cluster_counts = dbl_data.groupby('dbl_cluster')['cell_id'].count()
    cluster_dbl_counts = dbl_data[dbl_data['is_dbl']].groupby('dbl_cluster')['cell_id'].count()
    dbl_prop_map = (cluster_dbl_counts / cluster_counts).fillna(0).to_dict()
    dbl_data['dbl_prop'] = dbl_data['dbl_cluster'].map(dbl_prop_map).astype(float)

    # Flag enriched clusters
    dbl_data['in_dbl_cl'] = dbl_data['dbl_prop'] > dbl_cl_prop

    return dbl_data


# ---------------------------------------------------------------------------
# Filter out doublets
# ---------------------------------------------------------------------------

def _adata_filter_out_doublets(all_hvg_mat, cells_df, dbl_data):
    """Remove doublets and cells in doublet-enriched clusters. Return AnnData."""
    ok_ids = set(
        dbl_data.loc[
            (dbl_data['is_dbl'] == False) & (dbl_data['in_dbl_cl'] == False),
            'cell_id'
        ].tolist()
    )

    # Build AnnData from normalised matrix, subset to ok cells
    adata = ad.AnnData(X=all_hvg_mat.T, obs=cells_df.reset_index(drop=True))
    keep_idx = adata.obs['cell_id'].isin(ok_ids).to_numpy()
    adata = adata[keep_idx, :].copy()
    adata.obs_names_make_unique()

    return adata


# ---------------------------------------------------------------------------
# Main integration function
# ---------------------------------------------------------------------------

def run_integration(hvg_h5, dbl_hvg_h5, qc_csv_files, metadata_vars, exclude_mito,
                    n_dims, dbl_res, dbl_cl_prop, theta, res_ls,
                    out_csv):
    """Two-pass integration identical to scprocess."""

    _log('=== INTEGRATION (two-pass) ===')
    _log(f"Input singlet HVG H5: {hvg_h5}")
    _log(f"Input doublet HVG H5: {dbl_hvg_h5 if dbl_hvg_h5 else 'none'}")
    _log(f"QC files discovered: {len(qc_csv_files)}")
    _log(f"Metadata vars: {metadata_vars if metadata_vars else ['sample_id']}")
    _log(f"Exclude mito from normalization: {exclude_mito}")

    # ------------------------------------------------------------------
    # 1. Load HVG matrices (singlet + doublet)
    # ------------------------------------------------------------------
    with _timed_step('Loading HVG matrices'):
        all_hvg_mat, bcs_passed, bcs_dbl = _get_hvg_mat(hvg_h5, dbl_hvg_h5)
    n_singlets = len(bcs_passed)
    n_doublets = len(bcs_dbl)
    _log(f'  Singlets: {n_singlets:,}, Doublets: {n_doublets:,}')
    _describe_sparse_matrix('  Combined HVG matrix', all_hvg_mat)

    # ------------------------------------------------------------------
    # 2. Build cells_df from QC metrics + metadata
    # ------------------------------------------------------------------
    with _timed_step('Loading cell metadata'):
        cells_df = _get_cells_df(
            qc_csv_files, bcs_passed, bcs_dbl,
            metadata_vars=metadata_vars
        )
    _log(f'  cells_df: {len(cells_df):,} cells')

    # ------------------------------------------------------------------
    # 3. Normalise (library-size CPM10k + log1p on sparse)
    # ------------------------------------------------------------------
    with _timed_step('Normalising HVG matrix'):
        all_hvg_mat = _normalize_hvg_mat(all_hvg_mat, cells_df, exclude_mito)

    # ------------------------------------------------------------------
    # 4. Pass 1: Integration with doublets (theta=0, dbl_res)
    # ------------------------------------------------------------------
    _log('Pass 1: Integration with doublets (theta=0)')
    adata_dbl = ad.AnnData(X=all_hvg_mat.T, obs=cells_df.reset_index(drop=True))
    adata_dbl.obs_names_make_unique()
    _log(f'  AnnData (with doublets): {adata_dbl.shape[0]:,} cells x {adata_dbl.shape[1]:,} genes')

    # Pass 1 uses sample_id as batch var and theta=0
    int_dbl = _do_one_integration(
        adata_dbl,
        batch_var='sample_id',
        n_dims=n_dims,
        res_ls=[str(dbl_res)],
        theta=0
    )

    del adata_dbl
    gc.collect()

    # ------------------------------------------------------------------
    # 5. Calculate doublet data (proportion per cluster)
    # ------------------------------------------------------------------
    with _timed_step('Calculating doublet cluster enrichment'):
        dbl_data = _calc_dbl_data(int_dbl, cells_df, dbl_res, dbl_cl_prop)

    n_in_dbl_cl = dbl_data['in_dbl_cl'].sum()
    n_is_dbl = dbl_data['is_dbl'].sum()
    n_removed = ((dbl_data['is_dbl']) | (dbl_data['in_dbl_cl'])).sum()
    _log(f'  Doublets: {n_is_dbl:,}')
    _log(f'  Cells in doublet-enriched clusters: {n_in_dbl_cl:,}')
    _log(f'  Total cells to remove before pass 2: {n_removed:,}')

    del int_dbl
    gc.collect()

    # ------------------------------------------------------------------
    # 6. Pass 2: Clean integration (remove doublets + enriched clusters)
    # ------------------------------------------------------------------
    _log('Pass 2: Clean integration (doublets removed)')
    with _timed_step('Filtering out doublets and enriched clusters'):
        adata = _adata_filter_out_doublets(all_hvg_mat, cells_df, dbl_data)
    _log(f'  AnnData (clean): {adata.shape[0]:,} cells x {adata.shape[1]:,} genes')

    del all_hvg_mat
    gc.collect()

    # Determine batch_var for pass 2.
    # scprocess uses a single batch key for clean-data integration, so when
    # multiple metadata vars are requested we collapse them into one combined
    # batch label for Harmony.
    if metadata_vars and len(metadata_vars) > 0:
        if len(metadata_vars) == 1:
            batch_var_pass2 = metadata_vars[0]
            adata.obs[batch_var_pass2] = adata.obs[batch_var_pass2].astype(str)
        else:
            batch_var_pass2 = '_harmony_batch'
            adata.obs[batch_var_pass2] = (
                adata.obs[metadata_vars].astype(str).agg('__'.join, axis=1)
            )
    else:
        batch_var_pass2 = 'sample_id'
        adata.obs[batch_var_pass2] = adata.obs[batch_var_pass2].astype(str)

    n_batch_levels_pass2 = adata.obs[batch_var_pass2].nunique()
    _log(
        f"  Pass 2 batch variable: '{batch_var_pass2}' with "
        f"{n_batch_levels_pass2} level(s)"
    )

    int_ok = _do_one_integration(
        adata,
        batch_var=batch_var_pass2,
        n_dims=n_dims,
        res_ls=res_ls,
        theta=theta
    )

    del adata
    gc.collect()

    # ------------------------------------------------------------------
    # 7. Join pass-2 results with pass-1 doublet data
    # ------------------------------------------------------------------
    with _timed_step('Joining pass 1 and pass 2 results'):
        int_df = int_ok.merge(dbl_data, on='cell_id', how='outer')

    # Carry sample-level metadata into the final integration table so the
    # report can infer available covariates directly from the data.
    meta_cols = ['sample_id'] + [col for col in metadata_vars if col != 'sample_id']
    cell_meta_df = cells_df[['cell_id'] + meta_cols].drop_duplicates('cell_id')
    int_df = int_df.merge(cell_meta_df, on='cell_id', how='left', suffixes=('', '__meta'))
    for col in meta_cols:
        meta_col = f'{col}__meta'
        if meta_col not in int_df.columns:
            continue
        if col in int_df.columns:
            int_df[col] = int_df[col].fillna(int_df[meta_col])
            int_df = int_df.drop(columns=meta_col)
        else:
            int_df = int_df.rename(columns={meta_col: col})

    # ------------------------------------------------------------------
    # 8. Save
    # ------------------------------------------------------------------
    _log(f'Saving integration table to {out_csv}')
    with _timed_step('Writing integration output'):
        with gzip.open(out_csv, 'wt') as fh:
            int_df.to_csv(fh, index=False)
    _log(f'  Written {len(int_df):,} rows')
    _log('=== INTEGRATION done ===')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Two-pass Harmony integration (scprocess-identical)')
    parser.add_argument('--hvg_h5',          type=str, required=True,
                        help='Singlet HVG count matrix H5')
    parser.add_argument('--dbl_hvg_h5',      type=str, default=None,
                        help='Doublet HVG count matrix H5')
    parser.add_argument('--qc_pattern',      type=str, default='qc_metrics_*.csv.gz',
                        help='Glob pattern for QC metrics CSVs')
    parser.add_argument('--metadata_vars',   type=str, default='',
                        help='Space-separated metadata columns already present in QC metrics for Harmony pass 2')
    parser.add_argument('--exclude_mito',    action='store_true',
                        help='Exclude mito counts from library size denominator')
    parser.add_argument('--n_dims',          type=int, default=30,
                        help='Number of PCA dimensions')
    parser.add_argument('--dbl_res',         type=float, default=2.0,
                        help='Leiden resolution for doublet detection pass')
    parser.add_argument('--dbl_cl_prop',     type=float, default=0.5,
                        help='Doublet proportion threshold for cluster exclusion')
    parser.add_argument('--theta',           type=float, default=2.0,
                        help='Harmony theta for clean integration pass')
    parser.add_argument('--leiden_res',      type=str, default='0.2 0.5 1.0',
                        help='Space-separated Leiden resolutions for clean pass')
    parser.add_argument('--out_csv',         type=str, default='integration_dt.csv.gz')
    args = parser.parse_args()

    # Parse list arguments
    metadata_vars = args.metadata_vars.split() if args.metadata_vars.strip() else []
    res_ls = args.leiden_res.split()

    # Glob QC files
    qc_files = sorted(glob.glob(args.qc_pattern))
    if not qc_files:
        sys.exit(f"ERROR: no QC files matched pattern '{args.qc_pattern}'")

    run_integration(
        hvg_h5          = args.hvg_h5,
        dbl_hvg_h5      = args.dbl_hvg_h5,
        qc_csv_files    = qc_files,
        metadata_vars   = metadata_vars,
        exclude_mito    = args.exclude_mito,
        n_dims          = args.n_dims,
        dbl_res         = args.dbl_res,
        dbl_cl_prop     = args.dbl_cl_prop,
        theta           = args.theta,
        res_ls          = res_ls,
        out_csv         = args.out_csv,
    )
