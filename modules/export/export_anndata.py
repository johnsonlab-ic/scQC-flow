#!/usr/bin/env python3

import argparse
import gzip
import os
import re

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csc_matrix, hstack


def _parse_gtf(gtf_path):
    opener = gzip.open if gtf_path.endswith('.gz') else open
    rows = []
    with opener(gtf_path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attrs = fields[8]
            match_id = re.search(r'gene_id "([^"]+)"', attrs)
            if match_id is None:
                continue
            match_symbol = re.search(r'gene_name "([^"]+)"', attrs)
            match_type = re.search(r'gene_type "([^"]+)"', attrs)
            if match_type is None:
                match_type = re.search(r'gene_biotype "([^"]+)"', attrs)
            gene_id = match_id.group(1)
            symbol = match_symbol.group(1) if match_symbol else gene_id
            gene_type = match_type.group(1) if match_type else 'unknown'
            rows.append({
                'gene_key': f'{symbol}_{gene_id}',
                'gene_id': gene_id,
                'symbol': symbol,
                'gene_type': gene_type,
            })
    return pd.DataFrame(rows).drop_duplicates('gene_key')


def _read_h5_csc(h5_path):
    with h5py.File(h5_path, 'r') as handle:
        indptr = handle['matrix/indptr'][:]
        indices = handle['matrix/indices'][:]
        data = handle['matrix/data'][:]
        features = handle['matrix/features/name'][:].astype(str)
        barcodes = handle['matrix/barcodes'][:].astype(str)
        shape = tuple(handle['matrix/shape'][:])
    matrix = csc_matrix((data, indices, indptr), shape=shape)
    return matrix, features, barcodes


def _sum_sua(matrix, features):
    row_names = np.asarray(features, dtype=str)
    suffixes = ['S', 'U', 'A']
    mats = []
    for suffix in suffixes:
        idx = np.where(np.char.endswith(row_names, f'_{suffix}'))[0]
        if idx.size == 0:
            raise ValueError('Expected stacked S/U/A rows in filtered H5 matrix')
        mats.append(matrix[idx, :])
    summed = mats[0] + mats[1] + mats[2]
    gene_keys = np.char.replace(row_names[np.where(np.char.endswith(row_names, '_S'))[0]], '_S', '')
    return summed.tocsc(), gene_keys


def _infer_sample_id(h5_path):
    sample_id = os.path.basename(h5_path)
    sample_id = re.sub(r'^filt_counts_', '', sample_id)
    sample_id = re.sub(r'\.h5$', '', sample_id)
    return sample_id


def _combine_duplicate_columns(df, suffix):
    dup_cols = [col for col in df.columns if col.endswith(suffix)]
    for dup_col in dup_cols:
        base_col = dup_col[:-len(suffix)]
        if base_col in df.columns:
            df[base_col] = df[base_col].where(~df[base_col].isna(), df[dup_col])
            df = df.drop(columns=dup_col)
        else:
            df = df.rename(columns={dup_col: base_col})
    return df


def _load_qc_metadata(qc_pattern):
    qc_paths = sorted([path for path in os.listdir('.') if re.fullmatch(qc_pattern.replace('*', '.*'), path)])
    if not qc_paths:
        raise FileNotFoundError(f'No QC files matched pattern: {qc_pattern}')
    qc_frames = [pd.read_csv(path) for path in qc_paths]
    qc_df = pd.concat(qc_frames, ignore_index=True)
    qc_df['sample_id'] = qc_df['sample_id'].astype(str)
    qc_df['cell_id'] = qc_df['cell_id'].astype(str)
    qc_df['cell_id'] = qc_df['sample_id'] + '_' + qc_df['cell_id']
    qc_df['n_counts'] = qc_df['sum']
    qc_df['n_features'] = qc_df['detected']
    qc_df['splice_fraction'] = (qc_df['total_spliced'] + 1.0) / (
        qc_df['total_spliced'] + qc_df['total_unspliced'] + 2.0
    )
    qc_df['mito_fraction'] = qc_df['mito_pct']
    return qc_df


def _load_annotation(annotation_csv):
    if annotation_csv == 'NO_FILE':
        return None
    ann_df = pd.read_csv(annotation_csv)
    rename_map = {
        'cluster': 'annotation_cluster',
        'label': 'annotation_label',
        'label_score': 'annotation_label_score',
        'n_markers': 'annotation_n_markers',
    }
    ann_df = ann_df.rename(columns=rename_map)
    ann_cols = [
        'cell_id',
        'annotation_cluster',
        'annotation_label',
        'annotation_label_score',
        'annotation_n_markers',
    ]
    return ann_df[ann_cols]


def _load_cell_metadata(integration_csv, qc_pattern, annotation_csv):
    int_df = pd.read_csv(integration_csv)
    int_df['cell_id'] = int_df['cell_id'].astype(str)
    int_df = int_df.drop_duplicates('cell_id').reset_index(drop=True)

    qc_df = _load_qc_metadata(qc_pattern)
    merged = int_df.merge(qc_df, on='cell_id', how='left', suffixes=('', '__qc'))
    merged = _combine_duplicate_columns(merged, '__qc')

    ann_df = _load_annotation(annotation_csv)
    if ann_df is not None:
        merged = merged.merge(ann_df, on='cell_id', how='left', suffixes=('', '__ann'))
        merged = _combine_duplicate_columns(merged, '__ann')

    if 'sample_id' not in merged.columns:
        raise KeyError('Merged export metadata is missing sample_id')
    if merged['sample_id'].isna().any():
        missing = merged.loc[merged['sample_id'].isna(), 'cell_id'].tolist()[:5]
        raise ValueError(f'sample_id missing for exported cells: {missing}')

    merged['sample_id'] = merged['sample_id'].astype(str)
    return merged


def _build_count_matrix(h5_pattern, obs_df):
    h5_paths = sorted([path for path in os.listdir('.') if re.fullmatch(h5_pattern.replace('*', '.*'), path)])
    if not h5_paths:
        raise FileNotFoundError(f'No filtered H5 files matched pattern: {h5_pattern}')

    h5_map = {_infer_sample_id(path): path for path in h5_paths}
    matrices = []
    ordered_barcodes = []
    gene_keys = None

    for sample_id, sample_obs in obs_df.groupby('sample_id', sort=False):
        if sample_id not in h5_map:
            raise FileNotFoundError(f'Missing filtered H5 for sample {sample_id}')

        matrix, features, barcodes = _read_h5_csc(h5_map[sample_id])
        summed, this_gene_keys = _sum_sua(matrix, features)
        if gene_keys is None:
            gene_keys = np.asarray(this_gene_keys, dtype=str)
        elif not np.array_equal(gene_keys, np.asarray(this_gene_keys, dtype=str)):
            raise ValueError('Gene keys differ across filtered H5 files')

        global_barcodes = np.asarray([f'{sample_id}_{barcode}' for barcode in barcodes], dtype=str)
        index_map = {barcode: idx for idx, barcode in enumerate(global_barcodes)}
        wanted = sample_obs['cell_id'].astype(str).tolist()
        missing = [barcode for barcode in wanted if barcode not in index_map]
        if missing:
            raise ValueError(f'Filtered H5 for sample {sample_id} is missing cells: {missing[:5]}')
        sel_idx = [index_map[barcode] for barcode in wanted]
        matrices.append(summed[:, sel_idx])
        ordered_barcodes.extend(wanted)

    if not matrices:
        raise ValueError('No cells selected for export')

    combined = hstack(matrices, format='csc')
    return combined, pd.Index(gene_keys, name='gene_key'), ordered_barcodes


def _build_var_df(gtf_df, gene_keys):
    var_df = pd.DataFrame({'gene_key': gene_keys.astype(str)})
    var_df = var_df.merge(gtf_df, on='gene_key', how='left')
    missing = var_df['gene_id'].isna()
    if missing.any():
        fallback = var_df.loc[missing, 'gene_key'].str.extract(r'^(?P<symbol>.*)_(?P<gene_id>ENSG[^_]+|ENSMUSG[^_]+)$')
        var_df.loc[missing, 'symbol'] = fallback['symbol'].fillna(var_df.loc[missing, 'gene_key'])
        var_df.loc[missing, 'gene_id'] = fallback['gene_id'].fillna(var_df.loc[missing, 'gene_key'])
        var_df.loc[missing, 'gene_type'] = 'unknown'
    var_df = var_df.set_index('gene_key')
    return var_df


def _sanitize_obs(obs_df):
    obs_df = obs_df.copy()
    for col in obs_df.columns:
        if pd.api.types.is_object_dtype(obs_df[col]):
            obs_df[col] = obs_df[col].astype('string')
    obs_df.index = obs_df['cell_id'].astype(str)
    return obs_df


def _embedding_matrix(obs_df, columns):
    if not columns:
        return None
    return obs_df[columns].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)


def _add_embeddings(adata, obs_df, include_doublet_umap):
    if include_doublet_umap and {'dbl_UMAP1', 'dbl_UMAP2'}.issubset(obs_df.columns):
        adata.obsm['X_umap'] = _embedding_matrix(obs_df, ['dbl_UMAP1', 'dbl_UMAP2'])
    elif {'UMAP1', 'UMAP2'}.issubset(obs_df.columns):
        adata.obsm['X_umap'] = _embedding_matrix(obs_df, ['UMAP1', 'UMAP2'])

    hmny_cols = sorted(
        [col for col in obs_df.columns if re.fullmatch(r'hmny_pca_\d+', col)],
        key=lambda value: int(value.rsplit('_', 1)[1]),
    )
    pca_cols = sorted(
        [col for col in obs_df.columns if re.fullmatch(r'pca_\d+', col)],
        key=lambda value: int(value.rsplit('_', 1)[1]),
    )
    if hmny_cols:
        adata.obsm['X_harmony'] = _embedding_matrix(obs_df, hmny_cols)
    if pca_cols:
        adata.obsm['X_pca'] = _embedding_matrix(obs_df, pca_cols)


def _make_adata(counts_mat, obs_df, var_df, include_doublet_umap):
    obs_df = _sanitize_obs(obs_df)
    adata = ad.AnnData(X=counts_mat.T.tocsr(), obs=obs_df, var=var_df.copy())
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    _add_embeddings(adata, obs_df, include_doublet_umap)
    return adata


def main():
    parser = argparse.ArgumentParser(description='Export scQC-flow results as AnnData objects')
    parser.add_argument('--h5_pattern', required=True)
    parser.add_argument('--qc_pattern', required=True)
    parser.add_argument('--integration_csv', required=True)
    parser.add_argument('--annotation_csv', required=True)
    parser.add_argument('--genome_gtf', required=True)
    parser.add_argument('--out_all', required=True)
    parser.add_argument('--out_clean', required=True)
    args = parser.parse_args()

    obs_df = _load_cell_metadata(args.integration_csv, args.qc_pattern, args.annotation_csv)
    gtf_df = _parse_gtf(args.genome_gtf)

    all_obs = obs_df.copy()
    all_counts, gene_keys, ordered_all = _build_count_matrix(args.h5_pattern, all_obs)
    all_obs = all_obs.set_index('cell_id').loc[ordered_all].reset_index()
    var_df = _build_var_df(gtf_df, gene_keys)
    all_adata = _make_adata(all_counts, all_obs, var_df, include_doublet_umap=True)
    all_adata.write_h5ad(args.out_all, compression='gzip')

    clean_mask = (~all_obs['is_dbl'].fillna(False).astype(bool)) & (~all_obs['in_dbl_cl'].fillna(False).astype(bool))
    clean_obs = all_obs.loc[clean_mask].copy().reset_index(drop=True)
    clean_counts, _, ordered_clean = _build_count_matrix(args.h5_pattern, clean_obs)
    clean_obs = clean_obs.set_index('cell_id').loc[ordered_clean].reset_index()
    clean_adata = _make_adata(clean_counts, clean_obs, var_df, include_doublet_umap=False)
    clean_adata.write_h5ad(args.out_clean, compression='gzip')


if __name__ == '__main__':
    main()