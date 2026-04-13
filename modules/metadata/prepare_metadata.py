#!/usr/bin/env python3

import argparse
import gzip
import warnings

import pandas as pd


def _read_metadata_csv_compat(metadata_csv):
    try:
        return pd.read_csv(metadata_csv)
    except pd.errors.ParserError as err:
        warnings.warn(
            "Metadata CSV could not be parsed with pandas' default engine; "
            "retrying with the Python engine and skipping malformed rows. "
            f"Original error: {err}",
            RuntimeWarning,
        )
        return pd.read_csv(
            metadata_csv,
            engine='python',
            on_bad_lines='warn',
        )


def prepare_metadata(metadata_csv, metadata_id_col, sample_ids, metadata_vars, out_csv):
    meta = _read_metadata_csv_compat(metadata_csv)
    if metadata_id_col not in meta.columns:
        raise KeyError(f"metadata_id_col '{metadata_id_col}' not found in metadata CSV")

    if metadata_id_col != 'sample_id' and 'sample_id' in meta.columns:
        existing = meta['sample_id'].astype(str)
        normalized = meta[metadata_id_col].astype(str)
        if not existing.equals(normalized):
            raise ValueError(
                "Metadata CSV already has a 'sample_id' column that does not match "
                f"'{metadata_id_col}'. Resolve the conflict before running the pipeline."
            )
        meta = meta.drop(columns=['sample_id'])

    meta = meta.rename(columns={metadata_id_col: 'sample_id'})
    meta['sample_id'] = meta['sample_id'].astype(str)

    if metadata_vars:
        missing_vars = [column for column in metadata_vars if column not in meta.columns]
        if missing_vars:
            raise KeyError(f"metadata_vars not found in metadata CSV: {missing_vars}")

    filtered = meta[meta['sample_id'].isin(sample_ids)].copy()
    duplicate_mask = filtered['sample_id'].duplicated(keep=False)
    if duplicate_mask.any():
        duplicate_samples = sorted(filtered.loc[duplicate_mask, 'sample_id'].unique().tolist())
        raise ValueError(
            "Metadata CSV contains duplicate rows for one or more run samples: "
            f"{duplicate_samples}"
        )

    missing_samples = sorted(set(sample_ids) - set(filtered['sample_id']))
    if missing_samples:
        raise ValueError(
            "Metadata CSV is missing one or more discovered samples: "
            f"{missing_samples}"
        )

    column_order = ['sample_id'] + [column for column in filtered.columns if column != 'sample_id']
    filtered = filtered[column_order].sort_values('sample_id').reset_index(drop=True)

    with gzip.open(out_csv, 'wt') as handle:
        filtered.to_csv(handle, index=False)


def main():
    parser = argparse.ArgumentParser(description='Normalize and validate sample metadata for scQC-flow')
    parser.add_argument('--metadata_csv', required=True, type=str)
    parser.add_argument('--metadata_id_col', required=True, type=str)
    parser.add_argument('--sample_ids', required=True, type=str,
                        help='Comma-separated discovered sample IDs')
    parser.add_argument('--metadata_vars', default='', type=str,
                        help='Optional space-separated metadata vars to validate')
    parser.add_argument('--out_csv', default='sample_metadata.csv.gz', type=str)
    args = parser.parse_args()

    sample_ids = [sample_id for sample_id in args.sample_ids.split(',') if sample_id]
    metadata_vars = args.metadata_vars.split() if args.metadata_vars.strip() else []
    prepare_metadata(args.metadata_csv, args.metadata_id_col, sample_ids, metadata_vars, args.out_csv)


if __name__ == '__main__':
    main()