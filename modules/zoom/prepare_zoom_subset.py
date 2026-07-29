#!/usr/bin/env python3

import argparse
import base64
import glob
import gzip
import json
import sys

import pandas as pd


def _resolve_cluster_column(df, cluster_res):
    candidates = [f"RNA_snn_res.{cluster_res}", f"leiden_{cluster_res}"]
    for column in candidates:
        if column in df.columns:
            return column
    raise KeyError(f"Cluster column not found for resolution '{cluster_res}'")


def _decode_spec(spec_b64):
    return json.loads(base64.b64decode(spec_b64).decode("utf-8"))


def _load_annotation_labels(annotation_labels_csv):
    if annotation_labels_csv == "NO_FILE":
        raise FileNotFoundError("Annotation cell labels are required for annotation_label zooms")
    labels_df = pd.read_csv(annotation_labels_csv)
    required = {"sample_id", "cell_id", "label"}
    missing = required.difference(labels_df.columns)
    if missing:
        raise KeyError(f"Annotation labels file is missing columns: {sorted(missing)}")
    labels_df["sample_id"] = labels_df["sample_id"].astype(str)
    labels_df["cell_id"] = labels_df["cell_id"].astype(str)
    return labels_df


def _load_annotation_method_labels(annotation_method_id, annotation_method_pattern):
    labels_path = f"annotation_cells_{annotation_method_id}.csv.gz"
    matched_paths = sorted(glob.glob(annotation_method_pattern))
    if labels_path not in matched_paths:
        raise FileNotFoundError(
            f"Annotation method labels for method '{annotation_method_id}' were not found; "
            f"expected '{labels_path}' within pattern '{annotation_method_pattern}'"
        )
    labels_df = pd.read_csv(labels_path)
    required = {"sample_id", "cell_id", "label"}
    missing = required.difference(labels_df.columns)
    if missing:
        raise KeyError(
            f"Annotation method labels file '{labels_path}' is missing columns: {sorted(missing)}"
        )
    labels_df["sample_id"] = labels_df["sample_id"].astype(str)
    labels_df["cell_id"] = labels_df["cell_id"].astype(str)
    return labels_df


def _load_qc(qc_pattern):
    qc_files = sorted(glob.glob(qc_pattern))
    if not qc_files:
        raise FileNotFoundError(f"No QC files matched pattern '{qc_pattern}'")
    qc_df = pd.concat((pd.read_csv(path) for path in qc_files), ignore_index=True)
    qc_df["sample_id"] = qc_df["sample_id"].astype(str)
    qc_df["cell_id"] = qc_df["cell_id"].astype(str)
    qc_df["global_cell_id"] = qc_df["sample_id"] + "_" + qc_df["cell_id"]
    return qc_df


def run_prepare_zoom_subset(integration_csv, annotation_labels_csv, annotation_method_pattern, qc_pattern,
                            zoom_spec_b64, out_qc_csv, out_selection_csv):
    spec = _decode_spec(zoom_spec_b64)
    zoom_name = spec["name"]
    zoom_source = spec["source"]
    sample_min_cells = int(spec.get("sample_min_cells", 10))

    print(f"=== PREPARE_ZOOM_SUBSET: {zoom_name} ===")
    print(f"Source: {zoom_source}")
    if zoom_source == "annotation_label":
        zoom_label = str(spec.get("label", zoom_name))
        print(f"Label: {zoom_label}")
    elif zoom_source == "annotation_method_label":
        zoom_label = str(spec.get("label", zoom_name))
        annotation_method_id = str(spec["annotation_method_id"])
        print(f"Annotation method: {annotation_method_id}")
        print(f"Label: {zoom_label}")
    else:
        zoom_values = [str(value) for value in spec["values"]]
        print(f"Values: {zoom_values}")
    print(f"Min cells per sample: {sample_min_cells}")

    int_df = pd.read_csv(integration_csv)
    int_df["sample_id"] = int_df["sample_id"].astype(str)
    int_df["cell_id"] = int_df["cell_id"].astype(str)
    int_df = int_df.loc[int_df["UMAP1"].notna()].copy()

    if zoom_source == "cluster":
        zoom_values = [str(value) for value in spec["values"]]
        cluster_res = str(spec["cluster_res"])
        cluster_col = _resolve_cluster_column(int_df, cluster_res)
        selected_df = int_df.loc[
            int_df[cluster_col].astype(str).isin(zoom_values),
            ["sample_id", "cell_id", cluster_col],
        ].copy()
        selected_df = selected_df.rename(columns={cluster_col: "zoom_value"})
    elif zoom_source == "annotation_label":
        zoom_label = str(spec.get("label", zoom_name))
        labels_df = _load_annotation_labels(annotation_labels_csv)
        labels_df = labels_df.loc[labels_df["label"].astype(str) == zoom_label].copy()
        selected_df = labels_df.loc[:, ["sample_id", "cell_id", "label"]].copy()
        selected_df = selected_df.rename(columns={"label": "zoom_value"})
        selected_df = selected_df.merge(
            int_df.loc[:, ["sample_id", "cell_id"]].drop_duplicates(),
            on=["sample_id", "cell_id"],
            how="inner",
        )
    elif zoom_source == "annotation_method_label":
        zoom_label = str(spec.get("label", zoom_name))
        annotation_method_id = str(spec["annotation_method_id"])
        # Cluster-LEVEL selection: assign each cell its cluster's MAJORITY predicted label at
        # annotation_cluster_res, then take every cell in clusters whose majority == zoom_label.
        # Downstream analyses act on the cluster-level call, so the zoom must too -- filtering the
        # raw per-cell label under-selects (e.g. neurons predicted individually as glia inside a
        # neuron-majority cluster get dropped).
        cluster_res = str(spec.get("annotation_cluster_res", "0.5"))
        cluster_col = _resolve_cluster_column(int_df, cluster_res)
        labels_df = _load_annotation_method_labels(annotation_method_id, annotation_method_pattern)
        merged = labels_df.loc[:, ["sample_id", "cell_id", "label"]].merge(
            int_df.loc[:, ["sample_id", "cell_id", cluster_col]].drop_duplicates(),
            on=["sample_id", "cell_id"],
            how="inner",
        )
        merged[cluster_col] = merged[cluster_col].astype(str)
        merged = merged.loc[merged[cluster_col].str.len() > 0].copy()
        majority = (
            merged.groupby(cluster_col)["label"]
            .agg(lambda s: s.astype(str).value_counts().idxmax())
        )
        selected_clusters = set(majority.index[majority.astype(str) == zoom_label])
        print(f"Cluster resolution: {cluster_res} ({cluster_col})")
        print(f"Clusters with majority label '{zoom_label}': {sorted(selected_clusters)}")
        selected_df = merged.loc[
            merged[cluster_col].isin(selected_clusters), ["sample_id", "cell_id"]
        ].copy()
        selected_df["zoom_value"] = zoom_label
    else:
        raise ValueError(f"Unsupported zoom source '{zoom_source}'")

    if selected_df.empty:
        raise ValueError(f"Zoom '{zoom_name}' selected zero cells before sample filtering")

    selected_df = selected_df.drop_duplicates(subset=["sample_id", "cell_id"]).copy()
    per_sample = (
        selected_df.groupby("sample_id", as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
    )
    keep_samples = set(per_sample.loc[per_sample["n_cells"] >= sample_min_cells, "sample_id"])
    selected_df = selected_df.loc[selected_df["sample_id"].isin(keep_samples)].copy()
    selected_df["zoom_name"] = zoom_name
    selected_df["zoom_source"] = zoom_source

    if selected_df.empty:
        raise ValueError(
            f"Zoom '{zoom_name}' selected zero cells after applying sample_min_cells={sample_min_cells}"
        )

    qc_df = _load_qc(qc_pattern)
    subset_qc_df = qc_df.loc[qc_df["global_cell_id"].isin(selected_df["cell_id"])].copy()
    subset_qc_df = subset_qc_df.drop(columns=["global_cell_id"])

    if subset_qc_df.empty:
        raise ValueError(f"Zoom '{zoom_name}' matched zero QC rows")

    with gzip.open(out_qc_csv, "wt") as handle:
        subset_qc_df.to_csv(handle, index=False)
    with gzip.open(out_selection_csv, "wt") as handle:
        selected_df.to_csv(handle, index=False)

    print(f"Selected cells after filtering: {len(selected_df):,}")
    print(f"Samples retained: {selected_df['sample_id'].nunique()}")
    print(f"QC rows written: {len(subset_qc_df):,}")
    print("=== PREPARE_ZOOM_SUBSET done ===")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare subset QC inputs for zoom workflows")
    parser.add_argument("--integration_csv", required=True)
    parser.add_argument("--annotation_labels_csv", required=True)
    parser.add_argument("--annotation_method_pattern", default="annotation_cells_*.csv.gz")
    parser.add_argument("--qc_pattern", default="qc_metrics_*.csv.gz")
    parser.add_argument("--zoom_spec_b64", required=True)
    parser.add_argument("--out_qc_csv", required=True)
    parser.add_argument("--out_selection_csv", required=True)
    args = parser.parse_args()

    try:
        run_prepare_zoom_subset(
            integration_csv=args.integration_csv,
            annotation_labels_csv=args.annotation_labels_csv,
            annotation_method_pattern=args.annotation_method_pattern,
            qc_pattern=args.qc_pattern,
            zoom_spec_b64=args.zoom_spec_b64,
            out_qc_csv=args.out_qc_csv,
            out_selection_csv=args.out_selection_csv,
        )
    except Exception as err:
        print(f"ERROR: {err}", file=sys.stderr)
        raise