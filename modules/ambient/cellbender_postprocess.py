#!/usr/bin/env python3
import argparse
import csv
import gzip
import h5py


def read_barcodes_from_h5(h5_path):
    with h5py.File(h5_path, "r") as h5:
        barcodes = h5["matrix/barcodes"][:]
    return [b.decode("utf-8") if isinstance(b, bytes) else str(b) for b in barcodes]


def read_knee_rows(knee_csv):
    with open(knee_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"knee CSV is empty: {knee_csv}")
    return rows


def as_int(value, default=0):
    if value is None or value == "":
        return default
    return int(float(value))


def write_lines(path, values):
    with open(path, "w", newline="") as fh:
        for value in values:
            fh.write(f"{value}\n")


def main():
    parser = argparse.ArgumentParser(description="Build CellBender barcode comparison outputs")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--knee_csv", required=True)
    parser.add_argument("--filtered_h5", required=True)
    parser.add_argument("--barcodes_out", required=True)
    parser.add_argument("--summary_out", required=True)
    parser.add_argument("--labels_out", required=True)
    args = parser.parse_args()

    knee_rows = read_knee_rows(args.knee_csv)
    row0 = knee_rows[0]

    expected_cells = as_int(row0.get("expected_cells"), 0)
    total_droplets_included = as_int(row0.get("total_droplets_included"), 0)
    low_count_threshold = as_int(row0.get("low_count_threshold"), 0)

    ranked_rows = sorted(knee_rows, key=lambda row: float(row.get("rank", "inf")))
    estimated_barcodes = [row["barcode"] for row in ranked_rows[:expected_cells]]

    called_barcodes = read_barcodes_from_h5(args.filtered_h5)

    estimated_set = set(estimated_barcodes)
    called_set = set(called_barcodes)
    overlap = estimated_set & called_set

    estimated_only = estimated_set - called_set
    called_only = called_set - estimated_set
    union_n = len(estimated_set | called_set)
    jaccard = (len(overlap) / union_n) if union_n else 0.0

    write_lines(args.barcodes_out, sorted(called_set))

    with open(args.summary_out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sample_id",
                "estimated_cells",
                "cellbender_called",
                "overlap_barcodes",
                "estimated_only",
                "cellbender_only",
                "jaccard_index",
                "called_over_estimated_ratio",
                "total_droplets_included",
                "low_count_threshold",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": args.sample_id,
                "estimated_cells": len(estimated_set),
                "cellbender_called": len(called_set),
                "overlap_barcodes": len(overlap),
                "estimated_only": len(estimated_only),
                "cellbender_only": len(called_only),
                "jaccard_index": round(jaccard, 6),
                "called_over_estimated_ratio": round(len(called_set) / max(len(estimated_set), 1), 6),
                "total_droplets_included": total_droplets_included,
                "low_count_threshold": low_count_threshold,
            }
        )

    all_barcodes = sorted(estimated_set | called_set)
    with gzip.open(args.labels_out, "wt", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample_id", "barcode", "in_estimated", "in_cellbender", "class"],
        )
        writer.writeheader()
        for barcode in all_barcodes:
            in_est = barcode in estimated_set
            in_cb = barcode in called_set
            if in_est and in_cb:
                barcode_class = "overlap"
            elif in_est:
                barcode_class = "estimated_only"
            else:
                barcode_class = "cellbender_only"
            writer.writerow(
                {
                    "sample_id": args.sample_id,
                    "barcode": barcode,
                    "in_estimated": in_est,
                    "in_cellbender": in_cb,
                    "class": barcode_class,
                }
            )


if __name__ == "__main__":
    main()
