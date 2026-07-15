#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


STAGE_ORDER = ["mapping", "ambient", "qc", "per_sample_annotation", "cellsweep",
               "hvg", "integration", "annotation", "zoom", "other"]

STAGE_LABELS = {
    "overview": "Overview",
    "mapping": "Mapping",
    "ambient": "Ambient",
    "qc": "QC",
    "per_sample_annotation": "Per-sample annotation",
    "cellsweep": "CellSweep",
    "hvg": "HVG",
    "integration": "Integration",
    "annotation": "Annotation",
    "zoom": "Zoom",
    "other": "Other",
}

KNOWN_REPORTS = {
    "mapping_report.html": {
        "title": "Mapping diagnostics",
        "stage": "mapping",
        "blurb": "Barcode knee plots and read-level mapping diagnostics across samples.",
    },
    "ambient_report.html": {
        "title": "Ambient RNA diagnostics",
        "stage": "ambient",
        "blurb": "decontX summaries, contamination estimates, and ambient cleanup diagnostics.",
    },
    "cellbender_report.html": {
        "title": "CellBender diagnostics",
        "stage": "ambient",
        "blurb": "Estimated-vs-called barcode comparisons and CellBender call diagnostics across samples.",
    },
    "qc_report.html": {
        "title": "Cell-level QC",
        "stage": "qc",
        "blurb": "Filtering metrics, threshold summaries, and QC visualizations.",
    },
    "hvg_report.html": {
        "title": "HVG selection",
        "stage": "hvg",
        "blurb": "Highly variable gene selection diagnostics and ambient-gene context.",
    },
    "integration_report.html": {
        "title": "Integration diagnostics",
        "stage": "integration",
        "blurb": "Harmony integration, clustering, metadata overlays, and QC-feature UMAPs.",
    },
    "annotation_report.html": {
        "title": "Annotation diagnostics",
        "stage": "annotation",
        "blurb": "Marker-panel matching, cluster labels, and marker expression views.",
    },
    "emptydrops_report.html": {
        "title": "Cell calling (EmptyDrops)",
        "stage": "ambient",
        "blurb": "DropletUtils emptyDrops cell calling: cells vs empty droplets per sample.",
    },
    "cell_calling_report.html": {
        "title": "Cell calling (GMM)",
        "stage": "ambient",
        "blurb": "Splice-aware GMM cell calling: nuclei/ambient/damaged populations per sample.",
    },
    "knee_report.html": {
        "title": "Mapping + cell calling (knee/shin)",
        "stage": "ambient",
        "blurb": "Read mapping diagnostics and classic barcode-rank cell calling (decontX-style selection), combined.",
    },
    "cellsweep_input_report.html": {
        "title": "CellSweep inputs",
        "stage": "cellsweep",
        "blurb": "Audit of exactly which barcodes feed CellSweep: ambient background, excluded droplets, and which QC-called cells are kept vs dropped before CellSweep runs.",
    },
    "cellsweep_report.html": {
        "title": "CellSweep decontamination",
        "stage": "cellsweep",
        "blurb": "Ambient-fraction (alpha_hat) estimates per cell and per cluster from the EM denoiser.",
    },
    "per_sample_annotation_report.html": {
        "title": "Per-sample annotation",
        "stage": "per_sample_annotation",
        "blurb": "Per-sample normalize → HVG → PCA → Leiden → XGBoost annotation, used as CellSweep's per-cell-type input labels.",
    },
    "index.html": {
        "title": "Run overview",
        "stage": "overview",
        "blurb": "Pipeline summary, effective configuration, and links to generated reports.",
    },
}

def pretty_method_name(method_name):
    tokens = method_name.split("_")
    pretty_tokens = []
    for token in tokens:
        lowered = token.lower()
        if lowered == "xgboost":
            pretty_tokens.append("XGBoost")
        elif lowered == "singler":
            pretty_tokens.append("SingleR")
        else:
            pretty_tokens.append(token.title())
    return " ".join(pretty_tokens)


def pretty_report_name(file_name):
    known = KNOWN_REPORTS.get(file_name)
    if known:
        return known["title"]
    if file_name.startswith("annotation_report_") and file_name.endswith(".html"):
        method_name = file_name[len("annotation_report_"):-len(".html")]
        return f"Annotation: {pretty_method_name(method_name)}"
    if file_name.startswith("zoom_") and file_name.endswith("_report.html"):
        zoom_name = file_name[len("zoom_"):-len("_report.html")]
        return f"Zoom: {zoom_name.replace('_', ' ').title()}"
    return file_name[:-5].replace("_", " ").title() if file_name.endswith(".html") else file_name


def report_stage(file_name):
    known = KNOWN_REPORTS.get(file_name)
    if known:
        return known["stage"]
    if file_name.startswith("annotation_report_") and file_name.endswith(".html"):
        return "annotation"
    if file_name.startswith("zoom_"):
        return "zoom"
    return "other"


def report_blurb(file_name):
    known = KNOWN_REPORTS.get(file_name)
    if known:
        return known["blurb"]
    if file_name.startswith("annotation_report_") and file_name.endswith(".html"):
        return "Per-method annotation diagnostics, prediction confidence summaries, and cluster-level label agreement views."
    if file_name.startswith("zoom_"):
        return "Subset rerun report for a configured zoom, including reintegration and top-marker summaries."
    return "Generated HTML report from this pipeline run."


def build_pages(html_paths):
    basenames = sorted({path.name for path in html_paths if path.suffix == ".html"})
    report_names = [name for name in basenames if name != "index.html"]
    report_names.sort(key=lambda name: (STAGE_ORDER.index(report_stage(name)), pretty_report_name(name).lower()))

    pages = [
        {
            "file": "index.html",
            "title": pretty_report_name("index.html"),
            "stage": "overview",
            "stage_label": STAGE_LABELS["overview"],
            "blurb": report_blurb("index.html"),
        }
    ]

    for name in report_names:
        stage = report_stage(name)
        pages.append(
            {
                "file": name,
                "title": pretty_report_name(name),
                "stage": stage,
                "stage_label": STAGE_LABELS.get(stage, "Report"),
                "blurb": report_blurb(name),
            }
        )

    return pages


def write_manifest(outdir, pages, payload):
    manifest = {
        "pipeline": payload.get("pipeline", "scQC-flow"),
        "run_name": payload.get("run_name"),
        "profile": payload.get("profile"),
        "generated_at": payload.get("generated_at"),
        "pages": pages,
    }
    (outdir / "site_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Assemble the HTML report site")
    parser.add_argument("html_files", nargs="+", help="Rendered HTML files to include in the site")
    parser.add_argument("--payload", required=True, help="JSON payload with run metadata")
    parser.add_argument("--outdir", required=True, help="Output directory for the assembled site")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    payload = json.loads(Path(args.payload).read_text(encoding="utf-8"))
    html_paths = [Path(path) for path in args.html_files]
    pages = build_pages(html_paths)
    page_files = {page["file"] for page in pages}

    for html_path in html_paths:
        if html_path.name not in page_files:
            continue
        (outdir / html_path.name).write_text(html_path.read_text(encoding="utf-8"), encoding="utf-8")

    write_manifest(outdir, pages, payload)


if __name__ == "__main__":
    main()