#!/usr/bin/env python3

import argparse
import html
import json
import re
import shutil
from pathlib import Path


STAGE_ORDER = ["mapping", "ambient", "qc", "hvg", "integration", "annotation", "zoom", "other"]

STAGE_LABELS = {
    "overview": "Overview",
    "mapping": "Mapping",
    "ambient": "Ambient",
    "qc": "QC",
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
    "index.html": {
        "title": "Run overview",
        "stage": "overview",
        "blurb": "Pipeline summary, effective configuration, and links to generated reports.",
    },
}

BODY_OPEN_RE = re.compile(r"<body([^>]*)>", re.IGNORECASE)
BODY_CLOSE_RE = re.compile(r"</body>", re.IGNORECASE)
HEAD_CLOSE_RE = re.compile(r"</head>", re.IGNORECASE)
CLASS_ATTR_RE = re.compile(r'class=("|\')(.*?)(\1)', re.IGNORECASE | re.DOTALL)
QUARTO_MARGIN_SIDEBAR_RE = re.compile(
    r'<div id="quarto-margin-sidebar"[^>]*>.*?</div>\s*(?=<main class="content" id="quarto-document-content">)',
    re.IGNORECASE | re.DOTALL,
)
QUARTO_CONTENT_RE = re.compile(
    r'<div id="quarto-content" class="[^"]*">',
    re.IGNORECASE,
)
EXISTING_SHELL_OPEN_RE = re.compile(
    r'<div class="report-site-shell">.*?<div class="report-site-content">',
    re.IGNORECASE | re.DOTALL,
)
EXISTING_SHELL_CLOSE_RE = re.compile(
    r'</div></main></div>\s*<script src="site.js"></script>\s*',
    re.IGNORECASE | re.DOTALL,
)
SITE_CSS_LINK_RE = re.compile(r'\s*<link rel="stylesheet" href="site.css">\s*', re.IGNORECASE)


def pretty_report_name(file_name):
    known = KNOWN_REPORTS.get(file_name)
    if known:
        return known["title"]
    if file_name.startswith("zoom_") and file_name.endswith("_report.html"):
        zoom_name = file_name[len("zoom_"):-len("_report.html")]
        return f"Zoom: {zoom_name.replace('_', ' ').title()}"
    return file_name[:-5].replace("_", " ").title() if file_name.endswith(".html") else file_name


def report_stage(file_name):
    known = KNOWN_REPORTS.get(file_name)
    if known:
        return known["stage"]
    if file_name.startswith("zoom_"):
        return "zoom"
    return "other"


def report_blurb(file_name):
    known = KNOWN_REPORTS.get(file_name)
    if known:
        return known["blurb"]
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


def group_pages(pages):
    groups = []

    overview_pages = [page for page in pages if page["stage"] == "overview"]
    core_pages = [page for page in pages if page["stage"] not in {"overview", "zoom", "other"}]
    zoom_pages = [page for page in pages if page["stage"] == "zoom"]
    other_pages = [page for page in pages if page["stage"] == "other"]

    if overview_pages:
        groups.append(("Overview", overview_pages))
    if core_pages:
        groups.append(("Reports", core_pages))
    if zoom_pages:
        groups.append(("Zooms", zoom_pages))
    if other_pages:
        groups.append(("Other", other_pages))

    return groups


def render_nav(groups, current_file):
    parts = []
    for heading, pages in groups:
        parts.append(f'<div class="report-site-nav-group"><div class="report-site-nav-heading">{html.escape(heading)}</div>')
        for page in pages:
            active = " active" if page["file"] == current_file else ""
            parts.append(
                "<a class=\"report-site-nav-link%s\" href=\"%s\">%s</a>"
                % (
                    active,
                    html.escape(page["file"]),
                    html.escape(page["title"]),
                )
            )
        parts.append("</div>")
    return "".join(parts)


def add_body_class(html_text, class_name):
    match = BODY_OPEN_RE.search(html_text)
    if not match:
        raise ValueError("Could not locate <body> tag in HTML document")

    body_tag = match.group(0)
    attrs = match.group(1)
    class_match = CLASS_ATTR_RE.search(attrs)
    if class_match:
        classes = class_match.group(2).split()
        if class_name not in classes:
            classes.append(class_name)
        replacement = 'class="%s"' % " ".join(classes)
        new_tag = body_tag.replace(class_match.group(0), replacement, 1)
    else:
        new_tag = f"<body{attrs} class=\"{class_name}\">"

    return html_text[:match.start()] + new_tag + html_text[match.end():]


def inject_assets(html_text):
    asset_block = '  <link rel="stylesheet" href="site.css">\n'
    if HEAD_CLOSE_RE.search(html_text):
        return HEAD_CLOSE_RE.sub(asset_block + "</head>", html_text, count=1)
    return asset_block + html_text


def unwrap_existing_shell(html_text):
    html_text = SITE_CSS_LINK_RE.sub("\n", html_text, count=1)
    html_text = EXISTING_SHELL_OPEN_RE.sub("", html_text, count=1)
    html_text = EXISTING_SHELL_CLOSE_RE.sub("", html_text, count=1)
    return html_text


def simplify_quarto_layout(html_text):
    html_text = unwrap_existing_shell(html_text)
    html_text = QUARTO_MARGIN_SIDEBAR_RE.sub("", html_text, count=1)
    html_text = QUARTO_CONTENT_RE.sub(
        '<div id="quarto-content" class="report-site-quarto-content">',
        html_text,
        count=1,
    )
    html_text = html_text.replace(
        '<main class="content" id="quarto-document-content">',
        '<main class="content report-site-quarto-document" id="quarto-document-content">',
        1,
    )
    return html_text


def inject_shell(html_text, page, pages, payload):
    run_name = payload.get("run_name") or "run"
    profile = payload.get("profile") or "default"
    generated_at = payload.get("generated_at") or ""

    groups = group_pages(pages)
    nav_html = render_nav(groups, page["file"])
    overview_link = '' if page["file"] == "index.html" else '<a class="report-site-overview-link" href="index.html">Overview</a>'

    shell_open = (
        '<div class="report-site-shell">'
        '<aside class="report-site-sidebar" id="report-site-sidebar">'
        '<div class="report-site-brand">'
        '<a class="report-site-brand-link" href="index.html">scQC-flow</a>'
        f'<div class="report-site-brand-meta">{html.escape(run_name)}</div>'
        f'<div class="report-site-brand-submeta">Profile: {html.escape(profile)}</div>'
        '</div>'
        '<nav class="report-site-nav">'
        f'{nav_html}'
        '</nav>'
        f'<div class="report-site-footer">Generated {html.escape(generated_at)}</div>'
        '</aside>'
        '<div class="report-site-overlay" data-report-site-overlay hidden></div>'
        '<main class="report-site-main">'
        '<div class="report-site-toolbar">'
        '<button class="report-site-toggle" type="button" data-report-site-toggle aria-controls="report-site-sidebar" aria-expanded="false">Reports</button>'
        f'<div class="report-site-toolbar-label">{html.escape(page["title"])}</div>'
        f'{overview_link}'
        '</div>'
        '<div class="report-site-content">'
    )
    shell_close = '</div></main></div>\n  <script src="site.js"></script>\n'

    wrapped = simplify_quarto_layout(html_text)
    wrapped = add_body_class(wrapped, "report-site-page")

    body_open_match = BODY_OPEN_RE.search(wrapped)
    if not body_open_match:
        raise ValueError("Could not locate <body> tag after class injection")

    wrapped = wrapped[:body_open_match.end()] + shell_open + wrapped[body_open_match.end():]

    if BODY_CLOSE_RE.search(wrapped):
        wrapped = BODY_CLOSE_RE.sub(shell_close + "</body>", wrapped, count=1)
    else:
        wrapped += shell_close

    return inject_assets(wrapped)


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
    parser = argparse.ArgumentParser(description="Assemble wrapped HTML report site")
    parser.add_argument("html_files", nargs="+", help="Rendered HTML files to include in the site")
    parser.add_argument("--payload", required=True, help="JSON payload with run metadata")
    parser.add_argument("--outdir", required=True, help="Output directory for the assembled site")
    parser.add_argument("--css", required=True, help="Shared site CSS asset")
    parser.add_argument("--js", required=True, help="Shared site JS asset")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    payload = json.loads(Path(args.payload).read_text(encoding="utf-8"))
    html_paths = [Path(path) for path in args.html_files]
    pages = build_pages(html_paths)
    page_map = {page["file"]: page for page in pages}

    shutil.copy2(args.css, outdir / "site.css")
    shutil.copy2(args.js, outdir / "site.js")

    for html_path in html_paths:
        page = page_map.get(html_path.name)
        if page is None:
            continue
        html_text = html_path.read_text(encoding="utf-8")
        wrapped = inject_shell(html_text, page, pages, payload)
        (outdir / html_path.name).write_text(wrapped, encoding="utf-8")

    write_manifest(outdir, pages, payload)


if __name__ == "__main__":
    main()