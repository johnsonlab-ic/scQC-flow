# Design Notes

## Caching Strategy
- `DOUBLET_DETECTION` is **cached** — reruns only when the input H5 changes (scDblFinder is slow)
- `APPLY_QC` is **not cached** — intentionally reruns whenever QC threshold params change
- This split allows fast threshold iteration without recomputing doublets

## QC Hard vs Soft Filters
- Hard filters (`qc_hard_*`) permanently remove cells before HVG/integration
- Soft filters (`qc_min/max_*`) flag cells in `qc_metrics.csv.gz` but retain them
- Hard-filtered cells are excluded from downstream entirely; soft-filtered retained for QC reporting

## Doublet Handling
- Doublets **retained** through HVG selection and integration (matching scprocess behaviour)
- Integration pass-1: singlet + doublet H5s are hstacked; Harmony(theta=0) + high-res Leiden identifies doublet-enriched clusters (threshold: `--integration_dbl_cl_prop`)
- Integration pass-2: doublet-enriched cluster cells removed; clean re-integration runs
- `scdbl_class` appended to `integration_dt.csv.gz` for post-hoc analysis
- Integration report visualises doublets on UMAP overlay

## Ambient Gene Exclusion
- `AMBIENT_DE` runs edgeR pseudobulk (empty droplets vs. cells) → `edger_dt.csv.gz`
- `HVG_SELECTION` uses `edger_dt.csv.gz` to exclude top ambient genes from HVG candidates
- Only applies when `--run_ambient true`

## Chemistry Resolution (workflows.nf)
- User provides `--chemistry` string (e.g. `3v4`)
- Lookup tables in `workflows.nf` resolve to: `af_chemistry`, `orientation`, `whitelist_filename`
- Whitelist read from `<cellrangerPath>/lib/python/cellranger/barcodes/` and decompressed inside SIMPLEAF_QUANT
- No auto-detection; user must provide correct chemistry

## Report Conventions (from /memories/repo/report_conventions.md)
- Workflow names: `HVG`, `INTEGRATION` (no `_WF` suffixes)
- Keep process `HVG_SELECTION` distinct from workflow `HVG` to avoid Nextflow symbol collisions
- Integration output has both `RNA_snn_res.<res>` (scprocess-style) and `leiden_<res>` aliases
- `INTEGRATION_REPORT` expects staged `qc_metrics_*.csv.gz` for cluster QC distributions
- Cell IDs in QC CSVs are bare barcodes — must prefix as `<sample_id>_<barcode>` before joining with integration
- `integration_report.qmd` infers metadata covariates directly from `integration_dt.csv.gz`
- Doublet UMAP sections conditional on doublets being present in integration output

## Container Strategy
- simpleaf processes: conda (`bioconda::simpleaf`)
- All other processes: Docker `ghcr.io/johnsonlab-ic/landmark-sc_image`
- `imperial` profile uses Singularity (same image)

## Validation Policy
- For manual checks, run `Rscript` and `python` commands inside `ghcr.io/johnsonlab-ic/landmark-sc_image`.
- Avoid host-shell validation for pipeline scripts to prevent environment drift from runtime.

## Sample Discovery
- Each subdir of `--raw_data_dir` is one sample — no explicit samplesheet
- `--raw_data_dir` accepts comma-separated list for multiple parent dirs
- `sampleId` derived from directory name

## Output Structure
- Per-sample outputs: `${outputDir}/<stage>/<sampleId>/`
- Reports: `${outputDir}/reports/`
- Pipeline info (timeline, trace, DAG): `${outputDir}/pipeline_info/`

## scprocess Reference
When asked "what does scprocess do" or to compare with scprocess — look in:
`/mnt/data/projects/scQC-flow/dev/scprocess/`

Key files:
- `rules/*.smk` — Snakemake rule definitions (workflow logic)
- `scripts/*.R` / `scripts/*.py` — implementation scripts
- `envs/*.yaml` — conda environments
- `profiles/` — compute profiles

Notable scprocess rules relevant to annotation:
- `rules/marker_genes.smk` → `scripts/marker_genes.R::calculate_marker_genes()` — pseudobulk edgeR with BiocParallel
- `rules/label_celltypes.smk` → `scripts/label_celltypes.py` (CelltypeList) + `scripts/label_celltypes.R` (XGBoost scprocess labeller)
- scprocess reads H5AD files per batch with parallel `bplapply`; builds pseudobulk per-sample before merging

## Known Limitations / Open Questions
- RNA-only scope; multiome (ATAC) not supported
- Single chemistry per run
- Adaptive QC (SampleQC) not implemented; `sample_qc.R` exists but is not wired in
- `run_integration.R` (R Harmony) exists but Python `run_integration.py` is primary
- `chemistry_detection.py` exists but not used in main flow
- `.claude/README.md` describes an older architecture (Cell Ranger, CellBender comparison, multiome) — the actual codebase has moved past this
