# scQC-flow — Copilot Reference

Nextflow DSL2 pipeline: FASTQ → ambient removal → QC → HVG → integration → annotation → export.

## Entry Points
- `main.nf` — parameter validation, workflow dispatch, sample discovery
- `workflows/workflows.nf` — all workflow definitions
- `nextflow.config` — all params + defaults
- `configs/profiles.config` — compute profiles

## Key Facts
- Samples discovered from subdirs of `--raw_data_dir`; no samplesheet
- Channel pattern: `tuple(sampleId, file)` joined via `.join()`
- Docker container: `ghcr.io/johnsonlab-ic/landmark-sc_image` (all non-simpleaf processes)
- simpleaf processes use conda (`bioconda::simpleaf`)
- Execution artifacts: `${outputDir}/pipeline_info/`
- All HTML reports collected into landing page via `REPORT_SITE`

## Testing Rule
- Any ad hoc `Rscript` or `python` validation/test run must be executed inside `ghcr.io/johnsonlab-ic/landmark-sc_image`.
- Do not validate pipeline R/Python scripts against host-installed environments.
- Preferred pattern:
	`docker run --rm -u $(id -u):$(id -g) -v "$PWD":/work -w /work ghcr.io/johnsonlab-ic/landmark-sc_image <command>`

## Workflow Stages (all optional except MAPPING)
| Stage | Flag | Requires |
|-------|------|----------|
| MAPPING | always | — |
| AMBIENT | `--run_ambient true` | MAPPING |
| QC | `--run_qc true` | AMBIENT |
| HVG | `--run_hvg true` | QC |
| INTEGRATION | `--run_integration true` | HVG + `--metadata_csv` |
| ANNOTATION | `--run_annotation true` | INTEGRATION + `--annotation_marker_csv` |
| ZOOMS | `--run_zooms true` | INTEGRATION |
| EXPORT | `--export anndata\|seurat\|both` | INTEGRATION |

## Reference Files
- [architecture.md](architecture.md) — workflows, processes, inputs/outputs
- [params.md](params.md) — all params with defaults
- [modules.md](modules.md) — module files and scripts
- [design_notes.md](design_notes.md) — design decisions, caching, known issues
