
# scQC-flow (Nextflow re-implementation of scprocess)

## Based on [scprocess](https://www.biorxiv.org/content/10.64898/2026.03.09.710141v1)

Full credit for the analytical design, methods, and workflow goes to the scprocess authors.
scQC-flow is a Nextflow re-implementation of scprocess, adapted for HPC execution.
Please cite the scprocess publication if you use this pipeline.

## Overview

Nextflow DSL2 pipeline for single-cell/single-nucleus RNA-seq mapping, QC, and reporting.

**Features:**
- FASTQ → count matrix mapping via simpleaf (alevin-fry)
- Ambient RNA removal via decontX (or CellBender)
- Per-cell QC with hard and soft thresholds
- Doublet detection (scDblFinder, 2-pass design)
- HVG selection (Seurat VST, batch-aware)
- Harmony integration → Leiden clustering → UMAP
- Pseudobulk marker gene annotation (edgeR one-vs-rest)
- Cell type labelling from user-supplied marker panel
- Zoom subsets (re-analysis of cluster subsets)
- Export to AnnData (.h5ad) or Seurat (.rds)
- Automated HTML reports collected into a landing page

---

## Workflow and Module Map

Top-level workflow dispatch is in `main.nf` and workflow definitions are in `workflows/workflows.nf`.

| Workflow | Main modules/process groups |
|----------|-----------------------------|
| `MAPPING` | `SIMPLEAF_INDEX`, `SIMPLEAF_QUANT`, `BARCODE_ESTIMATION`, `MAPPING_REPORT` |
| `AMBIENT` | `DECONTX` or `CELLBENDER`, `AMBIENT_REPORT` or `CELLBENDER_REPORT`, `AMBIENT_DE_PROC` |
| `QC` | `DOUBLET_DETECTION`, `APPLY_QC`, `QC_REPORT` |
| `HVG` | `HVG_SELECTION`, `HVG_REPORT` |
| `INTEGRATION` | `RUN_INTEGRATION`, `INTEGRATION_REPORT` |
| `ANNOTATION` | `RUN_ANNOTATION_MARKERS`, `PREP_REPORT_INPUTS`, `ANNOTATION_REPORT` |
| `ZOOMS` | `PREPARE_ZOOM_SUBSET`, `ZOOM_AMBIENT_DE`, `ZOOM_HVG_SELECTION`, `RUN_ZOOM_INTEGRATION`, `RUN_ZOOM_MARKERS`, `ZOOM_REPORT` |
| Site/export | `REPORT_SITE`, `EXPORT_SCANPY`, `EXPORT_SEURAT`, `EXPORT_SCANPY_ZOOM`, `EXPORT_SEURAT_ZOOM` |

---

## Quick Start

### Requirements
- Nextflow (>=25.04)
- Docker (all non-mapping processes) or Singularity (HPC)
- conda (for simpleaf mapping processes)
- Cell Ranger installation (for barcode whitelists only)

### Basic Usage

```bash
nextflow run main.nf -profile offline \
  -w ./my_project/work \
  --raw_data_dir /path/to/fastq_parent_dir \
  --cellrangerPath /path/to/cellranger \
  --genome_fasta /path/to/genome.fa \
  --genome_gtf /path/to/genes.gtf \
  --metadata_csv /path/to/metadata.csv \
  --outputDir ./my_project/outputs
```

**Where:**
- `-w ./my_project/work` — Nextflow work directory (intermediate files; safe to delete after run)
- `--outputDir ./my_project/outputs` — Final results and reports

Output structure:
```
my_project/outputs/
├── mapping/           # simpleaf quant outputs + H5 matrices
├── ambient/           # decontX-filtered H5 + barcodes
├── qc/                # per-cell QC metrics CSVs
├── hvg/               # HVG count matrices + stats
├── integration/       # integration_dt.csv.gz (UMAP + clusters)
├── annotation/        # marker DE, logCPMs, cell labels
├── export/            # .h5ad / .rds exports
├── reports/           # all HTML reports + landing page
└── pipeline_info/     # timeline, trace, DAG
```

### Run Examples

**Mapping + full QC + integration + annotation:**
```bash
nextflow run main.nf -profile offline \
  --raw_data_dir /path/to/raw \
  --cellrangerPath /path/to/cellranger \
  --genome_fasta /path/to/genome.fa \
  --genome_gtf /path/to/genes.gtf \
  --metadata_csv /path/to/metadata.csv \
  --metadata_vars 'brainregion condition' \
  --run_annotation true \
  --annotation_marker_csv /path/to/markers.csv \
  --outputDir ./outputs
```

**Skip mapping (pre-mapped H5s already exist):**  
Not currently supported as a standalone mode — mapping is always run. Pre-existing mapping outputs can be reused via Nextflow caching if the work directory is preserved.

**Export to AnnData:**
```bash
nextflow run main.nf -profile offline \
  [... required params ...] \
  --export anndata \
  --outputDir ./outputs
```

---

## Input Discovery

- Each subdirectory of `--raw_data_dir` is treated as one sample; the folder name becomes the sample ID.
- No samplesheet CSV is needed.

No samplesheet CSV is needed — the pipeline discovers samples automatically from the directory structure.

---

## Execution Artifacts

All execution artifacts are collected in `${outputDir}/pipeline_info/`:
- **execution_timeline_*.html** — Task execution timeline
- **execution_report_*.html** — CPU, memory, retry summary
- **execution_trace_*.txt** — Detailed per-task trace
- **pipeline_dag_*.html** — Visual DAG of workflow dependencies

---

## Parameters

### Required
| Parameter | Description |
|-----------|-------------|
| `--raw_data_dir` | Parent dir of FASTQ sample folders (each subdir = one sample) |
| `--cellrangerPath` | Cell Ranger installation (for barcode whitelists) |
| `--genome_fasta` | Reference genome FASTA |
| `--genome_gtf` | Reference genome GTF |

### Workflow Control
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--chemistry` | `3v4` | 10x chemistry: `3v2`, `3v3`, `3v4`, `3LT`, `5v1`, `5v2`, `5v3`, `multiome` |
| `--outputDir` | `results` | Output directory |
| `--run_ambient` | `true` | Run ambient RNA removal |
| `--ambient_method` | `decontx` | `decontx` or `cellbender` |
| `--run_qc` | `true` | Run cell-level QC (requires `--run_ambient`) |
| `--run_hvg` | `true` | Run HVG selection (requires `--run_qc`) |
| `--run_integration` | `true` | Run Harmony integration (requires `--run_hvg` + `--metadata_csv`) |
| `--run_annotation` | `false` | Run pseudobulk annotation (requires `--run_integration` + `--annotation_marker_csv`) |
| `--export` | `none` | Export format: `none`, `anndata`, `seurat`, or `both` |

### Zoom Configuration
`ZOOMS` is controlled via nested `params.zoom` config rather than a `--run_zooms` flag.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `zoom.enabled` | `false` | Enable zoom workflows |
| `zoom.items` | `[]` | List of zoom specs (name, source, values/label, and optional overrides) |

### Metadata
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metadata_csv` | — | Sample metadata CSV (required for integration) |
| `--metadata_id_col` | `sample_id` | Column mapping to sample IDs |
| `--metadata_vars` | — | Space-separated columns for Harmony batch correction |
| `--annotation_marker_csv` | — | Long-format marker CSV (required for annotation); columns: `label`, `symbol` |

### QC Thresholds
| Parameter | Description |
|-----------|-------------|
| `--qc_hard_min_counts` | Hard min counts (cells removed before downstream) |
| `--qc_hard_min_feats` | Hard min features |
| `--qc_hard_max_mito` | Hard max mitochondrial % |
| `--qc_min_counts` | Soft min counts (flagged, retained) |
| `--qc_min_feats` | Soft min features |
| `--qc_max_mito` | Soft max mito % |
| `--qc_min_mito` | Soft min mito % |
| `--qc_max_splice` | Soft max splice % |
| `--qc_min_splice` | Soft min splice % |

### Algorithm Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hvg_n_hvgs` | `2000` | Number of HVGs to select |
| `--exclude_mito` | `false` | Exclude mito genes from HVG/integration |
| `--integration_n_dims` | `30` | PCA dimensions |
| `--integration_theta` | `2.0` | Harmony theta |
| `--integration_leiden_res` | `'0.2 0.5 1.0'` | Space-separated Leiden resolutions |
| `--integration_cluster_seed` | `42` | Random seed |
| `--integration_dbl_res` | `2.0` | Pass-1 doublet Leiden resolution |
| `--integration_dbl_cl_prop` | `0.5` | Doublet-enriched cluster threshold |
| `--annotation_sel_res` | `0.2` | Leiden resolution for annotation |
| `--annotation_min_cl_size` | — | Min cluster size for annotation |
| `--annotation_min_cells` | `10` | Min cells per pseudobulk sample |
| `--annotation_top_n` | `10` | Top N markers per cluster |
| `--export_write_combined` | `true` | Write combined export objects |

