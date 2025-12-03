# scQC-flow Pipeline Overview

## Purpose

scQC-flow is a Nextflow DSL2 pipeline for automated single-cell/single-nucleus RNA-seq quality control and reporting. It processes Cell Ranger outputs through multiple QC steps and produces annotated Seurat objects along with HTML reports.

---

## Pipeline Steps

### 1. Input Parsing
- Reads a CSV file (`--mapping_dirs`) containing sample names and paths to Cell Ranger output directories.
- Each row maps a `samplename` to its corresponding `path` (Cell Ranger `outs/` parent directory).

### 2. CellBender (Optional)
- **When enabled** (`--cellbender`): Runs CellBender `remove-background` to remove ambient RNA contamination.
- Supports CPU or GPU execution (`--gpu`).
- Extracts "Estimated Number of Cells" from Cell Ranger's `metrics_summary.csv` and passes it to CellBender as `--expected-cells`.
- Converts CellBender's output HDF5 to a Seurat-compatible H5 format using `ptrepack`.
- **Downstream effect**: When CellBender is enabled, DropletQC and scDblFinder operate on CellBender outputs (barcodes/H5), ensuring consistent cell sets.

### 3. DropletQC (Nuclear Fraction)
- Computes nuclear fraction for each cell using the `DropletQC` R package.
- Inputs:
  - BAM file (`possorted_genome_bam.bam`) and its index (`.bai`)
  - Barcodes file (Cell Ranger or CellBender, depending on `--cellbender`)
- Outputs: `*_dropletqc_metrics.csv` with per-cell nuclear fraction values.

### 4. scDblFinder (Doublet Detection)
- Detects doublets using the `scDblFinder` R/Bioconductor package.
- Inputs:
  - H5 counts file (Cell Ranger or CellBender-converted H5, depending on `--cellbender`)
- Outputs: `*_scdbl_metrics.csv` with per-cell doublet scores and classifications.

### 5. Seurat Object Creation
- Creates Seurat objects from the H5 counts file (Cell Ranger or CellBender output).
- Merges DropletQC and scDblFinder metrics into cell metadata.
- Applies QC filtering based on:
  - `--max_mito` (maximum mitochondrial %)
  - `--min_nuclear` (minimum nuclear fraction)
- Outputs:
  - `*_preQC.rds` — Seurat object before QC filtering
  - `*_postQC.rds` — Seurat object after QC filtering

### 6. Quarto Reporting (Optional)
- **Per-sample reports** (`--report`, default `true`): Generates HTML QC reports for each sample showing QC metrics, filtering thresholds, and cell distributions.
- **Combined book** (`--book`): Combines all per-sample reports into a single Quarto book for easy cross-sample comparison.

---

## Workflow Diagram

```
                    ┌─────────────────┐
                    │  mapping_dirs   │
                    │     (CSV)       │
                    └────────┬────────┘
                             │
              ┌──────────────┴──────────────┐
              │                             │
      [--cellbender]                 [no cellbender]
              │                             │
              ▼                             │
      ┌───────────────┐                     │
      │  CELLBENDER   │                     │
      │  (CPU/GPU)    │                     │
      └───────┬───────┘                     │
              │                             │
              ▼                             │
      ┌───────────────┐                     │
      │ H5 CONVERT    │                     │
      │ (ptrepack)    │                     │
      └───────┬───────┘                     │
              │                             │
              ▼                             ▼
      ┌───────────────┐             ┌───────────────┐
      │  DROPLETQC    │             │  DROPLETQC    │
      │ (CB barcodes) │             │ (CR barcodes) │
      └───────┬───────┘             └───────┬───────┘
              │                             │
              ▼                             ▼
      ┌───────────────┐             ┌───────────────┐
      │   SCDBL       │             │   SCDBL       │
      │ (CB H5)       │             │ (CR H5)       │
      └───────┬───────┘             └───────┬───────┘
              │                             │
              └──────────────┬──────────────┘
                             │
                             ▼
                    ┌─────────────────┐
                    │  CREATE_SEURAT  │
                    │ (pre/post QC)   │
                    └────────┬────────┘
                             │
                             ▼
                    ┌─────────────────┐
                    │ GENERATE_REPORTS│
                    │  (per-sample)   │
                    └────────┬────────┘
                             │
                    [--book] │
                             ▼
                    ┌─────────────────┐
                    │ COMBINE_REPORTS │
                    │  (Quarto book)  │
                    └─────────────────┘
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mapping_dirs` | `personal/mapping_dirs.csv` | CSV with `samplename` and `path` columns |
| `--outputDir` | `results` | Output directory |
| `--cellbender` | `false` | Run CellBender ambient RNA removal |
| `--gpu` | `false` | Use GPU for CellBender (requires `--cellbender`) |
| `--report` | `true` | Generate per-sample HTML reports |
| `--book` | `false` | Combine reports into a Quarto book |
| `--max_mito` | `10.0` | Max mitochondrial % threshold |
| `--min_nuclear` | `0.4` | Min nuclear fraction threshold |
| `--metadata` | `null` | Optional metadata CSV |
| `--help` | `false` | Show help message |
