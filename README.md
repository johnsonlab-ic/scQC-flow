# scQC-flow: Single-Cell QC and Reporting Pipeline

## Overview

**scQC-flow** is a modular Nextflow pipeline for single-cell RNA-seq quality control and reporting. It automates nuclear fraction analysis (DropletQC), doublet detection (scDblFinder), Seurat object creation (pre/post QC), and generates per-sample and combined HTML reports using Quarto.

## Workflow Steps

1. **Input Mapping**: Reads a CSV listing sample names and mapping directories (Cell Ranger outputs).
2. **DropletQC**: Computes nuclear fraction metrics for each sample.
3. **scDblFinder**: Detects doublets and annotates cells.
4. **Seurat Object Creation**: Builds Seurat objects (pre/post QC) with DropletQC and doublet metadata.
5. **Reporting**:
	 - Per-sample HTML reports (Quarto, customizable templates).
	 - Combined multi-sample report/book.

## Inputs

- `mapping_dirs` (CSV): Two columns, `samplename` and `path` (to mapping directory).
- Cell Ranger output folders (with `outs/filtered_feature_bc_matrix`).
- Optional: Custom report templates (Quarto QMD).

## Outputs

- Per-sample QC metrics (`*_dropletqc_metrics.csv`, `*_scdbl_metrics.csv`).
- Seurat RDS files (pre/post QC).
- Per-sample HTML reports (`*_qc_report.html`).
- Combined HTML report/book (`combined_qc_book/_book/index.html`).

## Usage

Basic run (offline profile, custom output folder):

```bash
nextflow run main.nf \
	--mapping_dirs personal/test_mapping_dirs.csv \
	--outputDir personal/outs \
	-profile offline \
	-resume
```

## Customization

- Templates: Edit Quarto QMD files in `modules/reports/` for custom report layouts.
- Parameters: Adjust QC thresholds and filtering in `make_seurat.R`.

## Requirements

- Nextflow
- Docker/Singularity (for containers)
- R (with Seurat, DropletQC, scDblFinder, Quarto installed in container)

## Citation

If you use scQC-flow, please cite the relevant tools (Nextflow, Seurat, DropletQC, scDblFinder, Quarto).
git merge --no-ff dropletQC -m "Merge dropletQC into main"git merge --no-ff dropletQC -m "Merge dropletQC into main"