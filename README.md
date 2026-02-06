
# scQC-flow

Nextflow pipeline for single-cell/single-nucleus RNA-seq quality control and reporting.

**Features:**
- Cell Ranger output processing
- Optional CellBender ambient RNA removal (CPU/GPU)
- Nuclear fraction analysis (DropletQC)
- Doublet detection (scDblFinder)
- CellBender vs Cell Ranger droplet calling comparison with knee-plots
- Seurat object creation with pre/post-QC filtering
- Automated HTML reports
- Support for 10x Multiome (GEX + ATAC)

---

## Quick Start

### Requirements
- Nextflow
- Docker or Singularity

### Basic Usage

```bash
nextflow run /path/to/scQC-flow \
  --mapping_dirs mapping_dirs.csv \
  --outputDir results
```

### With CellBender (GPU)

```bash
nextflow run /path/to/scQC-flow \
  --mapping_dirs mapping_dirs.csv \
  --outputDir results \
  --cellbender \
  --gpu \
  -profile imperial
```

### Multiome Data

```bash
nextflow run /path/to/scQC-flow \
  --mapping_dirs mapping_dirs.csv \
  --outputDir results \
  --multiome \
  --cellbender
```

---

## Input Format

CSV file with sample names and Cell Ranger output paths:

```csv
samplename,path
sample1,/path/to/sample1_mapped
sample2,/path/to/sample2_mapped
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mapping_dirs` | Required | CSV with samplename and path columns |
| `--outputDir` | `results` | Output directory |
| `--cellbender` | `false` | Run CellBender ambient RNA removal |
| `--gpu` | `false` | Use GPU for CellBender |
| `--multiome` | `false` | Process 10x Multiome data |
| `--report` | `true` | Generate HTML QC reports |
| `--book` | `false` | Combine reports into a book |
| `--max_mito` | `20.0` | Max mitochondrial % threshold |
| `--min_nuclear` | `0.4` | Min nuclear fraction threshold |
| `--metadata` | `null` | Optional metadata CSV |

---

## Outputs

| File | Description |
|------|-------------|
| `{sample}_dropletqc_metrics.csv` | Nuclear fraction metrics |
| `{sample}_scdbl_metrics.csv` | Doublet detection scores |
| `{sample}_cellbender_comparison_metrics.csv` | Droplet counts: raw vs CellRanger vs CellBender |
| `{sample}_cellbender_comparison_kneeplot.png` | Knee-plot showing droplet calling differences |
| `{sample}_seurat_object.rds` | Seurat object (pre-QC) |
| `{sample}_seurat_object_postqc.rds` | Seurat object (post-QC) |
| `{sample}_qc_report.html` | Per-sample QC report |
| `combined_qc_book/_book/index.html` | Combined report (if `--book` enabled) |

---

## Pipeline Steps

1. **Input**: Read sample mappings from CSV
2. **CellBender** (optional): Remove ambient RNA
3. **DropletQC**: Compute nuclear fraction per cell
4. **scDblFinder**: Detect doublets
5. **CellBender Comparison**: Analyze droplet calling differences (raw vs Cell Ranger vs CellBender)
6. **Seurat**: Create objects with QC metrics, apply filtering
7. **Reporting**: Generate HTML reports

---

## Multiome Workflow

When using `--multiome`:
- Automatically extracts Gene Expression and ATAC modalities
- Runs CellBender on GEX only (if enabled)
- Creates multiome-aware Seurat objects with both assays
- Preserves ATAC fragment and peak files
- Generates ATAC-specific reports

---

## Notes

- **CellBender**: Slow on CPU; use `--gpu` for reasonable runtime
- **BAM Index**: Pipeline expects `.bai` files alongside BAM files in Cell Ranger output
- **Profiles**: Use `-profile imperial` for Imperial HPC or configure for your system in `nextflow.config`
- **CellBender Comparison**: When CellBender is enabled, generates metrics CSV and knee-plot showing differences in droplet calling between Cell Ranger and CellBender

---

## Support

Issues? See the [GitHub repository](https://github.com/johnsonlab-ic/scQC-flow).

