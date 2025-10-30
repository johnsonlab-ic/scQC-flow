
# 🧬 scQC-flow

<img src="https://img.shields.io/badge/Nextflow-v22.10.0+-green.svg" alt="Nextflow Version">
<img src="https://img.shields.io/badge/Containers-Docker%2FSingularity-orange.svg" alt="Container Support">

**A Nextflow pipeline for single-cell RNA-seq quality control and reporting**

scQC-flow automates nuclear fraction analysis (DropletQC), doublet detection (scDblFinder), Seurat object creation (pre/post QC), and generates per-sample and combined HTML reports using Quarto.

---

## 🚀 Quick Start

### Installation

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
```

### Running the Pipeline (from the website)

```bash
nextflow run johnsonlab-ic/scQC-flow \
	--mapping_dirs <path/to/mapping_dirs.csv> \
	--outputDir <output_directory> \
	-profile <profile> \
	-resume
```

Replace `<path/to/mapping_dirs.csv>`, `<output_directory>`, and `<profile>` as needed. For example:

```bash
nextflow run johnsonlab-ic/scQC-flow \
	--mapping_dirs personal/test_mapping_dirs.csv \
	--outputDir personal/outs \
	-profile offline \
	-resume
```

---

## 📁 Input Files & Parameters

### Input and Output Paths

| Parameter         | Description                                 | Example/Default |
|-------------------|---------------------------------------------|-----------------|
| `--mapping_dirs`  | CSV with `samplename` and `path` columns    | `personal/mapping_dirs.csv` |
| `--outputDir`     | Output directory for results                | `results`       |

### Pipeline Options (optional)

The pipeline accepts the following optional parameters. Defaults shown in parentheses.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--cellbender` | Run CellBender ambient RNA removal for all samples. When enabled, downstream QC and Seurat use CellBender outputs. | `false` |
| `--gpu` | Use GPU acceleration for CellBender (requires `--cellbender` and a CUDA-enabled environment). | `false` |
| `--report` | Generate per-sample Quarto QC reports (HTML). | `true` |
| `--book` | Combine all per-sample reports into a single Quarto book. | `false` |
| `--max_mito` | Maximum mitochondrial percentage threshold used during Seurat QC (percent). | `10.0` |
| `--min_nuclear` | Minimum nuclear fraction threshold used during Seurat QC (0-1). | `0.4` |
| `--metadata` | Optional path to a CSV metadata file (will be passed into Seurat/reporting). | `null` |
| `--help` | Show the pipeline help message and exit. | `false` |

Notes:
- `--cellbender` will cause the workflow to run CellBender (CPU or GPU) and convert the CellBender HDF5 output to a Seurat-compatible H5. DropletQC and scDblFinder will use CellBender barcodes/H5 when this flag is set.
- `--gpu` only affects CellBender execution. Ensure your environment has a compatible GPU, drivers, and that the container image supports GPU access (or run on a node with GPU passthrough and use the appropriate container runtime flags).
- Common Nextflow flags such as `-profile <profile>` and `-resume` are also supported as usual.

### Required Input Format

The mapping CSV should look like:

```csv
samplename,path
LM0094_SN1,/path/to/LM0094_SN1_snRNA_mapped
LM0094_SN2,/path/to/LM0094_SN2_snRNA_mapped
```

---

## 📋 Overview

The pipeline performs these key steps:

1. **Input Mapping**: Reads a CSV listing sample names and mapping directories (Cell Ranger outputs)
2. **DropletQC**: Computes nuclear fraction metrics for each sample
3. **scDblFinder**: Detects doublets and annotates cells
4. **Seurat Object Creation**: Builds Seurat objects (pre/post QC) with DropletQC and doublet metadata
5. **Reporting**: Generates per-sample and combined HTML reports using Quarto

---

## 📦 Outputs

| Output File/Folder                        | Description                                 |
|-------------------------------------------|---------------------------------------------|
| `*_dropletqc_metrics.csv`                 | Per-sample DropletQC metrics                |
| `*_scdbl_metrics.csv`                     | Per-sample scDblFinder metrics              |
| `*_qc_report.html`                        | Per-sample HTML QC reports                  |
| `combined_qc_book/_book/index.html`       | Combined multi-sample HTML report/book      |
| `*_preQC.rds`, `*_postQC.rds`             | Seurat RDS files (pre/post QC)              |

---

## ⚙️ Requirements

- Nextflow
- Docker or Singularity (for containers)

---

## ⚠️ Notes & Warnings

> **System Requirements**: This pipeline is computationally intensive and best run on HPC systems. It is optimized for the Imperial College HPC system but can be adapted to other systems with sufficient resources.

---

## 📚 Repository Information

<img src="https://img.shields.io/badge/GitHub-scQC--flow-lightgrey?logo=github" alt="GitHub Repo">

This pipeline is maintained in a public repository:
- [johnsonlab-ic/scQC-flow](https://github.com/johnsonlab-ic/scQC-flow)

---

## 📞 Support

For questions or issues, please:
- Open an issue on the [GitHub repository](https://github.com/johnsonlab-ic/scQC-flow/issues)
- Contact the Johnson Lab at Imperial College London

