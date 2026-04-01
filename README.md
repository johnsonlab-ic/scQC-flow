
# scQC-flow

Nextflow pipeline for single-cell/single-nucleus RNA-seq quality control and reporting.

**Features:**
- Cell Ranger mapping (optional, integrated)
- Optional CellBender ambient RNA removal (CPU/GPU)
- Nuclear fraction analysis (DropletQC)
- Doublet detection (scDblFinder)
- CellBender vs Cell Ranger droplet calling comparison with knee-plots
- Seurat object creation with pre/post-QC filtering
- Automated HTML reports

---

## Quick Start

### Requirements
- Nextflow
- Docker or Singularity

### Configuration

All run parameters live in the `params {}` block at the top of `nextflow.config`. Edit that file to set your inputs and options, then run:

```bash
nextflow run main.nf -profile offline
```

Or override any param from the CLI:

```bash
nextflow run main.nf -profile offline --run_mode mapping \
  --raw_data_dir /path/to/raw_fastq_dirs \
  --cellrangerPath /path/to/cellranger \
  --transcriptome /path/to/refdata-gex-GRCh38-2024-A \
  --outputDir results
```

### Run Mode: QC Only (pre-mapped data)

```bash
nextflow run main.nf -profile offline \
  --run_mode qc \
  --mapped_data_dir /path/to/mapped_dirs \
  --outputDir results
```

### Run Mode: Mapping + QC

```bash
nextflow run main.nf -profile offline \
  --run_mode both \
  --raw_data_dir /path/to/raw_fastq_dirs \
  --cellrangerPath /path/to/cellranger \
  --transcriptome /path/to/refdata-gex-GRCh38-2024-A \
  --outputDir results
```

---

## Input Discovery

- For `run_mode = mapping` or `both`: set `--raw_data_dir`. Each subdirectory is treated as one sample (name = folder name).
- For `run_mode = qc`: set `--mapped_data_dir`. Each subdirectory is treated as one mapped sample (name = folder name).

No samplesheet CSV is needed — the pipeline discovers samples automatically from the directory structure.

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_mode` | `qc` | `mapping`, `qc`, or `both` |
| `--raw_data_dir` | `null` | Required for `mapping`/`both`; parent dir of FASTQ sample folders |
| `--mapped_data_dir` | `null` | Required for `qc`; parent dir of mapped sample folders |
| `--mapper` | `cellranger` | Mapper to run for mapping stage |
| `--cellrangerPath` | `null` | Cell Ranger binary path or install directory |
| `--transcriptome` | `null` | Cell Ranger transcriptome reference directory |
| `--outputDir` | `results` | Output directory |
| `--cellbender` | `false` | Run CellBender ambient RNA removal |
| `--gpu` | `false` | Use GPU for CellBender |
| `--report` | `true` | Generate HTML QC reports |
| `--book` | `false` | Combine reports into a book |
| `--max_mito` | `20.0` | Max mitochondrial % threshold |
| `--min_nuclear` | `0.4` | Min nuclear fraction threshold |
| `--metadata` | `null` | Optional metadata CSV |

