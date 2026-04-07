
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
- Nextflow (>=25.04)
- Docker or Singularity

### Basic Usage

The recommended way to run the pipeline is to specify the output directory and work directory separately:

```bash
nextflow run main.nf -profile offline \
  -w ./my_project/work \
  --run_mode qc \
  --mapped_data_dir /path/to/mapped_dirs \
  --outputDir ./my_project/outputs
```

**Where:**
- `-w ./my_project/work` — Nextflow work directory (temporary files, intermediate results)
- `--outputDir ./my_project/outputs` — Final results and reports

This keeps your project organized with a clean separation:
```
my_project/
├── work/              # Nextflow work directory (can be safely deleted after run)
└── outputs/           # Final results
    ├── mapping/
    ├── cellbender/
    ├── dropletqc/
    ├── scdblfinder/
    ├── seurat/
    ├── reports/
    └── pipeline_info/ # Execution timeline, report, trace, DAG
```

All execution artifacts (timeline, report, trace, DAG) are collected into `pipeline_info/` for easier debugging.

### Run Mode Examples

**QC Only (pre-mapped data):**
```bash
MY_PROJECT="/data/my_project"
mkdir -p "$MY_PROJECT/work" "$MY_PROJECT/outputs"

nextflow run main.nf -profile offline \
  -w "$MY_PROJECT/work" \
  --run_mode qc \
  --mapped_data_dir /path/to/mapped_dirs \
  --outputDir "$MY_PROJECT/outputs"
```

**Mapping + QC:**
```bash
MY_PROJECT="/data/my_project"
mkdir -p "$MY_PROJECT/work" "$MY_PROJECT/outputs"

nextflow run main.nf -profile offline \
  -w "$MY_PROJECT/work" \
  --run_mode both \
  --raw_data_dir /path/to/raw_fastq_dirs \
  --cellrangerPath /path/to/cellranger \
  --transcriptome /path/to/refdata-gex-GRCh38-2024-A \
  --outputDir "$MY_PROJECT/outputs"
```

**With CellBender and GPU:**
```bash
nextflow run main.nf -profile offline \
  -w "$MY_PROJECT/work" \
  --run_mode qc \
  --mapped_data_dir /path/to/mapped_dirs \
  --cellbender true \
  --gpu true \
  --outputDir "$MY_PROJECT/outputs"
```

---

## Input Discovery

- For `run_mode = mapping` or `both`: set `--raw_data_dir`. Each subdirectory is treated as one sample (name = folder name).
- For `run_mode = qc`: set `--mapped_data_dir`. Each subdirectory is treated as one mapped sample (name = folder name).

No samplesheet CSV is needed — the pipeline discovers samples automatically from the directory structure.

---

## Execution Artifacts

The pipeline automatically generates execution reports collected in `${outputDir}/pipeline_info/`:

- **execution_timeline_*.html** — Interactive timeline showing task execution order and duration
- **execution_report_*.html** — Comprehensive execution summary (CPU, memory, retries, etc.)
- **execution_trace_*.txt** — Detailed task execution trace for debugging
- **pipeline_dag_*.html** — Visual DAG showing workflow dependencies

These are timestamped to prevent collisions across multiple runs. They are useful for:
- Understanding bottlenecks and resource usage
- Debugging failed tasks
- Auditing pipeline execution history

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

