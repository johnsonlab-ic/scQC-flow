# Architecture

## Pipeline Flow
```
FASTQ → MAPPING → AMBIENT → QC → HVG → INTEGRATION → ANNOTATION → ZOOMS → EXPORT
                                                                  → REPORT_SITE
```

## Workflows & Processes

### MAPPING
Processes: `SIMPLEAF_INDEX` → `SIMPLEAF_QUANT` → `BARCODE_ESTIMATION` → `MAPPING_REPORT`  
Outputs: `h5_files` (sampleId, af_counts_mat.h5), `knee_data`, `ambient_params` (.env tuples), `report`  
Chemistry: user string → af_chemistry + orientation + whitelist_filename (lookup in workflows.nf)

### SAMPLE_METADATA
Process: `PREPARE_SAMPLE_METADATA`  
Outputs: `sample_metadata.csv.gz` (or NO_FILE placeholder)

### AMBIENT
Method selected via `--ambient_method` (default: `decontx`)  
**decontX path:** `DECONTX` → `AMBIENT_REPORT` + `STAGE_RAW_H5` + `AMBIENT_DE` (pseudobulk empty DE)  
**cellbender path:** `CELLBENDER` → `CELLBENDER_REPORT`  
Outputs: `h5_files` (filt_counts_<id>.h5), `barcodes`, `report`, `de_table` (edger_dt.csv.gz), `pb_empties`

### QC (2-pass caching design)
1. `DOUBLET_DETECTION` — scDblFinder; **cached** (reruns only if H5 changes)  
2. `APPLY_QC` — metrics + hard/soft filtering; **not cached** (reruns on threshold changes)  
3. `QC_REPORT`  
Outputs: `qc_metrics` (qc_metrics_<id>.csv.gz with all QC flags), `dbl_results`

### HVG (multi-sample, 2-pass in Python)
Process: `HVG_SELECTION` → `HVG_REPORT`  
Pass 1: per-sample VST stats; Pass 2: extract HVG-only, hstack combined matrix  
Input: collected H5s + QC CSVs + GTF + edger_dt.csv  
Outputs: `hvg_counts.h5` (singlets), `dbl_hvg_counts.h5` (doublets), `hvg_stats.csv.gz`

### INTEGRATION (2-pass Harmony in Python)
Process: `RUN_INTEGRATION` → `INTEGRATION_REPORT`  
Pass 1: singlet+doublet hstack, CPM10k+log1p, PCA → Harmony(theta=0) → Leiden(high res) → flag doublet-enriched clusters  
Pass 2: remove doublet clusters, re-PCA → Harmony(user_theta, user_vars) → Leiden(user_res_list) → UMAP  
Cluster names: cl01, cl02… (descending size)  
Outputs: `integration_dt.csv.gz` (embeddings + cluster cols + scdbl_class)

### ANNOTATION
Process: `RUN_ANNOTATION_MARKERS` → `ANNOTATION_REPORT`  
Method: edgeR pseudobulk one-vs-rest on selected resolution (`--annotation_sel_res`)  
Outputs: `annotation_markers.csv.gz`, `annotation_logcpms.csv.gz`, `annotation_marker_panel.csv.gz`, `annotation_cell_labels.csv.gz`

### ZOOMS
Processes: `PREPARE_ZOOM_SUBSET` → `ZOOM_HVG_SELECTION` → `RUN_ZOOM_INTEGRATION` → `RUN_ZOOM_MARKERS` → `ZOOM_REPORT`  
Subset by cluster ID or annotation label; reruns full HVG + integration on subset

### EXPORT
Processes: `EXPORT_SCANPY`, `EXPORT_SEURAT`  
Outputs: per-sample `*_all.h5ad|rds` + `*_clean.h5ad|rds`; optional combined objects

### REPORT_SITE
Process: `REPORT_SITE` — aggregates all HTML into landing page via `build_report_site.py`

## Process Labels (resource allocation)
| Label | CPUs (offline) | RAM |
|-------|---------------|-----|
| `process_single` | 1 | 6 GB |
| `process_low` | 2 | 12 GB |
| `process_medium` | 9 | 36 GB |
| `process_high` | 60 | 120 GB |
| `process_simpleaf` | 16 | 40 GB |
| `process_reports` | 60 | 80 GB |
| `process_cellbender` | 64 | 200 GB |
| `process_gpu` | 10 | 120 GB + GPU |

## Compute Profiles
`offline` (local Docker), `imperial` (PBS + Singularity), `slurm`, `dsi`
