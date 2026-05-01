# Modules

## modules/mapping/
| File | Purpose |
|------|---------|
| `mapping.nf` | SIMPLEAF_INDEX, SIMPLEAF_QUANT, BARCODE_ESTIMATION |
| `barcode_estimation.R` | Knee detection (DropletUtils), S+U+A H5 creation, ambient params (.env output) |
| `chemistry_detection.py` | (unused in main flow) Auto-detect chemistry from FASTQ headers |

## modules/ambient/
| File | Purpose |
|------|---------|
| `ambient.nf` | DECONTX, CELLBENDER |
| `decontx.R` | decontXcounts; outputs filtered H5, cell barcodes, QC metrics |
| `cellbender_postprocess.py` | Post-process CellBender H5; outputs filtered H5, barcode labels |

## modules/ambient_de/
| File | Purpose |
|------|---------|
| `ambient_de.nf` | STAGE_RAW_H5, AMBIENT_DE |
| `ambient_de.R` | edgeR pseudobulk DE (empty droplets vs cells); outputs edger_dt.csv.gz, pb_empties.rds |

## modules/qc/
| File | Purpose |
|------|---------|
| `qc.nf` | DOUBLET_DETECTION, APPLY_QC, QC_REPORT |
| `doublet_detection.R` | scDblFinder (rRNA-excluded); cached |
| `apply_qc.R` | Per-cell metrics + hard/soft threshold filtering; not cached |
| `qc_plots.R` | Plotting helpers for qc_report.qmd |
| `sample_qc.R` | (SampleQC helpers — adaptive QC, currently unused in main flow) |

## modules/hvg/
| File | Purpose |
|------|---------|
| `hvg.nf` | HVG_SELECTION, HVG_REPORT |
| `hvg_selection.py` | Seurat VST HVG selection; 2-pass per-sample; batch_key=sample_id |
| `hvg_plots.R` | Variance plots, top gene list plots |

## modules/integration/
| File | Purpose |
|------|---------|
| `integration.nf` | RUN_INTEGRATION, INTEGRATION_REPORT |
| `run_integration.py` | 2-pass Harmony (Python); primary integration script |
| `run_integration.R` | R Harmony wrapper (legacy/alternative; not primary) |
| `integration_plots.R` | UMAP and cluster plotting helpers |

## modules/annotation/
| File | Purpose |
|------|---------|
| `annotation.nf` | RUN_ANNOTATION_MARKERS, ANNOTATION_REPORT |
| `marker_genes.R` | edgeR one-vs-rest pseudobulk marker DE |
| `annotation_utils.R` | Heatmap generation, helper functions |

## modules/zoom/
| File | Purpose |
|------|---------|
| `zoom.nf` | PREPARE_ZOOM_SUBSET, ZOOM_HVG_SELECTION, RUN_ZOOM_INTEGRATION, RUN_ZOOM_MARKERS, ZOOM_REPORT |
| `prepare_zoom_subset.py` | Subset QC metrics by cluster ID or annotation label |
| `zoom_markers.R` | Recompute top markers on zoomed subset |

## modules/metadata/
| File | Purpose |
|------|---------|
| `metadata.nf` | PREPARE_SAMPLE_METADATA |
| `prepare_metadata.py` | Validate and normalize metadata CSV; keyed by metadata_id_col |

## modules/export/
| File | Purpose |
|------|---------|
| `export.nf` | EXPORT_SCANPY, EXPORT_SEURAT |
| `export_anndata.py` | Export per-sample + combined .h5ad |
| `export_seurat.R` | Export per-sample + combined .rds |
| `export_utils.R` | Object construction and metadata joining helpers |

## modules/reports/
| File | Purpose |
|------|---------|
| `reports.nf` | All report processes + REPORT_SITE |
| `mapping_report.qmd` | Knee plots, simpleaf QC |
| `ambient_report.qmd` | Barcode rank plots, contamination metrics |
| `cellbender_report.qmd` | CellBender summary |
| `qc_report.qmd` | Per-sample QC scatter plots, hard/soft filter counts |
| `hvg_report.qmd` | HVG stats, top genes, ambient gene exclusion |
| `integration_report.qmd` | UMAP, cluster plots, doublet overlay |
| `annotation_report.qmd` | Heatmaps, marker plots |
| `zoom_report.qmd` | Per-zoom UMAP + markers |
| `build_report_site.py` | Aggregates all HTML into landing page |
