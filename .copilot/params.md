# Parameters

## Required
| Param | Purpose |
|-------|---------|
| `--raw_data_dir` | Parent dir; each subdir = one sample (comma-separated for multiple) |
| `--cellrangerPath` | Cell Ranger install (for barcode whitelists) |
| `--genome_fasta` | Reference FASTA |
| `--genome_gtf` | Reference GTF |

## Workflow Control
| Param | Default | Notes |
|-------|---------|-------|
| `--chemistry` | `3v4` | 3v2, 3v3, 3v4, 3LT, 5v1, 5v2, 5v3, multiome |
| `--outputDir` | `results` | |
| `--run_ambient` | `true` | |
| `--ambient_method` | `decontx` | `decontx` or `cellbender` |
| `--run_qc` | `true` | requires `--run_ambient` |
| `--run_hvg` | `true` | requires `--run_qc` |
| `--run_integration` | `true` | requires `--run_hvg` + `--metadata_csv` |
| `--run_annotation` | `false` | requires `--run_integration` + `--annotation_marker_csv` |
| `--run_zooms` | `false` | requires `--run_integration` |
| `--export` | `none` | `none`, `anndata`, `seurat`, `both` |
| `--publish_mapping_simpleaf` | `false` | publish large simpleaf dirs |

## Metadata
| Param | Default | Notes |
|-------|---------|-------|
| `--metadata_csv` | — | Required for integration |
| `--metadata_id_col` | `sample_id` | Column mapping to sample IDs |
| `--metadata_vars` | — | Space-separated columns for Harmony batch correction |
| `--annotation_marker_csv` | — | Long-format marker CSV for annotation |

## Algorithm Parameters
| Param | Default | Purpose |
|-------|---------|---------|
| `--hvg_n_hvgs` | `2000` | HVGs to select |
| `--exclude_mito` | `false` | Exclude mito genes from HVG/integration |
| `--integration_n_dims` | `30` | PCA dimensions |
| `--integration_theta` | `2.0` | Harmony theta |
| `--integration_leiden_res` | `'0.2 0.5 1.0'` | Space-separated Leiden resolutions |
| `--integration_cluster_seed` | `42` | Random seed |
| `--integration_dbl_res` | `2.0` | Pass-1 doublet Leiden resolution |
| `--integration_dbl_cl_prop` | `0.5` | Doublet-enriched cluster threshold |
| `--annotation_sel_res` | `0.2` | Resolution used for annotation |
| `--annotation_min_cl_size` | — | Min cluster size |
| `--annotation_min_cells` | `10` | Min cells per pseudobulk group |
| `--annotation_top_n` | `10` | Top N markers per cluster |
| `--annotation_fdr_cut` | — | FDR cutoff for marker DE |
| `--export_write_combined` | `true` | Write combined export objects |

## QC Thresholds
| Param | Purpose |
|-------|---------|
| `--qc_hard_min_counts` | Hard minimum counts (drops cells before downstream) |
| `--qc_hard_min_feats` | Hard minimum features |
| `--qc_hard_max_mito` | Hard max mitochondrial % |
| `--qc_min_counts` | Soft min counts (flagged only) |
| `--qc_min_feats` | Soft min features |
| `--qc_max_mito` | Soft max mito % |
| `--qc_min_mito` | Soft min mito % |
| `--qc_max_splice` | Soft max splice % |
| `--qc_min_splice` | Soft min splice % |

## CellBender Parameters (only if `--ambient_method cellbender`)
| Param | Default |
|-------|---------|
| `--cellbender_env_name` | `cellbender` |
| `--cellbender_expected_cells` | `0` (auto) |
| `--cellbender_total_droplets_included` | `0` (auto) |
| `--cellbender_low_count_threshold` | `0` (auto) |
| `--cellbender_learning_rate` | — |
| `--cellbender_posterior_batch_size` | `0` |
| `--cellbender_mig_device` | — (NVIDIA MIG UUID) |
