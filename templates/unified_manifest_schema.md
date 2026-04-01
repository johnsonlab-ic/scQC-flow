# Unified Manifest Schema

This schema separates run-level configuration from sample rows.

## Files

- `templates/unified_run_config.yaml`: global run settings
- `templates/unified_samples_rna.csv`: one row per RNA sample

## Design rules

1. Keep global defaults in YAML.
2. Keep sample/library identity and input paths in CSV.
3. Allow sample-level overrides only for exceptional cases.

## Run config keys (`unified_run_config.yaml`)

Required top-level keys:

- `run_id`
- `output_dir`
- `execution`
- `assay`
- `mapping`
- `qc`
- `inputs`

Optional top-level keys:

- `cellbender`
- `reporting`
- `metadata`

### Execution block

- `run_mode` (`mapping`, `qc`, or `both`)
- `run_reporting` (bool)

### Assay block

- `type`: `rna`

### Mapping block

- `mapper`: `cellranger` or `alevin_fry`
- `cpus`, `memory_gb`, `time_hours`
- `cellranger` sub-block for Cell Ranger settings
- `alevin_fry` sub-block for alevin-fry settings

### CellBender block

- `enabled` (bool)
- `mode`: `cpu` or `gpu`
- key numeric controls such as `expected_cells`, `fpr`, `epochs`

### QC block

- `max_mito`
- `min_nuclear`
- module toggles (`run_dropletqc`, `run_scdblfinder`)
- optional Seurat cutoffs

### Metadata block

- `path`
- `join_key`
- `vars` (list of covariates)

### Inputs block

- `rna_samplesheet`: path to RNA CSV or `null`

Validation rule:

- if `run_mode=mapping` or `run_mode=both`, `rna_samplesheet` must be set.

## RNA sample CSV columns (`unified_samples_rna.csv`)

Required:

- `sample_id`
- `sample_name`
- `fastq_path`

Optional:

- `condition`, `brainregion`, `donor`, `batch`
- `mapper_override` (`cellranger` or `alevin_fry`)
- `cellbender_override` (`true` or `false`)
- `metadata_id`

## Precedence

From highest to lowest:

1. Sample-level override columns in CSV
2. Run-level YAML values
3. Pipeline defaults in `main.nf`

## Suggested pipeline parameter mapping

- `qc.max_mito` -> `--max_mito`
- `qc.min_nuclear` -> `--min_nuclear`
- `metadata.path` -> `--metadata`
- `cellbender.enabled` -> `--cellbender`
- `cellbender.mode=gpu` -> `--gpu`
