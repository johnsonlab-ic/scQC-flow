# Annotation Helpers

This directory contains user-facing helper scripts for preparing reference inputs used by the annotation workflows.

Files:

- `train_xgboost_reference.R`: train a multiclass XGBoost model from a labelled `SingleCellExperiment` reference.
- `prepare_singler_reference.R`: clean and standardise a `SingleCellExperiment` object for use as a `SingleR` reference.

Typical flow:

1. Prepare or clean a reference with `prepare_singler_reference.R`.
2. Point a `singler` annotation-method config entry at the resulting `reference_rds` and `reference_label_col`.
3. Optionally train an XGBoost model from the same reference with `train_xgboost_reference.R`.
4. Point an `xgboost` annotation-method config entry at the resulting `model_rds` and `class_csv`.

Example config entry:

```groovy
params.annotation_methods = [
    [
        id: 'singler_kamath',
        engine: 'singler',
        reference_name: 'kamath',
        reference_rds: '/path/to/kamath_reference.rds',
        reference_label_col: 'broad_label',
        cluster_res: '0.2',
    ],
    [
        id: 'xgboost_kamath',
        engine: 'xgboost',
        reference_name: 'kamath',
        model_rds: '/path/to/xgboost_obj_kamath_broad.rds',
        class_csv: '/path/to/allowed_cls_kamath_broad.csv',
        cluster_res: '0.2',
    ],
]
```