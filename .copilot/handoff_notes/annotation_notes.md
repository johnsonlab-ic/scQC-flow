Sure, here is a copy-paste summary you can give Agent 2:

I want to run a SingleR annotation test on my own dataset in two stages.

Current status:
- The reference build script is working and creates:
  kamath_singleR_reference.rds
- Query data is in AnnData format at:
  /workspace/data/mini_landmark_SN/outputs/export/anndata
- For testing, use one clean sample first:
  sample_LM0001_SN1_snRNAseq_clean.h5ad

What I want right now:
- Build/run the annotation test on one sample only
- Use SingleR at cluster level (not per-cell)
- Use broad_label first (from the Kamath reference)

Useful context:
- In the .h5ad, counts are available in layer: counts
- The sample has cluster-like columns in obs (e.g. RNA_snn_res.* and annotation_cluster)
- broad_label exists in the Kamath reference and should be used as labels for this first test

Goal for this test:
- A simple interactive annotation script flow (not over-engineered) that:
  1) loads the Kamath reference RDS
  2) loads one clean .h5ad query
  3) aligns genes
  4) runs SingleR cluster-level with broad_label
  5) writes basic outputs (cluster labels/scores, per-cell mapped labels)


Everything you need to focus on (except the input data, the .h5ad files) is in ${WORKSPACE}/dev/annotation . Zellkonverter is installed if needed.