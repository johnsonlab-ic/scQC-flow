# scQC-flow — GPU handover / large-cohort integration

Context + instructions for running the **all-region (105-sample, ~2M-cell) integration on a
GPU system**. Target system needs **Nextflow + Docker (with NVIDIA GPU support)**.

---

## 1. Why we want a GPU

scQC-flow processes single-nucleus RNA-seq: mapping → cell-calling → CellSweep
decontamination → per-sample annotation → **pooled integration** → clustering/annotation →
per-cell-type "zoom" re-integrations. Everything except the **pooled integration** scales
linearly per-sample and is fine on CPU.

The bottleneck is the pooled integration at cohort scale (105 samples ≈ 2M cells):
`scale → PCA → Harmony → neighbors → Leiden → UMAP`. Harmony itself is cheap; the killers are
**neighbors + Leiden + UMAP on the full N** (largely single-threaded on CPU).

**What we established (see `dev/cellcalling/INVESTIGATION_SUMMARY.md`):**
- Brute-force CPU scanpy integration works fine up to ~876k cells (one region).
- **Seurat v5 sketch (BPCells) under-delivered** on CPU: >1h on 876k, 61 GB RAM, LeverageScore
  broken (fell back to Uniform), sketch stages not cheap. Not recommended.
- **scprocess** (trusted reference pipeline, handles 100s of samples;
  `dev/_archive/scprocess/scripts/integration.py`) scales purely by **GPU via
  `rapids_singlecell`** — full-data, standard Harmony, no sketching/batching. PAGA is used only
  as UMAP `init_pos`, not as a scaling method. Its CPU fallback = plain scanpy = what scQC-flow
  already does.

**Conclusion:** the proven route for 100s of samples is **GPU (rapids_singlecell), full-data**,
exactly as scprocess does. That's what this handover is for.

---

## 2. The one code change needed on the GPU box

`modules/integration/run_integration.py` currently imports `scanpy as sc` (CPU only). Add a
`--use_gpu` branch that swaps in `rapids_singlecell` — a near drop-in, and scprocess already has
the exact pattern to copy from (`dev/_archive/scprocess/scripts/integration.py`, `_do_one_integration`):

```python
# top of main(), after arg parsing
if use_gpu:
    import cupy as cp, rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(managed_memory=False, pool_allocator=True, devices=0)
    cp.cuda.set_allocator(rmm_cupy_allocator)
    import rapids_singlecell as sc          # replaces `import scanpy as sc`
else:
    import scanpy as sc

# in the integration body:
if use_gpu:
    sc.get.anndata_to_GPU(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=n_dims)
if use_gpu:
    sc.pp.harmony_integrate(adata, key=batch_var, max_iter_harmony=5, dtype=cp.float32, theta=theta)
else:
    import scanpy.external as sce; sce.pp.harmony_integrate(adata, key=batch_var, theta=theta)
sc.pp.neighbors(adata, n_pcs=n_dims, use_rep='X_pca_harmony')
for res in res_ls: sc.tl.leiden(adata, key_added=f"RNA_snn_res.{res}", resolution=float(res))
sc.tl.umap(adata, init_pos='paga' if use_paga else 'spectral')
if use_gpu:
    sc.get.anndata_to_CPU(adata)
```
`rapids_singlecell` mirrors the scanpy API for all of scale/pca/harmony/neighbors/leiden/umap, so
the rest of the function is unchanged. Wire a `--use_gpu` flag through `RUN_INTEGRATION` in
`modules/integration/integration.nf` and a `--integration_use_gpu` param in `nextflow.config`.

---

## 3. Container (IMPORTANT)

The current image `ghcr.io/johnsonlab-ic/landmark-sc_image` has **CPU scanpy only — NO rapids**.
On the GPU box you need a GPU-enabled image with **`rapids_singlecell` + `cupy` + `rmm`**
(RAPIDS ≥ 24.x, matching the CUDA version). Options:
- Build a GPU integration image FROM the rapidsai base (`rapidsai/base` or `rapidsai/notebooks`)
  + `pip install rapids-singlecell`, and use it only for the `RUN_INTEGRATION` process
  (`container` directive per-process; other processes keep the CPU image).
- scprocess's conda env spec (`dev/_archive/scprocess/envs/`) lists the working package set —
  use it as the reference for versions.
Nextflow must launch that process with GPU access: `docker.runOptions = '--gpus all'` (or
`accelerator` directive / `--nv` for Singularity).

---

## 4. Hardware guidance

- ~2M cells × 30 PCs harmonized embedding + kNN graph. A single **NVIDIA A100 (40 or 80 GB)** or
  similar is the comfortable target; 24 GB cards may OOM on the neighbors/UMAP step at 2M — if so,
  reduce to per-region integration (~700k each) or subsample for the graph.
- RMM pool allocator (above) is important for Harmony's iterative allocs.

---

## 5. How to run (once container + code change are in place)

```bash
nextflow run main.nf -profile offline -c <run.config> -w <workdir> \
  --cell_caller gmm \                 # recommended (see §7)
  --ambient_method cellsweep \
  --integration_use_gpu true \        # the new flag
  --raw_data_dir <dir with all 105 sample subdirs> \
  --cellrangerPath <cellranger install> \
  --genome_fasta <genome.fa> --genome_gtf <genes.gtf.gz> \
  --metadata_csv metadataFull-14.csv --metadata_id_col experimentid \
  --metadata_vars "brainregion condition brainbank ageatdeath gender" \
  --outputDir <out>
```
`metadata_vars` includes `brainregion` so Harmony corrects across the 3 regions into **one joint
representation**. `run.config` reference: `personal/scqc_flow_cs_runs/run.config`.

---

## 6. Data to carry to the new system

- Raw FASTQs: `mini-landmark/data/raw_data/raw_data_{SN,CG,FC}/` (105 sample subdirs) — OR the
  per-sample mapping/CellSweep outputs if you want to skip re-mapping (mapping is the bulk of
  compute + disk; ~1.4 TB workdir per region).
- Reference: genome FASTA + GTF (`resources/downloads/refdata-gex-GRCh38-2024-A/`).
- Cell Ranger barcode whitelists (`resources/downloads/cellranger-10.0.0/`).
- Metadata: `mini-landmark/metadata/metadataFull-14.csv` (105 rows; id col = `experimentid`).
- The repo (this pipeline) + this doc + `dev/cellcalling/INVESTIGATION_SUMMARY.md`.
- Disk: budget ~4–5 TB workdir for a from-scratch 105-sample run.

---

## 7. Science context to carry over (decisions already made)

Full detail: `dev/cellcalling/INVESTIGATION_SUMMARY.md`. Short version:
- **Endpoint:** dopaminergic (DA) neurons in dissected SNpc **PD** samples (TH/SLC6A3/SLC18A2/DDC).
- **Recommended caller: GMM + CellSweep.** GMM uniquely recovers the low-UMI/degenerating DA
  (depth-controlled: 3–5× DA-marker detection over the ambient floor); knee discards them.
- **CAVEAT:** the pipeline's flat `alpha_hat < 0.5` filter (`modules/second_pass/cellsweep_to_h5.py`)
  removes much of that low-UMI DA. For DA sensitivity, keep α̂ as a covariate or relax the cut —
  don't blindly flat-filter. (Baseline runs can keep the default; revisit for the DA analysis.)
- CPU marker scores are depth-biased — for low-UMI populations use depth-controlled DETECTION
  rates, not raw CPM.
- Rare cells are recovered by the per-cell-type **zoom** re-integrations, so global-integration
  rare-cell loss is acceptable at first pass.

---

## 8. If no GPU after all — CPU fallback

Brute-force scanpy (the current default) works to ~876k; for 2M on CPU expect it to be slow but
possibly acceptable on a many-core/large-RAM box (neighbors/UMAP are multithreaded in scanpy,
unlike the single-threaded Seurat path). If too slow, add a **lean python sketch** (subsample →
cluster/UMAP on subset → kNN-project the rest) — NOT Seurat sketch. Per-region integration
(~700k each) instead of one joint 2M embedding is the other CPU-friendly fallback.
