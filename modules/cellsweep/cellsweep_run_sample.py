#!/usr/bin/env python3
"""
Per-sample CellSweep decontamination test (prototype).

Builds the CellSweep input from the pipeline outputs and runs denoise_count_matrix:
  - cells   = post-QC singlet nuclei (integration cells) with their Leiden cluster as
              `celltype` (CellSweep needs coherent groups, not biological names)
  - empties = a subsample of is_empty barcodes (raw counts), is_empty=True -> ambient profile
Gene space = summed S+U+A (Ensembl IDs). Outputs the denoised AnnData (alpha_hat per cell)
+ a per-cell alpha_hat table.

Usage:
  cellsweep_run_sample.py --sample_id <sid> --filt_h5 <filt_counts.h5> --raw_h5 <af_counts_mat.h5> \
      --empty_barcodes <empty_barcodes.csv> --integration_dt <integration_dt.csv.gz> \
      --out_dir <dir> [--n_empties 30000] [--celltype_col RNA_snn_res.0.5]
"""
import argparse, re, os
import numpy as np, pandas as pd, scipy.sparse as sp
import h5py, anndata as ad
import cellsweep


def load_sua(h5_path):
    """Read stacked S/U/A 10x H5 -> (cells x genes CSR, barcodes, gene_ids), summed per gene."""
    with h5py.File(h5_path, "r") as f:
        g = f["matrix"]
        data = g["data"][:]; indices = g["indices"][:]; indptr = g["indptr"][:]
        shape = g["shape"][:]
        bc = np.array([b.decode().split("-")[0] for b in g["barcodes"][:]])
        feat = np.array([x.decode() for x in g["features"]["name"][:]])
    m = sp.csc_matrix((data, indices, indptr), shape=tuple(shape))  # features x cells
    base = np.array([re.sub(r"_(S|U|A)$", "", x) for x in feat])
    genes, inv = np.unique(base, return_inverse=True)
    S = sp.csr_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))), shape=(len(genes), len(base)))
    summed = (S @ m).tocsr()              # genes x cells
    return summed.T.tocsr(), bc, genes    # cells x genes


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample_id", required=True)
    p.add_argument("--filt_h5", required=True)
    p.add_argument("--raw_h5", required=True)
    p.add_argument("--empty_barcodes", required=True)
    p.add_argument("--integration_dt", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--n_empties", type=int, default=30000)
    p.add_argument("--celltype_col", default="RNA_snn_res.0.5")
    p.add_argument("--keep_h5ad", action="store_true",
                   help="keep the (large) decontaminated h5ad; default drops it, keeps alpha_hat CSV")
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    sid = a.sample_id
    print(f"== CellSweep: {sid} ==", flush=True)

    # --- cells: post-QC singlet nuclei from integration, celltype = Leiden cluster ---
    integ = pd.read_csv(a.integration_dt)
    integ = integ[integ["sample_id"] == sid].copy()
    if "is_dbl" in integ.columns:
        integ = integ[integ["is_dbl"] != True]
    integ["barcode"] = integ["cell_id"].str.replace(f"^{re.escape(sid)}_", "", regex=True)
    integ = integ[integ[a.celltype_col].notna() & (integ[a.celltype_col].astype(str) != "")]
    ct = dict(zip(integ["barcode"], integ[a.celltype_col].astype(str)))
    print(f"  integration singlet cells: {len(ct)}", flush=True)

    X, bc, genes = load_sua(a.filt_h5)
    keep = np.array([b in ct for b in bc])
    cells = ad.AnnData(X=X[keep], obs=pd.DataFrame({"barcode": bc[keep]}),
                       var=pd.DataFrame(index=genes))
    cells.obs["celltype"] = [ct[b] for b in cells.obs["barcode"]]
    cells.obs["is_empty"] = False

    # --- empties: subsample of is_empty barcodes from raw counts ---
    empty_bc = set(pd.read_csv(a.empty_barcodes, header=None)[0].astype(str))
    Xr, bcr, genesr = load_sua(a.raw_h5)
    is_emp = np.array([b in empty_bc for b in bcr])
    emp_idx = np.where(is_emp)[0]
    rng = np.random.default_rng(42)
    if len(emp_idx) > a.n_empties:
        emp_idx = rng.choice(emp_idx, a.n_empties, replace=False)
    empties = ad.AnnData(X=Xr[emp_idx], obs=pd.DataFrame({"barcode": bcr[emp_idx]}),
                         var=pd.DataFrame(index=genesr))
    empties.obs["celltype"] = "empty"
    empties.obs["is_empty"] = True
    print(f"  empties used: {empties.n_obs} (of {is_emp.sum()})", flush=True)

    adata = ad.concat([cells, empties], join="inner")
    adata.obs_names = (adata.obs["barcode"].values + "_" +
                       np.where(adata.obs["is_empty"], "E", "C"))
    adata.obs["celltype"] = adata.obs["celltype"].astype("category")
    print(f"  combined adata: {adata.n_obs} x {adata.n_vars}; "
          f"celltypes={adata.obs.loc[~adata.obs.is_empty,'celltype'].nunique()}", flush=True)

    out_h5 = os.path.join(a.out_dir, f"{sid}_cellsweep.h5ad")
    res = cellsweep.denoise_count_matrix(adata, adata_out=out_h5, round_X=True,
                                         threads=8, verbose=1, random_state=42)

    cells_res = res[~res.obs["is_empty"].values]
    al = cells_res.obs["alpha_hat"].values
    print(f"  alpha_hat: median={np.median(al):.3f} mean={np.mean(al):.3f} "
          f">0.5={np.mean(al>0.5)*100:.1f}%", flush=True)
    cells_res.obs[["barcode", "celltype", "alpha_hat"]].to_csv(
        os.path.join(a.out_dir, f"{sid}_alpha_hat.csv.gz"), index=False)
    if not a.keep_h5ad:
        os.remove(out_h5)
        print(f"  dropped {out_h5} (slim); kept {sid}_alpha_hat.csv.gz", flush=True)
    else:
        print(f"  wrote {out_h5} and {sid}_alpha_hat.csv.gz", flush=True)


if __name__ == "__main__":
    main()
