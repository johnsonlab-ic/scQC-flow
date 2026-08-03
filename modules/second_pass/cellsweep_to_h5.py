#!/usr/bin/env python3
# cellsweep_to_h5.py
#
# Second-pass adapter: turn a per-sample CellSweep-decontaminated AnnData
# (<sid>_cellsweep.h5ad) back into the exact per-sample contract the existing
# downstream (HVG -> integration -> annotation) already consumes:
#
#   filt_counts_<sid>.h5     stacked S/U/A 10x H5 (genes x cells, CSC). CellSweep
#                            gives one corrected value per gene, so the corrected
#                            count goes in the _S row and _U/_A rows are zero.
#                            This keeps sum_sua_counts()/_sum_sua() (which REQUIRE
#                            all three S/U/A blocks) working unchanged.
#   qc_metrics_<sid>.csv.gz  minimal QC table flagging every kept cell as a
#                            QC-passing singlet (cells are already the pass-1
#                            QC singlets; CellSweep dropped doublets), with the
#                            `sum`/`detected` columns integration validates.
#
# Cells kept = is_empty==False AND alpha_hat < --alpha_max (drop high-ambient).

import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import anndata as ad


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample_id", required=True)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--pass1_qc", required=True,
                   help="pass-1 qc_metrics_<sid>.csv.gz — carries metadata columns "
                        "(brainregion/condition/...) that integration reads from qc")
    p.add_argument("--alpha_max", type=float, default=0.5)
    p.add_argument("--out_h5", required=True)
    p.add_argument("--out_qc", required=True)
    a = p.parse_args()
    sid = a.sample_id

    adata = ad.read_h5ad(a.h5ad)
    adata = adata[adata.obs["is_empty"].values == False].copy()
    n_all = adata.n_obs
    keep = adata.obs["alpha_hat"].values < a.alpha_max
    adata = adata[keep].copy()
    print(f"{sid}: cells {n_all} -> {adata.n_obs} kept (alpha_hat < {a.alpha_max})", flush=True)
    if adata.n_obs == 0:
        raise SystemExit(f"{sid}: no cells left after alpha filter")

    # deterministic gene order (same reference across samples -> identical order)
    order = np.argsort(adata.var_names.values)
    genes = adata.var_names.values[order]
    Xgc = adata.X.T.tocsr()[order, :].tocsc()        # genes x cells, corrected counts
    barcodes = adata.obs["barcode"].astype(str).values

    ngenes, ncells = Xgc.shape
    # stacked rows: [g_S ...] then [g_U ...] then [g_A ...]; counts only in _S block
    zero_block = sp.csc_matrix((ngenes, ncells), dtype=Xgc.dtype)
    stacked = sp.vstack([Xgc, zero_block, zero_block], format="csc")
    features = np.concatenate([[f"{g}_S" for g in genes],
                               [f"{g}_U" for g in genes],
                               [f"{g}_A" for g in genes]])

    with h5py.File(a.out_h5, "w") as f:
        f.create_dataset("matrix/indptr",        data=stacked.indptr,  compression="gzip")
        f.create_dataset("matrix/indices",       data=stacked.indices, compression="gzip")
        f.create_dataset("matrix/data",          data=stacked.data,    compression="gzip")
        f.create_dataset("matrix/shape",         data=np.array(stacked.shape))
        f.create_dataset("matrix/features/name", data=np.array(features, dtype="S"), compression="gzip")
        f.create_dataset("matrix/barcodes",      data=np.array(barcodes, dtype="S"), compression="gzip")

    per_cell_sum = np.asarray(Xgc.sum(axis=0)).ravel()
    per_cell_det = np.asarray((Xgc > 0).sum(axis=0)).ravel()

    # Carry over pass-1 QC metadata (brainregion/condition/... that integration
    # reads from qc) for the kept cells; override the QC-state columns to reflect
    # the corrected, alpha-filtered singlet set.
    qc1 = pd.read_csv(a.pass1_qc)
    qc1["cell_id"] = qc1["cell_id"].astype(str)
    qc = qc1.set_index("cell_id").reindex(barcodes)
    if qc.isna().all(axis=1).any():
        missing = int(qc.isna().all(axis=1).sum())
        raise SystemExit(f"{sid}: {missing} kept cells absent from pass-1 qc")
    qc = qc.reset_index().rename(columns={"index": "cell_id"})
    qc["sample_id"]   = sid
    qc["sum"]         = per_cell_sum.astype(np.int64)
    qc["detected"]    = per_cell_det.astype(np.int64)
    qc["scdbl_class"] = "singlet"
    qc["scdbl_score"] = 0.0
    qc["is_singlet"]  = True
    qc["keep"]        = True
    qc["alpha_hat"]   = adata.obs["alpha_hat"].values
    qc.to_csv(a.out_qc, index=False)
    print(f"{sid}: wrote {a.out_h5} ({stacked.shape[0]}x{stacked.shape[1]}) + {a.out_qc}", flush=True)


if __name__ == "__main__":
    main()
