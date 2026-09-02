#!/usr/bin/env python3
"""Per-sample pseudobulk of CALLED cells (pre-QC, pre-ambient) for the cell-calling report.

Reads filt_counts_<sid>.h5 (stacked S+U+A), sums S+U+A per gene, then sums across all called
cells -> one pseudobulk vector per sample. Samples are read in parallel.

  called_cells_pseudobulk.py <in_dir> <out_dir> [n_cores]

Outputs (in <out_dir>):
  pseudobulk_called_cells.csv.gz   genes x samples
  sample_ncells.csv                sample_id, n_called_cells
"""
import glob, re, os, sys
import numpy as np, pandas as pd, h5py
from scipy.sparse import csc_matrix
from multiprocessing import Pool

IN_DIR  = sys.argv[1]
OUT_DIR = sys.argv[2]
N_CORES = int(sys.argv[3]) if len(sys.argv) > 3 else max(1, (os.cpu_count() or 2) - 1)


def one_sample(h5):
    sid = re.sub(r"^filt_counts_", "", os.path.basename(h5)).replace(".h5", "")
    with h5py.File(h5, "r") as f:
        mat = csc_matrix((f["matrix/data"][:], f["matrix/indices"][:], f["matrix/indptr"][:]),
                         shape=tuple(f["matrix/shape"][:]))
        names = f["matrix/features/name"][:].astype(str)
    idx = {t: [i for i, n in enumerate(names) if n.endswith(t)] for t in ("_S", "_U", "_A")}
    genes = list(dict.fromkeys(re.sub(r"_[SUA]$", "", n) for n in names))
    summed = mat[idx["_S"], :] + mat[idx["_U"], :] + mat[idx["_A"], :]   # genes x cells
    vec = np.asarray(summed.sum(axis=1)).ravel().astype(np.int64)
    return sid, np.array(genes), vec, mat.shape[1]


def main():
    files = sorted(glob.glob(f"{IN_DIR}/filt_counts_*.h5"))
    if not files:
        sys.exit(f"ERROR: no filt_counts_*.h5 in {IN_DIR}")
    workers = max(1, min(N_CORES, len(files)))
    print(f"{len(files)} samples, {workers} workers", flush=True)
    with Pool(workers) as pool:
        results = pool.map(one_sample, files)

    gene_ref = None
    cols, ncells = {}, []
    for sid, genes, vec, nc in results:
        if gene_ref is None:
            gene_ref = genes
        elif not np.array_equal(gene_ref, genes):
            sys.exit(f"ERROR: gene lists differ (sample {sid})")
        cols[sid] = vec
        ncells.append({"sample_id": sid, "n_called_cells": int(nc)})
        print(f"  {sid}: {nc} cells", flush=True)

    pb = pd.DataFrame(cols, index=gene_ref)
    pb.index.name = "gene"
    pb.to_csv(f"{OUT_DIR}/pseudobulk_called_cells.csv.gz")
    pd.DataFrame(ncells).sort_values("sample_id").to_csv(f"{OUT_DIR}/sample_ncells.csv", index=False)
    print(f"wrote pseudobulk {pb.shape[0]} genes x {pb.shape[1]} samples", flush=True)


if __name__ == "__main__":
    main()
