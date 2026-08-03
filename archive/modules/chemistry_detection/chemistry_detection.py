#!/usr/bin/env python3
"""
Chemistry detection for alevin-fry mapping.

Given a Cell Ranger install and R1 FASTQ files:
  1. Extract barcodes from R1, overlap against every Cell Ranger whitelist.
  2. If the best whitelist maps to a single chemistry, use it directly.
  3. If ambiguous (e.g. 3v2/5v1/5v2 share the same whitelist),
     subsample FASTQs, test-map with fw and rc orientations, pick the winner.

Outputs:
  - whitelist.txt              Selected whitelist (decompressed)
    - Prints CHEMISTRY/ORIENTATION and chemistry stats as KEY=VALUE lines
        for Nextflow env capture.
"""

import argparse
import gzip
import json
import os
import pathlib
import random
import shutil
import subprocess
import sys
import warnings


# ---------------------------------------------------------------------------
# Cell Ranger whitelist mapping (universal across CR versions)
# ---------------------------------------------------------------------------

# Each whitelist filename -> list of chemistries that use it
WHITELIST_TO_CHEMISTRIES = {
    "737K-august-2016.txt":           ["3v2", "5v1", "5v2"],
    "3M-february-2018_TRU.txt.gz":   ["3v3"],
    "3M-5pgex-jan-2023.txt.gz":      ["5v3"],
    "3M-3pgex-may-2023_TRU.txt.gz":  ["3v4"],
    "9K-LT-march-2021.txt.gz":       ["3LT"],
    "737K-arc-v1.txt.gz":            ["multiome"],
}

# Chemistry -> simpleaf chemistry string
CHEMISTRY_TO_AF = {
    "3v2": "10xv2", "5v1": "10xv2", "5v2": "10xv2",
    "3v3": "10xv3", "5v3": "10xv3", "3v4": "10xv3",
    "3LT": "10xv3", "multiome": "10xv3",
}

# Chemistry -> default expected orientation (when unambiguous)
CHEMISTRY_TO_ORI = {
    "3v2": "fw", "5v1": "rc", "5v2": "rc",
    "3v3": "fw", "5v3": "rc", "3v4": "fw",
    "3LT": "fw", "multiome": "fw",
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def detect_chemistry(run, r1_files, cellranger_home, af_index_dir, t2g_f,
                     threads, sample_size=100_000):

    barcode_dir = _resolve_barcode_dir(cellranger_home)

    print(f"[{run}] Checking barcode whitelist overlap ...", file=sys.stderr)
    overlap = _get_whitelist_overlap(r1_files, barcode_dir, sample_size)

    # Find best-matching whitelist
    best_wl = max(overlap, key=lambda x: x["overlap"])
    max_overlap = best_wl["overlap"]
    if max_overlap < 0.7:
        warnings.warn(
            f"[{run}] Maximum whitelist overlap is {max_overlap:.1%}; "
            "chemistry guess may be incorrect"
        )

    # All whitelists that tie for best overlap
    tied = [x for x in overlap if x["overlap"] == max_overlap]
    # All chemistries associated with the best whitelist(s)
    chem_opts = set()
    for t in tied:
        chem_opts.update(t["chemistries"])

    whitelist_path = tied[0]["path"]

    if len(chem_opts) == 1:
        # Unambiguous
        sample_chem    = list(chem_opts)[0]
        af_chemistry   = CHEMISTRY_TO_AF[sample_chem]
        orientation    = CHEMISTRY_TO_ORI[sample_chem]
        n_cells_fw     = ""
        n_cells_rc     = ""

    else:
        # Ambiguous — need orientation test
        print(f"[{run}] Ambiguous chemistry {chem_opts}; testing fw/rc ...", file=sys.stderr)
        af_chemistry = "10xv2" if chem_opts == {"3v2", "5v1", "5v2"} else "10xv3"

        # Prepare a decompressed whitelist for test mapping
        test_wl = _decompress_whitelist(whitelist_path, "test_whitelist.txt")

        tmp_dir = f"tmp_chem_detect_{run}"
        os.makedirs(tmp_dir, exist_ok=True)

        sub_r1, sub_r2 = _subset_fastqs(tmp_dir, r1_files,
                                          _r2_from_r1(r1_files), sample_size)

        for ori in ["fw", "rc"]:
            _run_simpleaf_quant(
                out_dir   = f"{tmp_dir}/{ori}_mapping",
                r1_files  = [sub_r1],
                r2_files  = [sub_r2],
                threads   = threads,
                index_dir = af_index_dir,
                chemistry = af_chemistry,
                ori       = ori,
                t2g_f     = t2g_f,
                wl_f      = test_wl,
            )

        orientation, n_cells_fw, n_cells_rc = _infer_read_orientation(tmp_dir)
        shutil.rmtree(tmp_dir)
        os.remove(test_wl)

        sample_chem = "3v2" if orientation == "fw" else "5v1/5v2"

    # Write the selected whitelist as whitelist.txt (decompressed)
    _decompress_whitelist(whitelist_path, "whitelist.txt")

    chem_stats = {
        "run":                      run,
        "selected_chemistry":       sample_chem,
        "selected_af_chemistry":    af_chemistry,
        "selected_ori":             orientation,
        "selected_whitelist":       pathlib.Path(whitelist_path).name,
        "selected_whitelist_overlap": float(max_overlap),
        "n_cells_fw":               n_cells_fw,
        "n_cells_rc":               n_cells_rc,
    }

    print(f"[{run}] Chemistry: {sample_chem} | AF chemistry: {af_chemistry} | orientation: {orientation}", file=sys.stderr)

    # Print key=value pairs for Nextflow env capture (stdout only).
    print(f"SELECTED_CHEMISTRY={sample_chem}")
    print(f"CHEMISTRY={af_chemistry}")
    print(f"ORIENTATION={orientation}")
    print(f"SELECTED_WHITELIST={pathlib.Path(whitelist_path).name}")
    print(f"SELECTED_WHITELIST_OVERLAP={float(max_overlap)}")
    print(f"N_CELLS_FW={n_cells_fw}")
    print(f"N_CELLS_RC={n_cells_rc}")

    return chem_stats


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_barcode_dir(cellranger_home):
    cr = pathlib.Path(cellranger_home)
    if (cr / "cellranger").exists():
        home = cr
    elif cr.name == "cellranger" and cr.is_file():
        home = cr.parent.parent
    else:
        home = cr

    barcode_dir = home / "lib" / "python" / "cellranger" / "barcodes"
    if not barcode_dir.is_dir():
        raise FileNotFoundError(
            f"Cell Ranger barcode directory not found: {barcode_dir}\n"
            f"Looked under: {home}"
        )
    return barcode_dir


def _decompress_whitelist(src, dst):
    src = str(src)
    if src.endswith(".gz"):
        with gzip.open(src, "rt") as fin, open(dst, "w") as fout:
            fout.write(fin.read())
    else:
        shutil.copy2(src, dst)
    return dst


def _get_whitelist_overlap(r1_files, barcode_dir, sample_size):
    random.seed(1234)
    sel_r1 = random.sample(r1_files, 1)[0]

    print(f"  Extracting barcodes from {sel_r1}", file=sys.stderr)
    cmd = f"seqkit head -n {sample_size} {sel_r1} | seqkit subseq -r 1:16"
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    barcodes = set(_extract_seqs(res.stdout))
    n_bcs = len(barcodes)
    print(f"  Unique barcodes sampled: {n_bcs}", file=sys.stderr)

    results = []
    for wl_file, chemistries in WHITELIST_TO_CHEMISTRIES.items():
        wl_path = barcode_dir / wl_file
        if not wl_path.exists():
            print(f"  Skipping {wl_file} (not found)", file=sys.stderr)
            continue

        opener = gzip.open if str(wl_path).endswith(".gz") else open
        mode = "rt" if str(wl_path).endswith(".gz") else "r"
        with opener(wl_path, mode) as fh:
            wl_set = {line.strip() for line in fh}

        matches = sum(1 for bc in barcodes if bc in wl_set)
        overlap_pct = matches / n_bcs if n_bcs > 0 else 0.0
        results.append({
            "file":        wl_file,
            "path":        wl_path,
            "chemistries": chemistries,
            "overlap":     overlap_pct,
        })
        print(f"  {wl_file} ({','.join(chemistries)}): {overlap_pct:.1%}", file=sys.stderr)

    if not results:
        raise RuntimeError("No whitelist files found in Cell Ranger barcode directory")

    return results


def _r2_from_r1(r1_files):
    r2 = []
    for f in r1_files:
        f2 = f.replace("_R1_", "_R2_")
        if not os.path.exists(f2):
            raise FileNotFoundError(f"Could not find R2 file for {f}. Expected: {f2}")
        r2.append(f2)
    return r2


def _subset_fastqs(out_dir, r1_files, r2_files, sample_size):
    random.seed(12346)
    idx = random.randrange(len(r1_files))
    r1_f, r2_f = r1_files[idx], r2_files[idx]

    sub_r1 = f"{out_dir}/downsampled_{pathlib.Path(r1_f).name}"
    sub_r2 = f"{out_dir}/downsampled_{pathlib.Path(r2_f).name}"
    subprocess.run(["seqkit", "head", "-n", str(sample_size), r1_f, "-o", sub_r1], check=True)
    subprocess.run(["seqkit", "head", "-n", str(sample_size), r2_f, "-o", sub_r2], check=True)
    return sub_r1, sub_r2


def _run_simpleaf_quant(out_dir, r1_files, r2_files, threads, index_dir,
                         chemistry, ori, t2g_f, wl_f):
    env = os.environ.copy()
    env["ALEVIN_FRY_HOME"] = os.path.join(os.getcwd(), "af_home")
    subprocess.run(["simpleaf", "set-paths"], env=env, check=True)
    subprocess.run([
        "simpleaf", "quant",
        "--reads1",       ",".join(r1_files),
        "--reads2",       ",".join(r2_files),
        "--threads",      str(threads),
        "--index",        str(index_dir),
        "--chemistry",    chemistry,
        "--resolution",   "cr-like",
        "--expected-ori", ori,
        "--t2g-map",      str(t2g_f),
        "--unfiltered-pl", str(wl_f),
        "--min-reads",    "1",
        "--output",       out_dir,
    ], env=env, check=True)


def _infer_read_orientation(af_res_dir):
    paths = {
        "fw": os.path.join(af_res_dir, "fw_mapping", "af_quant", "quant.json"),
        "rc": os.path.join(af_res_dir, "rc_mapping", "af_quant", "quant.json"),
    }
    counts = {}
    for ori, p in paths.items():
        with open(p) as fh:
            counts[ori] = json.load(fh).get("num_quantified_cells", 0)

    if counts["fw"] > counts["rc"]:
        return "fw", counts["fw"], counts["rc"]
    elif counts["rc"] > counts["fw"]:
        return "rc", counts["fw"], counts["rc"]
    else:
        raise ValueError("Ambiguous orientation: equal cell counts for fw and rc")


def _extract_seqs(fastq_text):
    seqs = []
    for entry in fastq_text.strip().split("\n@"):
        lines = entry.split("\n")
        if len(lines) > 1:
            seqs.append(lines[1])
    return seqs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Detect 10x chemistry from R1 FASTQ barcode overlap"
    )
    parser.add_argument("run", help="Sample / run ID")
    parser.add_argument("--r1_files", nargs="+", required=True)
    parser.add_argument("--cellranger_home", required=True,
                        help="Cell Ranger install directory")
    parser.add_argument("--af_index_dir", required=True,
                        help="simpleaf index directory (for orientation test)")
    parser.add_argument("--t2g_f", required=True,
                        help="Transcript-to-gene map")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sample_size", type=int, default=100_000)

    args = parser.parse_args()
    detect_chemistry(
        run             = args.run,
        r1_files        = args.r1_files,
        cellranger_home = args.cellranger_home,
        af_index_dir    = args.af_index_dir,
        t2g_f           = args.t2g_f,
        threads         = args.threads,
        sample_size     = args.sample_size,
    )
