#!/usr/bin/env python
"""
Build the ChromFound input h5ad for two cell states (e.g. MEF vs mES).

ChromFound is coordinate-based and universe-agnostic: it encodes each open
chromatin region (OCR) by chromosome + sin/cos of its hg38 start/end + a
continuous accessibility scalar. So we feed OUR union peak set directly — no
foreign region universe to snap to. The model is hg38-trained, so a mm10 peak
set is lifted to hg38 here (per-region liftOver via pyliftover); pass a hg38
peak set with --no-lift to skip.

Output h5ad matches ChromFound's schema (src/data/dataset_ds.py):
  var: features, #Chromosome (str "chrN"), hg38_Start (str), hg38_End (str)
  obs: celltype  (one row per state)
  X  : accessibility per state per OCR  (row s = atpm_<state>)

Only chromosomes in ChromFound's vocab (chr1..22, X, Y) are kept.

Usage (chromfound env):
  python chromfound_build_input.py \
      --union union_named.bed --atpm atpm_union.tsv --states MEF mES \
      --chain /yutiancheng/yuhao/eCR/refs/mm10ToHg38.over.chain.gz \
      --out chromfound_input.MEF_mES.h5ad
"""
from __future__ import annotations

import argparse

import anndata as ad
import numpy as np
import pandas as pd

VOCAB_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def lift_region(lo, chrom, start, end):
    """liftOver one interval; return (chrom,start,end) in target or None if it
    doesn't map cleanly (unmapped end, chromosome switch, or inverted)."""
    a = lo.convert_coordinate(chrom, start)
    b = lo.convert_coordinate(chrom, end)
    if not a or not b:
        return None
    c1, s1 = a[0][0], a[0][1]
    c2, e1 = b[0][0], b[0][1]
    if c1 != c2 or c1 not in VOCAB_CHROMS:
        return None
    s, e = (s1, e1) if s1 < e1 else (e1, s1)
    if e <= s:
        return None
    return c1, s, e


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--union", required=True, help="union peaks BED (chrom,start,end[,name])")
    ap.add_argument("--atpm", required=True, help="atpm_union.tsv with atpm_<state> columns")
    ap.add_argument("--states", nargs="+", required=True, help="state names, e.g. MEF mES")
    ap.add_argument("--chain", default=None, help="liftOver chain (source->hg38); omit with --no-lift")
    ap.add_argument("--no-lift", action="store_true", help="peaks are already hg38")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    peaks = pd.read_csv(args.union, sep="\t", header=None,
                        usecols=[0, 1, 2], names=["chrom", "start", "end"])
    atpm = pd.read_csv(args.atpm, sep="\t")
    akey = atpm["chrom"].astype(str) + ":" + atpm["start"].astype(str) + "-" + atpm["end"].astype(str)
    amaps = {}
    for st in args.states:
        col = f"atpm_{st}"
        if col not in atpm.columns:
            raise SystemExit(f"{col} not in {args.atpm} (have {list(atpm.columns)})")
        amaps[st] = dict(zip(akey, atpm[col]))

    # per-region liftOver to hg38 (or pass-through), keeping vocab chromosomes only
    if args.no_lift:
        rows = [(c, int(s), int(e)) for c, s, e in peaks.itertuples(index=False)
                if c in VOCAB_CHROMS]
        keep_src = [(c, int(s), int(e)) for c, s, e in peaks.itertuples(index=False)
                    if c in VOCAB_CHROMS]
    else:
        from pyliftover import LiftOver
        if not args.chain:
            raise SystemExit("--chain required unless --no-lift")
        lo = LiftOver(args.chain)
        rows, keep_src = [], []
        for c, s, e in peaks.itertuples(index=False):
            r = lift_region(lo, c, int(s), int(e))
            if r is not None:
                rows.append(r)               # hg38 coords
                keep_src.append((c, int(s), int(e)))   # original key for aTPM lookup
    if not rows:
        raise SystemExit("no regions survived liftOver to hg38 vocab chromosomes")

    hg = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    var = pd.DataFrame({
        "features": [f"{c}:{s}-{e}" for c, s, e in rows],
        "#Chromosome": hg["chrom"].astype(str),   # vocab lookup is by name string
        "hg38_Start": hg["start"].astype(np.int64),   # must be numeric (torch.tensor)
        "hg38_End": hg["end"].astype(np.int64),
    })
    var.index = var["features"]

    # X: state x OCR accessibility, looked up by the ORIGINAL (pre-lift) region key
    X = np.zeros((len(args.states), len(rows)), dtype=np.float32)
    for i, st in enumerate(args.states):
        m = amaps[st]
        X[i] = [m.get(f"{c}:{s}-{e}", 0.0) for c, s, e in keep_src]

    obs = pd.DataFrame({"celltype": args.states}, index=list(args.states))
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.write_h5ad(args.out)
    print(f"wrote {args.out}: {len(args.states)} states x {len(rows)} OCRs "
          f"(kept {len(rows)}/{len(peaks)} union peaks after lift+vocab filter)")


if __name__ == "__main__":
    main()
