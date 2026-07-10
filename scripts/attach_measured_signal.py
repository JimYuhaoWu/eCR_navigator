#!/usr/bin/env python
"""
Attach a MEASURED per-region accessibility `signal` to an existing embedding
artifact — the model-independent way to give the contract's `direction` column a
value for models that have NO accessibility-prediction head of their own
(ATACformer, GET, ChromBERT, ChromFound: their driver_score is an *unsigned*
embedding-shift magnitude).

Unlike EpiAgent — whose `signal` comes from its own Signal-Reconstruction head, so
its direction is model-native — the signal added here is an EXTERNAL measurement
(e.g. GET's aTPM accessibility channel, or a MACS peak score). The model contributes
only the magnitude; the direction is measured data bolted on. That is the *trusted*
end of the direction caveat (see docs/region_weight_contract.md) — real data, not a
value synthesized from the embedding — but it is still not the model's own output, so
label it as measured in the pipeline doc.

Given an artifact and a per-region intensity table (BED/TSV, same assembly) for the
SAME cell state, it maps each artifact region to the max-overlapping intensity value
and rewrites the .npz with a `signal` array. navigate.py then differences the two
states' signals into `direction`. Do this once per state artifact, with that state's
intensity column, then run navigate.py as usual.

  python attach_measured_signal.py --artifact atac.kidney.hg38.npz \
      --intensity get_human_kp/atpm_union.tsv --value-col atpm_kidney \
      --out atac.kidney.hg38.sig.npz

numpy-only; runs in the light navigator env (no torch). Reuses the shared writer so
the artifact dtypes stay contract-correct.
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from embedding_artifact import write_embedding_artifact


def read_intensity(path, value_col, has_header):
    """(chrom,start,end,value) columns -> {chrom: (starts, ends, vals)} sorted by start."""
    chroms, starts, ends, vals = [], [], [], []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t") if has_header else None
        if has_header:
            try:
                vi = int(value_col)
            except ValueError:
                if value_col not in header:
                    raise SystemExit(f"--value-col {value_col!r} not in header {header}")
                vi = header.index(value_col)
        else:
            vi = int(value_col)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chroms.append(f[0]); starts.append(int(f[1])); ends.append(int(f[2]))
            vals.append(float(f[vi]))
    chroms = np.array(chroms); starts = np.array(starts, np.int64)
    ends = np.array(ends, np.int64); vals = np.array(vals, np.float64)
    by_chrom = {}
    for c in np.unique(chroms):
        m = chroms == c
        order = np.argsort(starts[m])
        by_chrom[c] = (starts[m][order], ends[m][order], vals[m][order])
    return by_chrom


def map_signal(chrom, start, end, by_chrom):
    """Per artifact region: max intensity value over all overlapping intervals.

    A region with NO overlapping interval is left as NaN ("unmeasured"), NOT 0.0 —
    a measured zero and a missing measurement are different, and collapsing them makes
    a region covered in one state but not the other diff to a spurious open/close
    direction. navigate.py leaves `direction` unset wherever a state's signal is NaN.
    """
    # bound the backward scan by the longest interval so no overlap is missed
    max_len = max((int((e - s).max()) for s, e, _ in by_chrom.values() if len(s)), default=0)
    out = np.full(len(chrom), np.nan, dtype=np.float64)
    for i in range(len(chrom)):
        c = chrom[i]
        rec = by_chrom.get(c)
        if rec is None:
            continue
        s_arr, e_arr, v_arr = rec
        qs, qe = start[i], end[i]
        # candidate intervals have start in [qs - max_len, qe); among those, keep end > qs
        lo = np.searchsorted(s_arr, qs - max_len, side="left")
        hi = np.searchsorted(s_arr, qe, side="left")
        if hi <= lo:
            continue
        ov = e_arr[lo:hi] > qs
        if ov.any():
            out[i] = v_arr[lo:hi][ov].max()
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--artifact", required=True, help="embedding artifact .npz (one state)")
    ap.add_argument("--intensity", required=True, help="per-region intensity BED/TSV, SAME assembly")
    ap.add_argument("--value-col", required=True, help="intensity column: header name or 0-based index")
    ap.add_argument("--no-header", action="store_true", help="intensity file has no header row")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    z = np.load(args.artifact, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    chrom, start, end = z["chrom"], z["start"], z["end"]

    by_chrom = read_intensity(args.intensity, args.value_col, not args.no_header)
    signal = map_signal(chrom, start, end, by_chrom)
    covered = int(np.isfinite(signal).sum())   # NaN = unmeasured, excluded from direction

    n, d = write_embedding_artifact(
        args.out, chrom, start, end, z["embedding"],
        model=meta["model"], cell_state=meta["cell_state"], assembly=meta["assembly"],
        source=meta.get("source", "") + "+attach_measured_signal.py", signal=signal)
    print(f"wrote {args.out}: {n} regions x {d} dims; signal from {args.value_col} "
          f"on {covered}/{n} regions ({100*covered/n:.1f}% overlapped an intensity interval)")


if __name__ == "__main__":
    main()
