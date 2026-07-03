#!/usr/bin/env python
"""
Convert ChromBERT get_region_emb output (hdf5) + its make_dataset table (tsv)
into the .npz embedding artifact eCR_navigator consumes (docs/embedding_artifact.md).

The hdf5 has `region` (N,4 int64: [chrom_id, start, end, build_region_index]) and
`emb` (N,768 float16). chrom_id is an internal index, so we recover the chrom
string by joining on build_region_index against the make_dataset tsv
(chrom, start, end, build_region_index, label).

Runs anywhere with numpy + h5py (no torch) — locally or on the mirror.
"""
from __future__ import annotations

import argparse
import json

import h5py
import numpy as np


def load_tsv_index(path: str) -> dict[int, tuple[str, int, int]]:
    """build_region_index -> (chrom, start, end) from the make_dataset tsv."""
    idx: dict[int, tuple[str, int, int]] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            bri = int(f[col["build_region_index"]])
            idx[bri] = (f[col["chrom"]], int(f[col["start"]]), int(f[col["end"]]))
    return idx


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hdf5", required=True)
    ap.add_argument("--dataset", required=True, help="make_dataset tsv")
    ap.add_argument("--genome", required=True, help="hg38 or mm10 (recorded as assembly)")
    ap.add_argument("--cell-state", required=True)
    ap.add_argument("--out", required=True, help="output .npz artifact")
    args = ap.parse_args()

    tsv = load_tsv_index(args.dataset)
    with h5py.File(args.hdf5, "r") as f:
        region = f["region"][:]           # (N,4) int64
        emb = f["emb"][:].astype(np.float32)  # (N,768)
    build_idx = region[:, 3]

    chrom, start, end = [], [], []
    for bri in build_idx:
        c, s, e = tsv[int(bri)]
        chrom.append(c); start.append(s); end.append(e)

    meta = json.dumps({
        "model": "chrombert",
        "cell_state": args.cell_state,
        "assembly": args.genome,
        "dim": int(emb.shape[1]),
        "source": "chrombert_get_region_emb -> hdf5_to_artifact.py",
    })
    np.savez_compressed(
        args.out,
        chrom=np.array(chrom), start=np.array(start, dtype=np.int64),
        end=np.array(end, dtype=np.int64), embedding=emb, meta=np.array(meta),
    )
    print("wrote %s : %d regions x %d dims (%s, %s)"
          % (args.out, emb.shape[0], emb.shape[1], args.cell_state, args.genome))


if __name__ == "__main__":
    main()
