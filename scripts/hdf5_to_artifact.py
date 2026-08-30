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

import h5py
import numpy as np

import json

from embedding_artifact import write_embedding_artifact


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
        # dtype BEFORE the cast: ChromBERT stores fp16, and that is the number
        # that explains its A800-vs-DCU divergence.
        emb_dtype = "fp16" if f["emb"].dtype == np.float16 else str(f["emb"].dtype)
        emb = f["emb"][:].astype(np.float32)  # (N,768)
        # This script runs on whatever machine does the conversion, NOT the GPU
        # that computed the embedding, so auto-capture would stamp the wrong
        # host. run_chrombert_region_emb.sh attaches the GPU's own provenance to
        # the hdf5 as an attribute; carry it through unchanged.
        prov = f.attrs.get("ecr_provenance")
    provenance = json.loads(prov) if prov else None
    if provenance is None:
        print("NOTE: hdf5 carries no ecr_provenance attribute (written before "
              "this was added, or produced outside run_chrombert_region_emb.sh); "
              "the GPU that produced it is unrecorded.")
    build_idx = region[:, 3]

    chrom, start, end = [], [], []
    for bri in build_idx:
        c, s, e = tsv[int(bri)]
        chrom.append(c); start.append(s); end.append(e)

    n, d = write_embedding_artifact(
        args.out, chrom, start, end, emb,
        model="chrombert", cell_state=args.cell_state, assembly=args.genome,
        source="chrombert_get_region_emb -> hdf5_to_artifact.py",
        dtype=emb_dtype, attn="eager", provenance=provenance)
    print("wrote %s : %d regions x %d dims (%s, %s)"
          % (args.out, n, d, args.cell_state, args.genome))


if __name__ == "__main__":
    main()
