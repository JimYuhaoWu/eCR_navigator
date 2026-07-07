#!/usr/bin/env python
"""
Build the per-peak motif-score matrix for GET (RegionMotif prep), genome-parameterized.

GET's input per region is `region_motif` = [282 motif-archetype scores | 1 ATAC] = 283.
This script produces the 282-motif part for an arbitrary peak set, for hg38 OR mm10
(the ONLY species difference is which Vierstra archetype-motif file is queried — so
human needs this exact same step, it's just precomputed in the demo zarr).

Motifs are pre-scanned genome-wide and shipped as a tabix-indexed BED per assembly:
  https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/{assembly}.archetype_motifs.v1.0.bed.gz
We query it REMOTELY (tabix -R) so only blocks overlapping the peaks are fetched — no
9.5 GB (mm10) / 12 GB (hg38) download.

Pipeline (matches get_model preprocess_utils.get_motif):
  tabix <remote motif URL> -R peaks.bed          # motif hits over peaks
  bedtools intersect peaks vs hits, groupby sum  # per (peak, archetype) summed score
  pivot + reindex to the canonical 282 order, fill missing = 0

Runs on the model mirror (needs tabix + bedtools; both in the `get` conda env).

Usage:
  python get_regionmotif_matrix.py --peaks MEF.mm10.bed --assembly mm10 \
      --motif-names motif_names.txt --out MEF.motif.npz
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
import os

import numpy as np
import pandas as pd

MOTIF_BASE = ("https://resources.altius.org/~jvierstra/projects/"
              "motif-clustering/releases/v1.0")


def motif_url(assembly: str) -> str:
    return f"{MOTIF_BASE}/{assembly}.archetype_motifs.v1.0.bed.gz"


def run(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True)


def build_matrix(peaks_bed: str, assembly: str, motif_names: list[str],
                 workdir: str) -> tuple[pd.DataFrame, np.ndarray]:
    """Return (peaks df[chrom,start,end], matrix [n_peaks x 282]) in canonical order."""
    peaks = pd.read_csv(peaks_bed, sep="\t", header=None,
                        usecols=[0, 1, 2], names=["chrom", "start", "end"])
    peaks3 = os.path.join(workdir, "peaks3.bed")
    peaks.to_csv(peaks3, sep="\t", header=False, index=False)

    # 1. remote tabix: motif hits overlapping the peaks
    hits = os.path.join(workdir, "hits.bed")
    run(f'tabix "{motif_url(assembly)}" -R "{peaks3}" > "{hits}"')

    # 2. intersect peaks x hits, keep (peak chrom/start/end, archetype col4, score col5),
    #    sum score per (peak, archetype)
    pm = os.path.join(workdir, "peak_motif.bed")
    run(f'bedtools intersect -a "{peaks3}" -b "{hits}" -wa -wb '
        f'| cut -f1,2,3,7,8 '
        f'| sort -k1,1 -k2,2n -k3,3n -k4,4 '
        f'| bedtools groupby -g 1-4 -c 5 -o sum > "{pm}"')

    if os.path.getsize(pm) == 0:
        mat = np.zeros((len(peaks), len(motif_names)), dtype=np.float32)
        return peaks, mat

    df = pd.read_csv(pm, sep="\t", header=None,
                     names=["chrom", "start", "end", "motif", "score"])
    wide = df.pivot_table(index=["chrom", "start", "end"], columns="motif",
                          values="score", fill_value=0.0)
    # reindex to canonical 282 order (fill motifs absent in this peak set with 0)
    wide = wide.reindex(columns=motif_names, fill_value=0.0)
    # align rows back to the input peak order (peaks with no hits -> all zeros)
    wide = wide.reindex(index=pd.MultiIndex.from_frame(peaks), fill_value=0.0)
    return peaks, wide.to_numpy(dtype=np.float32)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--peaks", required=True, help="peak BED (>=3 col), assembly-matched")
    ap.add_argument("--assembly", required=True, choices=["hg38", "mm10"])
    ap.add_argument("--motif-names", required=True,
                    help="text file, one of the 282 canonical motif names per line")
    ap.add_argument("--out", required=True, help="output .npz (chrom,start,end,motif[NxM],motif_names)")
    args = ap.parse_args()

    motif_names = [l.strip() for l in open(args.motif_names) if l.strip()]
    with tempfile.TemporaryDirectory() as wd:
        peaks, mat = build_matrix(args.peaks, args.assembly, motif_names, wd)
    np.savez_compressed(
        args.out, chrom=peaks["chrom"].to_numpy(),
        start=peaks["start"].to_numpy(), end=peaks["end"].to_numpy(),
        motif=mat, motif_names=np.array(motif_names), assembly=np.array(args.assembly),
    )
    nz = int((mat.sum(axis=1) > 0).sum())
    print(f"wrote {args.out}: {mat.shape[0]} peaks x {mat.shape[1]} motifs "
          f"({nz} peaks with >=1 motif hit), assembly={args.assembly}")


if __name__ == "__main__":
    main()
