#!/usr/bin/env python
"""
Build the EpiAgent input (a tokenized "cell sentence") for ONE cell state from a
peak BED. Runs on the EpiAgent env (needs epiagent + pyliftover + bedtools).

EpiAgent is a single-cell, hg38-only model, so two bridges apply to our bulk
peaks:
  * pseudobulk: the whole state is one "cell" whose accessible cCREs are those its
    peaks overlap (barcode = the state name).
  * species: a mouse (mm10) peak set is lifted to hg38 first (EpiAgent's cCRE
    vocabulary is hg38); pass --chain. A human (hg38) peak set skips this with
    --no-lift. The unmapped fraction is reported (cross-species loss).

liftOver uses pyliftover (pure Python) because the UCSC liftOver binary needs a
newer glibc than the model-zoo mirror provides.

Pipeline: [pyliftover mm10->hg38] -> tag barcode -> bedtools intersect vs cCRE.bed
  -> construct_cell_by_ccre_matrix -> global_TFIDF -> tokenization -> .h5ad
The .h5ad (one pseudobulk cell x 1,355,445 cCREs, with obs['cell_sentences'] and
cCRE coordinates as var_names) feeds epiagent_embed_regions.py.

Usage:
  # mouse (reference track): lift mm10 -> hg38
  python epiagent_build_input.py --peaks MEF.e7_peaks.bed --cell-state MEF \
      --ccre-bed refs/cCRE.bed --ccre-freq refs/cCRE_document_frequency.npy \
      --chain refs/mm10ToHg38.over.chain.gz --out epiagent.MEF.h5ad
  # human (primary track): peaks already hg38
  python epiagent_build_input.py --peaks HS_state.bed --cell-state STATE \
      --ccre-bed refs/cCRE.bed --ccre-freq refs/cCRE_document_frequency.npy \
      --no-lift --out epiagent.STATE.h5ad
"""
from __future__ import annotations

import argparse
import os
import subprocess
import tempfile

import numpy as np
import scanpy as sc  # noqa: F401  (AnnData IO parity with EpiAgent demos)
from epiagent.preprocessing import construct_cell_by_ccre_matrix, global_TFIDF
from epiagent.tokenization import tokenization


def write_frags(peaks: str, barcode: str, out_bed: str, chain: str | None) -> tuple[int, int]:
    """Write a (chrom, start, end, barcode) BED in hg38. With --chain, lift each
    interval mm10->hg38 (both endpoints must map to the same chrom, start<end);
    without, pass peaks through as-is. Returns (n_in, n_kept)."""
    lo = None
    if chain is not None:
        from pyliftover import LiftOver
        lo = LiftOver(chain)
    n_in = n_out = 0
    with open(peaks) as fin, open(out_bed, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            n_in += 1
            if lo is None:
                fout.write(f"{chrom}\t{start}\t{end}\t{barcode}\n")
                n_out += 1
                continue
            a = lo.convert_coordinate(chrom, start)
            b = lo.convert_coordinate(chrom, end)
            if not a or not b:
                continue
            ca, pa = a[0][0], a[0][1]
            cb, pb = b[0][0], b[0][1]
            if ca != cb or pa >= pb:
                continue
            fout.write(f"{ca}\t{pa}\t{pb}\t{barcode}\n")
            n_out += 1
    return n_in, n_out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--peaks", required=True, help="one state's peak BED")
    ap.add_argument("--cell-state", required=True, help="e.g. MEF or mES (the pseudobulk barcode)")
    ap.add_argument("--ccre-bed", required=True, help="EpiAgent hg38 cCRE.bed (1,355,445 rows)")
    ap.add_argument("--ccre-freq", required=True, help="cCRE_document_frequency.npy")
    ap.add_argument("--chain", default=None, help="mm10ToHg38.over.chain.gz (mouse); omit with --no-lift")
    ap.add_argument("--no-lift", action="store_true", help="peaks are already hg38 (human) — skip liftOver")
    ap.add_argument("--out", required=True, help="output tokenized .h5ad")
    args = ap.parse_args()

    if not args.no_lift and args.chain is None:
        raise SystemExit("mouse input needs --chain (mm10ToHg38); for hg38 peaks pass --no-lift")

    tmp = tempfile.mkdtemp(prefix=f"epi_{args.cell_state}_")
    frags = os.path.join(tmp, "frags.hg38.bed")
    inter = os.path.join(tmp, "intersect.bed")

    n_in, n_keep = write_frags(args.peaks, args.cell_state, frags,
                               None if args.no_lift else args.chain)
    if args.no_lift:
        print(f">> {args.cell_state}: {n_keep} hg38 peaks (no liftOver)")
    else:
        print(f">> liftOver mm10->hg38: {n_keep}/{n_in} peaks mapped "
              f"({100.0 * n_keep / max(n_in, 1):.1f}%), {n_in - n_keep} unmapped")

    subprocess.run(
        f"bedtools intersect -a {frags} -b {args.ccre_bed} -wa -wb > {inter}",
        shell=True, check=True,
    )

    adata = construct_cell_by_ccre_matrix(inter, args.ccre_bed)
    n_ccre = int((adata.X > 0).sum())
    print(f">> pseudobulk cell '{args.cell_state}': {adata.n_obs} cell x {adata.n_vars} cCREs; "
          f"{n_ccre} cCREs accessible")

    freq = np.load(args.ccre_freq)
    adata = global_TFIDF(adata, freq)
    tokenization(adata)
    adata.write(args.out)
    print(f">> wrote {args.out} ({n_ccre} accessible cCREs; embed keeps top 8190 by TF-IDF)")


if __name__ == "__main__":
    main()
