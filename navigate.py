#!/usr/bin/env python
"""
eCR_navigator entrypoint — zero-shot driver scoring from ChromBERT embeddings.

Consumes two mirror-side embedding artifacts (one per cell state, e.g. MEF and
mES; see docs/embedding_artifact.md), scores each shared region by its TRN
embedding shift, and writes the region-weight contract TSV that eCR_predictor
consumes (docs/region_weight_contract.md).

    python navigate.py \
        --emb-a chrombert.MEF.mm10.npz \
        --emb-b chrombert.mES.mm10.npz \
        --out driver_weights.mm10.tsv
"""
from __future__ import annotations

import argparse

from ecr_navigator.features import embedding_shift, load_artifact
from ecr_navigator.model import driver_scores
from ecr_navigator.weights import write_region_weights


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--emb-a", required=True, help="cell-state A artifact (.npz), e.g. MEF")
    ap.add_argument("--emb-b", required=True, help="cell-state B artifact (.npz), e.g. mES")
    ap.add_argument("--out", required=True, help="output region-weight contract TSV")
    ap.add_argument("--norm", default="rank", choices=["rank", "minmax"],
                    help="shift -> driver_score normalization (default: rank)")
    args = ap.parse_args()

    a = load_artifact(args.emb_a)
    b = load_artifact(args.emb_b)
    chrom, start, end, shift = embedding_shift(a, b)
    rows = driver_scores(chrom, start, end, shift, method=args.norm)
    write_region_weights(rows, args.out)
    print("wrote %s : %d regions (%s vs %s, norm=%s)"
          % (args.out, len(rows), a.meta.get("cell_state"),
             b.meta.get("cell_state"), args.norm))


if __name__ == "__main__":
    main()
