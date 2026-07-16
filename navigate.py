#!/usr/bin/env python
"""
eCR_navigator entrypoint — endpoint embeddings to the eCR_predictor run bundle.

Consumes two mirror-side embedding artifacts (one per cell state, e.g. MEF and mES; see
docs/embedding_artifact.md) and scores each shared region by its embedding shift.

Two output modes:

  --out FILE      the region-weight contract TSV alone (docs/region_weight_contract.md).

  --bundle DIR    the full run bundle (docs/run_bundle_contract.md): weights.tsv +
                  nominations.tsv + manifest.json. This runs the nomination policy --
                  Gate 1 (admissibility, needs --matrix) and Gate 2 (which score to trust,
                  needs --anchors/--signed) -- and records the verdict in the manifest, so
                  eCR_predictor never has to know which model won.

    # scores only (unchanged)
    python navigate.py --emb-a chrombert.MEF.mm10.npz --emb-b chrombert.mES.mm10.npz \
        --out driver_weights.mm10.tsv

    # full bundle
    python navigate.py --emb-a get.fib.hg38.npz --emb-b get.iN.hg38.npz \
        --bundle bundles/in_gse299923 \
        --matrix endpoints.matrix.tsv --state-a fib --state-b iN \
        --anchors neural.promoter.bed --signed signed_delta.tsv \
        --transition transition.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

from ecr_navigator.features import embedding_shift, load_artifact, signed_delta
from ecr_navigator.model import attach_direction, driver_scores
from ecr_navigator.nominate import DEFAULT_TOP_FRAC, nominate, write_bundle
from ecr_navigator.weights import write_region_weights

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts"))

# assembly -> species. Multi-species is first-class (CLAUDE.md): never hardcode mouse.
SPECIES = {"hg38": "human", "mm10": "mouse"}

# How each model's direction signal is sourced -- governs how much to trust the sign.
# See docs/direction.md; the three tiers are a property of the model, not of the run.
DIRECTION_PROVENANCE = {
    "GET": "input-measured", "ChromFound": "input-measured",
    "EpiAgent": "predicted-model-native",
    "ATACformer": "external-attach", "ChromBERT": "external-attach",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--emb-a", required=True, help="cell-state A artifact (.npz), e.g. MEF")
    ap.add_argument("--emb-b", required=True, help="cell-state B artifact (.npz), e.g. mES")
    ap.add_argument("--out", help="write the region-weight contract TSV here")
    ap.add_argument("--bundle", help="write the full run bundle to this directory")
    ap.add_argument("--norm", default="rank", choices=["rank", "minmax"],
                    help="shift -> driver_score normalization (default: rank)")
    ap.add_argument("--direction", choices=["off", "auto", "on"], default="auto",
                    help="add signed direction [-1,1] from the artifacts' scalar "
                         "signal: 'auto' adds it when both carry one (default), "
                         "'on' requires it, 'off' skips it")
    ap.add_argument("--direction-norm", default="maxabs",
                    choices=["maxabs", "raw", "signed-rank"],
                    help="signal-delta -> direction magnitude scaling (sign kept)")
    # --- bundle mode: the nomination policy ---
    ap.add_argument("--matrix", help="per-replicate accessibility matrix TSV (Gate 1)")
    ap.add_argument("--state-a", help="substring identifying start-state matrix columns")
    ap.add_argument("--state-b", help="substring identifying end-state matrix columns")
    ap.add_argument("--anchors", help="target-cell master-TF loci BED (Gate 2 positives)")
    ap.add_argument("--signed", help="signed-Delta track TSV (Gate 2 baseline scorer)")
    ap.add_argument("--all-regions", action="store_true",
                    help="Gate 2: do NOT restrict to opening regions (default opening-only)")
    ap.add_argument("--top-frac", type=float, default=DEFAULT_TOP_FRAC,
                    help="nomination confidence band (default: %(default)s)")
    ap.add_argument("--transition", help="JSON with transition metadata (label, geo, ...) "
                                         "merged into the manifest's transition block")
    args = ap.parse_args()

    if not args.out and not args.bundle:
        raise SystemExit("nothing to write: pass --out and/or --bundle")

    a = load_artifact(args.emb_a)
    b = load_artifact(args.emb_b)
    asm_a, asm_b = a.meta.get("assembly"), b.meta.get("assembly")
    if asm_a != asm_b:
        # the contract calls assembly load-bearing: mismatched coordinates silently
        # corrupt both the Tier-2 weighting and downstream site selection.
        raise SystemExit(f"assembly mismatch: {args.emb_a} is {asm_a}, {args.emb_b} is {asm_b}")

    chrom, start, end, shift = embedding_shift(a, b)
    rows = driver_scores(chrom, start, end, shift, method=args.norm)

    n_dir = 0
    if args.direction != "off":
        delta = signed_delta(a, b)
        if delta is None:
            if args.direction == "on":
                raise SystemExit("--direction on: an artifact has no scalar signal")
        else:
            attach_direction(rows, delta, method=args.direction_norm)
            n_dir = sum(1 for r in rows if r.direction is not None)

    if args.out:
        write_region_weights(rows, args.out)
        print("wrote %s : %d regions (%s vs %s, norm=%s, direction=%d, unmeasured=%d)"
              % (args.out, len(rows), a.meta.get("cell_state"), b.meta.get("cell_state"),
                 args.norm, n_dir, len(rows) - n_dir))

    if args.bundle:
        _write_bundle(args, a, b, rows, asm_a, n_dir)


def _gates(args, rows):
    """Run the two nomination gates. Either may be None (its inputs were not supplied) --
    nominate() turns that into a refusal, which is a valid bundle."""
    from preflight import (admissibility, gate2_inputs, load_matrix,  # noqa: E402
                           select_score, top_k_fold)

    gate1 = None
    if args.matrix:
        if not (args.state_a and args.state_b):
            raise SystemExit("--matrix requires --state-a and --state-b")
        adm = admissibility(*load_matrix(args.matrix, args.state_a, args.state_b))
        gate1 = {"admit": bool(adm.admit), "pc1_frac": round(adm.pc1_frac, 4),
                 "coherence_margin": round(adm.coherence_margin, 4),
                 "n_a": adm.n_a, "n_b": adm.n_b, "reasons": adm.reasons}

    gate2 = None
    if args.anchors and args.signed:
        chrom = np.array([r.chrom for r in rows])
        start = np.array([r.start for r in rows])
        end = np.array([r.end for r in rows])
        score = np.array([r.driver_score for r in rows])
        labels, signed = gate2_inputs(chrom, start, end, args.anchors, args.signed)
        sel = select_score(score, signed, labels, opening_only=not args.all_regions)
        gate2 = {"primary": sel.primary,
                 "anchor_set": os.path.basename(args.anchors),
                 "n_anchors": int(labels.sum()),
                 "delta_auroc": round(sel.delta_auroc, 4),
                 "delta_auroc_ci": [round(sel.delta_ci[0], 4), round(sel.delta_ci[1], 4)],
                 "lr_p": sel.perm_p, "lr_coef": round(sel.driver_coef, 4),
                 "fold_enrichment_top5pct": round(top_k_fold(score, labels), 4),
                 "reason": sel.reason}
        if labels.sum() < 20:
            print(f"WARNING: only {int(labels.sum())} anchors in-universe (<20 is "
                  f"underpowered) — the Gate-2 verdict is weakly supported")
    elif args.anchors or args.signed:
        raise SystemExit("Gate 2 needs BOTH --anchors and --signed")

    return gate1, gate2


def _write_bundle(args, a, b, rows, assembly, n_dir):
    gate1, gate2 = _gates(args, rows)
    noms, block = nominate(rows, gate1, gate2, top_frac=args.top_frac,
                           score_norm=args.norm)

    model = a.meta.get("model")
    transition = {"id": os.path.basename(str(args.bundle).rstrip("/\\")),
                  "species": SPECIES.get(assembly, "unknown"), "assembly": assembly,
                  "state_a": a.meta.get("cell_state"), "state_b": b.meta.get("cell_state")}
    if args.transition:
        with open(args.transition) as fh:
            transition.update(json.load(fh))

    manifest = {
        "run_id": transition["id"],
        "transition": transition,
        "region_universe": {"n_regions": len(rows), "source": a.meta.get("source")},
        "gate1": gate1,
        "gate2": gate2,
        "nomination": block,
        "direction": {
            "provenance": DIRECTION_PROVENANCE.get(model),
            "source": f"{model} signal ({b.meta.get('cell_state')} - "
                      f"{a.meta.get('cell_state')})" if n_dir else None,
            "norm": args.direction_norm if n_dir else None,
        },
        "models_run": [model],
    }
    write_bundle(args.bundle, manifest, rows, noms)

    verdict = (f"REFUSED — {block['refusal_reason']}" if block["refused"]
               else f"{block['n_nominated']} nominations from {block['score_source']} "
                    f"(top {args.top_frac:.0%})")
    print(f"wrote bundle {args.bundle} : {len(rows)} regions, {assembly}\n  {verdict}")


if __name__ == "__main__":
    main()
