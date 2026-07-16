#!/usr/bin/env python
"""
eCR_navigator entrypoint — endpoint embeddings to the eCR_predictor run bundle.

Consumes two mirror-side embedding artifacts (one per cell state, e.g. MEF and mES; see
docs/embedding_artifact.md) and scores each shared region by its embedding shift.

Two input modes:

  --emb-a/--emb-b  two embedding artifacts; scores them here.
  --contract FILE  an ALREADY-scored region-weight contract TSV, skipping the embed step.
                   The GPU mirrors are not persistent (see CLAUDE.md), so their .npz
                   artifacts are routinely gone while the contract they produced survives.
                   A contract carries no metadata, so --assembly and --model are required.

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

    # full bundle from artifacts
    python navigate.py --emb-a get.fib.hg38.npz --emb-b get.iN.hg38.npz \
        --bundle bundles/in_gse299923 \
        --matrix endpoints.matrix.tsv --state-a fib --state-b iN \
        --anchors neural.promoter.bed --signed signed_delta.tsv \
        --transition transition.json

    # full bundle from an archived contract (no GPU needed)
    python navigate.py --contract get_in.tsv --assembly hg38 --model GET \
        --bundle bundles/in_gse299923 \
        --matrix endpoints.matrix.tsv --state-a fib --state-b iN \
        --anchors neural.promoter.bed --signed signed_delta.tsv
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
from ecr_navigator.weights import read_region_weights, write_region_weights

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
    ap.add_argument("--emb-a", help="cell-state A artifact (.npz), e.g. MEF")
    ap.add_argument("--emb-b", help="cell-state B artifact (.npz), e.g. mES")
    ap.add_argument("--contract", help="an already-scored contract TSV, instead of "
                                       "--emb-a/--emb-b (mirror artifacts are ephemeral)")
    ap.add_argument("--assembly", help="assembly for --contract mode (artifacts self-report)")
    ap.add_argument("--model", help="model name for --contract mode (artifacts self-report)")
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
    if bool(args.contract) == bool(args.emb_a or args.emb_b):
        raise SystemExit("pass EITHER --emb-a/--emb-b OR --contract")

    if args.contract:
        rows, meta = _from_contract(args)
    else:
        rows, meta = _from_artifacts(args)
    n_dir = sum(1 for r in rows if r.direction is not None)

    if args.out:
        write_region_weights(rows, args.out)
        print("wrote %s : %d regions (%s vs %s, norm=%s, direction=%d, unmeasured=%d)"
              % (args.out, len(rows), meta["state_a"], meta["state_b"],
                 args.norm, n_dir, len(rows) - n_dir))

    if args.bundle:
        _write_bundle(args, rows, meta, n_dir)


def _from_artifacts(args):
    """Score two endpoint artifacts. They self-report assembly/model/cell states, so the
    manifest metadata cannot drift from the run that produced the scores."""
    if not (args.emb_a and args.emb_b):
        raise SystemExit("--emb-a and --emb-b must be given together")
    a = load_artifact(args.emb_a)
    b = load_artifact(args.emb_b)
    asm_a, asm_b = a.meta.get("assembly"), b.meta.get("assembly")
    if asm_a != asm_b:
        # the contract calls assembly load-bearing: mismatched coordinates silently
        # corrupt both the Tier-2 weighting and downstream site selection.
        raise SystemExit(f"assembly mismatch: {args.emb_a} is {asm_a}, {args.emb_b} is {asm_b}")

    chrom, start, end, shift = embedding_shift(a, b)
    rows = driver_scores(chrom, start, end, shift, method=args.norm)

    if args.direction != "off":
        delta = signed_delta(a, b)
        if delta is None:
            if args.direction == "on":
                raise SystemExit("--direction on: an artifact has no scalar signal")
        else:
            attach_direction(rows, delta, method=args.direction_norm)

    return rows, {"assembly": asm_a, "model": a.meta.get("model"),
                  "state_a": a.meta.get("cell_state"), "state_b": b.meta.get("cell_state"),
                  "source": a.meta.get("source"),
                  "direction_norm": args.direction_norm}


def _from_contract(args):
    """Load an already-scored contract TSV. Its driver_score is taken as-is (it was
    normalized when it was produced -- declare how with --norm), and its `direction` column
    is the measured signed-delta, already normalized. A contract carries no metadata, so
    assembly and model must be supplied."""
    if not (args.assembly and args.model):
        raise SystemExit("--contract requires --assembly and --model (a contract TSV "
                         "carries no metadata; artifacts self-report theirs)")
    rows = read_region_weights(args.contract)
    return rows, {"assembly": args.assembly, "model": args.model,
                  "state_a": args.state_a, "state_b": args.state_b,
                  "source": f"contract {os.path.basename(args.contract)}",
                  "direction_norm": None}


def _gates(args, rows):
    """Run the two nomination gates. Either may be None (its inputs were not supplied) --
    nominate() turns that into a refusal, which is a valid bundle."""
    from preflight import (GATE1_UNIVERSE_N, admissibility,  # noqa: E402
                           gate2_inputs, load_matrix, select_score, top_k_fold)

    gate1 = None
    if args.matrix:
        if not (args.state_a and args.state_b):
            raise SystemExit("--matrix requires --state-a and --state-b")
        adm = admissibility(*load_matrix(args.matrix, args.state_a, args.state_b))
        gate1 = {"admit": bool(adm.admit), "pc1_frac": round(adm.pc1_frac, 4),
                 "coherence_margin": round(adm.coherence_margin, 4),
                 "n_a": adm.n_a, "n_b": adm.n_b,
                 "universe_n": adm.universe_n,
                 "universe_truncated": adm.universe_truncated,
                 "reasons": adm.reasons}
        if adm.universe_truncated:
            print(f"WARNING: Gate-1 universe is only {adm.universe_n} regions (< the "
                  f"standard {GATE1_UNIVERSE_N}) — PC1/coherence are not comparable to "
                  f"other transitions")

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


def _write_bundle(args, rows, meta, n_dir):
    gate1, gate2 = _gates(args, rows)
    noms, block = nominate(rows, gate1, gate2, top_frac=args.top_frac,
                           score_norm=args.norm)

    assembly, model = meta["assembly"], meta["model"]
    transition = {"id": os.path.basename(str(args.bundle).rstrip("/\\")),
                  "species": SPECIES.get(assembly, "unknown"), "assembly": assembly,
                  "state_a": meta["state_a"], "state_b": meta["state_b"]}
    if args.transition:
        with open(args.transition) as fh:
            transition.update(json.load(fh))

    manifest = {
        "run_id": transition["id"],
        "transition": transition,
        "region_universe": {"n_regions": len(rows), "source": meta["source"]},
        "gate1": gate1,
        "gate2": gate2,
        "nomination": block,
        "direction": {
            "provenance": DIRECTION_PROVENANCE.get(model) if n_dir else None,
            "source": (f"{model} signal ({meta['state_b']} - {meta['state_a']})"
                       if n_dir else None),
            "norm": meta["direction_norm"] if n_dir else None,
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
