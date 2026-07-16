#!/usr/bin/env python
"""Claim 2B: is the MEASURED `direction` trustworthy where it picks an effector domain?

`direction` (signed accessibility Delta) chooses the ED for a nominated region: + -> activator
(VP64), - -> repressor (KRAB). On the real iN bundle 42% of nominations have |direction|<0.05,
so the sign is weakly supported for a large fraction. This asks whether the sign is CORRECT on
regions whose direction is known from biology, and at what |direction| it stops beating chance.
See docs/claim2_plan.md 2B.

Non-circular via known biology, two-sided:
  DESTINATION master-TF loci  -> must OPEN  (expected sign +1)   [--dest]
  SOURCE-cell master-TF loci  -> must CLOSE (expected sign -1)   [--source]
Two-sided is the whole design: a one-sided "do destination anchors open?" is beaten by a
trivial everything-opens predictor (a transition can be >50% opening). The source set has no
such escape, and BALANCED accuracy across both sets scores that trivial predictor at 0.5.

Metrics (all vs the marginal opening rate, never 50%):
  - sign-accuracy per anchor set, gene-clustered bootstrap CI (anchors of one gene are not
    independent -- promoter and neighborhood overlap, and a gene contributes many regions);
  - balanced accuracy across both sets (the headline; immune to the opening-rate baseline);
  - accuracy STRATIFIED by |direction| -> the |direction| where accuracy meets the base rate
    is the empirical ambiguity threshold (read off, not picked).

Pure numpy; reuses eval_driver_claim1's loaders. Endpoint-only: uses only the two states'
measured accessibility, already in the contract.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_driver_claim1 import load_contract  # noqa: E402


def load_bed_genes(path):
    """BED -> (chrom, start, end, gene). gene = col 4 if present, else the row's own id, so
    every anchor clusters as its own unit when the BED has no gene column."""
    chrom, start, end, gene = [], [], [], []
    with open(path) as fh:
        for i, ln in enumerate(fh):
            ln = ln.rstrip("\n")
            if not ln or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.split("\t")
            chrom.append(f[0]); start.append(int(f[1])); end.append(int(f[2]))
            gene.append(f[3] if len(f) > 3 else f"row{i}")
    return np.array(chrom), np.array(start), np.array(end), np.array(gene)


def assign_gene(chrom, start, end, bc, bs, be, bgene):
    """For each query region, the gene of an overlapping BED interval (or None). Half-open
    overlap. If several overlap, the first by start wins -- clustering only needs a stable
    label, and anchor windows of different genes rarely overlap.

    Direct per-chrom scan rather than a start-sorted sweep: anchor BEDs are tiny (tens of
    windows), and a sweep that stops at the first end<=start would MISS a longer earlier
    interval when widths vary (neighborhood = gene +/-50kb) -- the same max-end pitfall
    eval_driver_claim1.overlap_labels avoids. A direct scan is obviously correct here."""
    out = np.array([None] * len(chrom), dtype=object)
    start = np.asarray(start, np.int64); end = np.asarray(end, np.int64)
    for c in np.unique(chrom):
        bi = np.where(bc == c)[0]
        if len(bi) == 0:
            continue
        o = np.argsort(bs[bi]); bi = bi[o]        # by start, so the first overlap is stable
        cs, ce, cg = bs[bi].astype(np.int64), be[bi].astype(np.int64), bgene[bi]
        for q in np.where(chrom == c)[0]:
            hit = np.where((cs < end[q]) & (ce > start[q]))[0]   # half-open overlap
            if len(hit):
                out[q] = cg[hit[0]]
    return out


def sign_accuracy(direction, expected_sign):
    """Fraction of regions whose measured sign matches `expected_sign` (+1 opens / -1 closes).
    Ignores exact 0.0 (measured-flat: no open/close call) -- reported separately as `n_flat`."""
    nz = direction != 0.0
    if nz.sum() == 0:
        return float("nan"), 0
    correct = np.sign(direction[nz]) == expected_sign
    return float(correct.mean()), int(nz.sum())


def gene_clustered_ci(direction, expected_sign, genes, n_boot=2000, seed=0):
    """Bootstrap sign-accuracy by RESAMPLING GENES (not regions), so overlapping/duplicated
    anchors of one gene count once. Returns (mean, lo, hi) at 95%."""
    rng = np.random.default_rng(seed)
    ug = np.unique(genes)
    if len(ug) < 2:
        acc, _ = sign_accuracy(direction, expected_sign)
        return acc, float("nan"), float("nan")
    by = {g: np.where(genes == g)[0] for g in ug}
    accs = []
    for _ in range(n_boot):
        pick = rng.choice(ug, size=len(ug), replace=True)
        idx = np.concatenate([by[g] for g in pick])
        a, n = sign_accuracy(direction[idx], expected_sign)
        if n > 0:
            accs.append(a)
    return (float(np.mean(accs)), float(np.percentile(accs, 2.5)),
            float(np.percentile(accs, 97.5)))


def stratify_by_magnitude(direction, expected_sign, edges):
    """Sign-accuracy within |direction| bins. Returns rows (lo, hi, n, accuracy). The bin
    where accuracy falls to the base rate is the ambiguity threshold."""
    mag = np.abs(direction)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (mag >= lo) & (mag < hi) & (direction != 0.0)
        if m.sum() == 0:
            rows.append((lo, hi, 0, float("nan"))); continue
        acc = float((np.sign(direction[m]) == expected_sign).mean())
        rows.append((lo, hi, int(m.sum()), acc))
    return rows


def evaluate(chrom, start, end, direction, dest_bed, source_bed, mag_edges=None):
    """Full 2B readout for one transition. `direction` is NaN where unmeasured; those regions
    are dropped (no sign to check)."""
    meas = np.isfinite(direction)
    # base rate over regions with an actual open/close call (exclude measured-flat 0.0), so it
    # shares sign_accuracy's denominator; it is the accuracy a trivial all-opens predictor hits.
    called = meas & (direction != 0.0)
    opening_rate = float((direction[called] > 0).mean()) if called.any() else float("nan")

    out = {"n_measured": int(meas.sum()), "opening_rate": opening_rate, "sets": {}}
    for name, bed, exp in (("destination", dest_bed, +1), ("source", source_bed, -1)):
        bc, bs, be, bg = load_bed_genes(bed)
        g = assign_gene(chrom, start, end, bc, bs, be, bg)
        sel = meas & (g != None)                           # noqa: E711 (object array)
        d, genes = direction[sel], g[sel]
        acc, n = sign_accuracy(d, exp)
        m, lo, hi = gene_clustered_ci(d, exp, genes)
        base = opening_rate if exp == +1 else (1.0 - opening_rate)
        out["sets"][name] = {
            "expected_sign": exp, "n_regions": int(sel.sum()), "n_signed": n,
            "n_genes": int(len(np.unique(genes))) if sel.sum() else 0,
            "accuracy": acc, "boot_mean": m, "ci_lo": lo, "ci_hi": hi,
            "base_rate": base, "beats_base": (lo == lo) and lo > base,   # lo==lo: not NaN
            "strata": stratify_by_magnitude(
                d, exp, mag_edges if mag_edges is not None
                else np.array([0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.01])),
        }
    # balanced accuracy is the headline BECAUSE it is immune to the opening-rate baseline --
    # but only if BOTH sets are populated. If one has no in-universe anchors, a mean would
    # silently report the other side alone, which is exactly the one-sided number the balanced
    # metric exists to avoid; report NaN instead.
    dest, src = out["sets"]["destination"]["accuracy"], out["sets"]["source"]["accuracy"]
    out["balanced_accuracy"] = (float(np.mean([dest, src]))
                                if np.isfinite(dest) and np.isfinite(src) else float("nan"))
    return out


def _fmt(out, title):
    L = ["=" * 70, title, "=" * 70,
         f"measured regions: {out['n_measured']}   opening rate (base): "
         f"{out['opening_rate']:.3f}",
         f"BALANCED accuracy (both sets): {out['balanced_accuracy']:.3f}   "
         f"(trivial one-sided predictor = 0.50)", ""]
    for name, s in out["sets"].items():
        exp = "opens (+)" if s["expected_sign"] == +1 else "closes (-)"
        ci = (f"[{s['ci_lo']:.3f},{s['ci_hi']:.3f}]"
              if s["ci_lo"] == s["ci_lo"] else "[n/a]")
        L += [f"{name.upper()}  expected {exp}   n={s['n_regions']} regions / "
              f"{s['n_genes']} genes",
              f"  sign-accuracy {s['accuracy']:.3f}  boot {s['boot_mean']:.3f} {ci}  "
              f"vs base {s['base_rate']:.3f}  -> {'BEATS' if s['beats_base'] else 'ties/loses'}",
              "  by |direction|:"]
        for lo, hi, n, acc in s["strata"]:
            bar = f"{acc:.3f}" if acc == acc else "  -  "
            L.append(f"    [{lo:.2f},{hi:.2f})  n={n:<5d} acc={bar}")
        L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contract", required=True, help="contract TSV with a `direction` column")
    ap.add_argument("--dest", required=True, help="destination master-TF anchor BED (must open)")
    ap.add_argument("--source", required=True, help="source-cell master-TF anchor BED (must close)")
    ap.add_argument("--title", default="Claim 2B")
    args = ap.parse_args()

    chrom, start, end, _score, direction = load_contract(args.contract)
    if not np.isfinite(direction).any():
        raise SystemExit("contract has no `direction` column (or all-empty) -- 2B needs it")
    out = evaluate(chrom, start, end, direction, args.dest, args.source)
    print(_fmt(out, args.title))


if __name__ == "__main__":
    main()
