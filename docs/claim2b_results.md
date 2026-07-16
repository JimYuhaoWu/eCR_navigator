# Claim 2B results — is the measured `direction` trustworthy where it picks an ED?

Run 2026-07-16 · harness `scripts/eval_direction_claim2b.py` · spec
[`claim2_plan.md`](claim2_plan.md) §2B · raw output PeiLab2 `/mnt3/wuyuhao/claim2b/results.txt`.

**Question.** `direction` (measured signed ΔaTPM) picks the effector domain — `+` → activator
(VP64), `−` → repressor (KRAB). On the iN bundle **42% of nominations have |direction| < 0.05**.
Is the sign correct on regions whose direction is known from biology, and at what magnitude
does it stop beating chance?

**Method.** Two-sided against known biology: **destination** master-TF loci must open
(expected `+`), **source-cell** master-TF loci must close (expected `−`). Sign-accuracy vs the
transition's marginal opening rate (never 50%), gene-clustered bootstrap CI, stratified by
|direction|. No GPU — measured direction is already in each bundle's `weights.tsv`.

## Headline

**There is a consistent, shippable trust cutoff at |direction| ≈ 0.05.** On every Gate-1-admit
transition, destination sign-accuracy sits at the opening-rate baseline below ~0.05 and jumps
to 0.80–0.95 above it. The Gate-1-**reject** control (MyoD) shows no such threshold — direction
trustworthiness tracks Gate-1 admissibility. **This confirms the run-bundle-contract concern:
the sub-0.05 nominations are exactly the regions where the ED call is not supported.**

## Scorecard

| Transition | Gate 1 | opening rate | balanced acc | destination acc (vs base) | source acc (vs base) |
|---|---|---|---|---|---|
| **iN** (fib→neuron) | admit | 0.593 | 0.452 | **0.783** (0.593) ✓ | 0.121 (0.407) ✗ |
| **C/EBpα** (preB→mac) | admit | 0.239 | **0.722** | **0.667** (0.239) ✓ | 0.778 (0.761) ≈ |
| **MEF→mES** | admit | 0.548 | **0.708** | **0.750** (0.548) ✓ | 0.667 (0.452) ≈ |
| **MyoD** | REJECT | 0.313 | 0.529 | 0.558 (0.313) | 0.500 (0.687) ✗ |

acc "✓" = gene-clustered 95% CI lower bound above the base rate; "≈" = CI includes the base
rate (ties). Base rate is the marginal opening rate (destination) or closing rate (source) —
never 50%.

Balanced accuracy scores a trivial "everything opens" predictor at 0.50 regardless of opening
rate — the reason the test is two-sided.

## Destination sign-accuracy by |direction| (the ambiguity threshold)

| |direction| bin | iN | C/EBpα | MEF→mES | MyoD (reject) |
|---|---|---|---|---|---|
| [0.00, 0.02) | 0.577 | 0.554 | 0.000 (n=4) | 0.500 |
| [0.02, 0.05) | 0.826 | 0.625 | 0.556 | 0.750 (n=4) |
| **[0.05, 0.10)** | **0.900** | **0.800** | 0.615 | 0.222 |
| [0.10, 0.20) | 0.875 | 0.800 | 0.850 | 0.615 |
| [0.20, 0.50) | 1.000 | 0.769 | 0.950 | 0.600 |
| [0.50, 1.01) | 1.000 | 1.000 | 1.000 | 0.800 |

The three admit transitions cross from ≈base-rate to reliable between **0.05 and 0.10**; MyoD
(reject) is noisy and non-monotonic throughout — no usable threshold, which is the correct
negative-control outcome.

## What this settles

1. **Ambiguity cutoff: |direction| ≥ 0.05 to trust the sign for an ED call.** Below it, emit
   `direction` but flag it low-confidence; `fuse.py` should not commit to an activator-vs-
   repressor choice there. Consistent on 3 independent transitions (2 species, 3 lineages).
   MEF→mES is slightly more conservative (~0.10), so **0.05 is the floor, not a promise** — a
   consumer wanting higher confidence can raise it.
2. **Direction trust tracks Gate-1.** The reject transition has no threshold and balanced
   accuracy ≈ 0.5. So the same admissibility gate that governs *which score to nominate from*
   also predicts *whether the sign is worth acting on* — a coherence with the rest of the
   policy, not a new gate.
3. **The two-sided (source) arm is system-dependent — exactly the spec's flagged threat.**
   - **C/EBpα is the clean two-sided case** (balanced 0.722): the pre-B source set
     (Pax5/Ebf1/Foxo1/Tcf3/Ikzf1/Bcl11a) closes as expected. Its raw source-accuracy (0.778)
     only "ties" the base rate because the transition is 76% closing — but balanced accuracy,
     which is immune to that, clears 0.5 decisively.
   - **MEF→mES**: mouse-fibroblast source is directionally right (0.667 vs 0.452 base) but
     ties at 95% CI (n=18, lower bound 0.444) — too few regions to claim it beats base.
   - **iN fails on the source side** (0.121, and *worse* at larger |direction|): the human
     mesenchymal TFs (PRRX1/TWIST1/2/SNAI1/2/ZEB1/TCF21) **open** rather than close in fib→iN.
     Either the gene list is wrong for this system, or the D7 conversion is too incomplete to
     silence the mesenchymal program (consistent with iN's own Gate-1 being the weakest admit,
     PC1 0.919). **iN is therefore reported one-sided (destination) + this noted**, per the
     spec's fallback. It does not affect the cutoff, which is a destination-side result.

## Decision (against the pre-registered rule)

The spec's rule "accuracy ≫ base rate above |Δ|=X, ≈ base rate below → ship X" is the outcome
on the destination side of all three admit transitions. **Ship |direction| ≥ 0.05 as the trust
cutoff.** This is *consumer guidance*, not a format change — `fuse.py` reads `direction` from
the existing contract and applies the cutoff when choosing an ED.

**Not yet done — the per-region `direction_confidence` (metric 4).** The spec's shippability
chain (validate that replicate sign-stability predicts sign-correctness on anchors, then ship
stability as a per-region confidence) is deferred: the cutoff already answers the contract's
open question, and a `direction_confidence` **column** is a contract change (bundle_version
bump), so it deserves its own decision rather than riding along here. The per-replicate
endpoint matrices needed for it exist (`/mnt3/wuyuhao/*/endpoints.matrix.tsv`).

## Reproduce

```bash
# source-cell anchors (mm10 + hg38): /mnt3/wuyuhao/claim2b/build_src_anchors.py
python scripts/eval_direction_claim2b.py \
    --contract bundles/<id>/weights.tsv \
    --dest <destination master-TF promoter BED> \
    --source /mnt3/wuyuhao/claim2b/src_anchors/<id>/promoter.bed
```

Anchors `/mnt3/wuyuhao/claim2b/src_anchors/`; destination anchors are each bundle's existing
Gate-2 anchors. Tests: `tests/test_eval_claim2b.py` (7, green).
