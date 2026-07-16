# Claim 2 plan — is the navigator's DIRECTION informative beyond signed-Δaccessibility?

Follows Claim 1 (*is `driver_score` informative?* — see [`claim1_results.md`](claim1_results.md),
[`claim1_human_progress.md`](claim1_human_progress.md)). Claim 1's positive result (GET
recovers master-TF loci vs a |Δaccessibility|-matched background, on clean strong
transitions) turned out to be **largely directional** — it mostly says *these regions
open*. So Claim 2 asks the honest follow-up:

> **Does the navigator's directional information add anything beyond the plain
> signed-Δaccessibility you already measure from the ATAC data?**

(The same test that decided [AlphaGenome's](alphagenome_pipeline.md) two-state case: if
signed-Δ from your own data dominates, the model's real contribution is
regulatory-region *prioritization*, not direction.)

## Two sub-claims

### 2A — `driver_score` vs a signed-Δ baseline  *(SCOPED FIRST — harness built)*
Promote signed-Δ from Claim 1's *confound* to a *competing scorer* and ask whether
`driver_score` **beats or adds over** it at recovering master-TF loci.

- **Models:** GET only — the sole model with a rankable Claim-1 signal. The others
  (ChromBERT/ATACformer null, ChromFound loci-null, EpiAgent too sparse) would only
  re-confirm their Claim-1 nulls.
- **Systems:** where Claim 1 was positive — mouse **MEF→mES** (GSE201577) and human
  **iN** (fibroblast→induced neuron).
- **Positives:** the master-TF loci / promoter BEDs already assembled (`neural_gt/` for
  iN; JGES / pluripotency panels for mouse).
- **Baseline scorer + confound:** the per-model signed-Δ track already on disk
  (`bedtools map` of the cCRE signed-Δ onto each contract — built for Claim 1).
- **Readouts** (one matched, opening-only sample = positives + |signed-Δ|-matched
  negatives, same construction as Claim 1):
  1. `AUROC(driver_score)` vs `AUROC(signed-Δ)` — head to head.
  2. **ΔAUROC** with a *paired* bootstrap CI (same resampled indices for both scorers).
  3. **Incremental LR test** — does `driver_score` improve a logistic model that
     already has signed-Δ? `LR = 2(ll_full − ll_reduced)`, permutation p-value, and the
     standardized `driver_score` coefficient.
  4. top-5% fold for both scorers.
- **Verdict logic:** `driver_score` AUROC CI above signed-Δ's **and** a significant
  incremental LR → the model adds regulatory info beyond "what opens" (2A supported).
  Overlapping CIs / null increment → `driver_score` is dominated by signed-Δ for this
  endpoint pair (the honest null).

### 2B — is the `direction` column correct, and where does it stop being trustworthy?  *(SPEC — not yet run)*

> **Reassigned 2026-07-16.** 2B used to be "does a model's *predicted* direction match the
> measured Δ?". That is **model QC for a prediction head, not a platform question** — under
> the endpoint-only principle the measured signed-Δ is always computable and *is* the design
> target, so we would never prefer a prediction to it. It moved to
> [`alphagenome_pipeline.md`](alphagenome_pipeline.md) §Validation, gated on that model
> existing. The 2B slot now holds the question that is actually load-bearing.

**The question.** `direction` picks the **effector domain** — `+` → activator (VP64), `−` →
repressor (KRAB). On the real iN bundle **42% of nominations have |direction| < 0.05** and 51
are exactly 0.0. So: *is the measured sign strong enough to bet an ED on, and at what
magnitude does it stop being trustworthy?* This closes the ambiguity-threshold decision left
open in [`run_bundle_contract.md`](run_bundle_contract.md).

**Why it is non-circular.** Same argument that licenses Gate 2: **known biology fixes the
expected direction independently of any measurement.**

| Anchor set | Expected sign | Rationale |
|---|---|---|
| **destination** master-TF loci | `direction > 0` (opens) | the target cell's identity TFs must switch on — already built (Gate-2 anchors) |
| **source-cell** master-TF loci | `direction < 0` (closes) | the starting cell's identity TFs must switch off — **new asset, must be curated** |

**Two-sided is the whole design.** A one-sided "do destination anchors open?" is beaten by a
trivial *everything opens* baseline — the iN universe is ~58% opening, so that is the floor.
The source set has no such escape: a model that says "open" everywhere scores 0 on it.

#### Metrics
1. **Sign-accuracy per anchor set** vs that transition's marginal opening rate (the base rate
   is the baseline, not 50%). Bootstrap CI.
2. **Balanced accuracy across both sets** — the headline; immune to the trivial baseline.
3. **Stratified by |direction|** (bins/deciles): sign-accuracy *within each bin*, against
   that bin's own base rate. **Where accuracy falls to the base rate is the empirical
   ambiguity threshold** — derived, not picked. Stratifying also handles the confound that
   anchors may sit at systematically larger |Δ| than background.
4. **Replicate sign-stability** — from the per-replicate endpoint matrices built for Gate 1:
   resample/leave-one-out the replicates, recompute Δ, and record the fraction of resamples
   agreeing on the sign.

#### The chain that makes it shippable
Sign-accuracy is only measurable **on anchors** (that is where truth is known), but the
predictor needs a confidence for **every** nominated region. Replicate sign-stability is
computable **everywhere**. So:

> validate on anchors that **stability predicts sign-correctness** → then ship *stability* as
> the per-region `direction_confidence`.

If that link fails, we report the |direction| cutoff alone and do not invent a per-region
confidence.

#### Scope
- **Transitions:** the three Gate-1-admit bundles (iN, C/EBPα, MEF→mES). Add **MyoD** as a
  negative control — direction should be *less* trustworthy on a Gate-1-reject.
- **No GPU.** Contracts, endpoint matrices and destination anchors are all on PeiLab2;
  source-cell anchors are the only new asset (`benchmark/build_anchors.py` already
  generalizes — curate the TF list, rerun it).

#### Statistical care
- **Anchors are not independent** — promoter and neighborhood of the same gene overlap, and
  a gene contributes many regions. **Cluster by gene** for CIs, or restrict to one promoter
  per gene.
- Report accuracy **against the base rate**, never against 50%.
- Bin by |direction| **before** comparing, so the anchor-vs-background magnitude difference
  cannot masquerade as accuracy.

#### Decision rule (fixed before running)
| Result | Read | Action |
|---|---|---|
| Accuracy ≫ base rate above \|Δ\|=X, ≈ base rate below | direction is trustworthy only above X | ship **X** as the ambiguity cutoff; below X emit `direction` but flag low-confidence, and let `fuse.py` refuse to pick an ED |
| Accuracy ≫ base rate at **all** magnitudes | even small Δ carries a real sign | no cutoff; record that the 42%-below-0.05 worry is unfounded |
| Accuracy ≈ base rate **everywhere** | the direction column is not trustworthy on known-direction regions | **escalate** — this would undercut the ED choice platform-wide, not just at small Δ |
| Stability predicts accuracy | confidence is derivable per region | ship `direction_confidence` in the contract (bundle_version bump) |

#### Deliverables
`scripts/eval_direction_claim2b.py` (reuses `eval_driver_claim1.py`'s loaders /
`overlap_labels`) · tests · `docs/claim2b_results.md` · a contract decision (cutoff and/or
`direction_confidence`).

#### Main threat to validity
**Source-cell master-TF curation.** Pre-B is strong — Pax5/Ebf1/Foxo1 are canonical and
C/EBPα repressing Pax5 is *the* textbook result — so C/EBPα is the cleanest two-sided case.
**Fibroblast identity TFs are fuzzier** (Twist2/Prrx1/Snai2…), which weakens the source set
for iN, MEF→mES and MyoD. If the fibroblast set is too soft, fall back to reporting C/EBPα
two-sided and the others one-sided-with-base-rate, and say so.

### Claim 3 — counterfactual direction *(parked; hard)*
*If we place an eCR at this region, which way does accessibility move, and how far?*
Causal, not correlational — the question measurement cannot answer and the one that would
genuinely de-risk design. Needs a sequence→function model (**AlphaGenome**) and a validation
strategy we do not yet have: the honest ground truth is **perturbation data** (regime 3 in
CLAUDE.md), which we would have to generate. **Parked deliberately — interesting but very
hard to prove.** Revisit when AlphaGenome is built and/or wet-lab perturbation data exists.

## Tooling

- **`scripts/eval_driver_claim2.py`** — built. Reuses Claim 1's primitives
  (`overlap_labels`, `matched_negative_indices`, `auroc`, `topk_fold_enrichment`,
  loaders) and adds `logistic_fit`, `incremental_lr_test`, `paired_delta_auroc`,
  `evaluate_claim2`. Pure numpy. CLI:
  ```bash
  python scripts/eval_driver_claim2.py --contract get.iN.hg38.tsv \
      --positives neural_gt/masterTF_promoters.bed \
      --signed get.iN.signedDelta.tsv        # opening-only by default
  ```
- **`tests/test_eval_claim2.py`** — built, 6 tests green. Plants a signed-Δ confound and
  a separable driver-specific signal; asserts the harness reports *null* when
  `driver = signed-Δ + noise` and *positive* only when a real increment is present
  (logistic fit + LR calibration checked too).

## Status & next step

- Harness + tests: **done, mirror-independent**.
- 2A run + recorded (2026-07-15): see [`claim2_results.md`](claim2_results.md) +
  [`claim2_results.tsv`](claim2_results.tsv). **Supported on human iN** (GET driver_score
  beats + adds over signed-Δ, decisively on promoters: opening-only ΔAUROC +0.170
  CI[+0.081,+0.264], incremental LR p=0.001; the increment survives the tougher
  all-regions baseline); **null on mouse MEF→mES** (underpowered / signed-Δ dominates).
- 2B **reassigned + specced (2026-07-16), not yet run**: now *is the measured `direction`
  trustworthy at the magnitude where it picks an ED, and where is the cutoff?* — two-sided
  against destination (opens) and source-cell (closes) master-TF anchors. No GPU; only the
  source-cell anchors are a new asset. The old 2B (predicted-vs-measured direction) is model
  QC and moved to [`alphagenome_pipeline.md`](alphagenome_pipeline.md) §Validation.
- **Claim 3** (counterfactual direction) parked — see §Claim 3 above.

## Honest expectation (recorded before the run — kept for the record)

Because Claim 1's signal was "largely directional," 2A may well show `driver_score`
does **not** beat signed-Δ except possibly as an *increment* on the clean GET
transitions. Either way it is a defensible result: it tells the platform whether the
foundation model's value is *direction* (probably not — you already measure that) or
*regulatory-region prioritization at matched magnitude+direction* (the Claim-1 GET
signal). This directly informs whether direction is worth model compute at all.

**Outcome vs expectation:** the increment *did* materialize on human iN — and more
strongly than expected (driver also beats signed-Δ head-to-head on promoters, not only
as an LR increment). Mouse matched the pessimistic prior. So the platform read is
system-dependent: the model's value is **regulatory-region prioritization** on strong
clean transitions, and negligible over signed-Δ on weaker/already-directional ones.
