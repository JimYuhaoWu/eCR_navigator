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

### 2B — is the `direction` COLUMN itself correct?  *(deferred)*
Only a real (non-circular) test for models whose direction is a **prediction**:
- **EpiAgent** (SR head) and future **AlphaGenome** (DNase head) — testable: does the
  predicted open/close agree with measured signed-Δ (sign-agreement %, AUROC of
  `direction` → measured-opening label), restricted to master-TF loci?
- **GET / ChromFound / ATACformer / ChromBERT** — their `direction` **IS** the measured
  aTPM-Δ (input-measured or external-attach, see [`direction.md`](direction.md)), so
  testing it against signed-Δ is **circular**. Skip / sanity-check only.
- Caveat: EpiAgent is sparse (8,190-cCRE cap → 4–6 loci overlaps in Claim 1), so 2B is
  likely underpowered — worth running mainly to stand up the harness for AlphaGenome.

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
- 2B (direction-column correctness): still **deferred** — only non-circular for
  prediction-head models (EpiAgent / AlphaGenome).

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
