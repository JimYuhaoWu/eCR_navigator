# Validation summary — is `driver_score` informative, and which score to trust?

**Canonical, current-state summary. Read this first; the deep per-phase docs (linked at the
bottom) hold the full trail, effect sizes, and honest negatives.** Two questions have been
answered on real transitions; one (2B) is deferred.

- **Claim 1 — is `driver_score` informative?** (recovers true drivers vs a
  change-magnitude-matched background) → **YES, for GET, on the master-TF *loci* reframe, on a
  strong clean transition.**
- **Claim 2A — does it add over a signed-Δaccessibility baseline?** → **YES on human iN
  (strong/clean), NO on mouse MEF→mES (signed-Δ dominates).**
- **Nomination policy** — a runnable per-transition preflight decides GET-vs-signed.
- **Claim 2B — is the `direction` *column* itself correct?** → **deferred** (circular for all
  current models; only testable on prediction-head models EpiAgent / AlphaGenome).

Method throughout: matched-background AUROC (control for |Δaccessibility|) + bootstrap CI +
top-k fold, on clean verified endpoints. Report AUROC/CI/fold, **not** permutation-p (n large).

---

## Claim 1 — informativeness

The informative target is **not where master TFs bind** (footprints fail to generalize across
cocktails — mouse OSKM→JGES, human) but the **cis-regulatory loci of the target-cell master-TF
genes** (promoters TSS±2kb + gene±50kb neighborhood), which must open.

| Model | mouse MEF→mES (clean GSE201577) | human fib→iN (GSE299923) |
|---|---|---|
| **GET** | ✅ master-TF loci **0.57–0.58** (robust, broad) | ✅ master-TF **promoters 0.668**, robust to opening-only |
| ChromFound | loci tail-only (top-5% 2–4×, null AUROC); OSKM opening-only 0.643 | null on loci; **positive on pioneer (Ascl1) binding 0.572** |
| ChromBERT | null | null / below-chance |
| ATACformer | not tested (mm10 liftOver too sparse) | null / below-chance |
| EpiAgent | not tested (too sparse) | too sparse (8,190-cCRE cap) |

Three lessons: (1) binding footprints **don't generalize** — use master-TF gene loci. (2) The
signal is **GET-specific** and needs a **strong clean transition** (weak human iCM GSE179011
was dropped; all models null there). (3) On mouse the GET signal is **largely directional**
(collapses to ~0.50 opening-only) → motivated Claim 2.

## Claim 2A — beyond signed-Δaccessibility (GET only)

signed-Δ promoted from Claim-1 *confound* to a *competing scorer*; paired ΔAUROC CI +
incremental logistic-LR (does driver improve `label ~ signed`?).

| System | positive | AUROC driver / signed | ΔAUROC [CI] | incr. LR | verdict |
|---|---|---|---|---|---|
| **human iN** | **promoter** | **0.664 / 0.494** | **+0.170 [+0.081,+0.264]** | p=0.001 | ✅ beats + adds |
| human iN | neighborhood | 0.557 / 0.499 | +0.057 [+0.012,+0.103] | p=0.002 | ✅ beats + adds |
| human iN | promoter (all-regions) | 0.668 / 0.599 | +0.069 [−0.009,+0.151] | p=0.001 | ⚠️ AUROC wash, LR sig |
| mouse | promoter | 0.579 / 0.507 | +0.071 [−0.045,+0.188] | p=0.13 | ➖ null (underpowered) |
| mouse | promoter (all-regions) | 0.582 / **0.633** | −0.051 | p=0.37 | ❌ signed-Δ wins |

**GET's value is regulatory-region prioritization at matched magnitude+direction**, real on a
strong clean transition, negligible over signed-Δ on weak/already-directional ones. Partly
overturns the "largely directional" read from Claim 1.

## Top-k confidence (the metric that matches the navigator's job)

Fold-enrichment of master-TF loci in the top-k of each model's ranking (its own base rate = 1×):

- **GET is the only model with usable top-k confidence, and it is front-loaded.** Human iN
  top-1% **~9–10×**, decaying to ~1.7× by top-10%; mouse ~2.9× at top-1%.
- **On mouse the signed-Δ top-k beats GET** (5–9× vs 1.4–2.9×) — master-TF promoters *are* the
  biggest openings there.
- **Every other model is null/mid-tail/too-sparse at the top**, both species (ChromFound = weak
  mid-tail top-2–5% mouse only; ChromBERT, ATACformer null; EpiAgent 0–5 in-universe positives).

**⇒ For nomination, trust GET's top ~1% on a strong clean transition; no other model or weak
transition earns a top-k.**

## Nomination policy — which score to trust (`scripts/preflight.py`)

No endpoint-only indicator predicts which score wins (**PC1 cleanliness fails** — mouse is
cleanest, PC1 0.933, yet GET loses; **rank-divergence from |signed-Δ| fails** — ~equal both
species). So don't *predict* "clean enough"; *measure* it per transition:

- **Gate 1 — admissibility** (endpoint-only): ≥2 reps/state, replicate-coherence margin ≥0.10,
  PC1 ≥0.80, all computed on a **fixed universe of the 50,000 most accessible regions**. The
  fixed universe is load-bearing: PC1/coherence rise as low-signal regions are dropped, so on
  raw universes (the v1 panel spans 63k–1.06M regions, 17×) the values are **not comparable**
  and iN — our best transition — scored 0.792 and *rejected*. At the fixed universe it scores
  0.919, and all six panel transitions match their expected verdict with 0.80 inside a real
  0.785→0.919 gap. Gate 1 is a **coarse screen for severe failures**, still **not** a reliable
  admit (mouse clears it with the cleanest endpoints and GET still loses) — Gate 2 decides.
- **Gate 2 — decision:** eCR design always targets a *known* cell type → its canonical
  master-TF loci (known biology) go into the Claim-2A harness as positives. GET is **PRIMARY**
  iff ΔAUROC CI excludes 0 **or** incremental-LR p<0.05 (driver coef>0), else **signed-Δ** is
  primary. Measured signed-Δ always attached for open/close **direction**.

Validated: human iN → **PRIMARY=GET**; mouse MEF→mES → **ADMIT but PRIMARY=signed-Δ** (the
instructive "clean but still directional" case). Thresholds first-pass, **calibrated on n=2**
transitions — tighten as the benchmark grows.

## `driver_score` is a magnitude, not a signed call

All current models' "direction" is the *measured* aTPM-Δ (input-measured or external-attach),
so pair every nominated region with its measured signed-Δ for the predictor. Whether a model
can *predict* direction (Claim 2B) is only non-circular for prediction-head models
(EpiAgent SR head, AlphaGenome DNase head) — deferred; see [`claim2_plan.md`](claim2_plan.md) §2B.

---

## Where the detail lives

| Need | Doc / data |
|---|---|
| Claim 1 mouse full trail (phases 1–2, JGES, loci, H3K27ac) | [`claim1_results.md`](claim1_results.md) · TSV [`claim1_results.mm10.tsv`](claim1_results.mm10.tsv), [`claim1_results.mtf.tsv`](claim1_results.mtf.tsv) |
| Claim 1 human (all 5 models, iN + dropped iCM) | [`claim1_human_progress.md`](claim1_human_progress.md) · TSV [`claim1_results.human.tsv`](claim1_results.human.tsv) |
| Claim 1 session handoff / reproduction paths | [`claim1_progress.md`](claim1_progress.md) |
| Cross-model magnitude (non-)consistency | [`cross_model_consistency.md`](cross_model_consistency.md) · TSV [`cross_model_consistency.mm10.tsv`](cross_model_consistency.mm10.tsv) |
| Claim 2 plan (2A scope, 2B deferred rationale) | [`claim2_plan.md`](claim2_plan.md) |
| Claim 2A results, top-k sweeps, per-model×species confidence, nomination policy | [`claim2_results.md`](claim2_results.md) · TSV [`claim2_results.tsv`](claim2_results.tsv) |
| Eval + preflight code | `scripts/eval_driver_claim1.py`, `scripts/eval_driver_claim2.py`, `scripts/preflight.py` (+ `tests/`) |
| Server-side artifacts | PeiLab2 `/mnt3/wuyuhao/{claim1_work,in_clean,mtf_loci,neural_gt,jges_gse199612}/` (see `server_mirrors.md`) |

## Open items

- **Claim 2B** (direction-column correctness) — deferred; needs a prediction-head model.
- **Benchmark v1 — BUILT AND RUN** (2026-07-16): 6 transitions × GET (+ChromBERT on the mouse
  bundles). **The Claim-1/2 finding generalized: GET wins on BOTH strong clean transitions —
  iN (0.668) *and* C/EBPα macrophage (0.640, ΔvsSigned +0.106 p=0.001)** — a second lineage and
  species, so it is not an iN artifact. GET fails on all three Gate-1-reject transitions, and
  **Gate-1 separated admit/reject cleanly on all six**, so the preflight thresholds are now
  calibrated on n=6 rather than n=2. Results: [`benchmark_v1_results.md`](benchmark_v1_results.md);
  design + v2 backlog: [`benchmark_spec.md`](benchmark_spec.md).
