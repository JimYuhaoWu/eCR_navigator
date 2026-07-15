# Benchmark spec — a frozen, model-agnostic transition panel (DRAFT / sketch)

> **Status: proposal for review.** Sketches the dataset that turns Claim-1/2 validation
> ([`validation_summary.md`](validation_summary.md)) from n=2 ad-hoc transitions into a frozen
> panel any current *or future / self-trained* model can be scored against with one command.
> Nothing built yet — this is the design to react to before we assemble data.

## Why

Today's conclusions rest on **two** transitions (human iN positive, mouse MEF→mES null, iCM
dropped). That is too thin to (a) trust the preflight thresholds (PC1≥0.80, coherence≥0.10,
"ΔAUROC CI>0 or LR p<0.05") and (b) fairly rank a new model. The pipeline is deliberately
**contract-based** (a model only has to emit `chrom,start,end,driver_score[,direction]`), so
the benchmark should exploit that: a new model plugs in by producing that contract per
transition — **no benchmark code changes** — and gets the full scorecard.

## Unit of the benchmark: a *transition entry*

Each transition is a self-contained, frozen bundle:

| Component | What | Source / recipe |
|---|---|---|
| **endpoints** | bulk (or pseudobulk) ATAC, **≥2 reps/state**, both states | GEO; keep per-replicate matrix (Gate-1 needs it) |
| **signed-Δ track** | `chrom,start,end,signedΔaccessibility` over the union peak set | the baseline scorer **and** the |Δ|-matching confound (built as in Claim 1) |
| **anchors** (ground truth) | target-cell master-TF **promoter** (TSS±2kb) + **neighborhood** (gene±50kb) BEDs | curated **known biology** — see recipe below (NON-circular: from literature, not this transition) |
| **binding GT** (optional) | in-study master-TF ChIP/CUT&Tag top-decile BED | secondary axis; only where the study provides it |
| **Gate-1 label** | ADMIT/REJECT + PC1 + coherence margin | `scripts/preflight.py` admissibility on the endpoint matrix |
| **expected verdict** | model-positive / signed-Δ-primary / Gate-1-reject | our prior, for tracking calibration (not used in scoring) |
| **provenance** | GEO acc, assembly, start/end cell types, cocktail, notes | registry row |

### Anchor-assembly recipe (the ground truth — must stay non-circular)
1. For the **destination** cell type, curate its canonical master TFs from literature (a frozen
   per-cell-type list, e.g. neuron = Ascl1/Neurod1/Myt1l/Brn2; pluripotency = Pou5f1/Sox2/
   Nanog/Klf4/Esrrb/…). Independent of any ATAC data.
2. Map each TF gene to refGene (per assembly) → **promoter** = TSS±2kb; **neighborhood** =
   gene span ±50kb. Freeze the BEDs with the TF list that produced them.
3. Require **≥~20 anchors in a model's universe** for a reliable Gate-2 verdict; if a model's
   universe is sparse, widen promoter→+neighborhood (and record it).

## Eligibility (what may enter the panel)
- **Direct reprogramming / transdifferentiation with a *defined destination* cell type** (so
  master TFs are known). Excludes vague/gradient transitions.
- Bulk/pseudobulk ATAC at **both** endpoints, ≥2 reps/state, one assembly at minimum
  (hg38 and/or mm10; note which).
- **Gate-1 REJECTs are kept as negative controls**, not discarded — they test that the
  benchmark (and any model) correctly refuses to nominate on a bad transition.

## Composition target (v1)
Balance so the scorecard can discriminate, not just confirm:
- **≥3 strong clean transitions** (expect model-positive) — e.g. fib→iN (have), MEF→mES (have,
  but note: model-*null* under Claim 2 — a clean-yet-directional case, keep it),
  + new: B-cell→macrophage (C/EBPα), fib→myotube (MyoD), fib→hepatocyte, fib→iNSC.
- **≥2 weak/partial** (expect signed-Δ-primary or Gate-1 reject) — e.g. iCM (GSE179011, dropped).
- **both species** represented; where possible a **same-destination cross-species pair** (tests
  species-invariance of a model's signal).

## Per-(model × transition) scorecard (generated)
Reuses `eval_driver_claim1.py`, `eval_driver_claim2.py`, `preflight.py` — no new stats code:
- **Claim-1 informativeness:** matched AUROC + bootstrap CI on promoter & neighborhood.
- **Claim-2A vs signed-Δ:** ΔAUROC (paired CI) + incremental-LR p + driver coef.
- **Top-k confidence:** fold enrichment at top-1/2/5%.
- **Gate-2 verdict:** PRIMARY = *model* or *signed-Δ*.
A model **passes** a transition iff Gate-2 → model (it beats signed-Δ).

### Aggregate (across the panel, per model)
- **Win rate:** fraction of *admissible* transitions where the model is PRIMARY.
- **Mean top-1% confidence** on its wins.
- **Preflight calibration:** do Gate-1 ADMIT + model-win transitions separate cleanly from the
  rest? (this is how the n=2 thresholds get retuned into defensible ones.)

## Layout (frozen, contract-based)
```
benchmark/
  transitions.tsv                     # registry: id, species, start, end, GEO, assembly, gate1, expected
  <transition_id>/
    endpoints.matrix.tsv              # per-replicate (Gate-1 input)
    signed_delta.tsv                  # baseline scorer + confound
    anchors/{promoter,neighborhood}.bed
    anchors/binding.bed               # optional
    gate1.json                        # PC1, coherence, ADMIT/REJECT
  models/<model>/<transition_id>.driver_scores.tsv   # ← the ONLY thing a new model must add
  scorecard.tsv                       # generated: model × transition × metrics
  VERSION                             # v1 freezes panel + anchor lists + thresholds
```
**A new / self-trained model integrates by:** running its existing embed→`navigate.py` on each
transition's endpoints, dropping the contract into `models/<model>/`, and running the scorecard
generator. Nothing else changes — that is the model-agnostic promise, enforced by the contract.

## Versioning
Freeze **v1** (panel + per-cell-type TF lists + anchor BEDs + preflight thresholds). Bump the
version when transitions or anchor lists change, so cross-model scores stay comparable within a
version.

## Open questions for review
1. **v1 size** — how many transitions to target before first freeze (I'd suggest 5–6: 3–4
   strong, 1–2 weak, ≥1 Gate-1 reject)?
2. **Binding as a second axis** — include in-study ChIP/CUT&Tag ground truth where available, or
   keep the benchmark loci-only for consistency?
3. **New-transition sourcing** — which reprogramming systems to prioritize collecting, given the
   "defined destination + bulk ATAC + reps" bar?
4. **Self-trained-model hook** — do we also want a *training* split (transitions reserved for
   fine-tuning a driver-supervised model) vs a held-out *test* split, or is v1 test-only?
