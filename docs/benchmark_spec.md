# Benchmark spec — a frozen, model-agnostic transition panel

> **Status: v1 scope SETTLED 2026-07-15 (see "v1 decisions" below); data not yet assembled.**
> Turns Claim-1/2 validation ([`validation_summary.md`](validation_summary.md)) from n=2 ad-hoc
> transitions into a frozen panel any current *or future / self-trained* model can be scored
> against with one command. Next action is the data-availability pass on the sourcing list.

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

## Composition target (v1) — SETTLED
Freeze **5–6 transitions**, balanced so the scorecard can discriminate, not just confirm:
- **have already:** fib→iN (strong, model-positive) · MEF→mES (clean-but-**null** under Claim 2,
  the instructive "clean yet directional" control — keep) · iCM GSE179011 (weak/partial —
  the weak-or-Gate-1-reject slot).
- **collect 2–3 new strong** from the sourcing priority list below (whichever clear the data
  bar), to reach 3–4 strong total.
- Panel must contain **≥1 weak/partial and ≥1 Gate-1 reject** (iCM is the current candidate;
  confirm which role it plays once its per-replicate matrix is scored), and **both species**.

### Sourcing priority list (settled — attempt all, take those with usable data)
Bar: **defined destination cell type + public bulk/pseudobulk ATAC at both endpoints, ≥2
reps/state.** Priority order by cleanliness of the master-TF definition:
1. **MyoD fib→myotube** — single master TF; cleanest possible positive control.
2. **C/EBPα B-cell→macrophage** — classic, well-defined; fast transdifferentiation.
3. **fib→hepatocyte** (Foxa/Hnf) — diversifies lineage (endoderm) beyond neuro/muscle/blood.
4. **fib→iNSC** (Sox2-led) — distinct destination identity, neuro-adjacent to iN.

Any of these that lack usable data → deferred to v2. If ≥2 same-destination across species turn
up (e.g. mouse+human MyoD), keep the **cross-species pair** (tests species-invariance).

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

## v1 decisions (settled 2026-07-15)
1. **Size:** freeze **5–6 transitions** (3–4 strong, ≥1 weak, ≥1 Gate-1 reject).
2. **Ground truth:** **master-TF loci primary** (scored on every transition); **binding
   (ChIP/CUT&Tag) optional secondary** where the study provides it — never a required column.
3. **Sourcing priority:** MyoD fib→myotube · C/EBPα B-cell→macrophage · fib→hepatocyte ·
   fib→iNSC (attempt all; take those clearing the data bar; rest → v2).
4. **Split:** **test-only for v1** — a clean held-out yardstick for zero-shot models. A reserved
   *train* split is added later, when a driver-supervised fine-tune corpus is actually being
   built (regime 1/3; see CLAUDE.md endpoint-only principle).

## Data-availability pass — DONE 2026-07-15

Searched GEO/literature for public bulk ATAC at both endpoints, ≥2 reps/state, on the four
sourcing candidates:

| Candidate | Master TF(s) | GEO (ATAC) | Species / asm | States assayed | Reps | Verdict |
|---|---|---|---|---|---|---|
| **C/EBPα pre-B→macrophage** | Cebpa | **GSE151748** | mouse / mm10 | 0,1,3,6,18,120 hpi (endpoints 0 h pre-B, 120 h macrophage) | **N=2**/tp | ✅ **PASS** — clean, dense time course (use WT arm, not the R35A Carm1 mutant) |
| **MyoD fib→iMPC** | Myod1 (+F/R/C small molecules) | **GSE186271** (ATAC sub-series, 8 samples) | mouse / mm10 | MEF, MyoD-d2, MyoD+FRC-d2, iMPC | **N=2**/group | ✅ **PASS** — endpoints MEF & iMPC. Note: iMPC (progenitor) is the *stable* end; pure myotube is unstable |
| fib→iNSC (Ptf1a) | Ptf1a | SRP136063 | mouse | MEF, miNSC10, control NSC | **~N=1** | ⚠️ **LIKELY FAIL** — appears single-replicate; neuro slot already covered by iN. Skip unless GEO shows reps |
| fib→hepatocyte (iHep) | Foxa3/Gata4/Hnf1a/Hnf4a | — | — | — | — | ❌ **FAIL (v2)** — no bulk ATAC with reps found; the field is RNA-seq / methylation / scATAC |

**Two new strong transitions clear the bar (C/EBPα, MyoD)** — enough to reach v1.

## Proposed v1 panel (5 transitions) — for confirmation
| # | Transition | Role | Species | Source |
|---|---|---|---|---|
| 1 | fib→iN (Ascl1) | strong / model-positive | human | GSE299923 (have) |
| 2 | **pre-B→macrophage (C/EBPα)** | strong (new) | mouse | GSE151748 |
| 3 | **MEF→iMPC (MyoD)** | strong (new) | mouse | GSE186271 |
| 4 | MEF→mES | clean-but-**null** control (Claim-2 signed-Δ-primary) | mouse | GSE201577 (have) |
| 5 | MEF→iCM | weak / partial (Gate-1 stress / reject candidate) | human | GSE179011 (have) |

Hits the target (3 strong, 1 null control, 1 weak; both species — 3 mouse / 2 human).
iNSC/iHep → v2. **Next after confirmation:** assemble each bundle (endpoints matrix + signed-Δ
+ target-cell master-TF anchor BEDs + Gate-1 label) per the layout above, starting with the two
new mouse sets (GSE151748, GSE186271) which are native mm10 and need no liftOver.
