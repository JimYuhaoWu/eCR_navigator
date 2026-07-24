# Next steps under limited driver-label data

**Context.** Fine-tuning experiments ([`finetune_results.md`](finetune_results.md)) show that
**supervised fine-tuning on the current driver-label scale does not help** — the head never
beats the **zero-shot `driver_score`** baseline and usually hurts. Note the zero-shot score is
itself validated and, on strong transitions, **beats** signed-Δ (Claim 2A: iN 0.664 vs 0.494 —
[`validation_summary.md`](validation_summary.md)); so the model *is* the asset, it just should
be used **zero-shot, not fine-tuned**. The corpus is tiny (tens of loci across ~15–26 genes, or
a couple of non-independent binding panels), and a foundation-model fine-tune (even a
low-capacity head) overfits and *loses* the zero-shot prior. So the question is: what *does*
move the needle when labels are scarce?

## The principle
**Don't add model capacity — add label quality/quantity, or improve the label-free signal.**
Capacity is not the bottleneck (the frozen prior already encodes regulatory logic); the
bottleneck is (a) too few, (b) proxy (binding/loci ≠ causal driver) labels.

## Recommended, in priority order

### 1. Treat the model as a fixed zero-shot prior; combine signals with ≤3 parameters
The only "supervised" object that can't overfit at this scale is a **1–3 parameter
calibration** over signals we already trust: **zero-shot `driver_score`**, **signed-Δ**, and
(orthogonal, measured, endpoint-only) **RNA-Δ** of the region's target gene. Fit/validate
leave-one-transition-out. Which signal leads is **transition-dependent** — `driver_score` on
strong/clean transitions (iN), signed-Δ on weak/already-directional ones (mouse) — which is
exactly what a small combiner is for. Note the repo already ships this as the **nomination
preflight** (`nominate.py` Gate-1/Gate-2 picks GET-vs-signed per transition); a learned
≤3-param combiner is the natural generalization. This is the immediate, safe deliverable for
eCR_predictor.

### 2. Improve the ZERO-SHOT readout (label-free) instead of supervising
The prior is the asset — sharpen how we read it, no labels required:
- **In-silico region perturbation / ablation importance** (GET supports this) — more causal
  than embedding-shift magnitude, and it's what the "driver vs passenger" question really asks.
- **Multi-model ensemble** — rank-average the 5 models; cheap, may beat any single one.
- **Better endpoint definition** — Claim 1 showed endpoint quality dominates the signal;
  cleaner reference states buy more than any head.
- **Add RNA-Δ** as an orthogonal measured channel (endpoint-only).

### 3. Invest in LABELS, two tracks (the real Regime-1/3 unlocks)
- **(a) Breadth — weak-supervision corpus.** Aggregate driver labels across MANY completed
  reprogramming systems (public ChIP/CUT&Tag for dozens of cocktails/cell types). Tens of
  transitions × their master TFs is the scale at which a head could learn a *transferable*
  signature — the thing that failed at n=1 transition. This is "the open data-assembly
  problem" (CLAUDE.md). Caveat: our transfer result hints even breadth may be dominated by
  accessibility, so treat it as a test, not a foregone win.
- **(b) Causality — perturbation labels (Regime 3).** CRISPRi/a screens, or the eCR
  **design-build-test-learn loop**, give *causal* driver labels for THIS system — few but
  high-quality, and not obtainable from binding/loci proxies. Each iteration converts a
  regime-2 problem into a regime-1 one. **Highest long-term value**, and it is where the
  platform is already heading.

### 4. If (and only if) a real corpus exists, adapt parameter-efficiently
LoRA/adapters (few params), frozen backbone, early-stop on a held-out transition. Pointless
at n=1 transition (it only delays overfit); revisit once track 3(a) lands.

## Headline
Fix the model as a **zero-shot prior**; ship a **≤3-parameter combiner** (signed-Δ +
driver_score + RNA) now; put the real investment into **perturbation labels via the eCR
DBTL loop** (causal, matches the platform's direction), with a **breadth corpus** as the
computational alternative. Do **not** fine-tune the backbone on the current labels.
