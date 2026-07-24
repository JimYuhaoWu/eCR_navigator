# Next steps under limited driver-label data

**Context.** Four independent experiments (see `finetune_results.md`, `claim2_results.md`)
show that **supervised fine-tuning on the current driver-label scale does not help** — the
head never beats zero-shot `driver_score` **or** signed-Δaccessibility, and usually hurts;
**measured signed-Δ is the strongest simple driver signal**. The corpus is tiny (tens of
loci across ~15–26 genes, or a couple of non-independent binding panels), and a
foundation-model fine-tune (even a low-capacity head) overfits and *loses* the zero-shot
prior. So the question is: what *does* move the needle when labels are scarce?

## The principle
**Don't add model capacity — add label quality/quantity, or improve the label-free signal.**
Capacity is not the bottleneck (the frozen prior already encodes regulatory logic); the
bottleneck is (a) too few, (b) proxy (binding/loci ≠ causal driver) labels.

## Recommended, in priority order

### 1. Treat the model as a fixed zero-shot prior; combine signals with ≤3 parameters
The only "supervised" object that can't overfit at this scale is a **1–3 parameter
calibration** over signals we already trust: **signed-Δ** (strongest), **zero-shot
`driver_score`**, and (orthogonal, measured, endpoint-only) **RNA-Δ** of the region's
target gene. Fit/validate leave-one-transition-out. Expect signed-Δ to dominate; the
combiner just formalizes "rank by opening, use the model to break ties." This is the
immediate, safe deliverable for eCR_predictor.

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
