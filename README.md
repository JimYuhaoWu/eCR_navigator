# eCR_navigator

Part of the **eCR** engineered-chromatin-regulator platform for MEF→iPSC
reprogramming design, alongside
[`eCR_mod_lib`](https://github.com/JimYuhaoWu/eCR_mod_lib) (module library) and
[`eCR_predictor`](https://github.com/JimYuhaoWu/eCR_predictor) (DBD → structure →
fusion → off-target pipeline).

## What this is

A **region-weighting and target-nomination front-end**. During reprogramming,
some genomic regions open/close as **drivers** of the cell-state transition while
others move as **passengers** (downstream consequences). eCR_navigator assigns
each region a **driver-importance score** and can nominate the regions worth
targeting with an engineered chromatin regulator.

## How it connects to the rest of the platform

1. **Off-target severity (Tier 2).** `eCR_predictor`'s off-target module scores a
   DBD's unintended binding across accessible regions as `Σ (site strength ×
   region weight)`. Today it uses a provisional weight from accessibility
   *dynamics* (|Δaccessibility| proxy). eCR_navigator replaces that with a real
   **driver score** — an off-target hit in a driver region is far more dangerous
   than one in a passenger region. It plugs in through a stable region-weight
   contract; no change to eCR_predictor's scoring code.
2. **Target nomination.** Driver regions can be fed as targets into
   `eCR_predictor` Step 1.

## Output contract (what eCR_predictor consumes)

A TSV of scored regions — `chrom  start  end  driver_score` (driver_score in
[0, 1], higher = more driver-like) — on the **same genome assembly as the peaks**
for that run (mm10 for the mouse work, hg38 for human). See
[docs/region_weight_contract.md](docs/region_weight_contract.md). This contract is
the stable interface; the model behind it can evolve freely.

## Status

**Zero-shot driver scoring working across multiple foundation models (2026-07).**
Each model runs in its own GPU mirror and emits a per-region embedding artifact
([docs/embedding_artifact.md](docs/embedding_artifact.md)); `navigate.py` diffs two
cell states into the driver-score contract. Integrated so far: ChromBERT, GET,
ATACformer, ChromFound — all behind the *same* artifact contract, so adding the next
model is a new embed script, not a new pipeline. Species: hg38 + mm10 are
first-class; hg38-only models are bridged by liftOver (see the per-model docs and
CLAUDE.md). A supervised/fine-tuned head is the planned next axis once time-course
data is in hand.
