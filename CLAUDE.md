# CLAUDE.md — eCR_navigator

## Coding principles

- **No features beyond what was asked.** No speculative abstractions, configurability, or error handling for impossible scenarios.
- **Surgical changes.** Touch only what the task requires. Don't improve adjacent code, comments, or formatting. Match existing style.
- **Surface tradeoffs before coding.** If multiple interpretations exist, present them — don't pick silently. If something is unclear, ask.
- **If you notice unrelated dead code or issues, mention them — don't silently fix them.**

## What this project is

Driver-vs-passenger **genomic region weighting** for MEF→iPSC reprogramming, and
the target-nomination front-end of the **eCR** platform. It produces per-region
**driver scores** that (a) weight off-target severity in `eCR_predictor` (Tier 2)
and (b) nominate target loci for engineered chromatin regulators.

## Sibling repos (all under github.com/JimYuhaoWu)

- `eCR_mod_lib` — module library (default branch **main**)
- `eCR_predictor` — prediction/design pipeline (default branch **master**)
- `eCR_navigator` — this repo (default branch **main**)

They sit as siblings; dev on Windows/Mac, test on the Linux HPCC.

## Inputs on hand

- **ATAC differential accessibility** — MEF vs mES MACS peaks, 5-col BED, **mm10**
  (e.g. `MEF.e7_peaks.bed`, `mES.e7_peaks.bed`). The union of these already drives
  eCR_predictor's off-target background; eCR_navigator adds the *importance* layer.
- **RNA differential expression** — MEF vs mES.
- **Public reprogramming time-course** data — to be collected (drivers change
  early/causally, passengers late/downstream — temporal ordering is signal).
- `atac_embed.pkl` — an ATAC embedding of unknown provenance; confirm before use.

## Open design decision — DECIDE THIS FIRST (new session)

How to score driver vs passenger:

- **Supervised on the reprogramming trajectory.** Use ATAC/RNA time-course:
  drivers change *early and causally*, passengers *late*. Needs temporal (or
  perturbation) data with enough resolution to order events.
- **Zero-shot / pretrained genomic LM** (Enformer / Borzoi-style). Predict a
  region's regulatory impact from sequence; no training data, but not
  reprogramming-specific.
- **Hybrid** — LM embeddings as features + a light supervised head on the
  time-course.

Pick based on what time-course data is actually in hand. Everything downstream
(`model.py`) sits behind the stable output contract, so this choice is swappable.

## Output contract (STABLE — eCR_predictor depends on it)

TSV: `chrom, start, end, driver_score` with `driver_score ∈ [0, 1]`, same assembly
as the peaks (**mm10**). See `docs/region_weight_contract.md` and
`ecr_navigator/weights.py`. eCR_predictor's off-target Tier-2 maps its union
regions onto this table via a `weight_fn(region)`.

## Planned layout (not yet built)

```
ecr_navigator/
  io.py         # load ATAC / RNA / peaks
  features.py   # region featurization (dynamics, LM embeddings, ...)
  model.py      # driver/passenger scorer (approach per the decision above)
  weights.py    # read/write the region-weight contract TSV  ← stub exists
navigate.py     # entrypoint
docs/
  region_weight_contract.md
```
