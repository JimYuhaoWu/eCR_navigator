# Evo2 scoping — sequence-intrinsic prior, not a two-state track

> **STATUS (2026-07-08): SCOPING ONLY — feasibility recorded, nothing built.**
> Evo2 is not installed on any mirror. This note captures what it is, why it plays a
> *different* role than the other models, and what it would take to build. Revisit
> after AlphaGenome.

Evo 2 ([Arc Institute + NVIDIA, 2025](https://arcinstitute.org/news/evo);
[ArcInstitute/evo2](https://github.com/ArcInstitute/evo2)) is a DNA **language model**
— StripedHyena-2 (convolution + linear-attention + state-space), **7B and 40B**
params, trained on **9.3T bp across all domains of life**, context up to **1 Mb**.
Outputs: per-position **embeddings** (intermediate layers via
`evo2_model(input_ids, return_embeddings=True, layer_names=[...])`), per-position
**likelihood/logits** (→ zero-shot variant/mutation effect by Δlog-likelihood), and
sequence **generation**.

## Why Evo2 is NOT a sixth two-state driver track

Every integrated model (ChromBERT, GET, ATACformer, ChromFound, EpiAgent) reads cell
state from an **accessibility / region-token input** and diffs two states into a
driver score. **Evo2 takes only raw DNA.** Two states of the same species share the
same genome, so Evo2's embedding/likelihood for a region is **identical in both
states** — the embedding-shift readout yields exactly zero signal. Evo2 alone cannot
separate a driver from a passenger *between two states*, because it never sees the
state.

## The role it *can* play — the sequence-intrinsic leg

`CLAUDE.md`'s open design decision names three legs: supervised-trajectory /
zero-shot genomic LM / hybrid. Evo2 is the premier **zero-shot genomic LM**, but it
contributes a **state-invariant** signal that *complements* the accessibility-shift
tracks rather than replacing them:

| Use | Evo2 output | Role in the navigator |
|---|---|---|
| **Constraint / importance prior** | sequence likelihood or embedding-derived score | a **multiplier** on the accessibility-shift driver score — a region that *changes* accessibility **and** sits in high-constraint regulatory sequence is a more credible driver than one in junk |
| **In-silico mutagenesis** | Δlikelihood under mutation | regulatory-element / motif importance within a region (most causal, still sequence-intrinsic) |
| **Embedding backbone** | per-region embeddings | features for a **supervised/hybrid head** once time-course data exists |

**Multi-species advantage.** Being sequence-only, Evo2 handles **mm10 and hg38
natively** — no liftOver, no coverage loss. It is the cleanest multi-species model
surveyed, and directly satisfies the first-class multi-species requirement.

## Blockers

1. **Not installed** — no env/weights on any mirror. Needs `vortex`/StripedHyena
   (+ Transformer Engine for 20B/40B) and ~14–28 GB of 7B weights staged to
   `/yutiancheng`.
2. **Hardware.** The model-zoo instance is a **V100-32GB (sm_70, no native bf16)** — a
   poor/likely-unworkable host. Per the repo: **7B** needs bf16 on a supported GPU
   (realistically an **A100**); **40B** needs **FP8 + Hopper (H100)**, multi-GPU. The
   realistic target is **Evo2-7B on an A100-80GB** (as used for ChromBERT/ChromFound).

## If/when we build it (next steps)

Unlike the drop-in `<model>_embed_regions.py` scripts, Evo2 changes how scores
*compose*, not just adds another `.npz`:

1. Provision an A100; install `evo2` (7B) to `/yutiancheng/yuhao/models/evo2/`.
2. Score the peak-union sequences (mm10 and hg38) → per-region constraint scalar +
   banked embeddings.
3. Add a **combiner** step in `ecr_navigator/model.py`: `driver_score *= f(evo2_prior)`
   (keep the raw accessibility-shift score available; the prior is an optional
   modulator, off by default so existing runs are unchanged).
4. Later: a supervised head on the banked embeddings when time-course arrives.
