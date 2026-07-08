# AlphaGenome pipeline — driver scores (native mm10 + hg38, directional)

> **STATUS (2026-07-08): SCOPING — feasibility confirmed, scoring pipeline NOT yet
> built.** This doc records the model survey, the make-or-break track check, and the
> plan. The embed/scoring scripts are the next session's work (see "Next steps").

AlphaGenome ([DeepMind, *Nature* 2025](https://www.nature.com/articles/s41586-025-10014-0))
is a sequence→function model: it predicts hundreds of genomic tracks (chromatin
accessibility, expression, TF/histone ChIP, contacts, splicing) at bp resolution from
up to 1 Mb of DNA, for **both human (hg38) and mouse (mm10) natively**. That makes it the
first candidate that can be the **primary mouse** driver track — no liftOver, and a
quantitative, directional accessibility readout that fills the contract's
`direction ∈ [-1, 1]` column (unlike EpiAgent/ATACformer, which are hg38-only and only
a sparse *reference* for mouse — see [`epiagent_pipeline.md`](epiagent_pipeline.md)).

## Runs locally from HF weights (not the API)

Weights are on Hugging Face, so we run inference **on the model-zoo GPU**, not via the
gRPC API:
- PyTorch port: [`gtca/alphagenome_pytorch`](https://huggingface.co/gtca/alphagenome_pytorch)
  — `model_all_folds.safetensors` (0.92 GB), **450M params**, fp32, `pip install
  alphagenome-pytorch`. Weights staged at `/yutiancheng/yuhao/models/AlphaGenome/`.
- The API (`alphagenome` client) is only used **once, off-mirror**, to fetch the track
  ontology (below) — the mirror can't reach Google endpoints.

## GPU: A100 / A800 / H800 — not V100

AlphaGenome is 450M params (weights tiny); the GPU tier is about the 1 Mb-context
activations. **A100 / A800 (sm_80) or H800 (Hopper) recommended**; run it on the A800
(80 GB) instance. **Not V100** — sm_70 has no native bf16 and only 32 GB. For scoring
peaks, 16–100 kb windows keep memory light, so full 1 Mb is not required.

## Make-or-break: mouse cell-state tracks exist (via DNase)

AlphaGenome predicts a **fixed track ontology**; a MEF-vs-mES differential needs mouse
fibroblast + ESC accessibility tracks to exist. Checked with `output_metadata` (dumped
by `scripts/alphagenome_dump_tracks.py`):

- Mouse **ATAC** — 18 tracks, mostly tissues; **no fibroblast/ESC** →
  [`alphagenome_mouse_atac_tracks.tsv`](alphagenome_mouse_atac_tracks.tsv).
- Mouse **DNase** — 67 tracks, **has both** →
  [`alphagenome_mouse_dnase_tracks.tsv`](alphagenome_mouse_dnase_tracks.tsv).

So we use **DNase** (same accessibility readout). Chosen proxies:

| state | track | mouse-DNase idx | ontology |
|---|---|---|---|
| MEF (start) | **NIH3T3** (mouse embryonic fibroblast line) | 12 | EFO:0001222 |
| mES (target) | **ES-E14** (E14 mESC line) | 25 | EFO:0007075 |

## Planned pipeline

For each region in eCR_predictor's mm10 peak union:
1. Extract a window (≥16 kb, centered) of **mm10** sequence, one-hot encode.
2. `model.predict(onehot, organism_index=1)` → `preds['dnase']`.
3. Read the **NIH3T3** and **ES-E14** channels; summarize accessibility over the region.
4. `direction = sign(ES-E14 − NIH3T3)`, `driver_score = |Δ|` (normalized) →
   `alphagenome_driver_scores.mm10.tsv` via the region-weight contract
   ([`region_weight_contract.md`](region_weight_contract.md)). Native mm10 — no lift.

Human transdifferentiation reuses the same path with `organism_index=0` and the human
DNase (or ATAC) ontology.

## Open question to resolve first (next session)

The port outputs the **full 384-channel DNase head** (human+mouse union), but the API
metadata lists mouse's 67 DNase tracks separately (idx 12 = NIH3T3, 25 = ES-E14 within
that 67). **We must pin where those sit in the port's 384 channels** — likely by
combining the human+mouse metadata ordering, then verifying with a test prediction at a
known NIH3T3/ES-E14-accessible locus. Until that mapping is confirmed, the channel
indices above are the *mouse-subset* positions, not the port output positions.

## Next steps
1. Finish `alphagenome-pytorch` + torch install on the A800; load weights; smoke-test
   `predict` on a small mm10 window.
2. Resolve the 384-channel organism/track mapping (above); confirm NIH3T3/ES-E14 columns.
3. Write `scripts/alphagenome_embed_regions.py` (or a direct scorer) + wire the
   differential into the contract TSV.
4. Run MEF→mES on the e7 peak union; report coverage + score distribution.
