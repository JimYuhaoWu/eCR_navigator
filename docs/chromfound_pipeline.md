# ChromFound pipeline — driver scores (hg38; mm10 via liftOver)

ChromFound (NeurIPS 2025, single-cell scATAC foundation model — Mamba + windowed
self-attention, dim 128) is used the same zero-shot way as the other models: embed
each region per cell state, diff the two states, `navigate.py` → driver scores.

Unlike ATACformer, ChromFound is **coordinate-based / universe-agnostic**: it encodes
each open chromatin region (OCR) by chromosome + sin/cos of its hg38 start/end + a
**continuous accessibility scalar**. So any peak set works — we feed our own union
directly, with no foreign region universe to snap to (that is why coverage is much
better than ATACformer). Runs on the `chromfound` conda env; **GPU is mandatory**
(Mamba selective-scan + FlashAttention, no CPU path). Weights at
`/yutiancheng/yuhao/models/chromFound/`, repo at `/root/ChromFound`.

## Role at a glance — reference for mouse, primary for human

Same standing as ATACformer: the model is hg38/human-trained, so it is a **primary**
driver result for human and a **reference** that corroborates the native-mm10 tracks
(GET / ChromBERT) for mouse.

| | Human (hg38) | Mouse (mm10) |
|---|---|---|
| **Bridge** | none — native | liftOver mm10→hg38 (and back) |
| **Coverage** | full | ~42% of the union embeds in hg38; the mm10 back-lift round-trips at **99.8%** (87,615 regions) — near-full, because it's our own peaks, not a foreign universe |
| **Role** | **primary** driver result | **reference** (human-trained model on lifted mouse) — but usably dense, unlike ATACformer |

Note the contrast with ATACformer: ChromFound is coordinate-based, so lifting our own
peaks out to hg38 and the scores back to mm10 loses almost nothing on the return trip
(99.8% vs ATACformer's ~12%). The "reference for mouse" call is about the model being
human-*trained*, not about sparse coverage.

## Steps

### 1. Build the input h5ad (both states)

`chromfound_build_input.py` lifts the union peaks mm10→hg38 (per-region, pyliftover),
keeps vocab chromosomes (chr1–22, X, Y), and writes an h5ad matching ChromFound's
schema — `var[features, #Chromosome (str), hg38_Start (int), hg38_End (int)]`,
`obs[celltype]`, `X` = per-state accessibility (aTPM), one row per state.

```bash
python chromfound_build_input.py \
    --union union_named.bed --atpm atpm_union.tsv --states MEF mES \
    --chain /yutiancheng/yuhao/eCR/refs/mm10ToHg38.over.chain.gz \
    --out chromfound_input.MEF_mES.h5ad
# human (hg38 peaks already): add --no-lift, drop --chain
```

Accessibility scale note: ChromFound trained on raw-ish ATAC counts; we pass the
[0,1] aTPM. Because the readout is the MEF-vs-mES *shift* and navigate.py rank-
normalizes, the global accessibility scale cancels — aTPM is fine here.

### 2. Per-OCR embeddings (GPU)

`chromfound_embed_regions.py` loads the model like the shipped `cell_embedding.py`
but **skips the mean-pool** so it keeps the encoder's per-OCR 128-d vectors (the cell
path pools them into one cell vector). One `.npz` per state, in the shared artifact
contract; OCR order is preserved (no reorder/pad when `max_length = n_OCR`), so row i
↔ var region i.

```bash
cd /root/ChromFound && python chromfound_embed_regions.py \
    --h5ad chromfound_input.MEF_mES.h5ad \
    --ckpt-dir /yutiancheng/yuhao/models/chromFound \
    --assembly hg38 --out-dir chromfound_emb
```

### 3. Diff into driver scores

```bash
python navigate.py --emb-a chromfound.MEF.hg38.npz --emb-b chromfound.mES.hg38.npz \
    --out chromfound_driver_scores.hg38.tsv
```

Stable contract `chrom,start,end,driver_score∈[0,1]`. To carry to mm10 for the mouse
pipeline, liftOver the TSV hg38→mm10 (round-trips at ~99.8% since these are the
original peaks).

## Validated (2026-07-07)

- Env/weights/GPU confirmed (A100-80GB); model loads from `model.pt` (`state["module"]`).
- Input built: 206,313 union peaks → **87,779** hg38 OCRs (42.5% kept after lift+vocab).
- Per-OCR embeddings `(87,779, 128)`, distinct per region; MEF-vs-mES shift nonzero
  across all regions (median 0.556).
- `navigate.py` → `chromfound_driver_scores.hg38.tsv` (87,779 regions, `[0,1]`).
- **mm10 back-lift**: 87,779 → **87,615 (99.8%)** → `chromfound_driver_scores.mm10.tsv`.
- **Sanity vs accessibility change**: driver_score strongly tracks |ΔaTPM| —
  corr **0.87**, mean score **0.92** for |ΔaTPM|>0.3 vs **0.42** for stable regions.
  (Tighter than GET, whose score is deliberately decoupled from raw ΔaTPM.)
- Artifacts staged at `/yutiancheng/yuhao/eCR/artifacts/` (hg38 + mm10 TSVs, both npz).
