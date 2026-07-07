# GET pipeline — driver scores for hg38 or mm10

GET (General Expression Transformer, Nature 2025) is used **zero-shot**: we embed
each peak in cell state A and state B, then `navigate.py` diffs the two embeddings
into a per-region driver score. The pipeline is **species-parameterized** — the only
difference between human and mouse is the `--assembly` flag (which selects the
hg38 vs mm10 Vierstra archetype-motif file). Nothing is hardcoded to mouse.

Runs on the GET mirror (`get` conda env). See connection details in
[`server_mirrors.md`](server_mirrors.md).

## Input GET expects per region

`region_motif` = `[282 motif-archetype scores | 1 aTPM]` = 283 dims.
- **282 motif scores** — max-normalized per column (RegionMotif, `motif_scaler=1`).
- **aTPM** (last column) — normalized per-peak ATAC accessibility in `[0,1]`; this
  is the cell-state channel that makes A vs B embeddings differ.

## Steps

### 1. aTPM (on the ATAC data server, where the bigWigs live) — genome-agnostic

```bash
scripts/compute_atpm.sh MEF MEF.bw MEF_peaks.bed  mES mES.bw mES_peaks.bed  out/
# -> out/atpm_union.tsv (chrom,start,end,atpm_MEF,atpm_mES) + union_named.bed
```

The union peak set (`union_named.bed`) is the shared region list both states are
scored on. Transfer `atpm_union.tsv` + `union_named.bed` to the GET mirror.

### 2 + 3. Motif matrix + embedding (on the GET mirror) — per state

`get_embed_regions.py` builds the motif matrix internally (via
`get_regionmotif_matrix.build_matrix`), appends aTPM, and extracts the encoder's
per-region 768-d output. Run once per cell state:

```bash
# mouse (mm10)
python get_embed_regions.py --peaks union_named.bed --assembly mm10 \
    --atpm-tsv atpm_union.tsv --state MEF \
    --motif-names motif_names.txt \
    --motif-file mm10.archetype_motifs.v1.0.bed.gz \
    --checkpoint .../pretrain_fetal_adult/checkpoint-799.pth \
    --out get.MEF.mm10.npz

# human (hg38) — identical, just the assembly + motif file change
python get_embed_regions.py --peaks union_named.bed --assembly hg38 \
    --atpm-tsv atpm_union.tsv --state fibroblast \
    --motif-names motif_names.txt \
    --motif-file hg38.archetype_motifs.v1.0.bed.gz \
    --checkpoint .../checkpoint-799.pth \
    --out get.fibroblast.hg38.npz
```

**`--motif-file`**: a local tabix-indexed `{assembly}.archetype_motifs.v1.0.bed.gz`
(with its `.tbi` alongside). Strongly recommended — omitting it falls back to the
remote URL, which is fine for a handful of peaks but far too slow for a full peak
set (per-peak HTTP range requests). Download once from
`resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/`
(mm10 ~9.5 GB, hg38 ~12 GB) to persistent storage.

### 4. Diff the two states into driver scores

```bash
python navigate.py --a get.MEF.mm10.npz --b get.mES.mm10.npz \
    --out get_driver_scores.mm10.tsv
```

Produces the stable contract `chrom,start,end,driver_score∈[0,1]` (see
[`region_weight_contract.md`](region_weight_contract.md)), consumed by
eCR_predictor.

## Validated

- 2-peak mm10 smoke test: checkpoint loads `missing=0 unexpected=0`, produces a
  `(2, 768)` artifact with distinct per-region vectors (2026-07-07).
