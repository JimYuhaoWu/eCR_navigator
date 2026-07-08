# EpiAgent pipeline — driver scores (hg38; mm10 via liftOver)

EpiAgent ([xy-chen16/EpiAgent](https://github.com/xy-chen16/EpiAgent), *Nat. Methods*
2025) is a transformer over cCRE **tokens**. A "cell" is its *set* of accessible
cCREs, drawn from a fixed universe of **1,355,445 hg38 cCREs**, ranked by TF-IDF into
a "cell sentence"; the state signal is *which* cCREs are accessible and how they
rank. We embed one state's accessible cCREs and read the encoder's per-cCRE
contextual output (skipping the CLS cell-pool) — then `navigate.py` diffs two states,
exactly as for GET/ChromBERT/ATACformer/ChromFound.

Runs on the `EpiAgent` conda env of the model-zoo mirror. Weights + refs live under
`/yutiancheng/yuhao/` (see [`server_mirrors.md`](server_mirrors.md)).

## Role at a glance — reference for mouse, primary for human

EpiAgent is native hg38 **single-cell** ATAC, so its standing depends on the species
*and* on the bulk→single-cell bridge:

| | Human (hg38) | Mouse (mm10) |
|---|---|---|
| **Species bridge** | none — native | liftOver mm10→hg38 (and back) |
| **Data bridge** | real scATAC → native | bulk peaks → one pseudobulk "cell" |
| **Coverage** | full | double liftOver + 8190-cCRE rank cap (sparse) |
| **Role** | **primary** driver result | **reference / corroboration** on conserved loci |
| **Primary track instead** | — | GET / ChromBERT (native mm10) |

So treat EpiAgent's output as a **primary result in human** (real scATAC, native
hg38) and as **reference** backing the native-mm10 tracks (GET/ChromBERT) in mouse.
The mouse run below is the *reference* deliverable, matching ATACformer's framing
([`atacformer_pipeline.md`](atacformer_pipeline.md)).

## Two bridges for bulk mouse data

EpiAgent wants single-cell, hg38 input; our MEF/mES inputs are bulk, mm10 peaks:

1. **Pseudobulk** — each state becomes ONE "cell" whose accessible cCREs are those
   its peaks overlap (barcode = the state name). A real single cell has ~2–8k
   accessible cCREs; a bulk pseudobulk has ~200k, so most are dropped at the rank cap
   (below). This is out-of-distribution for EpiAgent — read the mouse scores as a
   coarse reference, not a calibrated single-cell result.
2. **Species** — mm10 peaks are lifted to hg38 (EpiAgent's cCRE vocabulary is hg38)
   before the cCRE intersect. Human hg38 peaks skip this (`--no-lift`).

## Build the input (one h5ad per state)

`epiagent_build_input.py`: liftOver → tag barcode → `bedtools intersect` vs
`cCRE.bed` → `construct_cell_by_ccre_matrix` → `global_TFIDF` → `tokenization`. The
`.h5ad` is one pseudobulk cell × 1,355,445 cCREs, carrying `obs['cell_sentences']`
and the cCRE coordinates as `var_names`. liftOver is `pyliftover` (the UCSC binary
needs a newer glibc than the mirror has).

```bash
cd /yutiancheng/yuhao/eCR; R=refs
# mouse (reference): lift mm10 -> hg38
python epiagent_build_input.py --peaks peaks/MEF.e7_peaks.bed --cell-state MEF \
    --ccre-bed $R/cCRE.bed --ccre-freq $R/cCRE_document_frequency.npy \
    --chain $R/mm10ToHg38.over.chain.gz --out artifacts/epiagent.MEF.h5ad
python epiagent_build_input.py --peaks peaks/mES.e7_peaks.bed --cell-state mES \
    --ccre-bed $R/cCRE.bed --ccre-freq $R/cCRE_document_frequency.npy \
    --chain $R/mm10ToHg38.over.chain.gz --out artifacts/epiagent.mES.h5ad
# human (primary): peaks already hg38 -> add --no-lift, drop --chain
```

Observed mm10→hg38 mapping on the e7 peaks: **MEF 28.9%** (36,202/125,467),
**mES 19.5%** (25,136/128,790) — narrow MACS peaks map poorly across species, and the
loss is *asymmetric*, which biases a MEF-vs-mES differential. Quantify it per run.

## Embed each state

`epiagent_embed_regions.py` loads a state's h5ad, feeds `[CLS] + sentence[:8190] +
[SEP]`, and keeps `transformer_outputs[0, 1:1+L, :]` — the per-cCRE contextual 512-d
rows. cCRE token `t` maps to universe row `t-4`, whose coordinate is `var_names[t-4]`,
so both states share coordinates and `navigate.py` aligns them with no extra work.

```bash
python epiagent_embed_regions.py --h5ad artifacts/epiagent.MEF.h5ad --state MEF \
    --ckpt /yutiancheng/yuhao/models/EpiAgent/pretrained_EpiAgent.pth \
    --assembly hg38 --out artifacts/epiagent.MEF.hg38.npz
python epiagent_embed_regions.py --h5ad artifacts/epiagent.mES.h5ad --state mES \
    --ckpt /yutiancheng/yuhao/models/EpiAgent/pretrained_EpiAgent.pth \
    --assembly hg38 --out artifacts/epiagent.mES.hg38.npz
```

`--use-flash-attn auto` picks by GPU: FlashAttention needs Ampere (A100 sm_80+); on
the V100 (sm_70) it uses flash-attn's non-flash fallback (same weights, slower).

## Diff into driver scores

```bash
python navigate.py --emb-a artifacts/epiagent.MEF.hg38.npz \
    --emb-b artifacts/epiagent.mES.hg38.npz --out artifacts/epiagent_driver_scores.hg38.tsv
```

Produces the stable `chrom,start,end,driver_score∈[0,1]` contract
([`region_weight_contract.md`](region_weight_contract.md)). To carry scores back to
mm10 for the mouse pipeline, liftOver the output TSV hg38→mm10 (chain in `eCR/refs`);
expect heavy additional loss (conserved core only), same as ATACformer.

## Coverage is the real limit (mouse)

EpiAgent's positional (rank) embedding caps a cell sentence at **8190 cCREs**, so of a
pseudobulk state's ~200k accessible cCREs only the top-8190 by TF-IDF are embedded
(deterministic; `is_random` is NOT used). `navigate.py` then scores the *intersection*
of the two states' top-8190 sets. Stacked on the asymmetric mm10→hg38 loss and the
hg38→mm10 back-lift, the mouse driver set is small and conservation-biased — corroborate
with it, don't rely on it. EpiAgent comes into its own for the **planned human
transdifferentiation** (real scATAC, native hg38, full coverage).

## Validated — real MEF→mES run (2026-07-08)

End-to-end with the real `pretrained_EpiAgent.pth` (5.8 GB; loads into
`use_flash_attn=False`, 0 missing keys) on a V100. The coverage funnel shows why the
mouse track is a sparse reference, not a primary result:

| step | MEF | mES |
|------|-----|-----|
| mm10 e7 peaks | 125,467 | 128,790 |
| → hg38 (liftOver) | 36,202 (28.9%) | 25,136 (19.5%) |
| → accessible cCREs | 227,354 | 191,734 |
| → embedded (rank cap 8190) | 8,190 | 8,190 |
| shared hg38 driver regions | 2,285 | |

`navigate.py` produced **2,285** shared hg38 driver regions
(`artifacts/epiagent_driver_scores.hg38.tsv`, `driver_score ∈ [0,1]`, rank-norm) — the
right *human*-shaped deliverable, sparse for mouse. The 8190-cCRE rank cap makes it
sparser than ATACformer's 11k. Report the hg38 track; treat the mm10 back-projection as
a conserved-subset sanity check only.

## Readout choice (per-cCRE embedding shift vs predicted accessibility)

This pipeline uses the **per-cCRE contextual embedding shift** (like the other
models). EpiAgent also exposes `signal_decoder(cell_embedding)` — a predicted
accessibility over *all* 1.35M cCREs from a single forward, which would give full
coverage and a signed open/close direction (the `direction` column of the contract).
That is a natural alternative for the human/primary case; it is not wired here to keep
EpiAgent consistent with the shared embedding-shift path.

**Caveat to validate (unlike ATACformer):** EpiAgent adds a **rank** embedding, so a
cCRE's contextual vector depends on its TF-IDF rank position, which differs between
states. A shared cCRE's MEF-vs-mES shift therefore mixes genuine context change with
rank-position change. Run the ATACformer-style context-sensitivity check (hold
neighbors, vary position) before trusting the magnitude of the shift.

## Provisioning notes

- Env fixes on the `EpiAgent` conda env (persist on `/yutiancheng`): `transformers==4.44.2`
  (4.57 breaks `epiagent.model` import under torch 2.0.1), `numpy==1.26.4`, plus
  `pyliftover`. bedtools = static binary at `/yutiancheng/yuhao/bin/bedtools`.
- Refs at `/yutiancheng/yuhao/eCR/refs/`: `cCRE.bed`, `cCRE_document_frequency.npy`,
  `mm10ToHg38.over.chain.gz`, `hg38ToMm10.over.chain.gz`.
- Checkpoint: `pretrained_EpiAgent.pth` (5.8 GB) is Google-Drive-only; the mirror
  can't reach Drive, so fetch on a workstation and scp to
  `/yutiancheng/yuhao/models/EpiAgent/`.
