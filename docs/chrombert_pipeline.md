# ChromBERT pipeline — driver scores (hg38 + mm10, both native)

ChromBERT (context-specific transcription-regulatory-network foundation model) is
used **zero-shot**: embed each peak's 1 kb locus as a TRN-context vector in cell
state A and state B, then `navigate.py` diffs the two embeddings into a per-region
driver score. It was the first model wired into the navigator; the run scripts live
in [`scripts/`](../scripts) (`run_chrombert_region_emb.sh` + `hdf5_to_artifact.py`).

## Role at a glance — primary for BOTH species

Unlike the hg38-only scATAC models, ChromBERT ships **both assemblies natively**
(hg38 = 6k regulators, mm10 = 5k), selected by a `-g` flag — **no liftOver**. So it
is a **primary** driver track for mouse *and* human, with full native coverage on
each. See [`server_mirrors.md`](server_mirrors.md) for mirror access and the cached
model/data layout.

## Input

One peak **BED per cell state** (3-col `chrom start end` is enough). `make_dataset`
overlaps the peaks onto ChromBERT's fixed **1 kb grid**; `get_region_emb` emits a
**768-d** TRN-context embedding per grid region. Coordinates are snapped to 1 kb, so
both states share region keys and `navigate.py` aligns them directly.

## Steps

```bash
# one call per cell state (runs on the ChromBERT mirror; port may move -> MIRROR_PORT)
MIRROR_PORT=35963 ./scripts/run_chrombert_region_emb.sh \
    --peaks kidney.MergedPeaks.bed --genome hg38 --cell-state kidney \
    --out artifacts/chrombert.kidney.hg38.npz
MIRROR_PORT=35963 ./scripts/run_chrombert_region_emb.sh \
    --peaks pancreas.MergedPeaks.bed --genome hg38 --cell-state pancreas \
    --out artifacts/chrombert.pancreas.hg38.npz
# mouse is identical with --genome mm10 (native; the CLI mask is passed explicitly
# per genome inside the script, since its auto-inference defaults to hg38's mask).
```

Each call does, on the mirror: `chrombert_make_dataset` → `chrombert_get_region_emb`
(GPU) → fetch hdf5+tsv → `hdf5_to_artifact.py` writes the `.npz` via the shared
[`embedding_artifact.py`](../scripts/embedding_artifact.py) writer.

### Diff into driver scores

```bash
python navigate.py --emb-a artifacts/chrombert.kidney.hg38.npz \
    --emb-b artifacts/chrombert.pancreas.hg38.npz \
    --out artifacts/chrombert_driver_scores.kidney_pancreas.hg38.tsv
```

Stable contract `chrom,start,end,driver_score∈[0,1]`
([`region_weight_contract.md`](region_weight_contract.md)).

## Validated — both species

**Mouse (mm10, native).** Region-embedding pipeline built and verified end-to-end on
mm10 MEF/mES peaks (smoke + subset runs → valid mm10 driver-weight TSV); the
genome-blind-mask bug was fixed by passing `mm10_5k_mask_matrix.tsv` explicitly. See
memory `ecr-extensions-roadmap` for the running notes.

**Human (hg38, native) — 2026-07-09.** Full kidney/pancreas ENCODE MergedPeaks run on
an A100-80GB, end-to-end with zero friction:

| step | kidney | pancreas |
|------|--------|----------|
| input peaks | 87,522 | 129,874 |
| → 1 kb grid regions embedded | 79,936 × 768 | 129,190 × 768 |

`navigate.py` intersects the two states → **42,305** shared driver regions.

Output `chrombert_driver_scores.kidney_pancreas.hg38.tsv` is a clean contract TSV
(`driver_score ∈ [0,1]`, rank-norm, all 24 chromosomes). The shared writer + strict
`allow_pickle=False` loader consumed human data unchanged — confirming the pipeline
is genuinely species-agnostic.

> **Caveat for this specific test pair:** the top-scoring regions concentrate on
> **chrY**, which almost certainly reflects a **donor-sex difference** between the
> kidney and pancreas samples (chrY accessibility swings hardest between a male and
> female sample), not tissue-identity biology. For a real driver analysis, sex-match
> the two states or drop chrX/chrY. This is a property of the test data, not the model.

Artifacts staged at `/yutiancheng/yuhao/eCR/artifacts/` (both npz + the driver TSV),
same as the mouse runs.
