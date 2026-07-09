# ATACformer pipeline — driver scores (hg38; mm10 via liftOver)

ATACformer (databio's Atacformer, `geniml.atacformer` + `gtars`) is a transformer
over ATAC **region tokens**. A "cell" is the *set* of accessible regions drawn from a
fixed universe of **890,704 regions**; the state signal is *which* regions are
present. We embed one state's accessible peak set as one bag-of-regions and read the
encoder's per-region output — then `navigate.py` diffs two states, exactly as for
GET/ChromBERT.

Runs on the `atacformer` conda env of the model-zoo mirror. Model + universe live at
`/yutiancheng/yuhao/models/atacformer/` (see [`server_mirrors.md`](server_mirrors.md)).

## Role at a glance — reference for mouse, primary for human

The model is native hg38, so its standing in eCR_navigator depends on the species:

| | Human (hg38) | Mouse (mm10) |
|---|---|---|
| **Bridge** | none — native | liftOver mm10→hg38 (and back) |
| **Coverage** | full | ~1% survives the double liftOver (conserved core) |
| **Role** | **primary** driver result | **reference / corroboration** on conserved loci |
| **Primary track instead** | — | GET / ChromBERT (native mm10) |

So: treat ATACformer's output as a **primary result when working in human**, and as a
**reference** that backs up the native-mm10 tracks (GET/ChromBERT) when working in mouse.

## Species — hg38 only, mm10 must be lifted

The pretrained model `databio/atacformer-base-hg38` and its 890k-region universe are
**hg38-only**, unlike ChromBERT (native mm10) or GET (motif-based). So for the mm10
MEF→mES data:

1. **liftOver mm10 → hg38** each state's peak BED (chains are on the instance under
   `eCR/refs`):
   ```bash
   liftOver MEF.mm10.bed mm10ToHg38.over.chain.gz MEF.hg38.bed MEF.unmapped
   liftOver mES.mm10.bed mm10ToHg38.over.chain.gz mES.hg38.bed mES.unmapped
   ```
   Document the unmapped fraction — lifted regions that fail to map are dropped.

Human transdifferentiation data (hg38) needs no bridge — pass its peaks directly.

## Embed each state

`atac_embed_regions.py` snaps a state's peaks to the universe (`tok.tokenize`),
encodes them, and runs the encoder in `--window` chunks to get per-region 192-d
embeddings. **Output coordinates are the universe regions** the peaks snap to, so
both states share coordinates and align with no extra work.

```bash
python atac_embed_regions.py --peaks MEF.hg38.bed --state MEF --assembly hg38 \
    --model-dir /yutiancheng/yuhao/models/atacformer --out atac.MEF.hg38.npz
python atac_embed_regions.py --peaks mES.hg38.bed --state mES --assembly hg38 \
    --model-dir /yutiancheng/yuhao/models/atacformer --out atac.mES.hg38.npz
```

## Diff into driver scores

```bash
python navigate.py --emb-a atac.MEF.hg38.npz --emb-b atac.mES.hg38.npz \
    --out atac_driver_scores.hg38.tsv
```

Produces the stable `chrom,start,end,driver_score∈[0,1]` contract (see
[`region_weight_contract.md`](region_weight_contract.md)). To carry scores back to
mm10 for the mouse pipeline, liftOver the output TSV hg38→mm10.

## Interpretation note

The state signal is region *presence*, not a per-region accessibility value (unlike
GET's aTPM channel). A region open in both states gets a driver score from how its
*contextual* embedding shifts under each state's differing global accessibility
landscape; regions exclusive to one state are driver-by-presence (already captured by
ATAC itself). Report ATACformer alongside GET, don't treat it as a drop-in equal.

**The shift is context-driven, not a positional artifact (validated 2026-07-07).**
The encoder is effectively permutation-invariant: moving a region within the token
list while holding its neighbors fixed leaves its embedding *identical* (L2 0.000).
The embedding changes only when the *set of other accessible regions* changes (same
position, different neighbors → L2 16.6, cos 0.28). So a shared region's MEF-vs-mES
shift reflects a genuine difference in the surrounding accessibility landscape, not
how many peaks happen to precede it. Region identity also contributes but less than
context (different region, same context → L2 7.6, cos 0.85).

## Validated (2026-07-07)

- Model loads from `databio/atacformer-base-hg38` (hidden **192**, 6 layers); 500-region
  hg38 smoke test → `(500, 192)` artifact, distinct per-region, `.npz` contract matches.
- **Context-sensitivity** confirmed (position-invariant; neighbor-driven — see note above).
- **End-to-end**: two overlapping states → `atac_embed_regions.py` ×2 → the real
  `navigate.py` → contract TSV (2,500 shared regions, `driver_score ∈ [0,1]`, rank-norm).
- Throughput ~ hundreds of regions/s on one A100 (GPU); a full ~100k-region state is minutes.
- **Real MEF/mES run done (2026-07-07):** MEF/mES e7 peaks lifted mm10→hg38, embedded,
  diffed → **11,059 shared hg38 driver regions** (`artifacts/atacformer_driver_scores.hg38.tsv`
  on the instance). Valid, and the right native deliverable for a *human* run.
- **Native-hg38 human run done (2026-07-09) — kidney vs pancreas** (the primary case):
  per-state peaks embedded directly, no liftOver. kidney 87,522 peaks → **66,573** universe
  regions (76%); pancreas 129,874 → **108,868** (84%); `navigate.py` → **40,753** shared
  driver regions (`atacformer_driver_scores.kidney_pancreas.hg38.tsv`, `[0,1]`, all 24
  chroms, top drivers autosomal chr18/chr19). This is the coverage ATACformer was built
  for — dense (40,753) vs the mouse liftOver's 1,334 — confirming "primary for human."
  Staged at `/yutiancheng/yuhao/eCR/artifacts/`.

## Cross-species coverage is the real limit (mouse)

Because Atacformer is hg38-only, the mouse track pays two liftOvers, and both are
lossy for regulatory regions (non-conserved enhancers don't map):

| step | MEF | mES |
|------|-----|-----|
| mm10 e7 peaks | 125,467 | 128,790 |
| → hg38 (liftOver) | 50,813 (40%) | 32,017 (25%) |
| shared hg38 driver regions | \multicolumn{2}{c}{11,059} | |
| → mm10 (liftOver back) | \multicolumn{2}{c}{**1,334 (12%)**} | |

Net: ~206k mm10 union peaks → **1,334 mm10 driver scores (<1%)**, and the survivors are
the *conserved core* — the regions least likely to be state-specific drivers. So for the
**mouse** MEF→mES task, ATACformer is a sparse, conservation-biased signal: use it only
as corroboration on conserved loci, and rely on **GET / ChromBERT (native mm10)** for the
primary mouse driver track. ATACformer comes into its own for the **planned human
transdifferentiation**, where it is native hg38 and keeps full coverage (the 11k-region
hg38 track above is the shape of that deliverable). Do not liftOver the mouse track back
to mm10 for production use — report the hg38 track and treat the mm10 back-projection as a
conserved-subset sanity check only.
