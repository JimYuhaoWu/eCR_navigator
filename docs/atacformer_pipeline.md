# ATACformer pipeline — driver scores (hg38; mm10 via liftOver)

ATACformer (databio's Atacformer, `geniml.atacformer` + `gtars`) is a transformer
over ATAC **region tokens**. A "cell" is the *set* of accessible regions drawn from a
fixed universe of **890,704 regions**; the state signal is *which* regions are
present. We embed one state's accessible peak set as one bag-of-regions and read the
encoder's per-region output — then `navigate.py` diffs two states, exactly as for
GET/ChromBERT.

Runs on the `atacformer` conda env of the model-zoo mirror. Model + universe live at
`/yutiancheng/yuhao/models/atacformer/` (see [`server_mirrors.md`](server_mirrors.md)).

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
python navigate.py --a atac.MEF.hg38.npz --b atac.mES.hg38.npz \
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
- **Pending:** the actual MEF/mES driver track (needs the mm10→hg38 liftOver of the peaks).
