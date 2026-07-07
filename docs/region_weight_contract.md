# Region-weight contract

The stable interface between **eCR_navigator** (produces) and **eCR_predictor**
(consumes, off-target Tier-2). The model behind the scores can change freely as
long as this format holds.

## Format

A tab-separated file, one row per scored region:

| Column | Type | Meaning |
|---|---|---|
| `chrom` | str | Chromosome, matching the genome assembly of the peaks (mm10 for the mouse work, hg38 for human) |
| `start` | int | 0-based start (BED convention) |
| `end` | int | Exclusive end |
| `driver_score` | float | Driver importance in **[0, 1]** — higher = more driver-like |

```
chrom	start	end	driver_score
chr1	3062469	3062724	0.81
chr1	3119611	3119740	0.12
```

## Contract rules

- **Assembly must match** the ATAC peaks and the genome FASTA used by
  eCR_predictor for that run (currently mm10 for the mouse pipeline; hg38 for human).
  Mismatched coordinates silently corrupt the weighting.
- `driver_score` is normalized to `[0, 1]`. eCR_predictor multiplies each
  off-target site's strength by the score of the region it falls in.
- Regions should cover (at least) the union of accessibility peaks eCR_predictor
  uses as its off-target background. Any eCR_predictor region with no overlap
  here falls back to a default weight (its Tier-1 dynamics weight).
- Overlap resolution (when an off-target region maps to several scored rows) is
  eCR_predictor's responsibility (e.g. max or mean); eCR_navigator just provides
  the scored regions.

## Consumption (eCR_predictor side)

`ecr_predictor.offtarget.score_offtarget(..., weight_fn=...)` takes a
`weight_fn(region) -> float`. The Tier-2 adapter builds that function from this
table (interval lookup), replacing the Tier-1 dynamics weight. No change to the
scoring loop itself.
