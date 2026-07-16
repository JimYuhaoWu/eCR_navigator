# Run-bundle contract — the unified eCR_navigator output

**Status: v1, frozen 2026-07-16.** The single deliverable of a navigator run, and the
boundary object between **eCR_navigator** (produces) and **eCR_predictor** (consumes).

A run emits **one directory** with three parts:

```
<run_id>/
  manifest.json      # what this run is, the verdict, and how much to trust it
  weights.tsv        # DENSE  — every scored region; eCR_predictor Tier-2 off-target weighting
  nominations.tsv    # SPARSE — the instruction: act on these loci, in this order
```

## Why three parts

The old output was a single per-model score table (`region_weight_contract.md`), which
served two consumers with opposite needs and answered neither well:

| Consumer | Needs | Gets |
|---|---|---|
| `eCR_predictor` off-target Tier-2 (`offtarget.py`) | **dense** — a weight for every background region | `weights.tsv` |
| eCR **target design** (`target.py` → `cli.py`) | **sparse, ranked, confident** — a few loci worth building against | `nominations.tsv` |

More importantly, the score table carried **no verdict**. It looked identical whether the
transition was fib→iN (GET wins, top-1% ~9–10× enriched) or MyoD (GET *anti*-informative,
AUROC 0.288). The nomination policy that decides this
([`validation_summary.md`](validation_summary.md)) lived in `scripts/preflight.py` and in
docs — never in the artifact. **The bundle moves the policy into the deliverable.**

## Two invariants

1. **The navigator resolves the score choice; the predictor never learns which model won.**
   Gate 2 picks GET or signed-Δ per transition; `nominations.tsv` carries a single
   `nomination_score`, and `manifest.json` records its provenance. A future self-trained
   model that wins Gate 2 changes nothing downstream — the same model-agnostic promise the
   benchmark enforces ([`benchmark_spec.md`](benchmark_spec.md)).
2. **Refusal is a first-class output.** A Gate-1 REJECT emits a **valid bundle with zero
   nominations** and a reason — not a confident-looking top-1% of noise. On the v1 panel
   this is the correct answer for 3 of 6 transitions (MyoD, iCM, ETV2).

---

## `manifest.json`

```json
{
  "bundle_version": "1",
  "run_id": "in_gse299923",
  "created": "2026-07-16",

  "transition": {
    "id": "in_gse299923",
    "label": "fibroblast -> induced neuron (Ascl1)",
    "species": "human",
    "assembly": "hg38",
    "state_a": "fibroblast",
    "state_b": "iN",
    "geo": "GSE299923"
  },

  "region_universe": {
    "n_regions": 329984,
    "source": "dataset peaks (fib union iN, MACS)"
  },

  "gate1": {
    "admit": true,
    "pc1_frac": 0.89,
    "coherence_margin": null,
    "n_a": 2,
    "n_b": 2,
    "reasons": []
  },

  "gate2": {
    "primary": "GET",
    "anchor_set": "promoter",
    "n_anchors": 85,
    "delta_auroc": 0.170,
    "delta_auroc_ci": [0.081, 0.264],
    "lr_p": 0.001,
    "lr_coef": 0.67,
    "fold_enrichment_top5pct": 2.46,
    "reason": "driver beats signed-Delta: Delta-AUROC CI[+0.081,+0.264] excludes 0"
  },

  "nomination": {
    "score_source": "GET",
    "score_norm": "rank",
    "top_frac": 0.01,
    "n_nominated": 3300,
    "refused": false,
    "refusal_reason": null
  },

  "direction": {
    "provenance": "input-measured",
    "source": "GET aTPM (atpm_iN - atpm_fib)",
    "norm": "raw"
  },

  "models_run": ["GET", "ChromBERT", "ChromFound", "ATACformer", "EpiAgent"],
  "provenance": { "benchmark_version": "v1", "scorecard": "benchmark/scorecard.tsv" }
}
```

### Field notes

| Field | Rule |
|---|---|
| `transition.assembly` | **Load-bearing.** Must match the peaks *and* the genome FASTA the predictor uses. Mismatched coordinates silently corrupt both weighting and site selection. |
| `gate1` | Verbatim from `scripts/preflight.py` `admissibility()`. `admit:false` ⇒ `nomination.refused:true`. A reliable **reject**, not a reliable admit. |
| `gate2.primary` | `"GET"` \| `"signed-Delta"` (or any future model id). PRIMARY = GET iff ΔAUROC CI excludes 0 **or** (LR p<0.05 **and** `lr_coef > 0`). The `coef > 0` clause is what keeps MyoD's significant-but-**negative** LR from reading as signal. |
| `gate2.*` stats | Everything **measured** on the anchors (including `fold_enrichment_top5pct`) lives here, not in `nomination` — that block holds only what the nominator *decided*. `null` when Gate 2 did not run. `reason` is `select_score`'s human-readable justification for the verdict. |
| `nomination.score_source` | Mirrors `gate2.primary`. The only place the winning model is named. |
| `nomination.score_norm` | `rank` \| `minmax` (`ecr_navigator/model.py`). **Load-bearing — see "rank is the orderable column" below.** |
| `nomination.top_frac` | The **validated confidence band** (0.01 — see below), not a design budget. |
| `direction.provenance` | One of the three tiers in [`direction.md`](direction.md): `input-measured` \| `predicted-model-native` \| `external-attach`. Governs how much to trust the sign. |
| `direction.norm` | `raw` \| `maxabs` \| `signed-rank`. For GET/ChromFound aTPM (already `[0,1]`) use **`raw`**, so `direction` stays a literal ΔaTPM. |

On refusal, every `nomination.*` field except `n_nominated` (0), `refused` (true), and
`refusal_reason` is `null`. `refusal_reason` states **why refused** and nothing more — the
supporting statistics stay in `gate1.reasons` / `gate2`, never duplicated into prose.

**Refusal is decided by Gate 1, not Gate 2.** Gate 2 only chooses *which* score to nominate
from, so `primary: "signed-Delta"` is a normal, non-refusing outcome (MEF→mES: admissible,
but signed-Δ already captures the signal). `nominate()` refuses iff Gate 1 is missing or
rejects, or Gate 2 is missing. It therefore requires **both** gates — deliberately stricter
than `scripts/preflight.py`, whose Gate 1 is optional because it is a diagnostic; nomination
is the production path, and an unverified endpoint pair must not produce targets.

## `nominations.tsv`

Tab-separated, **ranked** (rank 1 = strongest), already cut to `top_frac` of the PRIMARY score.

| Column | Type | Meaning |
|---|---|---|
| `chrom` | str | Assembly per `manifest.transition.assembly` |
| `start` | int | 0-based, BED convention |
| `end` | int | Exclusive |
| `rank` | int | 1-based rank by `nomination_score`, descending |
| `nomination_score` | float | `[0,1]`, from the Gate-2 PRIMARY scorer. **Not** necessarily a `driver_score` — if PRIMARY is `signed-Delta` this is the normalized \|signed-Δ\|. |
| `direction` | float | `[-1,1]` measured signed-Δ; `+1` = should open, `−1` = should close. **Empty = unmeasured** (never collapse to 0.0 — that means measured-flat). |

```
chrom	start	end	rank	nomination_score	direction
chr2	54963892	54964201	1	0.9998	0.7412
chr8	27324671	27324893	2	0.9997	0.6688
```

A **refusal** is the header line and nothing else. That is a valid bundle.

### `rank` is the orderable column, not `nomination_score`
With the default `score_norm: rank`, `driver_score` **is a percentile**, so a `top_frac` band
spans `[1-top_frac, 1.0]` *by construction*. In the iN fixture all 3,300 nominations score
0.99–1.0, with 17 exact ties at 1.0. **Consumers must order and cut by `rank`, never
threshold on `nomination_score`** — within a band it is near-constant and carries no
prioritization signal. `nomination_score` is retained for provenance and for `minmax` runs.

### `direction` is often small — check it before choosing an ED
In the iN fixture `|direction|` has median **0.064**, and **42% of nominations (1,382/3,300)
fall below 0.05**; 51 are exactly 0.0 (measured-flat). This is *expected* and corroborates
Claim 2A — GET's value is prioritizing regulatory regions at **matched** magnitude, so it
deliberately does not just pick the biggest openings ([`validation_summary.md`](validation_summary.md)).
But it means the open/close call is **weakly supported for a large fraction of nominations**,
and `direction` picks the ED (activator vs repressor). Consumers should treat small
`|direction|` as low-confidence rather than a confident call. A principled ambiguity
threshold is **not** set in v1 — see open decisions.

### Why `nomination_score` and not `driver_score`
`weights.tsv` keeps `driver_score` so `offtarget.py` and `read_region_weights()` stay
byte-compatible. `nominations.tsv` is a **new** consumer (`target.py`), and its score may not
be a driver score at all when Gate 2 selects signed-Δ. Naming it `driver_score` there would
be a lie in exactly the case the policy exists to catch.

## `weights.tsv`

**Unchanged** — the existing 4/5-column region-weight contract
([`region_weight_contract.md`](region_weight_contract.md)), covering the full region
universe. Emitted for **every** run, including refusals: a Gate-1 reject still yields usable
off-target weights, because Tier-2 weighting only needs relative accessibility importance,
not a trustworthy driver ranking.

## Downstream: how the predictor consumes this

```
weights.tsv     -> offtarget.py       (Tier-2 weight_fn; already built)
nominations.tsv -> target.py          (NEW: region -> target site)
                     |  extract region sequence (pyfaidx, manifest assembly)
                     |  scan library DBD motifs   [reuse scan.py PWM core]
                     |  rank by on-target strength x specificity  [reuse offtarget.py]
                     v
                   sites.tsv -> cli.py --sequence <site>   -> refine.py -> fuse.py
```

`direction` selects the **ED, not the DBD** (`+` → activator e.g. VP64; `−` → repressor
e.g. KRAB) and flows through to `fuse.py`. Site selection is **library-first**: the best PWM
hit per library DBD, filtered by off-target burden. It lives in **eCR_predictor**, which
already owns the genome FASTA (`offtarget.py` pyfaidx), the PWM core (`scan.py`), and the
module library — so the navigator takes no `eCR_mod_lib` dependency and the region contract
stays coordinates-only.

## Open decisions (deliberately not settled in v1)

- **`top_frac` is a confidence band, not a design budget.** The validated claim is that GET's
  top ~1% is enriched ~9–10× on a strong clean transition. For iN that is ~3,300 regions —
  far more than anyone will build eCRs against. The bundle therefore ships the *ranked*
  band and lets the predictor's own `--top-n` set the budget. Whether the navigator should
  instead hard-cap `n_nominated` is unresolved.
- **Anchor overlap is not a column.** Flagging which nominated regions are known master-TF
  loci would help the predictor prioritize, and Gate 2 already computes the overlap — but it
  risks circularity (those loci are the Gate-2 positives). Left out of v1.
- **No `--assume-primary` escape hatch.** Gate 2 requires curated destination master-TF
  anchors. No anchors ⇒ no Gate 2 ⇒ no nominations, only `weights.tsv`. Falling back to an
  unmeasured top-1% reintroduces exactly the MyoD failure.
- **No direction-ambiguity threshold.** 42% of iN nominations have `|direction| < 0.05`, so
  the ED call is weakly supported there. Whether to emit a `direction_confidence`, or drop
  ambiguous rows, needs a defensible cutoff — and we have no validation of the *direction*
  column itself (that is Claim 2B, deferred: circular for every current model). Left to the
  consumer for now.

## Producing a bundle

One command (`navigate.py`) runs scoring, both gates, and nomination:

```bash
python navigate.py \
    --emb-a get.fib.hg38.npz --emb-b get.iN.hg38.npz \
    --bundle bundles/in_gse299923 \
    --matrix endpoints.matrix.tsv --state-a fib --state-b iN   `# Gate 1` \
    --anchors neural.promoter.bed --signed signed_delta.tsv    `# Gate 2` \
    --transition transition.json --direction-norm raw
```

`--out FILE` (the bare contract TSV) still works and is unchanged; the two modes compose.
Most manifest fields are derived from the artifacts' own metadata — `assembly`, `state_a`/
`state_b`, `model`, `region_universe.source` — so they cannot drift from the run that
produced them. `--transition` supplies only what the artifacts cannot know (`label`, `geo`).
**Assembly is enforced:** two artifacts from different assemblies abort the run rather than
silently emitting corrupt coordinates.

Omitting a gate's inputs is not an error — it produces a **refusal** bundle, since
`nominate()` requires both gates.

## Versioning

`bundle_version` is bumped when field semantics change. It is stamped by `write_bundle()`,
not by the caller: the writer determines the on-disk format, so it is the only thing that
can honestly declare the version.

## Fixtures — `examples/run_bundle/`

Real bundles for the `target.py` handoff, generated from the v1 benchmark artifacts on
PeiLab2 (not mocks). Phase 2 (`nominate.py`) will regenerate these from the pipeline itself.

| Fixture | Role | Contents |
|---|---|---|
| `in_gse299923/` | **GET-wins** (Gate-1 admit, Gate-2 PRIMARY=GET) | 3,300 nominations, hg38, real scores + measured ΔaTPM |
| `myod_gse186271/` | **Refusal** (Gate-1 REJECT, GET anti-informative) | 0 nominations — header only, with `refusal_reason` |

Source: `in_clean/get_in.tsv` + `in_clean/atpm_union.tsv`;
`benchmark/myod_gse186271/{gate1,transition}.json`; Gate-1/Gate-2 statistics carried from
`benchmark/scorecard.tsv`. `weights.tsv` is not duplicated into the fixtures — the full
329,983-region iN contract is the existing `get_in.tsv` on PeiLab2.
