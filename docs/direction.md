# The `direction` channel — signed open/close across models

`driver_score ∈ [0,1]` (embedding-shift **magnitude**) says *how much* a region
reorganizes between two states; it is unsigned by construction (a norm can't tell
"opens" from "closes"). The optional `direction ∈ [-1,1]` column adds the **sign**
(+1 = should open, −1 = should close). The two are separate systems joined by row
order (see [`region_weight_contract.md`](region_weight_contract.md)).

Mechanically, direction always comes from a per-region **scalar accessibility signal**
attached to each state's embedding artifact (`signal`, see
[`embedding_artifact.md`](embedding_artifact.md)); `navigate.py` differences the two
states (`signal_b − signal_a`) into `direction`. What differs per model is **where that
signal comes from**, and that provenance is exactly what determines how much to trust it.

## Three provenance tiers

| Tier | Models | `signal` source | Trust |
|---|---|---|---|
| **input-measured** | GET, ChromFound *(both validated 2026-07-10; identical sign-split on the shared union)* | the **measured accessibility the model conditioned on** (aTPM / continuous accessibility), emitted straight from the input, aligned to the embedding by construction | highest — real measured data, native resolution, no overlap mapping |
| **predicted-model-native** | EpiAgent | the model's **own accessibility head** (`sigmoid(signal_decoder(cell_embedding))`, a predicted P(accessible)) | principled (a head *trained* to reconstruct accessibility) but **zero-shot, unvalidated** against measured Δ |
| **external-attach** | ATACformer, ChromBERT *(both validated 2026-07-10, 100% aTPM coverage)* | a measured accessibility track (GET's aTPM) **overlapped onto** the model's regions post hoc by `scripts/attach_measured_signal.py` — the model itself has no accessibility readout | real data, but bolted on by coordinate overlap (coverage gaps → unmeasured) |

`driver_score` (magnitude) is **unaffected** by all of this — it is always the
model's own embedding shift. Only the sign is sourced as above.

### Not a tier: embedding-synthesized direction
A direction *synthesized from the embedding itself* (e.g. projecting the shift onto some
axis) for a model with no accessibility readout would be a **navigator-side invention,
not validated** — we do **not** do this. Every direction we emit traces to measured or
predicted-accessibility data, per the table.

## Unmeasured vs measured-zero (important)

A region can be **measured and ~closed** (accessibility ≈ 0 → a real "flat", `direction`
near 0) or **unmeasured** (no accessibility value at all → `direction` left **empty**, the
region is dropped from the open/close tally). Collapsing these fabricates a sign — see the
Finding-1 fix. How each tier handles it:

- **external-attach** (ATACformer, ChromBERT): a region with no overlapping aTPM interval
  is genuinely **unmeasured** → `map_signal` returns `NaN` → empty `direction`.
- **GET** (input-measured): GET's region set is the motif matrix, *independent* of the
  aTPM table, so a region missing from the table is **unmeasured** → emitted as `NaN`.
  (Model input stays 0.0-filled; only the direction signal is `NaN`.)
- **ChromFound** (input-measured): the region set **is** the peak union the accessibility
  was built from, so a per-state absence is a **measured-low** (≈closed) `0.0`, not
  unmeasured — defensible for peak-union accessibility. ChromFound therefore does not flag
  truly-unmeasured regions; if that is ever needed, have `chromfound_build_input.py` emit a
  measured-mask layer and mask the signal in the embed step.
- **EpiAgent** (predicted): the SR head is defined for every emitted cCRE, so there is no
  unmeasured case.

## Normalization (`--direction-norm`)

`navigate.py` scales the signed delta to `[-1,1]` (the **sign is always kept**):
- `raw` — clamp Δ to `[-1,1]`; use when the signal is already a probability / `[0,1]`
  accessibility (EpiAgent SR, GET/ChromFound aTPM) so `direction = ΔaTPM` stays
  interpretable.
- `maxabs` — Δ / max|Δ|; any signal scale, preserves relative magnitudes (default).
- `signed-rank` — sign(Δ) · percentile-rank(|Δ|); robust to a cross-state scale mismatch,
  and the safest choice when the two states' accessibility may not be perfectly
  co-normalized.

**Cross-state scale confound (input-measured tier):** direction is `atpm_B − atpm_A`. If the
two states' accessibility were normalized independently, the Δ carries a batch/scale
artifact. Difference values that went through the **same** normalization (ideally the same
union quantification), and prefer `signed-rank` if unsure — the sign survives a monotonic
scale mismatch even when the magnitude does not.

**Channel coupling (input-measured tier):** because GET and ChromFound *condition on*
accessibility, their embedding shift (`driver_score`) already partly encodes the
accessibility change, so `direction` is less independent new information than for a
token-only model like EpiAgent. Still a useful explicit sign — just don't treat the two
channels as independent evidence for these models.

## Species

Direction is computed in the run's native assembly (hg38 for the scATAC models, mm10-native
for ChromBERT/GET), then the 5-column contract TSV is liftOver'd hg38→mm10 for the mouse
pipeline exactly like `driver_score` — direction is per-region, so lifting preserves it
(minus dropped rows).

## Producing it

- **GET / ChromFound** — the embed script emits `signal` from the aTPM/accessibility input;
  run `navigate.py --direction auto` (picks it up when both states carry `signal`).
- **EpiAgent** — the embed script emits `signal` from the SR head; same `navigate.py` call.
- **ATACformer / ChromBERT** — no accessibility input, so attach one first:
  `scripts/attach_measured_signal.py --artifact <state>.npz --intensity atpm_union.tsv
  --value-col atpm_<state> --out <state>.sig.npz` per state, then `navigate.py`.
