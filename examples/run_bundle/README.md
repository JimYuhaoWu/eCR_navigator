# Run-bundle fixtures

Boundary fixtures for **eCR_predictor**, per [`../../docs/run_bundle_contract.md`](../../docs/run_bundle_contract.md).
Everything here is written by the real `ecr_navigator` writers (`write_region_weights`,
`write_nominations`), so a consumer testing against these is testing against the actual
producer — not against its own idea of the format.

| Fixture | What it is |
|---|---|
| `in_gse299923/` | **GET-wins bundle** (Gate-1 admit, Gate-2 PRIMARY=GET). Real hg38 data. |
| `myod_gse186271/` | **Refusal bundle** (Gate-1 reject). Real mm10 data. |
| `_encoding_edge_cases/` | **Synthetic.** Pins format details real data doesn't reach. |

## The two real bundles

Both `weights.tsv` are **slices**, not full universes — the real ones are 329,983 (iN) and
232,788 (MyoD) rows and live on PeiLab2. The slices are real rows straight from the writer:

| | rows | nominations | measured-flat `0.0` | adjacent pairs |
|---|---|---|---|---|
| `in_gse299923` | 9,300 | all 3,300 present | 141 | 1,846 |
| `myod_gse186271` | 9,000 | **0 (refused)** | 122 | 183 |

The iN slice contains **every** nominated region plus each one's coordinate neighbours plus
random background, so `nominations.tsv` is a strict subset of `weights.tsv` — the mapping a
Tier-2 adapter has to perform.

**The refusal bundle is the point.** `myod_gse186271/` ships **0 nominations and 9,000
weights**, which is the contract's claim made concrete: a Gate-1 reject still yields usable
off-target weights, because Tier-2 weighting needs relative accessibility importance, not a
trustworthy driver ranking. On a refused transition `weights.tsv` is the *only* useful output.

## Two facts about real bundles worth knowing

**1. Our regions never overlap each other.** The region universe is a *merged* union, so
adjacent-but-disjoint is the densest it gets: the full iN contract has **0 overlapping
neighbour pairs and 124,485 adjacent pairs** (gap ≤200bp). A consumer's multi-overlap
resolution therefore fires when *its* region spans several of **our adjacent rows** — never
because our rows overlap one another. Both slices preserve adjacency runs.

**2. No real bundle has ever contained an unmeasured region.** Every GET contract in the
benchmark (iN, MyoD, C/EBPα, ETV2) has a measured direction for **every** row — 0 empties,
though all four do contain measured-flat `0.0` rows. The unmeasured case is permitted by the
contract and explained in [`../../docs/direction.md`](../../docs/direction.md) (GET's region
set is the motif matrix, independent of the aTPM table, so a region absent from that table
would be unmeasured), but in practice it has not occurred. It is more likely to arise for the
**external-attach** tier (ATACformer/ChromBERT), where a region with no overlapping aTPM
interval is genuinely unmeasured.

That is why `_encoding_edge_cases/` exists: the empty-vs-`0.0` distinction is
contract-critical and a consumer must be able to test its reader against our writer's
encoding, but manufacturing a fake unmeasured row inside a real bundle would misrepresent the
data. So the edge cases are separate and explicitly synthetic.

## `_encoding_edge_cases/`

**Synthetic — not a transition, not biology.** Writer output, for encoding tests only.

- **`weights.tsv`** — the 5-column form. Row 2 has a measured `0.0` (accessible but flat) and
  row 3 has an **empty** `direction` (unmeasured). Collapsing those two is the failure this
  fixture exists to catch: an unmeasured region must never be reported as an open/close call.
  Also includes adjacent rows and a second chromosome.
- **`weights.no_direction.tsv`** — the **4-column** form, emitted by any model that resolves
  importance but not direction. This is not exotic: **every ChromBERT contract in the
  benchmark is 4-column.** A consumer must handle a missing `direction` column entirely, which
  is distinct from an empty field in a 5-column file.
