# Pre-registration — navigator nominations as a prospective ATAC readout

**Status:** draft 2026-11-17, **not frozen.** Binding only once Section 7 is
complete and committed. Written before any construct ATAC data exists.

**Claim under test.** eCR_navigator's region nominations predict *where an eCR
construct actually opens chromatin* — prospectively, on data that did not exist
when the nominations were made. Every published navigator result so far is
retrospective (matched-background AUROC against known master-TF loci). This is
the first test that can fail for the right reason.

**System.** Mouse MEF → mES. Construct: Nanog (or its extracted 72-aa
homeodomain, 90–161) fused to an optimised SS18-SNH effector domain.

---

## 1. Why both rankings are registered

On MEF → mES the navigator's own Gate 2 selects **signed-Δ, not GET**:
ΔAUROC +0.080, CI [−0.037, +0.198], LR p=0.107. The 870 nominations in
`output/bundles/mef_mes_gse201577/` are ranked by **|direction|**, not by
`driver_score`. A readout validating that bundle validates **the baseline** and
says nothing about the model claim.

*Measured channel, verified:* `direction` **is** ΔaTPM, mES − MEF. The bundle
itself only records `direction.source = "GET signal (mESCS - MEFS)"` with
`provenance: "input-measured"`, which does not name the channel — but
`datasets/mef_mes_clean_GSE201577/atpm_union.tsv` (86,956 regions, columns
`atpm_MEF` / `atpm_mES`) reproduces it: across all 86,956 regions,
`atpm_mES − atpm_MEF` equals `direction` to within the 4-decimal rounding step,
with no exceptions. Sign convention matches the manifest string.

An intermediate draft of this section claimed no aTPM artifact existed on the
server and that the manifest did not support |ΔaTPM|. **Both were wrong** — the
artifact is on disk and the identity holds exactly. The only true part was the
narrow point that the manifest does not itself name the channel, which this
paragraph now closes.

So both rankings are frozen, with equal status, before the ATAC.

**The two lists are nearly disjoint.** Top-870 by GET vs top-870 by signed-Δ
over the same 86,956 regions share **39 regions — Jaccard 0.023**. They are not
two orderings of one candidate pool; they point at almost entirely different
chromatin. One ATAC experiment can therefore separate them cleanly, which is
what makes this worth pre-registering rather than reporting post hoc.

Either outcome is reportable:

- GET wins prospectively where it lost retrospectively → the anchor-set test was
  underpowered and the model carries information the baseline misses.
- signed-Δ wins again → confirms what the docs already say, on stronger evidence
  than an anchor AUROC.

### How the GET ranking was produced

Not by re-running nomination. `nominate.py` deliberately has **no
`--assume-primary` fallback** — defaulting to an unmeasured GET top-1% would
reintroduce exactly the MyoD failure — and that guard should not be subverted to
manufacture a bundle.

Instead the ranking is a deterministic cut of the already-frozen `weights.tsv`,
which carries `driver_score` for all 86,956 regions: k = ceil(0.01 × 86,956) =
**870**, descending `driver_score`, stable sort so ties keep input order (the
same rule `nominate.py` uses — the stable sort is at `nominate.py:91`, and k at
`nominate.py:95`). Kept score range 1.0000 … 0.9900.

**The cutoff is tied, and the tie-break is arbitrary.** Exactly 8 regions in the
universe carry `driver_score == 0.9900`; **4 are kept and 4 are dropped**, decided
purely by row order in `weights.tsv`. This is recorded here rather than left
implicit, because a frozen list whose boundary was set by file order is a fact a
reader of the pre-registration is entitled to.

One further nuance: `nominate.py`'s k is `max(1, ceil(len(scored) × top_frac))`,
and on the **`signed-Delta` branch only** (`nominate.py:117-129`) `scored` has
already dropped regions with a missing `direction`; the model branch (`:131`)
keeps all rows. This bundle is signed-Δ primary, so the drop applies — but
nothing is missing (0 of 86,956), so `len(scored)` equals the universe and the
arithmetic coincides. k is properly the count of *measured* regions on that
branch, not a general property of `nominate.py`.

**It is an analysis artifact, not a nomination.** The signed-Δ bundle remains the
production nomination for this transition and the run-bundle contract is
untouched.

---

## 2. Region sets (fixed in advance)

The construct can only act where its DBD binds, so the test is a three-way
partition, not nominated-vs-rest:

| set | prediction | what it tests |
|---|---|---|
| **A** nominated AND Nanog-bound | opens | the joint claim |
| **B** Nanog-bound, NOT nominated | opens less than A | nomination adds information beyond binding |
| **C** nominated, NOT Nanog-bound | does **not** open | opening requires the DBD, not the ED alone |
| **D** control regions | no change | specificity |

**Set C is load-bearing.** Without it a positive result cannot be distinguished
from "the navigator rediscovered Nanog binding sites."

Each of A, B, C is built twice — once per registered ranking.

### "Nanog-bound": primary = motif, secondary = CUT&RUN (decided 2026-11-17)

Primary definition is **Nanog motif occurrence** within the region; a construct
CUT&RUN in MEF is a pre-registered secondary arm, and if generated it must be
produced **blind to the nomination lists**.

> **REVISION NEEDED before freezing.** The motif route is weaker here than it
> looked when this was decided, and a stronger already-existing alternative was
> found on the server afterwards:
>
> - The motif matrix
>   (`refs/vierstra_pkg_mm10_hg38/MEF_mES.motif.mm10.npz`, 206,313 × 282)
>   **contains no motif named Nanog** — a case-insensitive search over all 282
>   archetype names returns nothing. **Nanog's PWM is clustered into the
>   composite `OCT4+SOX2` archetype**, which is verifiable in one call:
>   `tabix motif/mm10.archetype_motifs.v1.0.bed.gz chr1:3000000-60000000`, whose
>   column 4 is the archetype and column 7 the constituent source PWM, returns
>   `OCT4+SOX2 ← NANOG_MOUSE.H11MO.0.A` (and `NANOG_HUMAN`, `PO5F1_*`,
>   `SOX2_MOUSE`).
>
>   **This is the substantive objection, and it is stronger than a missing
>   motif.** A motif-positive region under `OCT4+SOX2` cannot distinguish Nanog
>   occupancy from Oct4/Sox2 occupancy — the archetype deliberately merges them.
>   For set C ("nominated, NOT Nanog-bound"), which is the load-bearing set, that
>   is disqualifying: the partition would separate nominations by OSN co-binding,
>   not by Nanog.
>
>   *Two intermediate drafts of this box were wrong and are retracted.* One
>   claimed `OCT4+SOX2` was merely "the only relevant entry"; the other claimed
>   that was an error and that Nanog "necessarily falls inside one of the 25
>   `HD/*` clusters", with the mapping "not determinable from the files on the
>   server". The mapping is in column 7 of the shipped BED, and it places Nanog in
>   `OCT4+SOX2`, not an `HD/*`.
>
> - **A de-novo Nanog scan is possible, but it would be a new measurement.**
>   `get_regionmotif_matrix.py` is not a PWM scanner — it `tabix`-queries the
>   pre-scanned Vierstra archetype BED and sums archetype scores into the 282
>   canonical columns; its own docstring says motifs are "pre-scanned genome-wide
>   and shipped as a tabix-indexed BED per assembly". And `refs/genome/` holds
>   only `hg38.fa` (+ `.fai`), no mm10 FASTA.
>
>   *Correction to an earlier draft, which concluded from those two facts that a
>   de-novo scan was "not available on this server" and would need "new inputs
>   and new code". That was wrong — it generalised from one directory and one
>   script without looking anywhere else.* All three ingredients are on PeiLab2:
>   **`fimo` (MEME Suite 5.4.1)** at `/usr/local/bin/fimo`; an **mm10 FASTA**
>   already symlinked into this project's own tree at
>   `eCR/eCR_predictor/mef_mes_testrun/refs/mm10.fa`; and a **Nanog PWM** at
>   `eCR/eCR_predictor/jaspar_cache/UN0805.1.jaspar` (header `>UN0805.1 NANOG`),
>   alongside `MA2339.1.jaspar`, the mouse Nanog matrix `module_library.tsv`
>   itself cites for `DBD_MU_NANOG`. The only missing step is a `jaspar2meme`
>   conversion of a ~0.5 KB file.
>
>   So feasibility is not the objection. The objection is that a fresh scan is a
>   **new measurement with free parameters** (p-value threshold, background
>   model, how hits aggregate to a region), chosen by us after the nominations
>   were frozen. That is exactly what a pre-registration exists to avoid, and it
>   is a weaker reason to decline the motif route than the one above.
> - Nanog's PWM is degenerate and Nanog co-binds with Oct4/Sox2, so a
>   motif-defined "Nanog-bound" set is intrinsically noisy for this TF
>   specifically — more so than for ASCL1 or C/EBPα.
> - **`ground_truth/chipatlas_OSKM_mm10/Nanog.consensus.bed` already exists** —
>   ChIP-Atlas Nanog consensus peaks (100,019 regions, 3-column BED, PSC class,
>   Q<1E-20, peaks reproducible across ≥3 experiments). It is **self-evidenced
>   mm10**, not merely mm10-by-directory-name: it carries GRCm38-only contigs
>   (`chr1_GL456210_random`, `chrUn_JH584304`, and ~30 other `*_random` /
>   `chrUn_*` scaffolds), no chr20–22, and chr1 max end 195,372,131 against
>   mm10's 195,471,971 (hg38 chr1 is
>   248,956,422). It needs no new data, is equally immune to being influenced by
>   the result, and — unlike the composite archetype — is Nanog-specific.
>
>   *Correction to an earlier draft:* this file was described as "already used
>   for Claim 1". It was an **intermediate**, not a label. Claim 1's positive set
>   was `master_tf.consensus.bed` (132,766 regions, the union of six TFs), and a
>   repo- and server-wide grep finds no script or log that reads
>   `Nanog.consensus.bed`. Using it here is a new use, and the prereg should say
>   so rather than borrow unearned provenance.
>
> Its one weakness: it is binding in the **destination** state (PSC), not in MEF,
> so the construct's actual footprint in MEF may differ — which is precisely what
> the CUT&RUN secondary arm measures.
>
> **Recommendation: primary = `Nanog.consensus.bed`, with the scanned motif kept
> as a supporting definition for a sensitivity analysis, and CUT&RUN unchanged as
> the secondary arm.** One line to flip; needs confirmation.

### Set D — control regions

Set: the housekeeping-REGIONS set — always-open **atlas peaks** (501 bp
fixed-width; 99–100% of them annotate to cCREs, which is why earlier drafts
called them cCREs, but the sets are peaks and the distinction seeded a wrong
universe claim below). These are already open, so the prediction is *no change*,
not *no signal*.

**The depletion this set is valued for is axis-specific, and that matters here.**
`hk_regions/RESULTS.md` T7 reports, over mm10 bundles:

| bundle | AUC on `driver_score` | AUC on \|direction\| |
|---|---|---|
| mef_mes_gse201577 | 0.458 | 0.361 |
| cebpa_gse204742 | **0.632** | **0.711** |
| myod_gse186271 (Gate-1 refused) | 0.461 | 0.339 |
| etv2_gse168636 (Gate-1 refused) | 0.486 | 0.406 |

So housekeeping regions are clearly depleted on **`direction`** (0.34–0.41) but
sit **near chance on `driver_score`** (0.46–0.49), and on C/EBPα — the one other
fully-admitted mm10 bundle, and the only mm10 bundle where GET won Gate 2 — the
relationship **inverts** (hk *enriched*, 0.632 / 0.711).

Consequence for this design: set D is a strong control for the **signed-Δ** arm
and a **materially weaker one for the GET arm**, where no depletion has been
demonstrated and one admitted bundle runs the other way. Earlier drafts said
flatly "already shown depleted from navigator nominations", which is true only of
the signed-Δ ranking. Register set D for both arms, but state the asymmetry in
the paper rather than reporting a GET-arm null as if it were a clean pass.

**No blocker — an earlier draft of this section had the assembly backwards.**
`eCR/hk_regions/sets/` **is** the native mm10 set, not a human one:
`call_hk.py`'s own docstring declares its input peaks `mm10` and it filters
against `ref/mm10-blacklist.bed.gz`. Spot-checked against reference rather than
trusting the docstring: no chr20–22 anywhere in the six sets, and across all six
the maxima are chr1 195,314,584 (mm10 chr1 = 195,471,971; hg38 = 248,956,422),
chrX 169,903,304 (mm10 chrX = 171,031,299) and chr19 61,309,625 — the last
exceeding hg38's entire chr19 (58.6 Mb). Unambiguously mm10.

*Correction to an earlier draft:* it gave "chr1–19 + chrX only" with maxima
195,131,870 / 169,320,712 / 61,141,097. Those figures are right for
`hk_regions.bed` and `hk_regions.full.bed` (20 chromosomes), but **`mid`,
`random` and `specific` each carry one chrY region** (21 chromosomes), so the
"only" was false for half the sets it was asserted over. The mm10 conclusion is
unaffected.

`RESULTS.md` / `PLAN.md` give the defining data as a **mouse** aging
single-cell ATAC atlas (EasySci-ATAC-seq, GSE288730), 63 CPM matrices over
**one shared 1,341,077-peak universe** (1,334,657 after the mm10 blacklist
removes 6,420), across 1,596 (tissue, age, cell-type) contexts.

*Correction to an earlier draft:* it gave the defining universe as "the mm10
SCREEN cCRE universe (313,838 regions)". That is a **different** universe, used
for the proxy-transition and driver-score tests (T6/T7/T10), not for calling the
sets. The sets come from the 1.33M atlas peak universe.

Available today:

| file | regions |
|---|---|
| `sets/hk_regions.bed` | 20,000 |
| `sets/hk_regions.full.bed` | 20,352 |
| `sets/hk_strict_regions.bed` | 6,297 |
| `sets/mid_regions.bed` | 20,000 |
| `sets/specific_regions.bed` | 20,000 |
| `sets/random_regions.bed` | 20,000 |

`sets_hg38/` is the **derived human** side: the `*_human_regions.bed` files are
natively called from the 16-organ ENCODE panel (column 4 carries hg38 SCREEN
cCRE accessions), and the plain-named files there are the mouse sets lifted
**mm10 → hg38** (chain `ref/mm10ToHg38.over.chain.gz`, source mm10 coordinates
preserved in column 4). `t8_work/partB_lifted_mouse.*` is likewise mouse regions
*in hg38 coordinates* — byte-identical to the corresponding `sets_hg38/` plain
files — not an mm10 file.

Worth knowing if that lift is ever reused in the other direction: it is lossy.
20,352 → 7,561 (hk.full), 6,297 → 2,700 (strict), 20,000 → 4,815 / 4,183 / 4,043
(mid / random / specific). The matching `*.unmapped` files live in
`hk_regions/t8_work/`, not in `sets_hg38/`.

So set D needs only a sha256, not a liftOver or a rebuild. Both remedies the
earlier draft proposed were wrong: lifting human→mm10 would degrade an
already-native set, and "rebuild natively with `aggregate_breadth.sh` +
`call_hk.py`" describes the pipeline that already produced it.

Use `mid` / `specific` / `random` as additional graded controls. Validation runs
to **T16**, not T1–T9 as earlier stated (`t1_internal` … `t16_universe_swap`);
note that `t11_S4_confounded_DO_NOT_CITE.json` is marked unusable.

---

## 3. Region universe — fixed now

The construct ATAC will have its own peak set. Quantifying over an
experiment-derived universe reintroduces the dependence that corrupted Gate 1
the first time: PC1 and coherence rise as low-signal regions drop, and the v1
panel's universes spanned 63k–1.06M regions (17×), which made the six PC1 values
incomparable and rejected iN, the best transition.

**Decision: the analysis universe is the 86,956-region mm10 set of the existing
`mef_mes_gse201577` bundle, fixed now.** Construct ATAC is quantified over it.
Peaks outside it are reported separately as description, never folded into the
primary test.

---

## 4. Endpoints and pass criteria (stated before the data)

**Primary endpoint.** Difference in accessibility change (construct vs control
transduction) between set A and set B, over the frozen universe, for each
registered ranking.

- Effect measure: AUROC of accessibility gain, A vs B, bootstrap CI clustered by
  region-to-gene assignment (the clustering the Claim-2B harness already uses, so
  neighbouring regions are not counted as independent). Note the harness's own
  statistic is sign-accuracy — the gene-clustering and `assign_gene` are
  reusable, the AUROC-on-clusters estimator is not yet written.
- **Pass:** CI lower bound above 0.5.
- **Fail:** CI includes 0.5 — reported as a failure, not re-analysed into a
  subgroup.

**Secondary endpoints.**

1. GET-ranked vs signed-Δ-ranked, side by side; the pre-specified comparison is
   which gives the larger A-vs-B separation, with its CI.
2. Set C: fraction opening, expected at or below the universe base rate.
3. Set D: no significant mean accessibility change.
4. Top-k enrichment at k = 1%, 5%, 10% of each ranking — top-k is the metric that
   matches the navigator's actual job.

**Not endpoints.** Transcriptional output, conversion efficiency, marker
expression. Those evaluate the construct, not the navigator.

---

## 5. What would falsify the claim

- A does not open more than B → nominations add nothing beyond DBD binding. The
  paper must be willing to report this.
- C opens as much as A → the ED acts without DBD targeting; the construct is not
  site-directed.
- D changes → the construct is non-specific and no nomination claim is
  interpretable.

---

## 6. Existing assets reused (checked on PeiLab2, 2026-11-17)

| Need | Reused from | Status |
|---|---|---|
| GET per-region scores | `output/bundles/mef_mes_gse201577/weights.tsv` | on disk, frozen |
| signed-Δ nominations | same bundle, `nominations.tsv` | on disk, frozen |
| Nanog occupancy | `ground_truth/chipatlas_OSKM_mm10/Nanog.consensus.bed` | on disk, **100,019 regions**, 3-col BED (see §2) |
| mm10 motif matrix | `refs/vierstra_pkg_mm10_hg38/MEF_mES.motif.mm10.npz` | on disk, 206,313 × 282 (the MEF∪mES peak union, **not** the 86,956 analysis universe of §3) — **no motif named Nanog**; Nanog is clustered into `OCT4+SOX2` per col. 7 of the archetype BED |
| de-novo PWM scan | `fimo` 5.4.1 + `mef_mes_testrun/refs/mm10.fa` + `jaspar_cache/UN0805.1.jaspar` | **feasible** (needs a `jaspar2meme` step), but declined — it is a new measurement with free parameters chosen after the freeze (§2) |
| gene clustering for a bootstrap CI | `scripts/eval_direction_claim2b.py:88` | **clustering** reusable; its statistic is **sign-accuracy**, not AUROC, so §4's CI needs a new statistic written |
| control region sets | `eCR/hk_regions/sets/` (native mm10) + companions | on disk — needs only a sha256 |
| mm10 gene annotation | `refs/refGene.mm10.txt.gz` | on disk |

Nothing in this design needs a GPU, and no step re-derives data that already
exists.

---

## 7. Freeze (complete before the first construct ATAC library is made)

- [x] Universe: `mef_mes_gse201577`, 86,956 regions
- [x] signed-Δ ranking: `nominations.tsv`
      sha256 `c9221ce3c67040ab62f1f0140d668d3abcefd5879aa5e9589344d23c61004b24`
- [x] GET ranking: `artifacts/prereg_nanog_readout/get_ranking.mef_mes_gse201577.tsv`
      sha256 `e446547cb4ca84f64cf8e128837596ea338c09815fca68f5d000e302b8827c0a`
      (server copy: `ecr_nav/eval_work/prereg_nanog_readout/`)
- [x] Source weights: `weights.tsv`
      sha256 `c0f7fbca816e540ede08c0eaf9d320a867f19c58493ab90005dd34de007bd50f`
- [ ] "Nanog-bound" definition confirmed (§2 revision)
- [ ] Reconcile the mouse Gate-2 numbers quoted in §1. The bundle manifest says
      ΔAUROC +0.080, CI [−0.037,+0.198], p=0.107, over `n_anchors: 69`,
      `anchor_set: mtf.promoter.bed`. **The manifest's number matches none of the
      three published mouse rows.** [`claim2_results.md`](claim2_results.md)
      carries promoter/opening +0.071 [−0.045,+0.188] p=0.13 at **n=51**
      (also in [`validation_summary.md`](validation_summary.md)),
      neighborhood/opening +0.023, and promoter/**all** −0.051. The n=51 row is
      the 69 promoter anchors after the opening-only filter — so these differ by
      **positive set, not by seed or run**, which an earlier draft of this item
      asserted. And §1's +0.080 is effectively a **fourth, unpublished number**,
      not one of two. A freeze document cannot cite it without saying which
      filter produced it.
- [ ] Control set D sha256 recorded (the mm10 set itself is already on disk —
      `eCR/hk_regions/sets/hk_regions.bed`, 20,000 regions)
- [ ] Frozen universe BED exported, with sha256
- [ ] Git commit hash of eCR_navigator at freeze
- [ ] Date + explicit confirmation that no construct ATAC exists yet

Until every line is checked, this is a design document, not a pre-registration.

---

## 8. Note on system choice

`fib_beta` (BJ / IMR-90 → pancreatic β) is **GET-primary** on both arms, and its
IMR-90 arm carries **ΔAUROC +0.217 [+0.140, +0.296]** — the largest margin of any
transition run (BJ→β is +0.182 [+0.103, +0.261]). It already carries
extracted-DBD designs with AF3 structures: 4 SS18-N70 openers and 4 SALL1-N12
closers. It is the
better system for testing GET specifically — but in-house validation is blocked
(no mature GSIS protocol, per `endpoint_selection/ENDPOINT_SELECTION.md`). That
is why the wet-lab arm runs on MEF → mES, and the paper should say so rather than
leave the choice looking arbitrary.
