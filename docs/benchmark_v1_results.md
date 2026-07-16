# Benchmark v1 — results (GET + native-mm10 ChromBERT)

First run of the frozen transition panel ([`benchmark_spec.md`](benchmark_spec.md)). Scores each
transition by the same scorecard: Gate-1 admissibility, Claim-1 matched AUROC (informative vs a
|Δaccessibility|-matched background), Claim-2A vs a signed-Δ baseline, and the Gate-2 preflight
verdict (PRIMARY = which score to nominate from). Machine-readable: `benchmark/scorecard.tsv`
(on PeiLab2). Models this round: **GET** (all transitions) + **ChromBERT** (native mm10, on the
new mouse bundles); the 3 hg38-only models are omitted on mouse (liftOver-sparse/null — see
`validation_summary.md`); iN/MEF→mES/iCM already carry all five from the Claim 1/2 work.

## Headline

**GET's `driver_score` recovers master-TF loci and beats the signed-Δ baseline on the strong,
clean transitions — now shown on TWO independent lineages (neuron + macrophage), not just iN —
and fails on every weak/Gate-1-reject transition. Gate-1 cleanly separates the two regimes.**

## Scorecard

Gate-1 columns were **recomputed 2026-07-16** on the fixed universe (see the correction
below); the values in the original v1 run are superseded.

| # | Transition | Species | Gate-1 (PC1 @50k) | GET Claim-1 AUROC (promoter) | GET vs signed-Δ (Claim-2A) | Gate-2 PRIMARY | ChromBERT |
|---|---|---|---|---|---|---|---|
| 1 | **fib→iN** (Ascl1) | human | ✅ admit **0.919** | **0.668** [0.606,0.727] | ΔAUROC **+0.170**, LR p=0.001 | **GET** ✅ | null |
| 2 | **preB→macrophage** (C/EBPα) | mouse | ✅ admit **0.963** | **0.640** [0.594,0.690] | ΔAUROC **+0.106** [+0.039,+0.177], LR p=0.001 | **GET** ✅ | null (0.47/0.52) |
| 3 | MEF→mES | mouse | ✅ admit **0.933** | 0.582 [0.501,0.660] | ΔAUROC +0.07, CI incl 0 | signed-Δ | (null, Claim 1) |
| 4 | MEF→iMPC (MyoD) | mouse | ❌ reject **0.785** | **0.288** [0.204,0.376] (anti) | ΔAUROC −0.24 (driver worse) | signed-Δ | null (0.44/0.51) |
| 5 | fib→iCM | human | ❌ reject **0.707** (coherence **−0.059**) | 0.465 (null) | — | signed-Δ | (null, Claim 1) |
| 6 | fib→endothelial (ETV2) | mouse | ❌ reject **0.561** | 0.605 [0.461,0.745] (n=19, underpowered) | ΔAUROC +0.11, CI incl 0, LR p=0.084 | signed-Δ | null (0.42/0.48) |

## Correction — Gate-1 recomputed on a fixed universe (2026-07-16)

Building the run bundles (`run_bundle_contract.md`) exposed a flaw in the original Gate-1
column. **PC1 and coherence are not universe-invariant** — both rise as low-signal regions
are dropped, since those add noise dimensions and no transition signal. The v1 panel's
universes span **63,562–1,061,709 regions (17×)**, so the six PC1 values were never
comparable to one another, and the threshold was partly measuring universe size.

Two concrete errors this hid:

- **iN's Gate 1 had never actually been run.** The recorded 0.89 was a QC estimate on a
  different region set (`claim2_results.md` said so: *"per-rep matrix not on disk"*). Built
  from the 6 endpoint samples and run for real, iN scores PC1 **0.792 over its full 1.06M
  cCRE universe — a REJECT.** That is the transition GET demonstrably wins on, so the gate
  was refusing our best result.
- **iN's replicate counts were wrong** in the first fixture (2/2). The real endpoints are
  **3/3** (`claim1_human_progress.md`: start = D0 ×3 WT/AAA/QA, end = D7 WT ×3).

**Fix:** Gate 1 now computes on a **fixed universe of the 50,000 most accessible regions**,
ranked on depth-normalized signal (so a deeply-sequenced sample can't decide what counts as
accessible). N=50k is forced by the panel — ETV2's 63,562 is the smallest universe.

**Result: the split now holds for the right reason.** All six match their expected verdict,
and `MIN_PC1 = 0.80` lands inside a real **0.785 → 0.919 gap** instead of among the values.
iCM additionally fails coherence (−0.059), which is honest: its "replicates" are the 5F/4F
cocktail arms, not true replicates.

**The earlier claim that "Gate-1 separated admit/reject cleanly on all six, so the thresholds
are calibrated on n=6" did not hold as stated** — that separation was partly universe
artifact. It holds now, on the fixed universe.

## What the panel shows

1. **The Claim-1/2 result generalizes — GET wins on two independent strong transitions.**
   iN (neuron, human) *and* C/EBPα (macrophage, mouse) both: Gate-1 admit, GET AUROC ~0.64–0.67
   with CI excluding 0.5, top-5% fold 3.5–3.7×, **and GET beats signed-Δ** (ΔAUROC CI excludes 0,
   incremental LR p=0.001). Different lineages, different labs, different species — the signal is
   not an iN artifact.
2. **Gate-1 separates the regimes — on the fixed universe.** Every Gate-1-**admit**
   transition (iN 0.919, C/EBPα 0.963, MEF→mES 0.933) has an informative or directional GET
   signal; every Gate-1-**reject** transition (MyoD 0.785, iCM 0.707, ETV2 0.561) has GET
   failing or underpowered, PRIMARY=signed-Δ. The gate earns its place *once its statistic is
   universe-invariant* — on the raw universes this split was partly an artifact (see the
   correction above). It remains a **coarse screen for severe failures, not the decision**:
   MEF→mES has the panel's cleanest endpoints and GET still loses there, which only Gate 2
   catches.
3. **Weak-transition failure modes differ, all correctly caught by Gate-2:**
   - **MyoD** — GET *anti*-informative (0.29): the incomplete iMPC conversion leaves master-TF
     loci among the *un*changed regions. Preflight not fooled by the significant-but-negative LR.
   - **ETV2** — GET *weakly positive but underpowered* (0.60, n=19) and doesn't beat signed-Δ.
     Cause: the endo_r1 replicate is low-quality (7,441 peaks vs ~44k), so the endothelial state
     is incoherent (PC1 0.561). A hint of signal (endothelial TFs do open), not enough to clear.
   - **iCM** — flat null (0.49).
4. **ChromBERT is null on every mouse bundle** (MyoD 0.44–0.51, C/EBPα 0.47–0.52), confirming the
   Claim-1 mouse null — a clean native-mm10 negative control. **GET is the informative model.**
5. **MEF→mES is the "clean but directional" case** — Gate-1 admit, GET informative (0.58), but
   signed-Δ already captures it (Claim-2A CI includes 0) → PRIMARY=signed-Δ. Strong ≠ model-wins;
   the model must beat the baseline, which Gate-2 measures.

## The panel as run bundles (2026-07-16)

The whole panel is now regenerated through the production entrypoint —
`navigate.py --contract --bundle` — as v1 run bundles ([`run_bundle_contract.md`](run_bundle_contract.md)),
on PeiLab2 at **`/mnt3/wuyuhao/bundles/`**. This is the end-to-end test of the nomination
policy: **no per-transition special-casing**, one command each, the verdict falls out.

| Bundle | regions | nominations | PRIMARY | Gate 1 (PC1 @50k) |
|---|---|---|---|---|
| `in_gse299923` | 329,983 | **3,300** | **GET** | admit 0.919 |
| `cebpa_gse151748` | 313,838 | **3,139** | **GET** | admit 0.963 |
| `mef_mes_gse201577` | 86,956 | **870** | **signed-Δ** | admit 0.933 |
| `myod_gse186271` | 232,788 | **0** | *(refused)* | REJECT 0.785 |
| `icm_gse179011` | 233,342 | **0** | *(refused)* | REJECT 0.707 (coherence −0.059) |
| `etv2_gse168636` | 63,562 | **0** | *(refused)* | REJECT 0.561 |

**2 nominate from GET, 1 from signed-Δ, 3 refuse — exactly the predicted split.** All six
pass the structural checks: three parts present, `nominations.tsv` a strict subset of
`weights.tsv`, ranks dense `1..k`, `bundle_version` stamped.

Three things this exercised for the first time:

1. **The signed-Δ nomination path.** MEF→mES is the "clean but directional" case — Gate 1
   admits it with the panel's second-cleanest endpoints, and Gate 2 still measures
   ΔAUROC +0.080 (CI [−0.037,+0.198], LR p=0.107) → **PRIMARY = signed-Δ**. It nominates 870
   regions ranked by |ΔaTPM| (rank-1 |direction| 0.996), *not* by driver score. Strong ≠
   model-wins, and the bundle records which score it actually used.
2. **Refusal at scale.** The three rejects ship **0 nominations and 63k–233k weights** each —
   the contract's claim that a Gate-1 reject still yields usable off-target weights, now true
   of the real artifacts rather than only the prose.
3. **`--contract` mode.** The GPU mirrors are not persistent and their `.npz` artifacts are
   gone, but every archived contract survived, so bundles rebuild with **no GPU**. Three of
   the six predate the `direction` column and are 4-column; their direction was re-attached
   from the same measured aTPM table the original run used (`--direction-norm raw`). All six
   came out **100% measured — 0 unmeasured regions**.

The pipeline reproduces the committed `examples/run_bundle/in_gse299923/nominations.tsv`
**byte-identically**, so the shipped fixture and the production path agree.

## Data-engineering notes (for reproduction)

- **Region universe:** dataset-specific peaks per transition (macs from CDesk, or my own
  bowtie2+macs3). **C/EBPα used a cCRE-fallback** (bigWig quantification over 313,838 mm10 cCREs)
  because its GEO deposit is bigWig-only and CDesk's peak/bigWig stage failed (persistent pyBigWig
  error); the strong result stands but a dataset-peaks redo is optional.
- **ETV2 peaks** were called with my own bowtie2+macs3 (CDesk aligned only ~1/4 samples here).
- **CDesk gotchas:** `ports=2` for paired-end (1=single); never run two CDesk jobs in parallel;
  its peak/bigWig stage is unreliable in this env — prefer direct macs3 on the BAMs.
- Anchors: master-TF loci per destination cell type (`benchmark/*/anchors/`), mm10 refGene.
- GET motif matrices built on PeiLab2 with the **local** mm10 Vierstra file (the mirror's Altius
  route is throttled). Embeds on the GET (A800) + ChromBERT (A100) mirrors.

## Status / open items

- **Complete.** All 6 transitions scored with GET; ChromBERT run on all 3 new mouse bundles
  (all null: MyoD 0.44–0.51, C/EBPα 0.47–0.52, ETV2 0.42–0.48).
- v1 composition ended **2 strong (iN, C/EBPα) + 1 directional control (MEF→mES) + 3 weak/reject
  (MyoD, iCM, ETV2)**. ETV2 was *intended* strong but came in weak (data quality: endo_r1 low
  peak count); a cleaner endothelial dataset (or better replicates) would restore a 3rd strong.
- Thresholds (PC1≥0.80 etc.) calibrated on **6** transitions, up from 2 — the admit/reject
  split holds on all six **once Gate 1 uses the fixed 50k universe** (see the correction
  above; on the raw universes it did not, and iN rejected).
- Optional follow-ups: redo C/EBPα on dataset peaks (used cCRE fallback); source a cleaner
  endothelial dataset for a 3rd strong; add the 3 hg38 models when a human bundle is added.
