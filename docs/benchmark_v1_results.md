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

| # | Transition | Species | Gate-1 | GET Claim-1 AUROC (promoter) | GET vs signed-Δ (Claim-2A) | Gate-2 PRIMARY | ChromBERT |
|---|---|---|---|---|---|---|---|
| 1 | **fib→iN** (Ascl1) | human | ✅ admit 0.89 | **0.668** [0.606,0.727] | ΔAUROC **+0.170**, LR p=0.001 | **GET** ✅ | null |
| 2 | **preB→macrophage** (C/EBPα) | mouse | ✅ admit 0.83 | **0.640** [0.594,0.690] | ΔAUROC **+0.106** [+0.039,+0.177], LR p=0.001 | **GET** ✅ | null (0.47/0.52) |
| 3 | MEF→mES | mouse | ✅ admit 0.98 | 0.582 [0.501,0.660] | ΔAUROC +0.07, CI incl 0 | signed-Δ | (null, Claim 1) |
| 4 | MEF→iMPC (MyoD) | mouse | ❌ reject 0.80 | **0.288** [0.204,0.376] (anti) | ΔAUROC −0.24 (driver worse) | signed-Δ | null (0.44/0.51) |
| 5 | fib→iCM | human | ❌ reject 0.68 | 0.465 (null) | — | signed-Δ | (null, Claim 1) |
| 6 | fib→endothelial (ETV2) | mouse | ❌ reject 0.53 | 0.605 [0.461,0.745] (n=19, underpowered) | ΔAUROC +0.11, CI incl 0, LR p=0.084 | signed-Δ | null (0.42/0.48) |

## What the panel shows

1. **The Claim-1/2 result generalizes — GET wins on two independent strong transitions.**
   iN (neuron, human) *and* C/EBPα (macrophage, mouse) both: Gate-1 admit, GET AUROC ~0.64–0.67
   with CI excluding 0.5, top-5% fold 3.5–3.7×, **and GET beats signed-Δ** (ΔAUROC CI excludes 0,
   incremental LR p=0.001). Different lineages, different labs, different species — the signal is
   not an iN artifact.
2. **Gate-1 separates the regimes.** Every Gate-1-**admit** transition has an informative or
   directional GET signal; every Gate-1-**reject** transition (MyoD 0.80, iCM 0.68, ETV2 0.53)
   has GET failing or underpowered, PRIMARY=signed-Δ. The admissibility gate earns its place.
3. **Weak-transition failure modes differ, all correctly caught by Gate-2:**
   - **MyoD** — GET *anti*-informative (0.29): the incomplete iMPC conversion leaves master-TF
     loci among the *un*changed regions. Preflight not fooled by the significant-but-negative LR.
   - **ETV2** — GET *weakly positive but underpowered* (0.60, n=19) and doesn't beat signed-Δ.
     Cause: the endo_r1 replicate is low-quality (7,441 peaks vs ~44k), so the endothelial state
     is incoherent (PC1 0.53). A hint of signal (endothelial TFs do open), not enough to clear.
   - **iCM** — flat null (0.49).
4. **ChromBERT is null on every mouse bundle** (MyoD 0.44–0.51, C/EBPα 0.47–0.52), confirming the
   Claim-1 mouse null — a clean native-mm10 negative control. **GET is the informative model.**
5. **MEF→mES is the "clean but directional" case** — Gate-1 admit, GET informative (0.58), but
   signed-Δ already captures it (Claim-2A CI includes 0) → PRIMARY=signed-Δ. Strong ≠ model-wins;
   the model must beat the baseline, which Gate-2 measures.

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
- Thresholds (PC1≥0.80 etc.) now calibrated on **6** transitions, up from 2 — the admit/reject
  split held on all six.
- Optional follow-ups: redo C/EBPα on dataset peaks (used cCRE fallback); source a cleaner
  endothelial dataset for a 3rd strong; add the 3 hg38 models when a human bundle is added.
