# PeiLab2 file-management rules

How eCR_navigator work is laid out and maintained on the shared PeiLab2 box
(`/mnt3/wuyuhao/`). These rules were adopted in the 2026-07-30 reorg; the live map of
what-is-where is the server's own `ECR_NAV_INDEX.md` (this doc is the *policy*, that is the
*inventory*). Access + per-model runtime notes: [`server_mirrors.md`](server_mirrors.md).

## 1. One tree, categorized

**Everything eCR_navigator lives under `/mnt3/wuyuhao/ecr_nav/`.** The top level of
`/mnt3/wuyuhao/` holds only `ecr_nav/` + `ECR_NAV_INDEX.md` for eCR; everything else there is
other people's projects and is **off-limits** (never move, rename, or delete
`AIEpiRegulation/`, `foldx/`, `htge/`, `eCR/` [repo checkouts], `netMHCpan*`, `t2t*`, … — if
in doubt, it is not ours).

Inside `ecr_nav/`, every item goes in exactly one category:

| dir | holds | example |
|---|---|---|
| `datasets/` | endpoint data + raw downloads, one dir per transition | `mef_mes_clean_GSE201577/`, `iN_raw_GSE299923/` |
| `ground_truth/` | anchors / ChIP / master-TF loci used to *score*, not the transition data | `neural_iN/`, `chipatlas_OSKM_mm10/` |
| `refs/` | genome-wide references shared across runs | `hg38_cCRE_motif/`, `vierstra_pkg_mm10_hg38/` |
| `eval_work/` | analysis code + its outputs, one dir per claim/study | `claim1_work/`, `claim2b/` |
| `benchmark/` | the frozen benchmark panel (bundles + build tooling) | `benchmark/cebpa_gse151748/` |
| `output/` | deliverables consumed downstream | `bundles/`, `fixtures/` |
| `code/` | shared, path-agnostic tooling not tied to one study | `get_regionmotif_matrix.py`, `nav_code/` |
| `archive/` | provenance stubs for **deleted** data (see §4) | `MEF_mESC_SRP524550/` |

If a new thing does not obviously fit a category, it usually belongs in `eval_work/<study>/` or
`code/` — do not invent a new top-level `ecr_nav/` category without updating this table.

## 2. Naming

- Dataset dirs carry the **public accession**: `<system>_<role>_<ACCESSION>` — e.g.
  `mef_mes_clean_GSE201577`, `iN_raw_GSE299923`, `iCM_clean_DROPPED`. Role ∈ {`raw`, `clean`,
  `qc`}. A dropped/superseded set keeps a status suffix so nobody reuses it by accident.
- Archive stubs are `<origname>_<ACCESSION>` so the deletion is traceable to a study.
- Descriptive over cryptic: `vierstra_pkg_mm10_hg38`, not `MEF_mESC_3`.

## 3. Real dirs + path discipline (no symlink view)

- Items are **real directories**, physically placed. The old "symlink view over dirs that live
  elsewhere" is retired — it drifts out of date and hides the true location. The only symlinks
  kept are **compat pointers for genuinely shared refs** (e.g.
  `refs/vierstra_motifs_mm10_hg38 -> vierstra_pkg_mm10_hg38/motif`), and they must be
  **relative** so the tree stays self-contained and movable.
- **New scripts take paths as arguments, or reference a single `ROOT=/mnt3/wuyuhao/ecr_nav`
  constant** — do not scatter absolute paths through the body. `navigate.py` (args-only) is the
  model to follow; a script that hardcodes twenty absolute paths is the thing this rule exists
  to prevent.
- `mv` within `/mnt3` is a metadata rename — **instant and reversible**, even for 100 GB dirs.
  Moving is cheap; deleting is the only irreversible act. Use `mv`, never `cp`+`rm`, within the
  filesystem.

### Moving a dir (or reorganizing)

1. **Snapshot** the current layout + symlink targets first (`find -printf '%y %p -> %l\n'`).
2. Build a **deterministic old→new path map** and rewrite references from it, **most-specific
   name first** (`MEF_mESC_3` before `MEF_mESC_2` before `MEF_mESC`) so prefixes don't collide.
   Apply with one `sed -f` pass over `*.sh *.py *.md *.json` only — never over data or `.log`.
   Inserting `/ecr_nav/<category>/` into every path makes the rewrite idempotent (a second pass
   is a no-op).
3. `sed` does **not** fix symlink *targets* — after moving, repoint any absolute symlink whose
   target moved (`find -xtype l` finds the ones you broke).
4. **Verify before trusting** (§5).

## 4. Deletion policy (public-or-keep, always leave a record)

Delete large data **only** when all of these hold:

1. **It is public and re-fetchable** — a GEO/SRA/ENA accession, or a scripted download. Verify
   the accession actually exists; "probably public" is not a yes.
2. **A record is written to `ecr_nav/archive/<name>_<ACCESSION>/`** containing: the run/sample
   manifest (`SraRunTable.csv` or GSM list), the processing script(s) that built our artifacts
   from it, and a `REDOWNLOAD.md` (accession + why-deleted + exact re-fetch commands).
3. **Anything derived from it that we still need is retained** elsewhere (e.g. the clean
   endpoints, or a derived confound track) — deletion must not orphan a live artifact.
4. **Deletion runs last**, after the moves are verified and a functional smoke test passes, and
   is **guarded**: the delete step aborts unless the provenance files from (2) already exist.

Never delete another project's data, and never delete anything for which you cannot write a
truthful `REDOWNLOAD.md`. If it is not public and not reproducible, it stays.

## 5. Verify after any move or delete

Run all four; each must be clean:

- **No broken symlinks:** `find /mnt3/wuyuhao/ecr_nav -xtype l` → empty.
- **No stale path refs:** `grep -rIl '/mnt3/wuyuhao/<oldname>' ecr_nav --include='*.sh' --include='*.py' --include='*.md' --include='*.json'` → empty.
- **Functional smoke test:** re-run one real command (a `navigate.py --contract` bundle build,
  or the Claim-2B harness) and confirm it **reproduces committed numbers**, not merely "runs".
- **Update the inventory:** edit `ECR_NAV_INDEX.md` on the server, and any repo doc that cites a
  server path — live docs get the new path, historical/handoff docs get a one-line redirect
  banner rather than a full rewrite.

## 6. Working on the box

- Access: `ssh -i ~/.ssh/ecr_navigator -p 2020 wuyuhao@172.16.78.234` (Python 3.8, numpy
  1.24.4). Filter the post-quantum SSH warning line out of captured output.
- `du -sh /mnt3/wuyuhao/*` over the **whole** tree is slow (other projects are large) — size
  only the dirs you care about.
- Keep reorg scripts (like `ecr_nav/code/reorg_ecr_nav.sh`) checked in on the server as
  provenance of how the layout got the way it is.
