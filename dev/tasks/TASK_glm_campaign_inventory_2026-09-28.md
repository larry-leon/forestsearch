# TASK — inventory the GLM simulation campaigns (ACTG 175 binary and continuous), read-only (2026-09-28)

**Repo:** `forestsearch`, branch `feature/glm-extension`.
**Kind:** read-only. No simulation is run and nothing is computed. R may be used only to read a committed or
untracked object's structure (`names()`, classes, `nrow()`, stored settings) in a session that writes nothing.
**Writes exactly two files:** this task document into `dev/tasks/`, and one report into `dev/reports/`. Nothing
else is added, modified, committed, moved or deleted — **untracked campaign files are read in place and left
exactly as they are.**

**Why:** the fs-glms-interpretable manuscript states that the operating characteristics of the generalized-linear
paths are not evaluated. Larry reports that binary simulations were run. The manuscript chat needs to know what
exists before anything is written into the paper.

**First action:** copy this file into `dev/tasks/` and commit it (that commit only).

## Gates (stop on failure)
- No tracked file modified before any step; untracked files are listed and left alone.
- The repository is `forestsearch` on `feature/glm-extension`.

## Step 1 — find the campaigns
Search the repo (tracked **and** untracked) for simulation campaigns on the ACTG 175 data or on any binary
(logistic / odds-ratio / risk-difference) or continuous (mean-difference) outcome path. Start from
`quarto/simulations/actg175/` (including `binary_020/` and `binary_020/mr_or_harm/`), any continuous or MD campaign
directory, `dev/tasks/` and `dev/reports/` for their task and report documents, and any `current_status.md` or
catalog those directories keep. List every campaign directory found, tracked or untracked.

## Step 2 — for each campaign, report
1. **Identity:** directory, campaign tag/stem, the task and report documents that launched and recorded it (paths).
2. **Status:** complete, partial or in progress, from its own record; replicates per cell run versus planned.
3. **Git status:** every results file tracked and committed, or untracked (list the untracked ones).
4. **Design:** outcome type and effect measure (OR, RD, MD); the data-generating design (planted region, its
   effect, prevalence, sample sizes, cells); the trial-wide effect if recorded; replicates per cell. File:line for
   each.
5. **Analysis settings:** identifiers run (FS, DINA, GRF), selection rule and ε, thresholds (c1, c2, p⋆), MR on/off,
   which bounds were computed (field, field-s, Bonferroni pair, IJ two-sided, unadjusted, oracle). File:line.
6. **Constructions current?** Whether the run used the current field / field-s constructions, per its own record
   (e.g. the pin or the template it names). Do not re-derive; quote the record.
7. **Results available:** which payloads or extracts hold per-cell summaries (coverage, bias, declaration rate,
   classification), with paths; whether a summary table or report of results already exists.
8. **Anything the record flags** as superseded, failed, or pending.

Write `not established from source` where the record doesn't settle a point. No recommendations.

## Step 3 — report
`dev/reports/REPORT_glm_campaign_inventory_2026-09-28.md` (new): one section per campaign with the eight points;
a short table at the top (campaign / outcome / status / committed? / cells × replicates / results path).
Copy it to `~/Downloads/glm_campaign_inventory_2026-09-28/`.

## Post-conditions (stop on failure)
- `git status` shows exactly two added files (this task in `dev/tasks/`, the report in `dev/reports/`) beyond the
  untracked files present at the start, which are unchanged; no modified tracked file.
- Nothing run, computed or written besides the two files.
- Commit the two files. **Do not push.**
