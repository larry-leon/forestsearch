# TASK — `gbsg_020` DGM: the trial-wide (ITT) hazard ratio per cell (read-only), narrowed (2026-09-28)

**Supersedes** `TASK_gbsg020_dgm_estimands_2026-09-24.md` (never run). Same purpose, narrowed: the
catalog `quarto/simulations/gbsg_020/current_status.md` §1 already records the prevalence pairing and
what the design holds fixed; only the trial-wide estimand is missing.

**Repo:** `forestsearch`, branch `feature/glm-extension`.
**Kind:** read-only. No simulation is re-run and no rate or estimate is computed. R may be used only to
read a committed object's structure and stored constants (`names()`, classes, the value of a stored truth
field) in a session that writes nothing.
**Writes exactly two files:** this task document into `dev/tasks/`, and one report into `dev/reports/`.

**First action:** copy this file into `dev/tasks/` and commit it.

## Gates (stop on failure)
- No tracked file is modified before any step; untracked files are listed and left alone.
- The repository is `forestsearch` on `feature/glm-extension`.
- The `gbsg_020` directory and at least one committed results bundle exist.

## Step 1 — Q1, the trial-wide estimand (the one that matters)
For each cell type — **harm** (planted region at HR 1.50 and 1.75), **attenuated-benefit** (planted region
at HR 1.00), and **uniform-benefit** (the `null` design at marginal HR 0.657 and 0.721) — and at each
prevalence where one applies, give the **trial-wide marginal Cox hazard ratio** the design targets or
records.
- Read the template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`, the DGM constructor it calls, and one
  committed bundle per cell type (its `truth` object and any stored design constant).
- For each value: the field or script line that carries it, and whether it is a **design constant** or a
  **per-replicate realized value** (if the latter, quote the record's own summary and say so).
- If no field or line carries a trial-wide value for a cell type, write the exact words
  **`not established from source`** for that cell type. Do not compute one.

## Step 2 — Q2 and Q3, confirmation only
The catalog's §1 states: complement marginal HR 0.657 at 12.4% prevalence and 0.721 at 31%; in the
planted-region design `k_inter` is calibrated to the target HR inside the region and `k_treat = 1`, so the
complement carries the base treatment effect unmodified. Confirm each against source with one file:line,
or record a disagreement. No more than that.

## Step 3 — Q4, only as needed
List the fields of the stored `truth` object, one line each with its scale (marginal Cox, patient-level,
other), **only** to the extent needed to state Q1's answers on one consistent scale.

## Step 4 — Report and bundle
- `dev/reports/REPORT_gbsg020_itt_estimand_2026-09-28.md` (new file): a table, rows = cell type × prevalence,
  columns = trial-wide marginal HR (or `not established from source`) / design constant or realized /
  source file:line; then the Step 2 confirmations; then Step 3's field list if written. No recommendations.
- Copy it to `~/Downloads/gbsg020_itt_estimand_2026-09-28/`.

## Post-conditions (stop on failure)
- `git status` shows exactly two added files (this task in `dev/tasks/`, the report in `dev/reports/`) and no
  modified tracked file.
- Nothing re-run, nothing computed, no payload written; any payload read keeps its checksum.
- Every row of the Q1 table carries a file:line or `not established from source`.
- Commit the two files. **Do not push.**
