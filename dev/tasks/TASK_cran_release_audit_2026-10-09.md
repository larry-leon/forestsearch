# TASK — CRAN release audit, Stage 0 (forestsearch 0.1.0 → next CRAN release)

- Date: 2026-10-09
- Owner: Larry (decides) / Chat (spec) / CC (executes)
- Repo: `larry-leon/forestsearch`, branch `feature/glm-extension` (verify; do not switch branches)
- Type: **read-only audit**. No edits to `R/`, `DESCRIPTION`, `NAMESPACE`, `man/`, `tests/`, `vignettes/`, `README*`, `NEWS.md`, `.Rbuildignore`. No package installs, reinstalls or removals. No `git fetch` / `pull` / `push`.
- R/ callout: **this task touches nothing in R/** (no move, no behaviour change, no method change).
- Authorization: Larry's kickoff of this task is his explicit request for the CRAN checks listed in Step 6 (an exception to the standing "no CRAN checks unless asked" policy, for this task only). Nothing is uploaded anywhere (no win-builder, no R-hub, no `submit_cran()`).

## Why

CRAN carries forestsearch 0.1.0 (published 2026-03-23, survival-only Title/Description, imports weightedsurv 0.1.0 from CRAN). The development line has moved well past it (GLM outcomes, MR/IJ/field inference, OC wrapper, threshold directives A–C, calibrated declaration, `pconsistency.digits`, admission alignment). Before any remediation is specified, the release gap must be measured from source, not from memory. This stage produces one report that the remediation task (Stage 1) will be written against.

## First action

1. Copy this file to `dev/tasks/TASK_cran_release_audit_2026-10-09.md` and commit it alone (explicit path).

## Preconditions (stop on failure)

- P1. Working tree clean apart from pre-existing untracked files (record `git status --porcelain` verbatim in the report; never stage untracked files you did not create).
- P2. Record branch, `HEAD` SHA, and `DESCRIPTION` `Version:`.
- P3. `pandoc` and `quarto` resolvable on `PATH` from the shell that runs the check (the 2026-09-18 vignette failure was a pandoc PATH issue). Record `pandoc --version` / `quarto --version` first lines. If missing, record and skip Step 6 vignette build (`--no-build-vignettes`, `--ignore-vignettes`), do not attempt to fix.

## Steps (each with a hard timeout; on timeout record "TIMEOUT" and continue to the next step)

### Step 1 — Release delta since 0.1.0 (≤ 5 min)
- Locate the commit corresponding to CRAN 0.1.0 (`git tag`; else the commit where `DESCRIPTION` last read `Version: 0.1.0` before the bump). Record SHA and how it was identified.
- `git diff --stat <0.1.0>..HEAD -- R/ DESCRIPTION NAMESPACE vignettes/ tests/ inst/ data/` summary counts.
- Exports: list exported symbols at 0.1.0 vs HEAD (from `NAMESPACE`). Report **added**, **removed**, **renamed** exports. Removed or renamed exports are user-breaking and must be listed prominently.
- Changed defaults of exported functions since 0.1.0 (compare formals of exported functions between the two versions by parsing source, not by installing): list each argument whose default changed.

### Step 2 — DESCRIPTION audit (≤ 5 min)
Report verbatim and flag each:
- `Title` (title case, ≤ 65 chars, no "package"/"R"/trailing period). Note it currently says "with Survival Endpoints" — record whether it still describes the package.
- `Description` (no "This package…", software/package names in single quotes, method references as `<doi:...>`).
- `Authors@R`, `Maintainer` email, `License` + `LICENSE` file, `URL`, `BugReports`, `Depends` R version.
- `Imports` / `Suggests` / `LinkingTo` / `VignetteBuilder` / `SystemRequirements` (quarto vignettes require `SystemRequirements: Quarto command line tools`).
- **`Remotes:` field present?** (CRAN-incompatible; weightedsurv is on CRAN as 0.1.0.)
- weightedsurv: version requirement stated; list every `weightedsurv::` call and every `importFrom(weightedsurv, …)`; for each, state whether it exists in CRAN weightedsurv 0.1.0 (check against the installed weightedsurv's namespace AND its version; if the installed version is not 0.1.0, flag it — a forestsearch release that needs an unreleased weightedsurv requires weightedsurv to be submitted first).
- Imports declared but unused / used but undeclared (scan `R/` for `pkg::` and `@importFrom`).

### Step 3 — Build contents and size (≤ 15 min)
- `R CMD build --no-build-vignettes --no-manual .` into a temp dir (not the repo). Record tarball size.
- List tarball contents; flag anything that should be in `.Rbuildignore`: `dev/`, `quarto/`, `_payloads/`, `manuscript/`, `claude/`, `.claude/`, `CLAUDE.md`, `docs/`, `pkgdown/`, `*.rds` outside `data/`/`inst/extdata`, `*.html` outside built vignettes, `REPORT_*`, `TASK_*`, `HANDOFF_*`, scratch scripts.
- Largest 20 files in the tarball with sizes.
- Delete the temp tarball afterwards.

### Step 4 — Documentation and code-policy scan (≤ 15 min, static only)
For each finding give file:line.
- Exported functions missing `@return` / `\value` in `man/*.Rd`.
- Exported functions with no `@examples`; examples wrapped in `\dontrun{}` (list each; CRAN expects `\donttest{}` unless the code genuinely cannot run).
- Parallelism in examples/tests/vignettes: any `plan(multisession|multicore|cluster)`, `workers =`, `detectCores()`, `mc.cores` that could use > 2 cores under check.
- Side effects in package code: `setwd`, `options(` / `par(` without `on.exit` restore, writes outside `tempdir()` (`saveRDS`, `write*`, `ggsave`, `dir.create`) in non-example code paths, `set.seed` inside exported functions, `installed.packages()`, `:::` to other packages, `T`/`F` literals, `print`/`cat` used for messaging in non-print methods, `Sys.setenv` without restore.
- Non-ASCII characters in `R/` (CRAN check flags).
- `utils::globalVariables` coverage (list the 2026-09-18 NOTE variables if still present).
- Roxygen: confirm `Roxygen: list(markdown = TRUE)`; confirm `RoxygenNote` vs installed roxygen2 version; do **not** re-document (report only whether `man/` is stale: run `roxygen2::roxygenise()` into a temp copy of the repo, diff `man/` and `NAMESPACE` against the repo copy, report the diff summary, delete the temp copy).

### Step 5 — Tests and vignettes inventory (≤ 5 min, static)
- Number of test files / `test_that` blocks; which use `skip_on_cran()`; any test that launches parallel workers or long bootstraps without `skip_on_cran()`.
- Vignettes: file names, engine (`%\VignetteEngine{…}`), whether precomputed, presence of large embedded outputs; `inst/doc` size if present. CRAN 0.1.0 shows installed size 6.1 Mb with `doc` 4.6 Mb (a NOTE on r-oldrel/4.4) — record the current equivalent.
- `NEWS.md` present? last heading; `cran-comments.md` present?

### Step 6 — Baseline CRAN check (authorized by kickoff; hard timeout 120 min total)
- Build with vignettes into a temp dir: `R CMD build .` (timeout 45 min). If it fails, record the error tail (last 60 lines) and run the check on a `--no-build-vignettes` tarball with `--ignore-vignettes` instead, stating so.
- `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE R CMD check --as-cran <tarball>` (timeout 75 min). Run with `_R_CHECK_LIMIT_CORES_=TRUE` set.
- Record the final status line and **every** ERROR / WARNING / NOTE verbatim with its check section. Compare against the 2026-09-18 reference (`0 errors | 1 warning | 2 notes`) from `dev/reports/STATUS_thresholds_workstream_2026-09-18.md`: identify which are the same findings.
- Record examples timing (`<pkg>.Rcheck/<pkg>-Ex.timings`): all examples > 5 s elapsed.
- Record tests total elapsed and vignette build elapsed.
- If `urlchecker` and `spelling` are already installed, run `urlchecker::url_check()` and `spelling::spell_check_package()` and record results; if not installed, record "not installed" (do not install).
- Delete temp build/check directories after extracting what the report needs (copy `00check.log` into the report directory as `dev/reports/cran_audit_00check_2026-10-09.log`).

## Deliverable

`dev/reports/REPORT_cran_release_audit_2026-10-09.md`, bullet form, sections in this order:

1. Summary — one bullet per blocker (anything CRAN would reject: ERROR, WARNING, `Remotes`, unreleased weightedsurv dependency, > 2 cores, writes outside tempdir, missing `\value`), then one bullet per NOTE-level item.
2. Repo state (P1–P3).
3. Release delta (Step 1), with breaking changes first.
4. DESCRIPTION (Step 2).
5. Build contents and size (Step 3).
6. Documentation and code policy (Step 4).
7. Tests and vignettes (Step 5).
8. CRAN check (Step 6), verbatim findings + comparison to 2026-09-18.
9. Open questions for Larry — only those a remediation spec cannot resolve from source.

Commit the report and the check log together (explicit paths). Do not propose follow-up tasks inside the report beyond section 9.

## Post-conditions (machine-checkable; stop on failure)

- PC1. `git diff --name-only <task-doc-commit>^ HEAD` lists only: `dev/tasks/TASK_cran_release_audit_2026-10-09.md`, `dev/reports/REPORT_cran_release_audit_2026-10-09.md`, `dev/reports/cran_audit_00check_2026-10-09.log`.
- PC2. `git status --porcelain` shows no tracked-file modifications.
- PC3. No temp build/check directories or tarballs left inside the repo.
- PC4. The report contains all nine section headings.
- PC5. The report's repo-state section names the HEAD SHA at which the audit was run.

## Out of scope

- Any edit that fixes a finding (Stage 1, specified after Larry reads this report).
- Version bump, NEWS, cran-comments, Title/Description rewrite.
- Deferred API items (`sg_objective` rename, `max.minutes` deprecation).
- win-builder / R-hub / mac-builder uploads and CRAN submission (Larry, later).
