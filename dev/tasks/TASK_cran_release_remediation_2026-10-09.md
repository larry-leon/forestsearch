# TASK — CRAN release remediation, Stage 1 (forestsearch → 0.4.0)

- Date: 2026-10-09
- Owner: Larry (decides) / Chat (spec) / CC (executes)
- Repo: `larry-leon/forestsearch`, branch `feature/glm-extension`. Work on this branch; do not create or switch branches.
- Written against: `dev/reports/REPORT_cran_release_audit_2026-10-09.md` (audit HEAD `541fef30`). Every file:line cited there must be re-read from current source before editing; if a cited site has moved or changed, use the current source and record the difference.
- Dispositions: decisions S1-A … S1-F below are taken at the defaults stated (Larry's kickoff of this file = approval of those defaults).
- Authorizations given by this kickoff, for this task only:
  - One full `R CMD build` (with vignettes) and one `R CMD check --as-cran` at the end (Step 9). The package's own test suite runs inside that check.
  - Installing CRAN `weightedsurv` 0.1.0 and `latentcor` into a **temporary library under the session scratchpad only** (never the user library), deleted at the end.
  - Downloading the CRAN 0.1.0 source tarball into the scratchpad (read-only use).
- Not authorized: `git fetch`/`pull`/`push`; installs into the user library; edits to `quarto/`, `dev/` (other than this task and its report), `vignettes/`, `data/`, test fixtures; any upload (win-builder, mac-builder, R-hub, CRAN).
- Gates are stop-on-failure. A STOP ends the task with the report written to the point reached and committed.

## R/ callout (each item classified)

| Item | Files | Class |
|---|---|---|
| 2. `@param` / `@return` additions | `R/mrct_inject_dgm.R`, `R/grf_helpers.R`, `R/plot_sg_results.R`, `R/compare_selection_rules.R`, `R/forestsearch_cross_validation.R`, `R/sim_aft_gbsg.R` | Documentation only — roxygen comment lines; no code line changes |
| 4. weightedsurv install messages | `R/plot_km_band_forestsearch.R`, `R/plot_sg_weighted_km.R` | Messaging only — two error strings + one roxygen line |
| 6. `save_leff_calibration()` default | `R/leff_calibration_io.R` | **Changes behaviour**: a call without `output_dir` now errors instead of writing to `./_data` |
| 7. grf thread cap under check | new `R/utils_check_threads.R` + R/ files with grf training calls | Changes behaviour **only when `_R_CHECK_LIMIT_CORES_` is set**; identical arguments to grf otherwise |
| — | — | Nothing changes the method |

## Decisions (defaults taken)

- **S1-A Version:** `0.4.0`. Branch handling: release from `feature/glm-extension` (no release branch); after CRAN acceptance Larry tags `v0.4.0`, and a later task bumps to `0.4.0.9000`.
- **S1-B Title** (≤ 65 chars): `Exploratory Subgroup Identification in Survival and GLM Outcomes` (64).
- **S1-C API breaks since 0.1.0** (4 de-exports, removed `showten_subgroups`, changed defaults): documented as breaking in NEWS; no shims, no re-exports.
- **S1-D `save_leff_calibration(output_dir = NULL)`:** error when `output_dir` is not supplied, with a message naming the old default (`"_data"`) so callers can pass it explicitly.
- **S1-E `run_null_calibration` CPU ratio:** diagnose, then cap grf threads at 2 only under `_R_CHECK_LIMIT_CORES_` (Step 7).
- **S1-F Deferred, recorded only (no edits):** `\dontrun` → `\donttest` conversion (145 Rd), missing `\examples` (39 exports + 22 methods), `detectCores()` defaults not exercised under check, `subgroup.search()` plan not restored, `cox_ahr_cde_analysis()` `output_dir = "plots"` (opt-in via `save_plots`), unconditional `set.seed` sites, spelling WORDLIST, the internal `plot_subgroup()` `par` reset, deferred API items (`sg_objective`, `max.minutes`).

## First action

1. Copy this file to `dev/tasks/TASK_cran_release_remediation_2026-10-09.md` and commit it alone (explicit path).

## Preconditions (stop on failure)

- P1. `git status --porcelain` shows no tracked-file modifications. Record it verbatim; never stage untracked files you did not create.
- P2. HEAD is `541fef30` or a descendant of it (`git merge-base --is-ancestor 541fef30 HEAD`). Record HEAD SHA.
- P3. `pandoc --version`, `quarto --version`, `pdflatex --version` all succeed in the shell CC runs in. Record first lines. **STOP if any fails** — Larry fixes PATH before this task runs; do not work around it.

## Steps (hard timeout per step; a timeout is a STOP)

### Step 1 — Confirm the 0.1.0 baseline from the CRAN tarball (≤ 10 min, read-only)
- Download `https://cran.r-project.org/src/contrib/forestsearch_0.1.0.tar.gz` (fall back to `.../src/contrib/Archive/forestsearch/forestsearch_0.1.0.tar.gz`) into the scratchpad and extract it there.
- Re-derive, against the **tarball** (not the `v0.1.0` tag): removed / added exports from `NAMESPACE`, removed arguments and changed defaults of exports present in both (parse formals; same method as the audit).
- Record in the report whether this matches the audit's tag-based §3 list. NEWS (Step 8) uses the tarball-derived list.

### Step 2 — Documentation fixes (≤ 30 min)
- Add `@param` for: `analysis_covariates` (`dgm`, `case`, `drop`), `mrct_truth_by_region` (`dgm`, `seed`), `plot_region_modifier` (`dgm`, `file`, `main`), `sim_region_metrics` (`sim`, `region_col`, `subgroup_index`, `pi_preserve`). Write each description from what the function body does with the argument (type, role, default meaning).
- Add `@return` to the 11 Rd topics in audit §1: the 4 above, `print_grf_details`, and the S3 methods `plot.fs_sg_plot`, `print.fs_sg_plot`, `print.forestsearch_comparison`, `print.fs_kfold`, `print.fs_tenfold`, `print.gbsg_dgm`. For print/plot methods state what is actually returned (e.g. "`x`, invisibly"), verified from the method body.
- Run `roxygen2::roxygenise()`.
- Gate 2: (a) `NAMESPACE` byte-identical to HEAD; (b) every changed `man/*.Rd` is one of the 11 topics (or the Rd file holding them); (c) in the six R/ files of this step, every added/removed diff line matches `^[+-]\s*#'` (no code line changed).

### Step 3 — DESCRIPTION (≤ 5 min)
- `Version: 0.4.0` (S1-A).
- `Title:` per S1-B.
- `Description:` delete the word "advanced"; remove the single quotes around `GRF` and `MRCT` (method names, not software). No other wording changes.
- Add `SystemRequirements: Quarto command line tools (https://github.com/quarto-dev/quarto-cli)`.
- `weightedsurv` → `weightedsurv (>= 0.1.0)` in `Imports`.
- Gate 3: `desc::desc_get_field` (or `read.dcf`) confirms each field; `nchar(Title) <= 65`; `tools::toTitleCase(Title) == Title`.

### Step 4 — weightedsurv messages (≤ 5 min)
- Replace the `devtools::install_github('larry-leon/weightedsurv')` advice at `R/plot_km_band_forestsearch.R:342`, `R/plot_sg_weighted_km.R:334` and the roxygen line at `R/plot_km_band_forestsearch.R:171` with `install.packages("weightedsurv")`. Re-roxygenise if an Rd changes.
- Gate 4: `grep -rn "install_github" R/ man/` returns nothing related to weightedsurv.

### Step 5 — Build exclusions (≤ 2 min)
- Add `^inst/patches$` to `.Rbuildignore`. Do not move or delete the file.

### Step 6 — `save_leff_calibration()` default (S1-D) (≤ 15 min)
- When `output_dir` is `NULL`, `stop()` with a message stating that an output directory must be supplied and that the former default was `"_data"`. No other change to the function; `latest_rds()` untouched.
- Grep the whole repo (`R/`, `tests/`, `vignettes/`, `quarto/`, `dev/`, `inst/`) for calls of `save_leff_calibration(` that omit `output_dir`; list every hit (file:line) in the report. Edit only hits under `tests/` (pass `output_dir = tempdir()` or a `withr::local_tempdir()`); do not edit any other caller.
- Add a test: calling without `output_dir` errors; with `output_dir = tempdir()` writes there.

### Step 7 — `run_null_calibration` CPU ratio (S1-E) (≤ 30 min)
- 7a Diagnose (no code change). Extract the `run_null_calibration` example from `man/run_null_calibration.Rd` into a scratch script. Run it twice with `Rscript` under `_R_CHECK_LIMIT_CORES_=TRUE OMP_THREAD_LIMIT=2`, timing with `system.time()`: as written, and with `use_grf = FALSE`. Record user, system, elapsed and `(user+system)/elapsed` for both.
- Gate 7a: if the ratio with `use_grf = FALSE` is ≤ 2.0 and with grf > 2.5, grf is the source → continue to 7b. Otherwise mark Step 7 **BLOCKED** in the report with the measurements, skip 7b–7c, and continue to Step 8 (not a STOP).
- 7b Fix. New file `R/utils_check_threads.R` with an internal helper `.fs_grf_num_threads(user = NULL)`: returns `user` if non-NULL; else `2L` when `_R_CHECK_LIMIT_CORES_` is set to anything other than `""`/`"false"`/`"FALSE"`; else `NULL`. At every grf training call in `R/` (`grf::causal_forest`, `grf::regression_forest`, `grf::probability_forest`, `grf::survival_forest`, `grf::causal_survival_forest`, and any other `grf::*forest` call found), pass `num.threads = .fs_grf_num_threads(<existing num.threads value if the call already forwards one, else NULL>)`. List each site (file:line) in the report. Where a call passes `NULL` to grf's `num.threads`, grf's own default applies — identical to the current behaviour.
- 7c Verify. Re-run the 7a script as written under the same environment: ratio ≤ 2.0. Add `tests/testthat/test-check-thread-cap.R`: helper returns `NULL` with the variable unset, `2L` when set to `"TRUE"`, and the user value when supplied (use `withr::local_envvar`).
- Gate 7: ratio ≤ 2.0; the new test passes.

### Step 8 — NEWS.md and cran-comments.md (≤ 30 min)
- NEWS: rename `# forestsearch (development version)` to `# forestsearch 0.4.0`. Directly under it insert `## Changes since the last CRAN release (0.1.0)` in bullets:
  - `### Breaking changes` — first, from Step 1's tarball-derived list: each removed export (now internal), `showten_subgroups` removed from `forestsearch()` and `subgroup.consistency()`, and each changed default with old → new value.
  - `### New features` — one bullet per area, each supported by an existing NEWS entry or source (cite nothing that is not there): GLM outcomes (binary OR, continuous MD, count IRR), DINA and GRF identifiers, MR/IJ/field post-selection inference, operating-characteristics functions (`fs_oc_*`), threshold rules (c2 = 0.8·c1 when only c1 is given; c2 > c1 an error, survival HR and binary OR only), calibrated declaration threshold, `pconsistency.digits`, plus the Step 6 and Step 7 changes.
  - Demote the existing intermediate headings (`# forestsearch 0.2.x` … `0.3.5`) to `## forestsearch <v> (development; not released on CRAN)`. Their bodies are unchanged.
- cran-comments.md: rewrite for this submission:
  - `## Submission` — update from 0.1.0 to 0.4.0; intermediate version numbers were development-only and never submitted; the API changes listed in NEWS.
  - `## Test environments` — the local check of Step 9 (platform, R version); placeholders `win-builder (devel, release): TBD`, `mac-builder: TBD`.
  - `## R CMD check results` — filled from Step 9 verbatim, each NOTE with a one-line explanation.
  - `## Reverse dependencies` — from `tools::package_dependencies("forestsearch", db = available.packages(repos = "https://cloud.r-project.org"), reverse = TRUE, which = "all")`; record the result.

### Step 9 — Release check (authorized; build ≤ 45 min, check ≤ 75 min)
- Create the temp library `<scratch>/rlib`; install from CRAN source into it only: `weightedsurv` (must report version 0.1.0 and `Repository: CRAN`) and `latentcor`. Never into the user library.
- Build: `R CMD build .` (with vignettes) in a scratch directory. If the vignette build fails: record the last 60 lines and **STOP**.
- Check, with environment `R_LIBS=<scratch>/rlib _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE _R_CHECK_LIMIT_CORES_=TRUE OMP_THREAD_LIMIT=2 _R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5 _R_CHECK_TEST_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5`: `R CMD check --as-cran <tarball>`. No `_R_CHECK_FORCE_SUGGESTS_` override.
- Record the status line and every ERROR / WARNING / NOTE verbatim; examples over 5 s (user, system, elapsed, ratio); tests elapsed; vignette build elapsed; installed size; tarball size.
- If the check fails in a way traceable to CRAN weightedsurv 0.1.0 (missing function or argument), record it and **STOP** — that means weightedsurv must be updated on CRAN first.
- If a test fails only under the Step 7 cap (e.g. a grf-based fixture), **STOP** and report; do not edit fixtures.
- Copy `00check.log` to `dev/reports/cran_release_00check_2026-10-09.log`; keep the built tarball out of the repo; delete the temp library and the scratch build/check directories.
- Gate 9: status line has **0 ERRORs and 0 WARNINGs**. NOTEs are allowed only if each is explained in cran-comments.md.

## Deliverable

`dev/reports/REPORT_cran_release_remediation_2026-10-09.md`, bullet form, sections:

1. Outcome — final check status line; one bullet per gate (PASS / BLOCKED / STOP).
2. Repo state (P1–P3) and HEAD at start and end.
3. Step 1 — tarball vs tag comparison.
4. Steps 2–6 — what changed, file:line, gate results; the list of `save_leff_calibration` callers relying on the old default.
5. Step 7 — measurements before/after, grf call sites changed.
6. Step 8 — NEWS and cran-comments summary.
7. Step 9 — verbatim check findings, timings, sizes.
8. Recorded, not actioned — the S1-F items (one bullet each, no tasks).

Commit in logical commits with explicit paths (suggested: docs; DESCRIPTION + messages + .Rbuildignore; save_leff_calibration; thread cap; NEWS + cran-comments; report + log).

## Post-conditions (machine-checkable; stop on failure)

- PC1. `git diff --name-only <task-doc-commit>^ HEAD` lists only: this task doc, the report, the check log, `DESCRIPTION`, `NEWS.md`, `cran-comments.md`, `.Rbuildignore`, the R/ files named in Steps 2, 4, 6, `R/utils_check_threads.R`, R/ files containing a `grf::` forest call, `man/*.Rd` regenerated by roxygen, and `tests/testthat/*` files added or edited in Steps 6–7.
- PC2. `NAMESPACE` byte-identical to `541fef30`.
- PC3. `git status --porcelain` shows no tracked-file modifications; no tarball, `.Rcheck` directory or temp library inside the repo.
- PC4. `Version: 0.4.0` in `DESCRIPTION`; `NEWS.md` first heading is `# forestsearch 0.4.0`.
- PC5. The committed check log's status line matches `Status: (OK|[0-9]+ NOTEs?)$` (no ERROR, no WARNING) — unless the task STOPPED, in which case the report's Outcome says so.
- PC6. User library unchanged: `packageDescription("weightedsurv")$RemoteSha` (user library) equals its value at task start (record both).

## Out of scope

- Any S1-F item.
- Uploads (win-builder, mac-builder, R-hub) and the CRAN submission — Stage 2, on Larry's go.
- Tagging, pushing, and the post-release version bump.
