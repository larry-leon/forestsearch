# TASK: forestsearch — the check that closes the `mrct-analytic` series (r2)

Date: 2026-11-01 (r2 supersedes the unlaunched r1: §1 items 5–6 added)
Repository: **forestsearch**, branch `mrct-analytic`, from `450b5492` (the end of `TASK_fs-mrct-article-r2_2026-11-01`). CC does not push; Larry merges from GitHub Desktop after this task.
Authorisation: Larry asked for the full suite and the CRAN-style check for this series (2026-11-01). This is the "before a release" case CLAUDE.md names.

**`R/` category:** fixes only, called out in §1; no new behaviour.

## 0. First action

Copy this file into `dev/tasks/` and commit it on its own, on `mrct-analytic`.

## 1. Small fixes first (one commit)

1. `R/mrct_simulation.R` about line 1210: the non-ASCII em dash in a roxygen example comment (flagged in `NOTES_fs-mrct-region-sims-record_2026-11-01`). Replace with ASCII.
2. `R/forestsearch_main.R` about line 992: the roxygen for `mr_inference_args$return_reselection` says default `FALSE`, while `forestsearch()` passes `TRUE` (`NOTES_fs-exact-betaHhat-mr-record_2026-11-01_r2`, decision 4). Make the documentation state what the code does; do not change the code.
3. Any other non-ASCII character in `R/` that the check below reports.
5. `_pkgdown.yml` reference index: every export this series added must be listed (`ahr_population`, `ahr_limit`, `ahr_design`, `ahr_solve`, `scale_dgm_hazard`, `mrct_select_case`, `mrct_case_study`), under the existing MRCT section or a new one beside it; `pkgdown::check_pkgdown()` must pass.
6. `plot.mrct_case_study()`: the Kaplan–Meier x-axis label is "time"; give the method a `time_unit` argument (default `NULL` → "time") and have the article pass "months".
4. `NEWS.md`: the 0.4.1.9000 entry gets a one-line heading sentence naming the series (the MRCT consistency analytic layer, DGM arguments, the per-replicate record, the exact β(Ĥ) option, `mr_estimates_table()` without a template, the case-study helpers, the article).

## 2. The full suite

`devtools::test()` on the installed build from this tree. Record counts per file for the files this series added or touched, and the totals. Stop on any failure that traces to this series' files; a failure elsewhere that predates the series (check by running that file at `19b703b` in a scratch worktree) is recorded, not fixed.

## 3. The check

`rcmdcheck::rcmdcheck(args = "--as-cran")`, the certification surface of CLAUDE.md (it builds the PDF manual). Record the full output in `dev/check/check_mrct-analytic_<date>.txt` (create the directory if absent; `dev/` is build-ignored). Stop on any ERROR or WARNING. NOTEs are recorded with a one-line disposition each (pre-existing at `19b703b`, or introduced by this series and why acceptable; if a NOTE introduced by the series is fixable in a line, fix it and re-run).

Points the series touches that the check will exercise: the five new Rd files of the analytic layer and the two helpers; the examples of `scale_dgm_hazard()`, `ahr_*()`, `mrct_select_case()`, `mrct_case_study()` (they must run within CRAN example time, or be wrapped in `\donttest{}` with a comment saying why); `mvtnorm` in Imports with the `importFrom`; the test fixtures' sizes; the `.Rbuildignore` coverage of `vignettes/articles`, `docs`, `dev` and `HANDOFF_*`.

## 4. The site article

`pkgdown::build_article("articles/mrct-consistency")` completes on this tree (the site itself is rebuilt by Larry's usual route after the merge). Record its time.

## 5. Notes file

`~/Downloads/NOTES_fs-mrct-check_2026-11-01_r2.md`: status, commits, the suite totals, the check's ERROR/WARNING/NOTE lines with dispositions, the article build time, and the branch's full commit list from `19b703b` for the merge.

## 6. Commits (on `mrct-analytic`)

Two: this task; the fixes of §1 with the check record. If §2 or §3 stops, the notes say where, and nothing after the stop is committed.
