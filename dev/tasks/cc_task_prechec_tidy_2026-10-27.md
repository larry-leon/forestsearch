# CC task — pre-check tidy on `feature/docs-refs-maxminutes` (A: roxygen stamp, B: stale note, C: dead `fit_causal_forest_glm()`)

**Date:** 2026-10-27. **Repo:** forestsearch. **Branch:** `feature/docs-refs-maxminutes`, HEAD `b36da9dc` (verify with `git rev-parse --abbrev-ref HEAD` and `git rev-parse HEAD`; stop if either differs).
**Authorised by:** Larry, chat of 2026-10-27. Items are independent; Larry's disposition per item is recorded here: **A: yes  B: yes  C: yes** (fill in `yes`/`skip` before kickoff; CC stops if any is blank).
**Task doc path:** copy from `~/Downloads/` to `dev/tasks/cc_task_prechec_tidy_2026-10-27.md`, commit first and alone.

## 0. Rules

Same as the morning task: no `R CMD check`, no vignette build, no full suite; never touch `quarto/applications/**`; gates stop-on-failure; CC never pushes/merges; verify from the working tree, not from this doc's line hints.

## 1. Items

### A. Accept the roxygen2 stamp
Change `DESCRIPTION` `Config/roxygen2/version: 8.1.0` → `8.1.1`. From this commit on, `devtools::document()` no longer needs the restore step. Commit: `DESCRIPTION: accept roxygen2 8.1.1 stamp`.

### B. Stale note in `R/fs_family_report.R`
Read the comment at L≈354 (the one CC flagged in `REPORT_docs_refs_maxminutes_2026-10-27.md`). Rewrite it to state the current behaviour (the VI screen is deprecated and inert; the report row is `inert`). Comment text only; no code. Commit: `docs: fix stale note in fs_family_report.R`.

### C. Delete `fit_causal_forest_glm()` — gated
1. `grep -rn "fit_causal_forest_glm" R/ tests/ vignettes/ inst/ quarto/ NAMESPACE man/` — list every hit in the report.
2. Gate: proceed **only if** (i) `NAMESPACE` has no `export(fit_causal_forest_glm)`, and (ii) the only hits in `R/` are its own definition and roxygen, and (iii) there are no hits in `tests/`, `vignettes/`, `inst/`, `quarto/`. A hit in `quarto/` is a reason to **skip** (legacy documents are never edited), recorded, not fixed.
3. If the gate holds: delete the function and its roxygen block from `R/grf_helpers.R`; if that leaves an `@importFrom` in that file that nothing else in the package uses, leave the `@importFrom` in place (NAMESPACE must stay unchanged in this task). Add a NEWS line under a "Removed"/"Internal" heading: internal helper `fit_causal_forest_glm()` removed (unused since the `vi.grf.min` deprecation). Commit: `cleanup: remove dead fit_causal_forest_glm()`.
4. If the gate fails: no edit; the report says why and C stays on the December queue.

## 2. Procedure and gates

1. Branch/HEAD check; task doc commit.
2. Apply A, B, C in that order, one commit each (skipped items: no commit).
3. `devtools::document()`. Gate: `git diff --stat -- NAMESPACE` empty; `git diff -- man/` shows only `man/fit_causal_forest_glm.Rd` deleted (if C ran) and nothing else; `DESCRIPTION` differs from `b36da9dc` only in the stamp line (if A ran). Commit `docs: regenerate man pages` if anything in `man/` changed.
4. `NOT_CRAN=true` re-run of `tests/testthat/test-max-minutes-deprecation.R` and `test-vi-deprecation.R` (`load_package = "source"`) and, if C ran, `devtools::install(upgrade = FALSE)` then `test-lean-fits.R` with `load_package = "installed"`. Gate: all green.
5. `git status --porcelain quarto/` empty.
6. Append a section `## Pre-check tidy (same day)` to `dev/reports/REPORT_docs_refs_maxminutes_2026-10-27.md` with: the disposition line, the SHA of each commit, the full grep listing from C.1, the gate result for C, test pass counts. Commit `dev: report — pre-check tidy`. Copy the report to `~/Downloads/` (it overwrites the morning copy; that is intended).

End of task.
