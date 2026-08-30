# CC task — correct the roxygen for `use_grf` (role misattributed to variable importance) — v3

**Date:** 2026-08-31 · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`
**v3 change:** G1 updated after the push-size-fix history rewrite — the expected-tip note now references `13d5987e` (the size-fix report commit; `4ca38eec` was replayed as `0e9d48fa` and no longer exists), and the untracked-by-design `comparison_continuous.rds` is added to the expected untracked paths. Nothing else differs from v2.
**Type:** documentation-only edit **touching `R/`** — roxygen comments and regenerated `man/` plus a NEWS.md bullet. **Zero behavior change: every added or removed line under `R/` must begin with `#'`.** No renders, nothing under `quarto/`.
**Report:** `dev/tasks/REPORT_roxygen_use_grf_2026-08-31.md`
**Compute:** negligible (`devtools::document()` only).

## Decision being implemented (Larry, 2026-08-31)

The roxygen misstates the roles of `use_grf` and `vi.grf.min`. Ground truth from source: `use_grf` gates **GRF candidate-cut generation** (Section 3A, `if (use_grf && ...)`, ~L2491 — `grf.subg.harm.survival()`/`.glm()` cuts entering the family via `get_FSdata()`); `vi.grf.min` gates **variable-importance ordering/screening of the cut columns** (Section 5, `if (!is.null(vi.grf.min))`, ~L2791), independent of `use_grf`. The `@param vi.grf.min` and `@param max_n_confounders` blocks are **already correct — do not touch them**. The fix is the `@param use_grf` text plus two minor same-conflation roxygen touches.

## Gates — stop semantics

On any STOP: write the report with what was gathered, banner `STOPPED AT <gate>: <reason>` first, commit only the report, leave everything else uncommitted, end. Do not ask questions.

- **G1 — repo state (v3).** Branch `feature/glm-extension`. **No tracked modifications:** `git status --porcelain --untracked-files=no` must be empty; any tracked modification ⇒ STOP. Untracked entries are permitted and expected: `?? dev/tasks/_baseline_html_2026-08-31/`, `?? dev/tasks/_diffs_2026-08-31/` (render-task scratch), and `?? quarto/applications/actg175/_payloads_2026-08-31/analysis_actg175_continuous_compare_all/comparison_continuous.rds` (excluded from tracking by the push-size-fix policy — on disk by design). Record any *other* untracked paths in the report and proceed — they are irrelevant to a docs-only edit. Record `git log -1 --oneline` (expected at or after `13d5987e`, the push-size-fix report commit; whether the branch has been pushed does not matter to this task).
- **G2 — exact-text match.** Each `OLD` block below must occur **exactly once** in its file at HEAD. Line numbers are hints from the 29 Aug snapshot; the text is the contract. Any mismatch or multiple match ⇒ STOP, quoting the current text found.
- **G3 — roxygen-only diff.** After all edits, `git diff -- R/` must show added/removed lines beginning with `#'` only (ignoring hunk headers). Any other changed line ⇒ STOP and revert nothing (leave for inspection).
- **G4 — document() outcome.** `devtools::document()` succeeds; `man/forestsearch.Rd` contains "candidate-cut generation" and no longer contains "Use GRF for variable importance". Failure ⇒ STOP.

## Edits (transplant-form; apply with exact string replacement)

### E1 — `R/forestsearch_main.R`, `@param use_grf` (~L255)

OLD (three lines, verbatim):
```
#' @param use_grf Logical. Use GRF for variable importance. Default FALSE.
#'   The default changed from \code{TRUE} in an earlier release, for the same
#'   fixed-candidate-family reason given under \code{use_lasso}.
```

NEW:
```
#' @param use_grf Logical.  Generate additional candidate cuts from a GRF
#'   fit: tree-derived cutpoints are extracted via
#'   \code{grf.subg.harm.survival()} / \code{grf.subg.harm.glm()} (or taken
#'   from \code{grf_res} / \code{grf_cuts} when supplied) and enter the
#'   candidate family through \code{\link{get_FSdata}} alongside the
#'   quartile (and any DINA) cuts.  This flag governs candidate-cut
#'   generation only: GRF variable-importance ordering of the cut columns
#'   is a separate step controlled solely by \code{vi.grf.min} and runs,
#'   or is skipped, regardless of \code{use_grf}.  Default FALSE.
#'   The default changed from \code{TRUE} in an earlier release, for the same
#'   fixed-candidate-family reason given under \code{use_lasso}.
```

### E2 — `R/forestsearch_main.R`, `max.minutes` details (~L690–691)

OLD (two lines, verbatim; the second is the target, the first anchors uniqueness):
```
#'   behavior. Search scope is governed by \code{maxk}, candidate-factor
#'   screening (\code{use_grf}, \code{use_lasso}, \code{conf_force}),
```

NEW:
```
#'   behavior. Search scope is governed by \code{maxk}, the candidate-factor
#'   front ends (\code{use_grf}, \code{use_lasso}, \code{conf_force}),
```

### E3 — `R/interpret_search_config.R`, `@param use_grf` (~L31)

OLD (one line, verbatim):
```
#' @param use_grf Logical.
```

NEW:
```
#' @param use_grf Logical.  Whether GRF candidate-cut generation is on;
#'   see \code{\link{forestsearch}}.
```

### E4 — `NEWS.md` bullet

Read the top of `NEWS.md`. If a development/unreleased heading exists (e.g. `# forestsearch (development version)`), add under it; otherwise create that heading above the topmost version heading and add under it:

```
* Documentation: corrected the `use_grf` `@param`, which described the
  variable-importance role belonging to `vi.grf.min`. `use_grf` governs GRF
  candidate-cut generation only; variable-importance ordering/screening of
  the cut columns is controlled solely by `vi.grf.min` (`NULL` default since
  0.3.0). No behavior change.
```

Quote the heading context before/after in the report.

## Steps

1. **Environment + G1.** Record hostname, branch, `git log -1 --oneline`, installed `packageVersion("forestsearch")`, and the untracked-paths listing per G1.
2. **G2 checks**, then apply E1–E4. Show `git diff --stat` and the full `git diff -- R/ NEWS.md` in the report (it is short).
3. **G3** roxygen-only verification on the `R/` diff.
4. `Rscript -e 'devtools::document()'`. Apply **G4**. Expected regenerated files: `man/forestsearch.Rd` and the Rd for `interpret_search_config`. If `document()` regenerates additional `man/` files, that is pre-existing drift from an un-regenerated earlier state: record their names and include them (they reflect current roxygen truth), noting this in the report.
5. **Commit 1:** the edited `R/` files + `NEWS.md` + all regenerated `man/` files — message: `docs: use_grf @param corrected (GRF candidate-cut generation; VI ordering belongs to vi.grf.min)`.
6. **Report** per the skeleton; **Commit 2:** the report — message: `docs: roxygen use_grf correction report`. Print `git log --oneline -3` and stop. **No push.**

## Report skeleton

```
# REPORT — roxygen use_grf correction — 2026-08-31
[banner if stopped]
## 0 Environment (incl. untracked paths per G1)
## 1 G2 exact-text matches (per edit: found at line N, once)
## 2 The diff (verbatim; G3 roxygen-only verification)
## 3 document() outcome (regenerated man/ files; G4 greps; any drift files)
## 4 NEWS.md context (heading before/after)
## 5 Commits
## 6 Notes (anything observed, no action taken)
```

## Hard rules
- Do not touch the `@param vi.grf.min` or `@param max_n_confounders` blocks — they are already correct.
- Do not edit message strings, code lines, tests, or anything under `quarto/`; roxygen comments, `NEWS.md`, and generated `man/` only.
- Do not commit, move, or delete the scratch dirs `dev/tasks/_baseline_html_2026-08-31/` and `dev/tasks/_diffs_2026-08-31/`, nor the untracked `comparison_continuous.rds`.
- Verify from source; no push; no follow-up work on anything noted in §6.
