# REPORT — roxygen use_grf correction — 2026-08-31

All gates passed. Documentation-only: three roxygen edits, one NEWS.md bullet, two regenerated `man/` files. Zero behavior change — every added/removed line under `R/` begins with `#'`. No push.

## 0 Environment (incl. untracked paths per G1)

- Host `pop-os`; branch `feature/glm-extension`; installed `packageVersion("forestsearch")` = `0.3.1`.
- `git log -1 --oneline` at start of execution: `26cfce21 tasks: add roxygen use_grf task (2026-08-31 v3)` — the task-doc commit made on arrival atop `13d5987e` (the push-size-fix report commit the v3 note expects).
- **G1 pass:** `git status --porcelain --untracked-files=no` → empty (0 lines). Untracked entries were exactly the three expected paths, no others:
```
?? dev/tasks/_baseline_html_2026-08-31/
?? dev/tasks/_diffs_2026-08-31/
?? quarto/applications/actg175/_payloads_2026-08-31/analysis_actg175_continuous_compare_all/comparison_continuous.rds
```

## 1 G2 exact-text matches

| edit | file | OLD found at | occurrences |
|---|---|---|---|
| E1 | `R/forestsearch_main.R` | L255 (`#' @param use_grf Logical. Use GRF for variable importance. Default FALSE.` + the two continuation lines, verbatim) | 1 (`grep -cF "Use GRF for variable importance"` = 1) |
| E2 | `R/forestsearch_main.R` | L690–691 (anchor pair `behavior. Search scope is governed by \code{maxk}, candidate-factor` / `screening (\code{use_grf}, ...)`) | 1 each |
| E3 | `R/interpret_search_config.R` | L31, the exact whole line `#' @param use_grf Logical.` | 1 |
| E4 | `NEWS.md` | no development heading existed — top of file was `# forestsearch 0.3.1` — so `# forestsearch (development version)` was created above it with the bullet under it | — |

Edits were applied by a count-verified exact-string replacement (each OLD asserted to occur exactly once before writing; any other count would have aborted → STOP).

## 2 The diff (verbatim; G3)

`git diff --stat` (pre-`document()`):
```
 NEWS.md                     |  8 ++++++++
 R/forestsearch_main.R       | 14 +++++++++++---
 R/interpret_search_config.R |  3 ++-
 3 files changed, 21 insertions(+), 4 deletions(-)
```

`git diff -- R/`:
```
diff --git a/R/forestsearch_main.R b/R/forestsearch_main.R
@@ -252,7 +252,15 @@
 #'   candidate family data-dependent, which multiplier resampling
 #'   (\code{mr_inference = TRUE}) requires be fixed.  Set \code{TRUE} to
 #'   restore the prognostic-lasso prefilter used in Leon et al. (2024).
-#' @param use_grf Logical. Use GRF for variable importance. Default FALSE.
+#' @param use_grf Logical.  Generate additional candidate cuts from a GRF
+#'   fit: tree-derived cutpoints are extracted via
+#'   \code{grf.subg.harm.survival()} / \code{grf.subg.harm.glm()} (or taken
+#'   from \code{grf_res} / \code{grf_cuts} when supplied) and enter the
+#'   candidate family through \code{\link{get_FSdata}} alongside the
+#'   quartile (and any DINA) cuts.  This flag governs candidate-cut
+#'   generation only: GRF variable-importance ordering of the cut columns
+#'   is a separate step controlled solely by \code{vi.grf.min} and runs,
+#'   or is skipped, regardless of \code{use_grf}.  Default FALSE.
 #'   The default changed from \code{TRUE} in an earlier release, for the same
 #'   fixed-candidate-family reason given under \code{use_lasso}.
 #' @param grf_res GRF results object (optional, for reuse).
@@ -687,8 +695,8 @@
 #'   deprecation in v0.3.0.} Previously intended as a wall-clock time
 #'   budget for the combination search, this argument is no longer
 #'   enforced in the parallelized search path and has no effect on
-#'   behavior. Search scope is governed by \code{maxk}, candidate-factor
-#'   screening (\code{use_grf}, \code{use_lasso}, \code{conf_force}),
+#'   behavior. Search scope is governed by \code{maxk}, the candidate-factor
+#'   front ends (\code{use_grf}, \code{use_lasso}, \code{conf_force}),
 #'   and the number of parallel workers. The default of 3 is retained
 #'   only for signature compatibility.
 #' @param minp Numeric. Minimum prevalence threshold. Default 0.025.
diff --git a/R/interpret_search_config.R b/R/interpret_search_config.R
@@ -28,7 +28,8 @@
 #'   log scale for ratio measures, identity for others).
 #' @param consistency_threshold Numeric. Resolved consistency threshold.
 #' @param use_lasso Logical.
-#' @param use_grf Logical.
+#' @param use_grf Logical.  Whether GRF candidate-cut generation is on;
+#'   see \code{\link{forestsearch}}.
 #' @param outcome.name Character. Name of outcome column.
 #' @param event.name Character. Name of event column.
 #' @param treat.name Character. Name of treatment column.
```

**G3 pass:** `git diff -U0 -- R/ | grep '^[+-]' | grep -v '^+++\|^---' | grep -v "^[+-]#'"` → empty. Every added/removed line under `R/` begins with `#'`. The `@param vi.grf.min` and `@param max_n_confounders` blocks were not touched.

## 3 document() outcome (G4)

`Rscript -e 'devtools::document()'` succeeded; wrote exactly `forestsearch.Rd` and `interpret_search_config.Rd` — **no additional `man/` files regenerated, so no pre-existing drift**.

G4 greps on `man/forestsearch.Rd`:
- `"Use GRF for variable importance"` → **0 hits** (gone, as required).
- `"candidate-cut generation"` → 0 hits **as a single-line grep only**: Rd line-filling wraps the phrase across L163–164 (`…This flag governs candidate-cut` / `generation only: GRF variable-importance ordering…`). The full E1 text is present verbatim in the `\item{use_grf}` block (Rd L158–168, quoted below from the man diff). G4's substance is satisfied; recorded here so the gate's literal single-line grep isn't misread as a failure.

`git diff -- man/` (complete):
```
-\item{use_grf}{Logical. Use GRF for variable importance. Default FALSE.
+\item{use_grf}{Logical.  Generate additional candidate cuts from a GRF
+fit: tree-derived cutpoints are extracted via
+\code{grf.subg.harm.survival()} / \code{grf.subg.harm.glm()} (or taken
+from \code{grf_res} / \code{grf_cuts} when supplied) and enter the
+candidate family through \code{\link{get_FSdata}} alongside the
+quartile (and any DINA) cuts.  This flag governs candidate-cut
+generation only: GRF variable-importance ordering of the cut columns
+is a separate step controlled solely by \code{vi.grf.min} and runs,
+or is skipped, regardless of \code{use_grf}.  Default FALSE.
 The default changed from \code{TRUE} in an earlier release, for the same
 fixed-candidate-family reason given under \code{use_lasso}.}
   (man/forestsearch.Rd @@ -155,7 +155,15 @@)

-behavior. Search scope is governed by \code{maxk}, candidate-factor
-screening (\code{use_grf}, \code{use_lasso}, \code{conf_force}),
+behavior. Search scope is governed by \code{maxk}, the candidate-factor
+front ends (\code{use_grf}, \code{use_lasso}, \code{conf_force}),
   (man/forestsearch.Rd @@ -634,8 +642,8 @@)

-\item{use_grf}{Logical.}
+\item{use_grf}{Logical.  Whether GRF candidate-cut generation is on;
+see \code{\link{forestsearch}}.}
   (man/interpret_search_config.Rd @@ -35,7 +35,8 @@)
```

## 4 NEWS.md context

Before: file began `# forestsearch 0.3.1` (no development heading). After — heading created and bullet added; version heading follows immediately:
```
# forestsearch (development version)

* Documentation: corrected the `use_grf` `@param`, which described the
  variable-importance role belonging to `vi.grf.min`. `use_grf` governs GRF
  candidate-cut generation only; variable-importance ordering/screening of
  the cut columns is controlled solely by `vi.grf.min` (`NULL` default since
  0.3.0). No behavior change.

# forestsearch 0.3.1

## The consistency screen treats `plan = "sequential"` as sequential
```

## 5 Commits

| sha | message |
|---|---|
| `26cfce21` | `tasks: add roxygen use_grf task (2026-08-31 v3)` |
| `ec116ffe` | `docs: use_grf @param corrected (GRF candidate-cut generation; VI ordering belongs to vi.grf.min)` — `R/forestsearch_main.R`, `R/interpret_search_config.R`, `NEWS.md`, `man/forestsearch.Rd`, `man/interpret_search_config.Rd` |
| (report commit — see final `git log`) | `docs: roxygen use_grf correction report` |

No push. Scratch dirs and the untracked `comparison_continuous.rds` untouched.

## 6 Notes (no action taken)

1. G4's positive check assumed the phrase would survive Rd line-filling on one line; it wrapped (L163–164). A future gate of this shape could grep with `-z` or on a whitespace-normalized stream.
2. `devtools::document()` printed only the two expected `Writing` lines — `man/` was fully in sync with roxygen before this task (no drift), and `NAMESPACE`/`DESCRIPTION` were untouched.
3. The installed 0.3.1 package's help still carries the old `use_grf` text until the package is reinstalled; the source tree and `man/` are corrected. Not acted on (docs-only task; no install was requested).
