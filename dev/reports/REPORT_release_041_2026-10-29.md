# REPORT — release 0.4.1 (2026-10-29)

Task: `dev/tasks/cc_task_release_041_2026-10-29.md`. Branch `release/0.4.1`
created from `master` at `624bd457fbdf6136bb956b259e041558f2327577`. Nothing
fetched, pulled, pushed or merged.

## 1. Provenance (GATE — passed)

- `HEAD` on `master` = `624bd457fbdf6136bb956b259e041558f2327577`; no tracked
  modifications; `git merge-base --is-ancestor 624bd457 HEAD` exit 0.
- `grep -c 'actions/workflows/R-CMD-check' README.md` → `0`;
  `grep -n '^Version:' DESCRIPTION` → `3:Version: 0.4.0.9000`.
- `pandoc --version` → `pandoc 3.1.3` (on `PATH`).

## Commits

| SHA | Message |
|---|---|
| `6d2eea1c` | Task: release 0.4.1 (as received) |
| `5e1e78f8` | Release 0.4.1: version and NEWS heading |
| `d814ec54` | Release 0.4.1: documentation fixes (ci_method default, article reference, suggest_thresholds example) |
| `a24b07f2` | README: lifecycle badge stable |
| `dfd3016f` | cran-comments for 0.4.1 |
| (this report) | Report and record (see `git log`) |

`devtools::document()` rewrote exactly `man/forestsearch.Rd`,
`man/mr_estimates_table.Rd`, `man/suggest_thresholds.Rd`; no `NAMESPACE`
change and no `DESCRIPTION` roxygen re-stamp (six files committed).
`grep -n "survival-post-selection vignette" R/*.R man/*.Rd` printed nothing.
New lines in `R/forestsearch_main.R` are 65–75 columns.

## Diffs

```
 DESCRIPTION                                 |   2 +-
 NEWS.md                                     |   2 +-
 R/forestsearch_main.R                       |   9 +-
 R/mr_estimates_table.R                      |   6 +-
 R/suggest_thresholds.R                      |   8 +-
 README.md                                   |   2 +-
 cran-comments.md                            |  57 ++++----
 dev/tasks/cc_task_release_041_2026-10-29.md | 199 ++++++++++++++++++++++++++++
 man/forestsearch.Rd                         |   9 +-
 man/mr_estimates_table.Rd                   |   6 +-
 man/suggest_thresholds.Rd                   |   8 +-
 11 files changed, 256 insertions(+), 52 deletions(-)
```

```diff
diff --git a/DESCRIPTION b/DESCRIPTION
--- a/DESCRIPTION
+++ b/DESCRIPTION
@@ -1,6 +1,6 @@
 Package: forestsearch
 Title: Exploratory Subgroup Identification in Survival and GLM Outcomes
-Version: 0.4.0.9000
+Version: 0.4.1
 Authors@R: person("Larry", "Leon", , "larry.leon.05@post.harvard.edu", role = c("aut", "cre"))
 Description: Implements statistical methods for exploratory subgroup 
     identification in clinical trials. Provides tools for identifying patient 
diff --git a/NEWS.md b/NEWS.md
--- a/NEWS.md
+++ b/NEWS.md
@@ -1,4 +1,4 @@
-# forestsearch (development version)
+# forestsearch 0.4.1
 
 ### Deprecated
 
diff --git a/R/forestsearch_main.R b/R/forestsearch_main.R
--- a/R/forestsearch_main.R
+++ b/R/forestsearch_main.R
@@ -1432,10 +1432,11 @@
 #' @examples
 #' \dontrun{
 #' ## Example 1 -- Survival (GBSG, hormone therapy): the analysis of
-#' ## León and Anderson (2026a, Sec. 5) and the survival-post-selection
-#' ## vignette. Forest search over the cut grid (no model-based front ends),
-#' ## effMaxSG selection in a 10% effect neighbourhood, multiplier-resampling
-#' ## inference for the reported subgroup computed inside the fit.
+#' ## León and Anderson (2026a, Sec. 5) and of the survival-post-selection
+#' ## article on the package website. Forest search over the cut grid (no
+#' ## model-based front ends), effMaxSG selection in a 10% effect
+#' ## neighbourhood, multiplier-resampling inference for the reported
+#' ## subgroup computed inside the fit.
 #' library(survival)
 #' df_gbsg <- within(gbsg, {
 #'   id          <- seq_len(nrow(gbsg))
diff --git a/R/mr_estimates_table.R b/R/mr_estimates_table.R
--- a/R/mr_estimates_table.R
+++ b/R/mr_estimates_table.R
@@ -28,8 +28,10 @@
 #' robust to row ordering.  The MR de-biased confidence interval is a
 #' leading-order approximation of the full-bootstrap infinitesimal-jackknife
 #' interval -- it comes from the IJ variance computed on the multiplier draws
-#' under the default `ci_method = "ij"`, and from the subgroup robust SE only
-#' under `ci_method = "wald"`.  Either way the bias-corrected cell may differ
+#' under `ci_method = "ij"` or `"field"` (the forest-search path defaults to
+#' `"field"`, which adds the field bounds; the DINA and GRF paths default to
+#' `"ij"`), and from the subgroup robust SE only under `ci_method = "wald"`.
+#' Either way the bias-corrected cell may differ
 #' slightly from the full bootstrap even when the naive effect matches exactly.
 #'
 #' @seealso [forestsearch()] for the `mr_inference` switch and the vocabulary
diff --git a/R/suggest_thresholds.R b/R/suggest_thresholds.R
--- a/R/suggest_thresholds.R
+++ b/R/suggest_thresholds.R
@@ -86,10 +86,10 @@
 #' d <- d_eff_survival(n_sg = 60, prop_cens = 0.45)
 #' suggest_thresholds(d_eff = d, N = 700)
 #'
-#' # Binary: needs higher thresholds
-#' d <- d_eff_binary(n_sg = 150, p_event = 0.30)
-#' suggest_thresholds(d_eff = d, N = 500,
-#'                    L_eff_C = 0.220, L_eff_alpha = 1.298)
+#' # Binary: needs higher thresholds (shown, not run -- about 10 s)
+#' # d <- d_eff_binary(n_sg = 150, p_event = 0.30)
+#' # suggest_thresholds(d_eff = d, N = 500,
+#' #                    L_eff_C = 0.220, L_eff_alpha = 1.298)
 #'
 #' # Use calibration object from run_null_calibration()
 #' # cal <- run_null_calibration(...)
diff --git a/README.md b/README.md
--- a/README.md
+++ b/README.md
@@ -2,7 +2,7 @@
 # forestsearch <img src="man/figures/logo.png" align="right" height="139" alt="" />
 
 <!-- badges: start -->
-[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
+[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
 [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
 <!-- badges: end -->
 
```

## 6. Install, URLs, spelling, reverse dependencies

- `devtools::install(upgrade = FALSE)`: `Version` `0.4.1`; `Built`
  `R 4.6.1; ; 2026-10-29 18:10:20 UTC; unix`.
- `urlchecker::url_check()` (installed) — **GATE passed**: no URL with an
  error status. Three entries, all redirects (recorded, not gated), all in
  `vignettes/articles/references.bib` (excluded from the tarball by
  `.Rbuildignore`). `url_check()` defaults to `fail = TRUE` and errors with
  "Found 3 invalid URLs" when any row is listed, redirects included; it was
  re-run with `fail = FALSE` to record the table
  (`url_check_2026-10-29.csv`, `url_check.rds`):

| URL | Status | New |
|---|---|---|
| http://www.bibsonomy.org/bibtex/220316799a942085b678c6105d9a74dc8/dblp | 200 | https://www.bibsonomy.org/bibtex/220316799a942085b678c6105d9a74dc8/dblp |
| http://www.jthoracdis.com/article/view/149 | 200 | https://jtd.amegroups.org/article/view/149 |
| http://www.R-project.org (two places) | 200 | https://www.r-project.org/ |

- `spelling`: **not installed**; no spelling output (`spelling_2026-10-29.txt`
  not produced).
- Reverse dependencies (`tools::package_dependencies("forestsearch", reverse
  = TRUE, which = "all", db = available.packages())`, CRAN cloud mirror):
  `$forestsearch character(0)` (`revdeps_2026-10-29.txt`).

## 7. Pre-build hygiene

`vignettes/forestsearch.html` (git-ignored by `.gitignore:15
vignettes/*.html`; dated 2026-10-28) moved to
`dev/reports/release_041_2026-10-29/pre_build_render/`. Not a tracked change.
`vignettes/articles/` left as is.

## 8. Final check (GATE — passed)

`NOT_CRAN=true`, pandoc on `PATH`, `timeout 60m`:
`R CMD build . && R CMD check --as-cran --run-donttest forestsearch_0.4.1.tar.gz`,
11:11:32 → 11:23:56 PDT (12 min 24 s), exit 0. Log:
`check_as_cran_release_2026-10-29.log`.

1. Tarball vignettes (GATE passed — no `.html`):
```
forestsearch/vignettes/
forestsearch/vignettes/forestsearch.qmd
```
2. `grep -n '^Status:' forestsearch.Rcheck/00check.log` → `70:Status: OK` (GATE passed).
3. Timing lines, verbatim:
```
* checking CRAN incoming feasibility ... [6s/14s] OK
* checking whether package ‘forestsearch’ can be installed ... [17s/17s] OK
* checking R code for possible problems ... [39s/39s] OK
* checking examples ... [41s/42s] OK
  Running ‘testthat.R’ [410s/394s]
 [410s/394s] OK
* checking re-building of vignette outputs ... [73s/72s] OK
* checking HTML version of manual ... [16s/16s] OK
Status: OK
```
   No "Examples ... took more than 5s" table (`grep -c 'took more than'` on
   the log → 0). Examples 41 s against 51 s in the morning run: the binary
   `suggest_thresholds()` call is no longer run.
4. `cran-comments.md` fields filled (Status OK, 0 errors | 0 warnings |
   0 notes); committed `dfd3016f`.
5. Tarball moved to `~/Downloads/forestsearch_0.4.1.tar.gz` (1368765 bytes);
   `forestsearch.Rcheck/` removed. Tarball not committed.

## Deviations and anchors

- Anchor moved: the `R/mr_estimates_table.R` text was at **L30–L32**, not
  L31–L33 (one line up). `R/forestsearch_main.R` L1434–L1436,
  `R/suggest_thresholds.R` L89–L92 and `README.md` L5 matched; NEWS and
  DESCRIPTION anchored by content.
- `urlchecker::url_check()` default errors on redirects; re-run with
  `fail = FALSE` (above). Gate judged on status: no error status.
- `spelling` not installed; §6 spelling step not done.
- `pre_build_render/forestsearch.html` is not git-ignored at its new location;
  it is a 1.4 MB render and was **left untracked**, not committed with the
  record directory.
- Flag, not edited: `cran-comments.md` says the examples take "about 50 s
  with `--run-donttest`"; this run measured 41 s. Text left as written by
  the task.
