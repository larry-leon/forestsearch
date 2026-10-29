# CC task — release 0.4.1: version, release-task edits, cran-comments, final check

**File:** `dev/tasks/cc_task_release_041_2026-10-29.md` · **Issued:** 2026-10-29 by chat, commissioned by Larry (plan: `HANDOFF_vignettes_release_0.4.1_2026-10-28_v2.md` §3 step 2, with the items recorded since: decision memo of 2026-10-28, robustness and example-timing tasks of 2026-10-29, the two checks of 2026-10-29, Larry's decisions of 2026-10-29 on the `suggest_thresholds` example and the lifecycle badge)
**Repo:** `forestsearch` · **Start on `master`** containing `624bd457` (README badge removal). CC creates **`release/0.4.1`** from it. CC never fetches, pulls, pushes or merges; Larry runs win-builder on the branch, merges in GitHub Desktop and submits on 10 November.
**Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone).
**Anchors:** line numbers are from `master` at `606a78ea` for files unchanged since (`R/mr_estimates_table.R`, `R/suggest_thresholds.R`, `R/forestsearch_main.R` around L1435, `README.md` L5, `cran-comments.md`); NEWS and DESCRIPTION anchors are by content. Verify every anchor by content and report any that moved.

## ⚠ CATEGORY — release metadata and documentation only

- `R/` edits are roxygen text only (three files; `man/` regenerates for exactly those three topics). No function body, formal or default changes. No behaviour change of any kind.
- `DESCRIPTION`: version only. `NEWS.md`: heading only. `README.md`: one badge line. `cran-comments.md`: rewritten for this submission.
- **Check policy for this task (the release gate):** CC runs the full `R CMD build` + `R CMD check --as-cran --run-donttest` once, as the last step (§8). No test suite or check before that.
- **Compute:** the one check, about 13 min on this machine (today's run: install 17 s, examples 51 s, tests 412 s, vignette 73 s), plus `urlchecker`/`spelling` seconds. Hard cap 60 min.
- **Unattended.** Gates stop the task; never ask.

## 1. Provenance and branch — GATE

1. `git rev-parse --abbrev-ref HEAD` = `master`; `git status --porcelain` has no tracked modifications (untracked dirt tolerated, never staged); `git merge-base --is-ancestor 624bd457 HEAD` exits 0. Record HEAD.
2. `grep -c 'actions/workflows/R-CMD-check' README.md` must print 0 (the badge is gone) and `grep -n '^Version:' DESCRIPTION` must print `Version: 0.4.0.9000`.
3. `git checkout -b release/0.4.1`; copy this document to `dev/tasks/cc_task_release_041_2026-10-29.md`; add by path; commit alone.
4. Record directory `dev/reports/release_041_2026-10-29/`; report `dev/reports/REPORT_release_041_2026-10-29.md`.
5. `pandoc --version` must succeed on `PATH` (installed 2026-10-29). GATE.

## 2. Version and NEWS heading

- `DESCRIPTION`: `Version: 0.4.0.9000` → `Version: 0.4.1`. Nothing else in the file (there is no `Date:` field; do not add one).
- `NEWS.md` line 1: `# forestsearch (development version)` → `# forestsearch 0.4.1`. Nothing else.

Commit both by path ("Release 0.4.1: version and NEWS heading").

## 3. Documentation-only `R/` fixes (roxygen text)

**(a) `R/mr_estimates_table.R` L31–33** — the `@details` sentence that says `"ij"` is the default. Replace

```
#' interval -- it comes from the IJ variance computed on the multiplier draws
#' under the default `ci_method = "ij"`, and from the subgroup robust SE only
#' under `ci_method = "wald"`.  Either way the bias-corrected cell may differ
```
with
```
#' interval -- it comes from the IJ variance computed on the multiplier draws
#' under `ci_method = "ij"` or `"field"` (the forest-search path defaults to
#' `"field"`, which adds the field bounds; the DINA and GRF paths default to
#' `"ij"`), and from the subgroup robust SE only under `ci_method = "wald"`.
#' Either way the bias-corrected cell may differ
```
(the following line `#' slightly from the full bootstrap ...` stays).

**(b) `R/forestsearch_main.R` L1434–L1436** — the example comment. Replace
```
#' ## Example 1 -- Survival (GBSG, hormone therapy): the analysis of
#' ## León and Anderson (2026a, Sec. 5) and the survival-post-selection
#' ## vignette. Forest search over the cut grid (no model-based front ends),
```
with
```
#' ## Example 1 -- Survival (GBSG, hormone therapy): the analysis of
#' ## León and Anderson (2026a, Sec. 5) and of the survival-post-selection
#' ## article on the package website. Forest search over the cut grid (no
#' ## model-based front ends),
```
and re-join the next line so the paragraph reads on (`effMaxSG selection ...`); keep lines under 80 columns. Afterwards `grep -n "survival-post-selection vignette" R/*.R man/*.Rd` must print nothing (after `document()`).

**(c) `R/suggest_thresholds.R` L89–L92** — Larry's decision: the binary call is shown, not run (it takes about 10 s, over CRAN's 5 s example limit). Replace
```
#' # Binary: needs higher thresholds
#' d <- d_eff_binary(n_sg = 150, p_event = 0.30)
#' suggest_thresholds(d_eff = d, N = 500,
#'                    L_eff_C = 0.220, L_eff_alpha = 1.298)
```
with
```
#' # Binary: needs higher thresholds (shown, not run -- about 10 s)
#' # d <- d_eff_binary(n_sg = 150, p_event = 0.30)
#' # suggest_thresholds(d_eff = d, N = 500,
#' #                    L_eff_C = 0.220, L_eff_alpha = 1.298)
```
The survival call above it stays live.

Then `devtools::document()`: `git status --porcelain` must show exactly the three `R/` files and `man/mr_estimates_table.Rd`, `man/forestsearch.Rd`, `man/suggest_thresholds.Rd` (plus a `DESCRIPTION` roxygen re-stamp only if it appears — accept it, it is the current roxygen2). No `NAMESPACE` change. Commit the six (or seven) files by path ("Release 0.4.1: documentation fixes (ci_method default, article reference, suggest_thresholds example)").

## 4. README lifecycle badge

`README.md` L5 (Larry's decision: `stable`):
```
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
```
→
```
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
```
Commit by path ("README: lifecycle badge stable").

## 5. `cran-comments.md`

Replace the whole file with the text below; the two `<<...>>` fields are filled in §8 from the check (CC), and the win-builder lines are left for Larry as written.

```
## Submission

This is an update from forestsearch 0.4.0 (the current CRAN version) to
0.4.1. NEWS.md lists the changes under "forestsearch 0.4.1":

* Three arguments of `forestsearch()` are deprecated (`vi.grf.min`,
  `max_n_confounders`, `max.minutes`): the formals are kept at the same
  positions and defaults, a non-default value warns once and has no effect.
* Two per-candidate fits call the fitting engines directly
  (`survival::coxph.fit()`, `stats::glm.fit()`); results are bit-identical
  to the formula-interface fits and are covered by identity tests. The
  package now requires survival (>= 3.5-8).
* Bug fixes: the log axis of `plot_subgroup_results_forestplot()` (#7); the
  weighted fit for `ps_adjust_method = "iptw"` with a continuous outcome.
* Documentation: a new getting-started vignette reproducing a published
  analysis; the longer post-selection analysis is now an article on the
  package website, so the package ships one vignette instead of two.

## Test environments

* Local: Pop!_OS 24.04 (x86_64-pc-linux-gnu), R 4.6.1, `R CMD check --as-cran
  --run-donttest` on the built tarball (vignette built): <<local status>>
* win-builder R-release: (Larry fills in)
* win-builder R-devel: (Larry fills in)

## R CMD check results

<<local status line, and each NOTE in full if any>>

Check time: the vignette rebuild takes about 75 s locally; the examples
about 50 s with `--run-donttest`; the test suite runs its long tests only
with `NOT_CRAN=true` (`skip_on_cran()`), so the tests take under a minute on
CRAN. CRAN's checks of 0.4.0 shipped two vignettes; 0.4.1 ships one, which
reproduces a published analysis at its settings.

## Reverse dependencies

There are no reverse dependencies (`tools::package_dependencies(reverse =
TRUE, which = "all")` against CRAN returns none).
```

Commit by path after §8 fills the fields ("cran-comments for 0.4.1").

## 6. `devtools::install(upgrade = FALSE)`; URL and spelling checks

- Install from the tree; record `packageDescription("forestsearch")[c("Version", "Built")]` (`0.4.1`).
- `urlchecker::url_check()` if `urlchecker` is installed (`requireNamespace`): **GATE** — no URL with an error status (404 or similar) may remain; redirects are recorded, not gated. If `urlchecker` is not installed, say so; the as-cran check's URL check (§8) stands in.
- `spelling::spell_check_package()` if `spelling` is installed: informational only — save the output as `spelling_2026-10-29.txt` in the record directory; do not edit anything for it (there is no `inst/WORDLIST`; Larry decides what, if anything, to fix).
- `tools::package_dependencies("forestsearch", reverse = TRUE, which = "all", db = available.packages())` — record the result (expected `character(0)`) for cran-comments.

## 7. Pre-build hygiene

- If a git-ignored render `vignettes/forestsearch.html` exists in the working tree, move it to `dev/reports/release_041_2026-10-29/pre_build_render/` (it is ignored, so this is not a tracked change; say so). `R CMD build` regenerates the vignette.
- The untracked renders under `vignettes/articles/` are excluded from the tarball by `.Rbuildignore` (`^vignettes/articles$`); leave them.

## 8. Final check — GATE (the release gate)

From the repo root, with `NOT_CRAN=true` and pandoc on `PATH`, inside one `timeout 60m`:

```sh
R CMD build . && R CMD check --as-cran --run-donttest forestsearch_0.4.1.tar.gz
```
console captured to `dev/reports/release_041_2026-10-29/check_as_cran_release_2026-10-29.log`.

1. `tar tzf forestsearch_0.4.1.tar.gz | grep '^forestsearch/vignettes/'` — must list the `.qmd` source and no `.html` (GATE: no stale render in the tarball).
2. `grep -n '^Status:' forestsearch.Rcheck/00check.log` — **GATE: `Status: OK`** (0 errors, 0 warnings, 0 notes). Anything else: stop and report every NOTE/WARNING in full; do not fix.
3. Record verbatim in the report: the Status line, the timing lines (`[..s/..s]`), the "Examples that took more than 5s" table if present (expected absent).
4. Fill the two `<<...>>` fields of `cran-comments.md` (local status `Status: OK`, 0 errors | 0 warnings | 0 notes) and commit it (§5).
5. Keep the tarball: move `forestsearch_0.4.1.tar.gz` to `~/Downloads/` (it is what `devtools::submit_cran()` rebuilds anyway; the copy is for Larry's record) and remove `forestsearch.Rcheck/`. Do not commit the tarball.

## 9. Report and post-conditions — GATE

Report `dev/reports/REPORT_release_041_2026-10-29.md`: provenance; every commit SHA; the diffs (`git diff --stat <base> HEAD`, and the full `git diff` of `DESCRIPTION`, `NEWS.md` line 1, `README.md`, the three roxygen hunks); the §6 outputs; the §8 figures verbatim; deviations; the closing block below. Commit the report and the record directory by path.

Closing block (run after the last commit; output pasted into the report in a final report-only commit):
```
git diff --name-only <base> HEAD | grep -v -E '^(dev/|NEWS\.md$|DESCRIPTION$|README\.md$|cran-comments\.md$|R/(mr_estimates_table|forestsearch_main|suggest_thresholds)\.R$|man/(mr_estimates_table|forestsearch|suggest_thresholds)\.Rd$)'
    # must print nothing
git diff <base> HEAD -- DESCRIPTION | grep -E '^[-+]Version'
    # must print exactly:  -Version: 0.4.0.9000   and   +Version: 0.4.1
git diff <base> HEAD -- NEWS.md | grep -E '^[-+]#'
    # must print exactly the two heading lines
git diff <base> HEAD -- R/ | grep -E '^[-+]' | grep -v -E "^[-+]#'" | grep -v -E '^(\+\+\+|---)'
    # must print nothing: every changed R/ line is a roxygen line
grep -rn "survival-post-selection vignette" R/ man/ ; grep -c 'actions/workflows/R-CMD-check' README.md
    # first must print nothing; second must print 0
grep -n '^Status:' dev/reports/release_041_2026-10-29/check_as_cran_release_2026-10-29.log
    # must print Status: OK
git status --porcelain | grep -v '^?? '
    # must print nothing
```
Copy the report, the check log and `cran-comments.md` to `~/Downloads`; print `ls -la` of the copies.

## 10. After this task (Larry, not CC)

`devtools::check_win_release()` and `devtools::check_win_devel()` on `release/0.4.1` (a few days before 10 November; fill the two win-builder lines of `cran-comments.md` from the result emails — a one-line edit each, direct or via CC); merge `release/0.4.1` into `master` in GitHub Desktop and push; on 10 November `devtools::submit_cran()` from `master`, confirm the email link. After acceptance: tag `v0.4.1`, bump to `0.4.1.9000` with a new `# forestsearch (development version)` NEWS heading.

## 11. Out of scope

Any function body; the `.gitignore` decision for `vignettes/articles/*.html` (handoff v1 §6 item 10 — tarball unaffected); the other open findings of that list (none blocks); the ACTG175 GLM article (after the release).
