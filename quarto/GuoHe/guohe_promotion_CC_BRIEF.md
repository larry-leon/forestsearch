---
title: "Brief: promote Guo & He functions into the forestsearch package"
bibliography: []
---

# Promote the Guo & He analysis functions from `quarto/GuoHe/` into `R/`

The Guo & He (2021) de-biasing method is becoming a core evaluated method, so its
implementation should be package API rather than sourced scratch scripts.

## Scope

**Promote these three files** from `quarto/GuoHe/` to `R/` (git mv, preserving
history):

| file | public functions | internal helpers |
|---|---|---|
| `guohe_algorithm3.R` | `guohe_algorithm3()`, `print.guohe_a3` | `.g3_*` |
| `guohe_adaptive_r.R` | `guohe_adaptive_r()`, `print.guohe_ar` | `.ar_*` |
| `fs_to_guohe.R` | `fs_to_guohe()`, `fs_export_maxeff_family()`, `fs_assert_membership()` | `.fs_gh_*` |

**Do NOT promote** `guohe_reproduction_sim.R` or `guohe_reproduction_run.R`.
Those are paper-supporting validation of the Section 5 reproduction, not package
API. They stay in `quarto/GuoHe/`.

**Keep all function names exactly as they are.** No renaming — the committed
analysis qmds call `guohe_algorithm3()` and `fs_to_guohe()` by name, and
`fs_to_guohe()` is the primary user-facing entry point.

## Work items

1. **Move the three files** to `R/`.

2. **Roxygen.** Add complete roxygen to the five public functions:
   `@title`, `@description`, `@param` for *every* argument, `@return`,
   `@examples`, `@export`. Mark all `.g3_*` / `.ar_*` / `.fs_gh_*` helpers
   `@noRd`. Register the two S3 print methods
   (`@export` on `print.guohe_a3` / `print.guohe_ar`; roxygen emits
   `S3method()` in NAMESPACE).
   - Roxygen markdown is ON for this package: write literal `<`, `>`, `&` in
     roxygen text; never Rd-escape manually. `@section` titles must be plain.

3. **`@references`** on `guohe_algorithm3()`, `guohe_adaptive_r()`, and
   `fs_to_guohe()` — verbatim, DOI verified:

   ```
   #' @references
   #' Guo, X. and He, X. (2021). Inference on Selected Subgroups in Clinical
   #' Trials. \emph{Journal of the American Statistical Association},
   #' \strong{116}(535), 1498--1506. \doi{10.1080/01621459.2020.1740096}
   ```

   Note in the description that this is an independent implementation of the
   published method, validated by reproducing the paper's Section 5 simulations.

4. **Examples must run fast under `R CMD check`.** This is the main real work.
   Production settings (`B = 2000`) are far too slow. For each exported
   function write a small self-contained example — tiny simulated data, `B`
   around 50, a handful of candidates — that exercises the API and returns in
   well under a second. Do NOT use `\dontrun{}` to dodge this; a runnable
   example is the point. `\donttest{}` only if genuinely unavoidable.
   - `fs_to_guohe()` / `fs_export_maxeff_family()` / `fs_assert_membership()`
     need a `forestsearch(sg_focus = "maxeff")` fit. If a real fit is too slow
     for an example, wrap in `\donttest{}` and keep it minimal.

5. **DESCRIPTION.** Verify these are in **Imports** (not Suggests):
   `survival`, `stats`, `utils`, `future`, `future.apply`. `future.apply` is the
   one most likely to be mis-declared — `guohe_algorithm3()` calls
   `future.apply::future_lapply()` on its main path, so it must be Imports.
   Add any that are missing.

6. **Update the analysis qmds** — remove the now-redundant `source()` lines:
   ```
   source(file.path(gh_dir, "guohe_algorithm3.R"))
   source(file.path(gh_dir, "fs_to_guohe.R"))
   ```
   in `quarto/GuoHe/analysis_gbsg_cox_maxeff.qmd` and
   `analysis_gbsg_frozen_family.qmd` (and anywhere else that sources them).
   Leaving them in would shadow the installed package functions with sourced
   copies, which then silently diverge. Check whether `gh_dir` is still needed
   for anything else (e.g. the LOO cache path) before removing it.

7. **`document()`**, then **`R CMD check`** (or `devtools::check()`) and report
   the full result. Target zero ERRORs/WARNINGs; report any NOTEs.

8. **Run the full testthat suite** — confirm no regressions vs the baseline.

## Verification before committing

- `R CMD check` clean (report NOTEs).
- All examples execute (not skipped).
- `library(forestsearch); args(fs_to_guohe)` works from a *fresh* R session
  against the *installed* package.
- One of the analysis qmds still renders its GH chunk with the `source()` lines
  removed — proving the functions resolve from the package.
- testthat: no new failures.

## Conventions

- Show the diff before committing. Show `git status --short` to confirm only
  intended files are staged (no stray scratch/`.rds`/`.html`).
- Conventional Commits on `feature/glm-extension`, e.g.
  `feat: promote Guo & He de-biasing functions to package API`.
  The qmd `source()` removals can ride in the same commit or a follow-up —
  your call, but do not leave the qmds sourcing files that have moved.
- No PR; the owner commits/pushes via GitHub Desktop.
- `devtools::install()` (not `load_all()`) before anything exercising
  `doFuture`/`future` workers — workers only see the installed package.

## Out of scope

Do not refactor the estimator logic, change defaults, or alter numerical
behavior. This is a promotion: same code, now documented, exported, and checked.
Any behavioral change must be raised separately.
