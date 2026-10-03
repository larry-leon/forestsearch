# TASK (forestsearch): band term in the region model of mrct_inject_dgm.R

Date: 2026-10-03
Owner: Larry · Spec: chat · Execution: Claude Code (CC)
Repo: forestsearch, branch `feature/glm-extension`. No gates. No R CMD check. CC never pushes.

## 0. First action

Copy `~/Downloads/TASK_fs-region-band-term.md` to `dev/tasks/` and commit it.

## 1. R/ change — called out

| File | Change | Classification |
|---|---|---|
| `R/mrct_inject_dgm.R` — `make_region_model()` | new optional argument `band = NULL`; when `band = list(var, cut, or)` is given, the logit gains the term `log(or) * 1{df[[var]] <= cut}` alongside the existing `s(x_pred)` and `x1_vars` terms; `alpha0` is solved as now for the target prevalence; the returned list records `band` | existing function, **new argument, default leaves behaviour unchanged** |
| `R/mrct_inject_dgm.R` — `inject_mrct_structure()` | accepts `region$band = list(cut, or)` (and optional `var`, default `x_pred`), passes it to `make_region_model()` for both the super-population and the seed rows, and records it under `dgm$mrct$region$band` | same |
| `tests/testthat/test-mrct_inject_dgm.R` | one added test on `survival::cgd`: `region = list(prevalence = 0.20, or_pred = 1, band = list(cut = 10, or = 20))` with `x3 = list(vars = "age", cuts = list(age = 10), loghr = log(5))` (flat spline); assert prevalence in [0.17, 0.23], that `flag_harm == (z_age <= 10)` on `df_super`, and that the band's share is larger in the region than outside it | new test |

No other file changes. `devtools::document()`, install, `devtools::test(filter = "mrct_inject_dgm")`, commit, bump the dev version.

## 2. Report

`REPORT_fs-region-band-term.md` beside the other REPORT_* files: the diff summary, test output, and the commit SHA (the project re-pins its local install to it).
