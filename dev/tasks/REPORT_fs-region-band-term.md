# REPORT: band term in the region model of mrct_inject_dgm.R

Date: 2026-10-03
Task: `dev/tasks/TASK_fs-region-band-term.md` (committed as received, `8edb4d0b`)
Branch: `feature/glm-extension`. Not pushed. `R CMD check` not run (per task).

## Commit SHA to pin

**`c98ed047662f9d5635e07610b342f95d2a29aa7b`** (`c98ed047`) on `feature/glm-extension`.

The local install (`R CMD INSTALL .`) was made from this tree: forestsearch 0.3.5.9002. This
report is committed separately, after that SHA.

## Outcome

| Step | Result |
|---|---|
| `devtools::document()` | Clean; `man/make_region_model.Rd` and `man/inject_mrct_structure.Rd` rewritten, `NAMESPACE` unchanged |
| `R CMD INSTALL .` | `* DONE (forestsearch)`; installed version 0.3.5.9002, `band` is a formal of `make_region_model` |
| `devtools::test(filter = "mrct_inject_dgm")` | `[ FAIL 0 | WARN 6 | SKIP 0 | PASS 19 ]` |
| Default-unchanged check (`band = NULL`) | Identical to the previous code (see below) |
| `R CMD check` | Not run |

## Diff summary

```
 DESCRIPTION                           |  2 +-
 R/mrct_inject_dgm.R                   | 34 +++++++++++++++++++++++-----
 man/inject_mrct_structure.Rd          |  4 +++-
 man/make_region_model.Rd              |  9 ++++++--
 tests/testthat/test-mrct_inject_dgm.R | 42 +++++++++++++++++++++++++++++++++++
 5 files changed, 82 insertions(+), 9 deletions(-)
```

`R/mrct_inject_dgm.R`, `make_region_model()`:

- New last argument `band = NULL`.
- With `band = list(var, cut, or)`, the linear predictor becomes
  `a0 + s-terms + log(or) * 1{df[[var]] <= cut}`. `var` is used on its natural scale; it does
  not pass through `s()`.
- `alpha0` is solved by the same `uniroot` call, now on the linear predictor that includes the
  band term, so the target prevalence still holds.
- The returned list gains `band` (NULL when not given). `coefs` and `ors` are unchanged: they
  still hold only the `x_pred` and `x1_vars` terms.
- Input check: `band` must be a list with `var`, `cut`, `or`, and `var` must be a column of
  `pop`; otherwise `stop()`.

`R/mrct_inject_dgm.R`, `inject_mrct_structure()`:

- Reads `region$band = list(cut, or)`; `var` defaults to `x_pred`.
- Passes it to `make_region_model()`. There is one region model, built on the
  super-population; its `draw()` is applied to both the super-population and the seed rows, so
  both get the band term.
- Records the resolved `list(var, cut, or)` under `dgm$mrct$region$band` (NULL when not given).
  It is also at `dgm$mrct$region_model$band`.

Also: two lines added to the file's header comment describing the band term; roxygen `@param`
text for `band` and `region`; `@return` of `make_region_model` lists `band`.

`DESCRIPTION`: version 0.3.5.9001 -> 0.3.5.9002. No dependency change.

## Existing behaviour unchanged

The previous `R/mrct_inject_dgm.R` (from `8edb4d0b`) and the new one were run on the existing
cgd test call (`or_pred = 10`, `x1_vars = "female"`, `or_x1 = 3`, `n_super = 5000`, `seed = 1`).
With `identical()`:

| Object | Identical |
|---|---|
| `dgm$df_super` | TRUE |
| `dgm$mrct$region$alpha0` | TRUE |
| `dgm$hazard_ratios` | TRUE |
| `dgm$model_params$gamma` | TRUE |

The only difference in the default case is one extra name, `band` (value NULL), in
`dgm$mrct$region` and in the list returned by `make_region_model()`.

## Added test

`tests/testthat/test-mrct_inject_dgm.R`, second `test_that` block, on `survival::cgd` (first
enrolment rows, as in the existing test):

- `region = list(prevalence = 0.20, or_pred = 1, band = list(cut = 10, or = 20))`
- `x3 = list(vars = "age", cuts = list(age = 10), loghr = log(5))`
- flat spline: `spline_spec = list(knot = 12, zeta = 25, log_hrs = rep(log(0.70), 3))`
- `x_pred = "age"`, `n_super = 5000`, `expand = "copula"`, `seed = 1`

The task did not give the flat spline's level; `log(0.70)` is my choice.

| Assertion | Value on this run |
|---|---|
| Region prevalence in [0.17, 0.23] | 0.2066 |
| `flag_harm == (z_age <= 10)` on `df_super` | TRUE for all 5000 rows |
| Band share larger in the region than outside | 0.8867 vs 0.3063 (overall 0.4262) |
| `dgm$mrct$region$band` equals `list(var = "age", cut = 10, or = 20)` | TRUE |
| `dgm$mrct$region_model$band` equals the above | TRUE |

The last two assertions go beyond the three the task named; they check the recording the task
asks for. Solved `alpha0` on this call: -3.3108.

## Test output

```
[ FAIL 0 | WARN 6 | SKIP 0 | PASS 19 ]
```

Before this task the file gave `[ FAIL 0 | WARN 3 | SKIP 0 | PASS 13 ]`. The new block adds 6
expectations and 3 warnings.

All six warnings are the same message from `survival::survreg.fit`, "Ran out of iterations and
did not converge", from the Weibull censoring-model fit in `prepare_censoring_model()`
(`R/generate_aft_dgm_helpers.R:312`) on the cgd seed. Each `inject_mrct_structure()` call gives
three: two inside `generate_aft_dgm_flex()` (`R/generate_aft_dgm_main.R:358` and `:393`) and one
from the module's own call (`R/mrct_inject_dgm.R:380`). This is the same non-converged
censoring fit noted in `REPORT_fs-add-mrct-inject-module.md`; it is not caused by the band term
and was not touched.

## Versions

- R 4.6.1 (2026-06-24)
- forestsearch 0.3.5.9002 (bumped from 0.3.5.9001)
- roxygen2 8.1.0
- latentcor: not installed on this machine, so both cgd tests use the normal-score correlation
  fallback in `expand_covariates()`.

## Side issues, not fixed

- The cgd censoring model does not converge (above). Any cgd-based DGM from this module takes
  its censoring parameters from a non-converged Weibull fit.
- `band` supports only the `<= cut` direction, as specified. A `>` band or a two-sided range
  would need a further argument.
- `or_pred = 1` with a `band` on `x_pred` still carries the `s(x_pred)` term with coefficient
  `log(1) = 0`; it has no effect, but `x_pred` cannot be dropped from the region model through
  `inject_mrct_structure()`, which always passes it.
