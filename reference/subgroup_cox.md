# Cox per-subgroup fit function

Returns a function `data -> c(HR, UB)` fitting `formula` by
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) and
extracting the first row of `summary(fit)$conf.int`: the hazard-ratio
point estimate (`"exp(coef)"`) and the upper 95% bound (`"upper .95"`).
Estimation failure or an empty coefficient table returns `c(NA, NA)`.
This reproduces the vignettes' `cox_sR()` exactly when called with
`survival::Surv(y_sim, event_sim) ~ treat_sim + survival::strata(grade)`.

## Usage

``` r
subgroup_cox(formula, lean = TRUE)
```

## Arguments

- formula:

  A model formula for
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html),
  evaluated against each subgroup's data.

- lean:

  Re-home the formula environment as described above (default `TRUE`).

## Value

A function of one argument (`data`) returning `c(hr, ub)`; carries the
formula in `attr(, "formula")`.

## Details

Custom fitters may be passed to
[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md)
directly: any function of a data frame returning a length-2 numeric
`(estimate, upper bound)` is accepted. Fitters run on parallel workers,
so they should reference non-base functions by `pkg::fun` or rely on
packages listed in `future_packages`.

## Formula environment (`lean`)

A formula captures the environment where it was created. In a script or
a knitr chunk that is the global environment, and `future`'s globals
machinery gives formulas special treatment: a fit whose formula points
at a populated environment imposes a multi-second per-future cost on
parallel runs (measured ~2-4 s x ~116 futures in the Phase 2 timing
diagnostic), independent of the serialized size. With `lean = TRUE`
(default) the formula is re-homed into a fresh, empty environment
parented on the survival namespace, cutting the costly chain while
keeping
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html),
[`survival::strata()`](https://rdrr.io/pkg/survival/man/strata.html),
and base/stats symbols resolvable. Consequence: the formula must
reference only data columns and functions reachable from the survival
namespace (use `pkg::` prefixes for anything else). Set `lean = FALSE`
to retain the calling environment – restoring arbitrary-symbol
resolution at that parallel cost. Independently of `lean`, the returned
closure is stripped of source references and given a minimal environment
holding only the formula (keep-source installs otherwise attach ~250 KB
of package-source baggage to every worker shipment).

## Examples

``` r
fit <- subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim)
```
