# Build an MRCT data-generating mechanism, keeping its warnings

A thin wrapper:
[`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md),
then
[`scale_dgm_hazard()`](https://larry-leon.github.io/forestsearch/reference/scale_dgm_hazard.md)
when `m != 1`. The build's warnings (usually the censoring model's fit)
are muffled and kept in the `"build_warnings"` attribute.

## Usage

``` r
mrct_build_dgm(data, ..., m = 1)
```

## Arguments

- data:

  The seed data: `seed_data` of
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md).

- ...:

  Passed to
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md).

- m:

  The hazard multiplier of
  [`scale_dgm_hazard()`](https://larry-leon.github.io/forestsearch/reference/scale_dgm_hazard.md).
  Default 1, no scaling.

## Value

The DGM, with the attribute `"build_warnings"`: a character vector,
`character(0)` for a clean build.

## See also

[`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md),
[`scale_dgm_hazard()`](https://larry-leon.github.io/forestsearch/reference/scale_dgm_hazard.md).

## Examples

``` r
set.seed(1)
n <- 400
d <- data.frame(age = rnorm(n, 50, 10), lab = rnorm(n, 300, 100),
                f1 = rbinom(n, 1, 0.4), treat = rbinom(n, 1, 0.5))
d$tte <- rexp(n, 0.05 * exp(-0.5 * d$treat))
d$event <- as.integer(d$tte < 40)
d$tte <- pmin(d$tte, 40)
dgm <- mrct_build_dgm(
  d, outcome_var = "tte", event_var = "event", treatment_var = "treat",
  continuous_vars = c("age", "lab"), factor_vars = "f1", x_pred = "lab",
  spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
  region = list(prevalence = 0.25, or_pred = 1, x1_vars = "f1", or_x1 = 3,
                loghr = -log(3), band = list(cut = 230, or = 1)),
  n_super = 2000, expand = "bootstrap", seed = 7L,
  continuous_vars_cens = c("age", "lab"), factor_vars_cens = "f1",
  m = 3)
length(attr(dgm, "build_warnings"))
#> [1] 9
```
