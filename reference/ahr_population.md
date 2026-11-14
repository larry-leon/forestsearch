# The population of the analytic layer

Builds, from a Weibull AFT DGM made by
[`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md),
the per-patient quantities the analytic layer works on: each
super-population patient's log hazard under control (\\\theta_0\\) and
log hazard ratio (\\\psi\\, so \\\theta_1 = \theta_0 + \psi\\), the H\*
(`flag_harm`) and region indicators, the Weibull `mu` and `tau`, and
each patient's entry window: entry at \\a_i + E\\, \\E \sim U(0, R)\\,
analysis at calendar time `analysis_time`.

## Usage

``` r
ahr_population(dgm, analysis_time, max_entry, entry_delay = NULL)

# S3 method for class 'ahr_population'
print(x, ...)
```

## Arguments

- dgm:

  An `aft_dgm_flex` object made by
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md).

- analysis_time:

  Calendar time of the analysis, in the DGM's time units.

- max_entry:

  Length R of each patient's uniform entry window.

- entry_delay:

  `NULL` (no delay), a non-negative numeric vector of length
  `nrow(dgm$df_super)` (each patient's delay \\a_i\\), or
  `list(var = <column>, delay = a)`: delay `a` for patients with
  `dgm$df_super[[var]] == 1`, none otherwise (e.g.
  `list(var = "z_region", delay = 6)`: the region enrols 6 months
  later). The analytic counterpart of
  `simulate_from_dgm(entry_delay = )`.

- x:

  An `ahr_population`.

- ...:

  Unused.

## Value

An object of class `ahr_population`: a list with `theta0`, `psi`, `band`
(H\*), `region`, `delay`, `mu`, `tau`, `x` (`analysis_time`), `R`
(`max_entry`), `N`, the distinct `delays`, and an environment holding
the quadrature grids and Laplace-sum caches (keyed by horizon).

## Details

The DGM's `theta_0`, `theta_1` and `loghr_po` columns are checked to be
consistent with its linear predictors and with each other (as the
mrct-consistency `ahr_pop()` did); the layer reads the hazards from them
and uses `dgm$model_params$mu` as it stands, so apply
[`scale_dgm_hazard()`](https://larry-leon.github.io/forestsearch/reference/scale_dgm_hazard.md)
before building the population.

## Limits

Weibull AFT DGMs from
[`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md)
only; entry uniform on \\(a_i, a_i + R)\\ for each entry group;
administrative censoring only (the random censoring model of the DGM,
e.g. `cens_adjust`, is not in the layer); correlations exact for nested
sets. Each distinct delay is one quadrature group with its own breaks,
so a delay vector with many distinct values is slow. Derivations: the
technical report (`docs/manuscript/techreport/analytic-layer.qmd`,
Chapter 6, of the mrct-consistency repository) and the theory notes
(`docs/theory/AHR_mrct_dgm.qmd`, `docs/theory/AHR_complete.qmd`).

## See also

[`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md),
[`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md),
[`ahr_solve()`](https://larry-leon.github.io/forestsearch/reference/ahr_solve.md),
[`scale_dgm_hazard()`](https://larry-leon.github.io/forestsearch/reference/scale_dgm_hazard.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- ahr_population(dgm, analysis_time = 40, max_entry = 11)
pop6 <- ahr_population(dgm, 40, 11, entry_delay = list(var = "z_region", delay = 6))
} # }
```
