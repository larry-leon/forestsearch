# Bagged DINA with IJ variance (data-frame interface)

`dina_bagged()` is the high-level user-facing wrapper around
[`dina_fit_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_fit_bagged.md).
It accepts a data frame plus column-name arguments and dispatches to
[`dina_fit_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_fit_bagged.md),
which fits the DINA estimator on `n_bags` bootstrap replicates and
returns the infinitesimal-jackknife variance per Wager, Hastie and Efron
(2014).

## Usage

``` r
dina_bagged(
  df,
  outcome,
  treatment,
  covariates,
  family = c("gaussian", "binomial", "poisson", "cox"),
  status = NULL,
  ...
)
```

## Arguments

- df:

  data frame (or tibble / `data.table`) containing the outcome,
  treatment, and covariate columns.

- outcome:

  character(1); the name of the response column in `df`. For
  `family = "cox"` this is the event/censoring time column.

- treatment:

  character(1); the name of the binary treatment indicator column in
  `df`. Values must be 0/1 or logical.

- covariates:

  character vector of length `>= 1`; the names of the covariate columns
  in `df`. All must be numeric.

- family:

  one of `"gaussian"`, `"binomial"`, `"poisson"`, `"cox"`.

- status:

  character(1) or `NULL`; for `family = "cox"` only, the name of the
  event indicator column in `df` (`1` = event, `0` = censored). Ignored
  for non-Cox families.

- ...:

  further arguments forwarded to
  [`dina_fit_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_fit_bagged.md).

## Value

An object of class `c("dina_bagged", "dina")`, as returned by
[`dina_fit_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_fit_bagged.md).
The `$call` component is overwritten with the original `dina_bagged()`
call.

## Details

All bagging-related arguments (`n_bags`, `cross_fitting_per_bag`,
`n_folds_per_bag`, `parallel`, `ij_finite_sample_correction`,
`project_psd`, `verbose`, ...) are forwarded via `...`. See
[`dina_fit_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_fit_bagged.md)
for full documentation.

For `parallel = "bags"`, set
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
to a parallel backend *before* calling `dina_bagged()`. Restore with
`future::plan("sequential"); gc()` afterwards.

## See also

[`dina_fit_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_fit_bagged.md)
for the underlying matrix-based bagged estimator;
[`dina()`](https://larry-leon.github.io/forestsearch/reference/dina.md)
for the single-pass DINA with sandwich variance.

## Examples

``` r
set.seed(1)
n <- 400
df_demo <- data.frame(
  y         = NA,
  w         = stats::rbinom(n, 1, 0.5),
  x1        = stats::runif(n, -1, 1),
  x2        = stats::runif(n, -1, 1),
  x3        = stats::runif(n, -1, 1)
)
tau_x      <- 0.5 + 0.8 * df_demo$x1 - 0.3 * df_demo$x2
df_demo$y  <- 1 + df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)

fit <- dina_bagged(
  df         = df_demo,
  outcome    = "y",
  treatment  = "w",
  covariates = c("x1", "x2", "x3"),
  family     = "gaussian",
  n_bags     = 30L,
  seed       = 1L
)
#> Warning: IJ variance had 1 negative eigenvalue(s); clipped to zero. Consider increasing `n_bags`.
coef(fit)
#> (Intercept)          x1          x2          x3 
#>   0.3929142   0.7464394  -0.1283292  -0.1705473 
sqrt(diag(vcov(fit)))
#> (Intercept)          x1          x2          x3 
#>  0.02663329  0.19603330  0.25479381  0.25317176 
```
