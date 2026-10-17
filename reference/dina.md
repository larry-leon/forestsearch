# Estimate heterogeneous treatment effects via DINA (data-frame interface)

`dina()` is the high-level user-facing wrapper around
[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md).
It accepts a data frame plus column-name arguments, matching the
conventions used elsewhere in the `forestsearch` package, and extracts
the covariate matrix, response, and treatment indicator before
dispatching to
[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md).

## Usage

``` r
dina(
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
  [`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md).

## Value

An object of class `"dina"`, as returned by
[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md).
The `$call` component is overwritten with the original `dina()` call so
that [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) show the high-level
invocation rather than the internal matrix-based dispatch.

## Details

All other arguments are forwarded to
[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md)
via `...`, including `propensity_method`, `baseline_method`,
`cross_fitting`, `n_folds`, `cens_type`, `cens_params`, `n_grid`, `eps`,
and `seed`. See
[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md)
for their full documentation.

Covariate columns must be numeric (or integer). Convert factors to dummy
variables (e.g., with
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html))
before calling if needed. Rows with any `NA` value in the referenced
columns are rejected with an error; use
[`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html) or impute
beforehand.

## See also

[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md)
for the underlying matrix-based estimator;
[`dina_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_bagged.md)
for the bagged version with infinitesimal-jackknife variance.

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

fit <- dina(
  df         = df_demo,
  outcome    = "y",
  treatment  = "w",
  covariates = c("x1", "x2", "x3"),
  family     = "gaussian",
  seed       = 1L
)
coef(fit)
#> (Intercept)          x1          x2          x3 
#>   0.4798796   0.7702880  -0.1455416  -0.2867361 
confint(fit)
#>                  2.5 %     97.5 %
#> (Intercept)  0.2763565 0.68340270
#> x1           0.4274624 1.11311351
#> x2          -0.5106812 0.21959790
#> x3          -0.6352202 0.06174792
```
