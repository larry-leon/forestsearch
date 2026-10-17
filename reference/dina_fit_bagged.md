# Bagged DINA with infinitesimal-jackknife variance

`dina_fit_bagged()` fits the DINA estimator on `n_bags` bootstrap
replicates of the data and aggregates the resulting per-bag coefficient
vectors into a bagged estimator. The variance is estimated by the
infinitesimal jackknife (IJ) of Wager, Hastie and Efron (2014), with the
finite-sample-corrected form \$\$\widehat{V}\_{IJ-U} =
\widehat{V}\_{IJ} - (n / B) \widehat{V}\_{boot},\$\$ where
\\\widehat{V}\_{boot}\\ is the empirical covariance of the per-bag
coefficients. Negative eigenvalues introduced by the Monte-Carlo
correction (an artifact of small `n_bags`) are optionally clipped to
zero via PSD projection so that downstream
[`confint()`](https://rdrr.io/r/stats/confint.html) is always
well-defined.

## Usage

``` r
dina_fit_bagged(
  X,
  Y,
  W,
  family = c("gaussian", "binomial", "poisson", "cox"),
  propensity_method = "logistic",
  baseline_method = NULL,
  n_bags = 100L,
  cross_fitting_per_bag = TRUE,
  n_folds_per_bag = 2L,
  cens_type = NULL,
  cens_params = list(),
  n_grid = 500L,
  ij_finite_sample_correction = TRUE,
  project_psd = TRUE,
  eps = 1e-06,
  seed = NULL,
  parallel = c("none", "bags"),
  verbose = FALSE
)
```

## Arguments

- X:

  a numeric covariate matrix or data frame with `n` rows and `d`
  columns. All columns are treated as numeric predictors. Column names,
  when present, are propagated to the coefficient names of the fitted
  object.

- Y:

  the response. For `family = "cox"`, a data frame (or matrix coercible
  to one) with columns `time` (event/censoring time) and `status` (1 =
  event, 0 = censored). For all other families, a numeric vector of
  length `n`.

- W:

  binary treatment indicator (0 = control, 1 = treatment), length `n`.

- family:

  character; one of `"gaussian"`, `"binomial"`, `"poisson"`, `"cox"`.

- propensity_method:

  how to estimate the propensity score \\e(x) = P(W = 1 \mid X = x)\\.
  One of:

  - the string `"logistic"` (the default): a logistic GLM on `X`;

  - a numeric scalar in `(0, 1)`: a known constant propensity (e.g.,
    `0.5` for a balanced randomized trial);

  - a function taking a covariate matrix and returning a numeric vector
    of fitted propensities, allowing the user to plug in random forests,
    boosting, etc.

- baseline_method:

  how to estimate the natural-parameter functions \\\eta_0(x)\\,
  \\\eta_1(x)\\ (or the marginal mean for the Gaussian family). One of:

  - `NULL` (the default): use `"glm"` for non-Cox families and `"cox"`
    for the Cox family;

  - `"glm"`: a GLM with the natural link for the family;

  - `"cox"`:
    [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
    (Cox family only);

  - a list with components `eta0` and `eta1`, each a function taking a
    covariate matrix and returning a vector of linear predictors. For
    Gaussian, a list with a single component `nu` may be supplied
    instead.

- n_bags:

  positive integer giving the number of bootstrap bags (default `100L`).
  For stable IJ variance estimates the literature recommends `n_bags` on
  the order of the sample size `n`; smaller values are fine for
  diagnostic use but produce noisier variance.

- cross_fitting_per_bag:

  logical; if `TRUE` (the default) each bag uses `n_folds_per_bag`-fold
  cross-fitting internally, matching the paper. If `FALSE`, each bag
  uses a single in-bag fit (faster but loses the Chernozhukov debiasing
  step).

- n_folds_per_bag:

  positive integer giving the number of cross-fitting folds within each
  bag (default `2L`). Ignored when `cross_fitting_per_bag = FALSE`.

- cens_type:

  for `family = "cox"`, the censoring mechanism. One of `NULL` (no
  censoring), `"fixed"`, `"uniform"`, `"exponential"`, or `"unknown"`
  (estimate the probability of *not* being censored from `X` by logistic
  regression). Ignored for non-Cox families.

- cens_params:

  named list of censoring parameters. Required entries depend on
  `cens_type`:

  - `"fixed"`, `"uniform"`: `T0` (the administrative censoring time).

  - `"exponential"`: `rate_c` (the exponential rate).

  - `"unknown"`, `NULL`: ignored.

- n_grid:

  positive integer giving the number of time-grid points used for
  numerical integration of the not-censoring probability for
  `cens_type %in% c("uniform", "exponential")` (default `500`).

- ij_finite_sample_correction:

  logical; if `TRUE` (the default) the IJ variance includes the
  Monte-Carlo correction \\(n / B) \widehat{V}\_{boot}\\. If `FALSE` the
  raw IJ formula is returned, which tends to overstate variance for
  small `n_bags`.

- project_psd:

  logical; if `TRUE` (the default) the corrected IJ variance matrix is
  projected to its nearest positive semi-definite matrix by clipping
  negative eigenvalues to zero. Strongly recommended unless you want to
  inspect raw variance values yourself.

- eps:

  small positive constant added to denominators for numerical stability
  (default `1e-6`).

- seed:

  optional integer seed for the random cross-fitting fold assignment.
  When supplied, the global RNG state is restored on exit.

- parallel:

  character, one of `"none"` (default) or `"bags"`. When `"bags"`, the
  per-bag loop is dispatched via
  `foreach::foreach() %dofuture% { ... }`, using whatever
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  is currently active. Call e.g.
  `future::plan(future::multisession, workers = N)` *before* this
  function to enable actual parallel execution; otherwise the parallel
  branch will execute serially with a warning.

- verbose:

  logical; if `TRUE`, print a one-line progress message every 10 bags.
  Default `FALSE`. Per-bag progress messages are only emitted in
  `parallel = "none"` mode; in `parallel = "bags"` mode a single
  start-of-loop summary is printed if `verbose = TRUE`.

## Value

An object of class `c("dina_bagged", "dina")`, a list containing all the
components of a `"dina"` fit (so that the standard S3 methods work via
inheritance) plus:

- vcov:

  the corrected, PSD-projected IJ variance matrix (or the raw IJ if
  `project_psd = FALSE`).

- vcov_ij_raw:

  the IJ variance without the finite-sample correction.

- vcov_ij_corrected:

  the IJ variance with the finite-sample correction, before PSD
  projection.

- vcov_boot:

  the empirical bootstrap covariance of the per-bag coefficients,
  \\\widehat{V}\_{boot}\\.

- bag_coefficients:

  `n_bags_used` by `(d + 1)` matrix of per-bag coefficient vectors.

- bag_inclusion:

  `n` by `n_bags_used` integer matrix of bag inclusion counts,
  \\N\_{bi}^\*\\.

- n_bags, n_bags_used:

  requested and successfully-fit bag counts.

- n_folds_per_bag, cross_fitting_per_bag:

  bag-level cross-fitting settings.

## Details

The implementation calls
[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md)
once per bag, with cross-fitting optionally enabled within each bag. All
other arguments (`propensity_method`, `baseline_method`, `cens_type`,
`cens_params`, `n_grid`, `eps`) are passed through unchanged.

## References

Wager, S., Hastie, T. and Efron, B. (2014). Confidence intervals for
random forests: the jackknife and the infinitesimal jackknife. *Journal
of Machine Learning Research* 15: 1625–1651.

Gao, Z. and Hastie, T. (2025). Estimating heterogeneous treatment
effects for general responses. *Biometrics* 81(4): ujaf162.
[doi:10.1093/biomtc/ujaf162](https://doi.org/10.1093/biomtc/ujaf162) .

## See also

[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md)
for the single-pass DINA estimator with sandwich variance.

## Examples

``` r
set.seed(1)
n <- 400
d <- 3
X <- matrix(stats::runif(n * d, -1, 1), n, d)
colnames(X) <- paste0("x", seq_len(d))
W <- stats::rbinom(n, 1, stats::plogis(0.5 * X[, 1]))
tau_x <- 0.5 + 0.8 * X[, 1] - 0.3 * X[, 2]
Y <- 1 + X[, 1] + W * tau_x + stats::rnorm(n)
fit <- dina_fit_bagged(X, Y, W, family = "gaussian",
                       n_bags = 30L, seed = 1L)
#> Warning: IJ variance had 1 negative eigenvalue(s); clipped to zero. Consider increasing `n_bags`.
coef(fit)
#> (Intercept)          x1          x2          x3 
#>   0.4956289   0.9512341  -0.5060264  -0.0614721 
sqrt(diag(vcov(fit)))
#> (Intercept)          x1          x2          x3 
#>  0.13755332  0.19385712  0.09391081  0.33511258 
```
