# Estimate heterogeneous treatment effects via DINA

`dina_fit()` implements the two-step DINA (DIfference in NAtural
parameters) estimator of Gao and Hastie (2025) for heterogeneous
treatment effects (HTE) with general response types: Gaussian (mean
difference), binomial (log odds ratio), Poisson (log mean ratio), and
Cox (log hazard ratio). The treatment effect is assumed to be linear in
the covariates, \$\$\tau(x) = \beta_0 + x^\top \beta\_{1:d},\$\$
following Proposition 1 of the paper. The procedure has two steps: (1)
estimate nuisance functions \\a(x)\\ and \\\nu(x)\\ on a training fold,
and (2) plug those estimators into a GLM (or Cox) likelihood on the
held-out fold to estimate \\\beta\\. Cross-fitting averages the two
estimators obtained by swapping the folds.

## Usage

``` r
dina_fit(
  X,
  Y,
  W,
  family = c("gaussian", "binomial", "poisson", "cox"),
  propensity_method = "logistic",
  baseline_method = NULL,
  cross_fitting = TRUE,
  n_folds = 2L,
  cens_type = NULL,
  cens_params = list(),
  n_grid = 500L,
  eps = 1e-06,
  seed = NULL
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

- cross_fitting:

  logical; if `TRUE` (the default) use `n_folds`-fold cross-fitting to
  separate nuisance estimation from target-parameter estimation. If
  `FALSE`, the same data are used for both steps.

- n_folds:

  positive integer giving the number of cross-fitting folds (default
  `2`, matching the paper).

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

- eps:

  small positive constant added to denominators for numerical stability
  (default `1e-6`).

- seed:

  optional integer seed for the random cross-fitting fold assignment.
  When supplied, the global RNG state is restored on exit.

## Value

An object of class `"dina"`, a list with components:

- coefficients:

  named numeric vector of estimated DINA coefficients \\\hat\beta\\,
  with `(Intercept)` first followed by one entry per covariate.

- vcov:

  estimated variance-covariance matrix of \\\hat\beta\\.

- family:

  the response family, as supplied.

- n, d:

  sample size and number of covariates.

- cross_fitting, n_folds:

  cross-fitting settings.

- call:

  the matched call.

## Details

The Gaussian family uses the simplification of Section 3.2: \\a(x) =
e(x)\\ and \\\nu(x) = E\[Y \mid X = x\]\\, so the conditional outcome
means \\\eta_0(x)\\, \\\eta_1(x)\\ are not required separately. For the
binomial and Poisson families, separate within-arm fits to \\\eta_0(x)\\
and \\\eta_1(x)\\ are used to construct \\a(x)\\ via the modified
propensity in Eq. (12) of the paper. For the Cox family, \\a(x)\\ uses
the not-censoring probabilities of Eq. (18); how those are estimated is
controlled by `cens_type`.

Variance estimation follows the standard sandwich form:
\\\widehat{V}(\hat\beta) = \widehat\Sigma_1^{-1} \widehat\Sigma_2
\widehat\Sigma_1^{-1} / n\\, where \\\widehat\Sigma_1\\ and
\\\widehat\Sigma_2\\ are the bread (negative score derivative) and meat
(outer product of scores) matrices. Under cross-fitting these are
averaged across folds before the sandwich is formed. For the Cox step 2
we use the partial-likelihood information-matrix variance returned by
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

## References

Gao, Z. and Hastie, T. (2025). Estimating heterogeneous treatment
effects for general responses. *Biometrics* 81(4): ujaf162.
[doi:10.1093/biomtc/ujaf162](https://doi.org/10.1093/biomtc/ujaf162) .

Robinson, P. M. (1988). Root-n-consistent semiparametric regression.
*Econometrica* 56(4): 931–954.

Nie, X. and Wager, S. (2021). Quasi-oracle estimation of heterogeneous
treatment effects. *Biometrika* 108(2): 299–319.

Dandl, S., Bender, A. and Hothorn, T. (2024). Heterogeneous treatment
effect estimation for observational data using model-based forests.
*Statistical Methods in Medical Research* 33(3): 392–413.

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
fit <- dina_fit(X, Y, W, family = "gaussian", seed = 1)
coef(fit)
#> (Intercept)          x1          x2          x3 
#>  0.52749233  0.90471230 -0.50526887 -0.06332475 
confint(fit)
#>                  2.5 %     97.5 %
#> (Intercept)  0.3245447  0.7304399
#> x1           0.5815606  1.2278640
#> x2          -0.8512792 -0.1592586
#> x3          -0.4234146  0.2967651
```
