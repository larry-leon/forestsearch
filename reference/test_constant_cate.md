# Test for Constant Conditional Average Treatment Effect

H0': tau(x) = tau for some tau and all x.

## Usage

``` r
test_constant_cate(
  Y,
  W,
  X,
  poly_order = 1L,
  covariate_select = c("all", "top_down", "bottom_up"),
  t_threshold = 2,
  regression = c("ols", "logistic", "poisson"),
  offset = NULL
)
```

## Arguments

- Y:

  Numeric outcome vector.

- W:

  Treatment indicator (0/1).

- X:

  Numeric covariate matrix (no factors; use `as.numeric(factor) - 1` to
  convert).

- poly_order:

  Integer. Default: 1.

- covariate_select:

  "all", "top_down", or "bottom_up".

- t_threshold:

  Numeric. Default: 2.0.

- regression:

  "ols", "logistic", or "poisson".

- offset:

  Numeric offset vector (Poisson only).

## Value

A list with components `test`, `chi_sq`, `df`, `p_value_chi`, `normal`
(normal-approximation statistic), `p_value_normal`, `K` (number of
covariates), `N0`, `N1`, `regression`, `covariates_selected`, `ate_diff`
(intercept-difference point estimate), and `diff_slope`, `V_slope`.
Returns `NULL` if per-arm fits fail or fewer than two covariates are
retained.
