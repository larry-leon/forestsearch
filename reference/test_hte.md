# Run Both Crump et al. Tests

Run Both Crump et al. Tests

## Usage

``` r
test_hte(
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

Object of class `hte_test`.
