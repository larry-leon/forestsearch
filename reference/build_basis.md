# Build Polynomial Basis Matrix

Build Polynomial Basis Matrix

## Usage

``` r
build_basis(X, poly_order = 1L, include_intercept = TRUE)
```

## Arguments

- X:

  Numeric matrix (N x d). Must be numeric (no factors).

- poly_order:

  Integer. 1 = linear, 2 = squares + interactions.

- include_intercept:

  Logical. Default: TRUE.

## Value

Numeric matrix (N x K).
