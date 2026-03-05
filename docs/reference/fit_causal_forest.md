# Fit causal survival forest

Wrapper function to fit GRF causal survival forest with appropriate
settings

## Usage

``` r
fit_causal_forest(X, Y, W, D, tau.rmst, RCT, seedit)
```

## Arguments

- X:

  Matrix. Covariate matrix

- Y:

  Numeric vector. Outcome variable

- W:

  Numeric vector. Treatment indicator

- D:

  Numeric vector. Event indicator

- tau.rmst:

  Numeric. Time horizon for RMST

- RCT:

  Logical. Is this RCT data?

- seedit:

  Integer. Random seed

## Value

Causal survival forest object
