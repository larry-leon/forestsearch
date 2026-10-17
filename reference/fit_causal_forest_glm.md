# Fit Causal Forest for GLM Outcomes

Fits a GRF causal forest for binary or continuous outcomes using
[`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html).
This is the GLM counterpart to
[`fit_causal_forest`](https://larry-leon.github.io/forestsearch/reference/fit_causal_forest.md),
which uses `causal_survival_forest`.

## Usage

``` r
fit_causal_forest_glm(X, Y, W, RCT, seedit, tune_grf = FALSE)
```

## Arguments

- X:

  Numeric matrix. Covariates.

- Y:

  Numeric vector. Outcome (binary 0/1 or continuous).

- W:

  Numeric vector. Treatment indicator (0/1).

- RCT:

  Logical. Is this RCT data?

- seedit:

  Integer. Random seed.

- tune_grf:

  Logical. If TRUE, enables cross-validated hyperparameter tuning via
  `tune.parameters = "all"`. Default FALSE preserves existing behavior.

## Value

Causal forest object.
