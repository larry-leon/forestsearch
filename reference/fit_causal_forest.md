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

## Examples

``` r
# \donttest{
library(survival)
library(grf)
df <- survival::gbsg
X <- as.matrix(df[, c("age", "meno", "size", "nodes", "pgr", "er")])
tau <- quantile(df$rfstime[df$status == 1], 0.6)
cs <- fit_causal_forest(X = X, Y = df$rfstime, W = df$hormon,
                        D = df$status, tau.rmst = tau,
                        RCT = FALSE, seedit = 42L)
# }
```
