# Fit policy trees up to specified depth

Fits policy trees of depths 1 through maxdepth and computes metrics

## Usage

``` r
fit_policy_trees(X, data, dr.scores, maxdepth, n.min)
```

## Arguments

- X:

  Matrix. Covariate matrix

- data:

  Data frame. Original data

- dr.scores:

  Matrix. Doubly robust scores

- maxdepth:

  Integer. Maximum tree depth (1-3)

- n.min:

  Integer. Minimum subgroup size

## Value

List with trees and combined values

## Examples

``` r
if (FALSE) { # \dontrun{
# fit_policy_trees() is called internally by grf.subg.harm.survival().
# See grf.subg.harm.survival() for the standard entry point.
} # }
```
