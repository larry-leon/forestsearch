# Bootstrap Ystar Matrix

Generates a bootstrap matrix for Ystar using parallel processing.

## Usage

``` r
bootstrap_ystar(df, nb_boots, seed = 8316951L)
```

## Arguments

- df:

  Data frame.

- nb_boots:

  Integer. Number of bootstrap samples.

- seed:

  Integer. Random seed for reproducibility. Default 8316951L. Must match
  the seed used in
  [`bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/bootstrap_results.md)
  to ensure bootstrap index alignment with the Ystar matrix.

## Value

Matrix of bootstrap samples (nb_boots x nrow(df)).

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- data.frame(id = 1:50, tte = rexp(50), event = rbinom(50, 1, 0.6),
                 treat = rep(0:1, 25))
future::plan(future::sequential)
ystar <- bootstrap_ystar(df, nb_boots = 10)
dim(ystar)
future::plan(future::sequential)
} # }
```
