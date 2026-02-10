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
  [`bootstrap_results`](https://github.com/larry-leon/forestsearch/reference/bootstrap_results.md)
  to ensure bootstrap index alignment with the Ystar matrix.

## Value

Matrix of bootstrap samples (nb_boots x nrow(df)).
