# Discretize Continuous Variable into Quantile-Based Categories

Discretize Continuous Variable into Quantile-Based Categories

## Usage

``` r
cut_numeric(x, probs = c(0.25, 0.5, 0.75))
```

## Arguments

- x:

  Numeric vector to discretize

- probs:

  Numeric vector of probabilities for quantile breaks. Default: c(0.25,
  0.5, 0.75) creates quartiles coded as 1, 2, 3, 4

## Value

Integer vector with category codes (1 = lowest, max = highest)
