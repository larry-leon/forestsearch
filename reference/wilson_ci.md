# Wilson Score Confidence Interval

Computes Wilson score confidence interval for a proportion, which has
better coverage properties than the normal approximation for small
samples and proportions near 0 or 1.

## Usage

``` r
wilson_ci(x, n, conf.level = 0.95)
```

## Arguments

- x:

  Integer. Number of successes.

- n:

  Integer. Number of trials.

- conf.level:

  Numeric. Confidence level (default 0.95).

## Value

Named numeric vector with elements: estimate, lower, upper.
