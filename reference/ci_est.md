# Confidence Interval for Estimate

Calculates confidence interval for an estimate, optionally on log(HR)
scale.

## Usage

``` r
ci_est(x, sd, alpha = 0.025, scale = "hr", est.loghr = TRUE)
```

## Arguments

- x:

  Numeric estimate.

- sd:

  Numeric standard deviation.

- alpha:

  Numeric significance level (default: 0.025).

- scale:

  Character. "hr" or "1/hr".

- est.loghr:

  Logical. Is estimate on log(HR) scale?

## Value

List with length, lower, upper, sd, and estimate.
