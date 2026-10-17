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

## Examples

``` r
wilson_ci(90, 100)
#>  estimate     lower     upper 
#> 0.9000000 0.8256343 0.9447709 
wilson_ci(5, 20, conf.level = 0.90)
#>  estimate     lower     upper 
#> 0.2500000 0.1273772 0.4322018 
```
