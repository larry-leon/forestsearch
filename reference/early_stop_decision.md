# Early Stopping Decision

Evaluates whether enough evidence exists to stop early based on
confidence interval for consistency proportion.

## Usage

``` r
early_stop_decision(
  n_success,
  n_total,
  threshold,
  conf.level = 0.95,
  min_samples = 20
)
```

## Arguments

- n_success:

  Integer. Number of splits meeting consistency.

- n_total:

  Integer. Total number of valid splits.

- threshold:

  Numeric. Target consistency threshold.

- conf.level:

  Numeric. Confidence level for decision (default 0.95).

- min_samples:

  Integer. Minimum samples before allowing early stop.

## Value

Character. One of "continue", "pass", or "fail".

## Examples

``` r
early_stop_decision(95, 100, threshold = 0.90)
#> [1] "continue"
early_stop_decision(60, 100, threshold = 0.90)
#> [1] "fail"
early_stop_decision(10, 15, threshold = 0.90)  # below min_samples
#> [1] "continue"
```
