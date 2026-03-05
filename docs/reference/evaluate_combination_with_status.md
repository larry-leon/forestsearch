# Evaluate a Single Factor Combination with Status Tracking

Tests whether a specific combination meets all criteria and returns a
status code indicating how far the evaluation progressed.

## Usage

``` r
evaluate_combination_with_status(
  covs.in,
  yy,
  dd,
  tt,
  zz,
  n.min,
  d0.min,
  d1.min,
  hr.threshold,
  minp,
  rmin,
  kk
)
```

## Arguments

- covs.in:

  Numeric vector. Factor selection indicators.

- yy:

  Numeric vector. Outcome values.

- dd:

  Numeric vector. Event indicators.

- tt:

  Numeric vector. Treatment indicators.

- zz:

  Matrix. Factor indicators.

- n.min:

  Integer. Minimum sample size.

- d0.min:

  Integer. Minimum control events.

- d1.min:

  Integer. Minimum treatment events.

- hr.threshold:

  Numeric. HR threshold.

- minp:

  Numeric. Minimum prevalence.

- rmin:

  Integer. Minimum size reduction.

- kk:

  Integer. Combination index.

## Value

List with:

- status:

  Integer status code: 0 = failed variance check, 1 = passed variance,
  failed prevalence, 2 = passed prevalence, failed redundancy, 3 =
  passed redundancy, failed events, 4 = passed events, failed sample
  size, 5 = passed sample size, failed Cox fit, 6 = passed Cox fit,
  failed HR threshold, 7 = passed all criteria (success)

- result:

  Result row if successful, NULL otherwise
