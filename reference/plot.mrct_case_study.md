# Plot an MRCT case study

The forest plot of the four sets' hazard ratios with Wald 95\\ (log
scale) and the Kaplan-Meier curves by arm in each set (one panel per
set, 2 x 2), the latter with weightedsurv.

## Usage

``` r
# S3 method for class 'mrct_case_study'
plot(x, which = c("forest", "km"), xlim = c(0.3, 2), time_unit = NULL, ...)
```

## Arguments

- x:

  An `mrct_case_study` object.

- which:

  Any of `"forest"`, `"km"`. Default both.

- xlim:

  Hazard-ratio range of the forest plot. Default `c(0.3, 2)`; intervals
  are clipped to it.

- time_unit:

  Character. The unit of the Kaplan-Meier time axis, used as its label
  (e.g. `"months"`); `NULL` (default) labels it `"time"`.

- ...:

  Unused.

## Value

`x`, invisibly.
