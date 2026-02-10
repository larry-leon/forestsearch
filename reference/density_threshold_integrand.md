# Vectorized Density for Integration

Wrapper around density_threshold_both for use with cubature integration.

## Usage

``` r
density_threshold_integrand(x, theta, prop_cens, n_sg, k_avg, k_ind)
```

## Arguments

- x:

  Numeric vector of length 2. Log hazard ratio estimates from the two
  split-halves.

- theta:

  Numeric. True hazard ratio in the subgroup.

- prop_cens:

  Numeric. Proportion censored (0-1). Default: 0.3

- n_sg:

  Integer. Subgroup sample size.

- k_avg:

  Numeric. Threshold for average log(HR) across splits. Typically
  log(hr.threshold).

- k_ind:

  Numeric. Threshold for individual split log(HR). Typically
  log(hr.consistency).

## Value

Numeric density value.
