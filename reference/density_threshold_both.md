# Bivariate Density for Split-Sample HR Threshold Detection

Computes the joint density for the two-split detection criterion where
both split-halves must exceed individual thresholds and their average
must exceed a consistency threshold.

## Usage

``` r
density_threshold_both(x, theta, prop_cens = 0.3, n_sg, k_avg, k_ind)
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

Numeric. Joint density value at x, or 0 if thresholds not met.

## Details

The detection criterion requires:

- Average of two splits: (x1 + x2)/2 \>= k_avg

- Individual splits: x1 \>= k_ind AND x2 \>= k_ind

Under the Leon et al. (2024, Section 2.1) approximation, each 50/50
split-half log(HR) estimator follows N(log(theta), 8/d) where d = n_sg
\* (1 - prop_cens) is the expected TOTAL number of events in the
subgroup. Each split holds d/2 events, so its variance is 4/(d/2) = 8/d.
