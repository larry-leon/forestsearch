# Find Minimum Sample Size for Target Detection Power

Determines the minimum subgroup sample size needed to achieve a target
detection probability for a given true hazard ratio.

## Usage

``` r
find_required_sample_size(
  theta,
  target_power = 0.8,
  prop_cens = 0.3,
  hr_threshold = 1.25,
  hr_consistency = 1,
  n_range = c(20L, 500L),
  tol = 1,
  verbose = TRUE
)
```

## Arguments

- theta:

  Numeric. True hazard ratio in subgroup.

- target_power:

  Numeric. Target detection probability (0-1). Default: 0.80

- prop_cens:

  Numeric. Proportion censored. Default: 0.3

- hr_threshold:

  Numeric. HR threshold. Default: 1.25

- hr_consistency:

  Numeric. HR consistency threshold. Default: 1.0

- n_range:

  Integer vector of length 2. Range of sample sizes to search. Default:
  c(20, 500)

- tol:

  Numeric. Tolerance for bisection search. Default: 1

- verbose:

  Logical. Print progress. Default: TRUE

## Value

A list with:

- n_sg_required:

  Minimum sample size (rounded up)

- achieved_power:

  Actual detection probability at n_sg_required

- theta:

  Input hazard ratio

- target_power:

  Input target power
