# Compute Detection Probability for Single Theta (Internal)

Compute Detection Probability for Single Theta (Internal)

## Usage

``` r
compute_detection_probability_single(
  theta,
  n_sg,
  prop_cens,
  k_avg,
  k_ind,
  method,
  n_mc,
  tol
)
```

## Arguments

- theta:

  Numeric. True hazard ratio in the subgroup. Can be a vector for
  computing detection probability across multiple HR values.

- n_sg:

  Integer. Subgroup sample size.

- prop_cens:

  Numeric. Proportion censored (0-1). Default: 0.3

- k_avg:

  Log of hr_threshold

- k_ind:

  Log of hr_consistency

- method:

  Character. Integration method: "cubature" (recommended for accuracy)
  or "monte_carlo" (faster for exploration). Default: "cubature"

- n_mc:

  Integer. Number of Monte Carlo samples if method = "monte_carlo".
  Default: 100000

- tol:

  Numeric. Relative tolerance for cubature integration. Default: 1e-4

## Value

Numeric probability
