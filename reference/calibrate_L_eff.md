# Calibrate L_eff from Simulation Results

Given per-subgroup approximation probabilities and simulated FPR values
at multiple sample sizes, fits the power-law model \\L\_{\text{eff}}(N)
= C \cdot (N / n\_{\min})^\alpha\\.

## Usage

``` r
calibrate_L_eff(N, P1, sim_fpr, n_min = 60)
```

## Arguments

- N:

  Integer vector. Sample sizes at which simulations were run.

- P1:

  Numeric vector. Per-subgroup detection probability from
  [`compute_detection_probability_glm`](https://larry-leon.github.io/forestsearch/reference/compute_detection_probability_glm.md)
  at each N.

- sim_fpr:

  Numeric vector. Simulated procedure-level FPR at each N.

- n_min:

  Integer. Minimum subgroup size used in ForestSearch. Default: 60.

## Value

A list of class `"leff_calibration"` with components:

- C:

  Numeric. Fitted constant.

- alpha:

  Numeric. Fitted power-law exponent.

- n_min:

  Integer. Minimum subgroup size.

- r_squared:

  Numeric. R-squared of the log-linear fit.

- data:

  Data frame with N, P1, sim_fpr, L_eff, and fitted values.

- fit:

  The `lm` object.

## Details

At each sample size, the implied \\L\_{\text{eff}}\\ is:
\$\$L\_{\text{eff}} = \frac{\log(1 - \text{FPR}\_{\text{sim}})}
{\log(1 - P_1)}\$\$

A log-linear model is then fit: \\\log(L\_{\text{eff}}) = \log(C) +
\alpha \log(N / n\_{\min})\\.

## Examples

``` r
# \donttest{
# From binary threshold calibration results
cal <- calibrate_L_eff(
  N       = c(200, 500, 700, 1000),
  P1      = c(0.152, 0.109, 0.091, 0.071),
  sim_fpr = c(0.185, 0.220, 0.270, 0.280),
  n_min   = 60
)
print(cal)
#> L_eff Calibration
#>   Model: L_eff = 0.454 * (N / 60)^0.795
#>   R-squared: 0.976
#>   Calibrated at N = 200, 500, 700, 1000
#> 
#>   Implied L_eff:
#>     N =  200: L_eff = 1.24 (fitted: 1.18)
#>     N =  500: L_eff = 2.15 (fitted: 2.45)
#>     N =  700: L_eff = 3.30 (fitted: 3.20)
#>     N = 1000: L_eff = 4.46 (fitted: 4.25)
# }
```
