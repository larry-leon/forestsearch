# Predict Corrected Null FPR

Uses a calibrated \\L\_{\text{eff}}\\ model to predict the
procedure-level FPR at arbitrary sample sizes and thresholds.

## Usage

``` r
predict_fpr_corrected(
  calibration,
  theta,
  d_eff,
  c1,
  c2,
  N,
  effect_scale = "ratio"
)
```

## Arguments

- calibration:

  An `"leff_calibration"` object from
  [`calibrate_L_eff`](https://larry-leon.github.io/forestsearch/reference/calibrate_L_eff.md).

- theta:

  Numeric. True treatment effect (1.0 for null).

- d_eff:

  Numeric. Effective information in the subgroup.

- c1:

  Numeric. Screening threshold.

- c2:

  Numeric. Consistency threshold.

- N:

  Integer. Total sample size.

- effect_scale:

  Character. `"ratio"` or `"difference"`.

## Value

A list with components:

- P1:

  Per-subgroup detection probability.

- L_eff:

  Effective number of candidates at this N.

- fpr_corrected:

  Corrected procedure-level FPR.

## Examples

``` r
# \donttest{
cal <- calibrate_L_eff(
  N = c(200, 500, 700, 1000),
  P1 = c(0.152, 0.109, 0.091, 0.071),
  sim_fpr = c(0.185, 0.220, 0.270, 0.280)
)

# Predict FPR at N = 800
d <- d_eff_binary(n_sg = round(0.30 * 800), p_event = 0.30)
predict_fpr_corrected(cal, theta = 1.0, d_eff = d,
                      c1 = 1.25, c2 = 1.25, N = 800)
#> $P1
#> [1] 0.083182
#> 
#> $L_eff
#> [1] 3.556482
#> 
#> $fpr_corrected
#> [1] 0.2657218
#> 
# }
```
