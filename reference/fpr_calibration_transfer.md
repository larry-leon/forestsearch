# Transfer Corrected FPR Prediction to a Nearby Threshold Pair

Uses the implied \\L\_{\text{eff}}\\ from a calibrated `fpr_calibration`
object to predict the procedure-level FPR at \\(c_1 + \delta_1, c_2 +
\delta_2)\\, then optionally runs a direct simulation at the new pair to
validate the prediction.

## Usage

``` r
fpr_calibration_transfer(
  cal_result,
  delta = c(0.05, 0.1, 0.2),
  delta_c1 = NULL,
  delta_c2 = NULL,
  nsims_check = 200L,
  ...
)
```

## Arguments

- cal_result:

  An object of class `"fpr_calibration"` produced by
  [`fpr_calibration`](https://larry-leon.github.io/forestsearch/reference/fpr_calibration.md).
  Must contain a finite `L_eff_implied`.

- delta:

  Numeric scalar. Added to both \\c_1\\ and \\c_2\\. Positive values
  tighten both thresholds; negative values relax them. Ignored if
  `delta_c1` or `delta_c2` are supplied.

- delta_c1:

  Numeric or `NULL`. Increment for \\c_1\\ only.

- delta_c2:

  Numeric or `NULL`. Increment for \\c_2\\ only.

- nsims_check:

  Integer. Simulations to run at the new pair for empirical validation.
  Default 200. Set to 0 to skip.

- ...:

  Additional arguments passed to `fpr_calibration` for the validation
  run (e.g. `n_workers`, `quiet`).

## Value

A data frame with one row per delta value containing: `c1_new`,
`c2_new`, `P1_new`, `fpr_pred` (predicted corrected FPR using
transferred `L_eff_implied`), `fpr_sim` (direct simulation, if
`nsims_check > 0`), `fpr_sim_lo`, `fpr_sim_hi` (95\\ (whether `fpr_pred`
falls within the simulation CI).

## Examples

``` r
if (FALSE) { # \dontrun{
# Calibrate at base pair
cal <- fpr_calibration(c1 = 1.70, c2 = 1.70, null_fn = my_dgm,
                        n_sample = 500, fs_params = fs_binary,
                        nsims = 200, n_min = 60,
                        outcome_type_approx = "binary", p_event = 0.30)

# Predict at base + 0.05, +0.10, +0.20 and validate
transfer_results <- fpr_calibration_transfer(
  cal, delta = c(0.05, 0.10, 0.20),
  nsims_check = 200, n_workers = 4)

print(transfer_results)
} # }
```
