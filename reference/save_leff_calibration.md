# Save L_eff Calibration Results

Writes a standardized `.rds` bundle containing the fitted L_eff
parameters and all supporting data. Called at the end of calibration
documents (e.g., `calibration_binary_leff.qmd`).

## Usage

``` r
save_leff_calibration(
  C,
  alpha,
  n_min,
  outcome_type,
  cal_data = NULL,
  dgm_description = "",
  n_sims_per_N = NA_integer_,
  output_dir = NULL,
  extra = list()
)
```

## Arguments

- C:

  Numeric. Fitted scale parameter.

- alpha:

  Numeric. Fitted power-law exponent.

- n_min:

  Integer. Reference minimum subgroup size.

- outcome_type:

  Character. One of "binary", "survival", "count", "continuous".

- cal_data:

  Data frame with columns N, sim_fpr, P1, L_eff (the raw calibration
  points).

- dgm_description:

  Character. Free-text description of the DGM used for calibration.

- n_sims_per_N:

  Integer. Simulations per sample size.

- output_dir:

  Character. Directory for .rds files. Required: an error is raised when
  it is `NULL` (the default). The former default was `"_data"`.

- extra:

  List. Any additional items to store.

## Value

Invisible path to the saved file.

## Examples

``` r
if (FALSE) { # \dontrun{
save_leff_calibration(
  C = 0.220, alpha = 1.298, n_min = 60,
  outcome_type = "binary",
  cal_data = cal_df,
  dgm_description = "Binary DGM, 4 confounders (bm1, bm2, age, ecog)",
  n_sims_per_N = 5000,
  output_dir = tempdir()
)
} # }
```
