# Load L_eff Calibration Results

Reads a calibration bundle saved by
[`save_leff_calibration()`](https://larry-leon.github.io/forestsearch/reference/save_leff_calibration.md).
Prints a summary and returns the full bundle.

## Usage

``` r
load_leff_calibration(outcome_type, output_dir = NULL, verbose = TRUE)
```

## Arguments

- outcome_type:

  Character. One of "binary", "survival", "count", "continuous".

- output_dir:

  Character. Directory containing .rds files. Default: `"_data"`.

- verbose:

  Logical. Print summary on load. Default: TRUE.

## Value

A list with components C, alpha, n_min, outcome_type, cal_data,
dgm_description, timestamp, etc.

## Examples

``` r
if (FALSE) { # \dontrun{
cal <- load_leff_calibration("binary")
C_prior <- cal$C
alpha_prior <- cal$alpha
} # }
```
