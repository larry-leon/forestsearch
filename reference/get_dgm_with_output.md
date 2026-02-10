# Create DGM with Output File Path

Wrapper function that creates a GBSG DGM and generates a standardized
output file path for saving results.

## Usage

``` r
get_dgm_with_output(
  model_harm,
  n,
  k_treat = 1,
  target_hr_harm = NULL,
  cens_type = "weibull",
  out_dir = NULL,
  file_prefix = "sim",
  file_suffix = "",
  include_hr_in_name = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- model_harm:

  Character. Model type ("alt" or "null")

- n:

  Integer. Planned sample size (for filename)

- k_treat:

  Numeric. Treatment effect multiplier

- target_hr_harm:

  Numeric. Target HR for harm subgroup (used for calibration when model
  = "alt")

- cens_type:

  Character. Censoring type

- out_dir:

  Character. Output directory path. If NULL, no file path is generated

- file_prefix:

  Character. Prefix for output filename

- file_suffix:

  Character. Suffix for output filename

- include_hr_in_name:

  Logical. Include achieved HR in filename. Default: FALSE

- verbose:

  Logical. Print diagnostic information. Default: FALSE

- ...:

  Additional arguments passed to
  [`create_gbsg_dgm`](https://github.com/larry-leon/forestsearch/reference/create_gbsg_dgm.md)

## Value

List with components:

- dgm:

  The gbsg_dgm object

- out_file:

  Character path to output file (NULL if out_dir is NULL)

- k_inter:

  The k_inter value used (either calibrated or default)
