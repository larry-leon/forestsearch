# Format Bootstrap Diagnostics Table with gt

Creates a publication-ready diagnostics table from bootstrap results.

## Usage

``` r
format_bootstrap_diagnostics_table(
  diagnostics,
  nb_boots,
  results,
  H_estimates = NULL,
  Hc_estimates = NULL
)
```

## Arguments

- diagnostics:

  List. Diagnostics information from summarize_bootstrap_results()

- nb_boots:

  Integer. Number of bootstrap iterations

- results:

  Data.table. Bootstrap results with bias-corrected estimates

- H_estimates:

  List. H subgroup estimates

- Hc_estimates:

  List. Hc subgroup estimates

## Value

A gt table object
