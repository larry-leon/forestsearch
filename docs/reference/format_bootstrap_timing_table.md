# Format Bootstrap Timing Table with gt

Creates a publication-ready timing summary table from bootstrap results.

## Usage

``` r
format_bootstrap_timing_table(timing_list, nb_boots, boot_success_rate)
```

## Arguments

- timing_list:

  List. Timing information from summarize_bootstrap_results()\$timing

- nb_boots:

  Integer. Number of bootstrap iterations

- boot_success_rate:

  Numeric. Proportion of successful bootstraps

## Value

A gt table object
