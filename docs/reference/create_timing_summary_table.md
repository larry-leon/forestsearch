# Create Timing Summary Table

Creates a data frame summarizing bootstrap timing information.

## Usage

``` r
create_timing_summary_table(
  overall_timing,
  iteration_stats,
  fs_stats,
  overhead_stats,
  nb_boots,
  boot_success_rate
)
```

## Arguments

- overall_timing:

  List. Overall timing statistics

- iteration_stats:

  List. Per-iteration timing statistics

- fs_stats:

  List. ForestSearch-specific timing statistics

- overhead_stats:

  List. Overhead timing statistics

- nb_boots:

  Integer. Number of bootstrap iterations

- boot_success_rate:

  Numeric. Proportion of successful bootstraps

## Value

Data frame with timing summary
