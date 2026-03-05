# Run GRF Analysis

Helper function to run standalone GRF analysis using
grf.subg.harm.survival().

## Usage

``` r
run_grf_analysis(
  data,
  confounders_name,
  params,
  dgm,
  cox_formula = NULL,
  cox_formula_adj = NULL,
  analysis_label = "GRF",
  verbose = FALSE,
  debug = FALSE
)
```

## Arguments

- data:

  Data frame with simulated trial data

- confounders_name:

  Character vector of confounder names

- params:

  List of GRF parameters (from grf_merged)

- dgm:

  DGM object for computing true HRs

- cox_formula:

  Cox formula for estimation

- cox_formula_adj:

  Adjusted Cox formula

- analysis_label:

  Character label for this analysis

- verbose:

  Print details

- debug:

  Print detailed debugging information

## Value

data.table with analysis estimates
