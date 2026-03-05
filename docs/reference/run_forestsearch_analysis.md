# Run ForestSearch Analysis

Helper function to run ForestSearch and extract estimates. Aligned with
forestsearch() parameters including use_twostage.

## Usage

``` r
run_forestsearch_analysis(
  data,
  confounders_name,
  params,
  dgm,
  cox_formula = NULL,
  cox_formula_adj = NULL,
  analysis_label = "FS",
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame with simulated trial data

- confounders_name:

  Character vector of confounder names

- params:

  List of ForestSearch parameters

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

## Value

data.table with analysis estimates
