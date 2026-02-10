# Extract Estimates from GRF Results

Aligned with new DGM output structure including AHR metrics. Correctly
handles grf.subg.harm.survival() output structure:

- sg.harm.id: subgroup definition string

- data: data frame with treat.recommend column (0 = harm, 1 =
  complement)

## Usage

``` r
extract_grf_estimates(
  df,
  grf_est,
  dgm,
  cox_formula = NULL,
  cox_formula_adj = NULL,
  analysis = "GRF",
  frac_tau = 1,
  verbose = FALSE,
  debug = FALSE
)
```

## Arguments

- df:

  Simulated data frame

- grf_est:

  GRF estimation result from grf.subg.harm.survival()

- dgm:

  DGM object

- cox_formula:

  Cox formula

- cox_formula_adj:

  Adjusted Cox formula

- analysis:

  Analysis label

- frac_tau:

  Fraction of tau used

- verbose:

  Print extraction details

- debug:

  Print detailed debugging information about GRF result structure

## Value

data.table with extracted estimates
