# Extract Estimates from ForestSearch Results

Extracts operating characteristics (HRs, classification metrics, etc.)
from ForestSearch analysis results. Aligned with new DGM output
structure.

## Usage

``` r
extract_fs_estimates(
  df,
  fs_res,
  dgm,
  cox_formula = NULL,
  cox_formula_adj = NULL,
  analysis = "FS",
  fs_full = NULL,
  verbose = FALSE
)
```

## Arguments

- df:

  Simulated data frame

- fs_res:

  ForestSearch result table (grp.consistency\$out_sg\$result, or NULL)

- dgm:

  DGM object containing true HRs (supports both old and new formats)

- cox_formula:

  Cox formula for estimation (optional)

- cox_formula_adj:

  Adjusted Cox formula (optional)

- analysis:

  Analysis label (e.g., "FS", "FSlg")

- fs_full:

  Full forestsearch result object (for df.est access)

- verbose:

  Logical. Print extraction details. Default: FALSE

## Value

data.table with extracted estimates including AHR metrics
