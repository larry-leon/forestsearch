# Summarize Factor Presence Across Bootstrap Subgroups

Analyzes how often each individual factor appears in identified
subgroups, extracting base factor names from full definitions and
identifying common specific definitions.

## Usage

``` r
summarize_factor_presence_robust(
  results,
  maxk = 2,
  threshold = 10,
  as_gt = TRUE
)
```

## Arguments

- results:

  Data.table or data.frame. Bootstrap results with M.1, M.2, etc.
  columns

- maxk:

  Integer. Maximum number of factors allowed

- threshold:

  Numeric. Percentage threshold for including specific definitions
  (default: 10)

- as_gt:

  Logical. Return gt tables (TRUE) or data.frames (FALSE)

## Value

List with base_factors and specific_factors data.frames or gt tables
