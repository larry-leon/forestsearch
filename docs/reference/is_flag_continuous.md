# Check if cut expression is for a continuous variable (OPTIMIZED)

Determines if a cut expression refers to a continuous variable. This
optimized version avoids redundant lookups by using word boundary
matching instead of partial string matching.

## Usage

``` r
is_flag_continuous(thiscut, confounders.name, df, cont.cutoff)
```

## Arguments

- thiscut:

  Character string of the cut expression.

- confounders.name:

  Character vector of confounder names.

- df:

  Data frame.

- cont.cutoff:

  Integer. Cutoff for continuous.

## Value

Logical; TRUE if continuous, FALSE otherwise.
