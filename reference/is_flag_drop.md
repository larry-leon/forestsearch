# Check if cut expression should be dropped

Determines if a cut expression should be dropped (e.g., variable has
\<=1 unique value).

## Usage

``` r
is_flag_drop(thiscut, confounders.name, df)
```

## Arguments

- thiscut:

  Character string of the cut expression.

- confounders.name:

  Character vector of confounder names.

- df:

  Data frame.

## Value

Logical; TRUE if should be dropped, FALSE otherwise.
