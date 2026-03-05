# Remove Near-Duplicate Subgroups

Removes subgroups with nearly identical statistics (HR, n, E, etc.) to
reduce redundancy in candidate list.

## Usage

``` r
remove_near_duplicate_subgroups(
  hr_subgroups,
  tolerance = 0.001,
  details = FALSE
)
```

## Arguments

- hr_subgroups:

  Data.table of subgroup results.

- tolerance:

  Numeric. Tolerance for numeric comparison (default 0.001).

- details:

  Logical. Print details about removed duplicates.

## Value

Data.table with near-duplicate rows removed.
