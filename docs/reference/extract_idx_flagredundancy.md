# Extract redundancy flag for subgroup combinations

Checks if adding each factor to a subgroup reduces the sample size by at
least `rmin`.

## Usage

``` r
extract_idx_flagredundancy(x, rmin)
```

## Arguments

- x:

  Matrix of subgroup factor indicators.

- rmin:

  Integer. Minimum required reduction in sample size.

## Value

List with `id.x` (membership vector) and `flag.redundant` (logical).
