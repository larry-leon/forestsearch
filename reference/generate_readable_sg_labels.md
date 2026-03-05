# Generate Readable Subgroup Labels from ForestSearch Object

Extracts human-readable subgroup labels that are also valid R
expressions for use with
[`plotKM.band_subgroups()`](https://rdrr.io/pkg/weightedsurv/man/plotKM.band_subgroups.html).
Attempts to extract the actual subgroup definition (e.g., "er \<= 0")
rather than column references.

## Usage

``` r
generate_readable_sg_labels(fs.est, verbose = FALSE)
```

## Arguments

- fs.est:

  A forestsearch object.

- verbose:

  Logical. Print diagnostic messages.

## Value

Character vector of length 2: c(harm_label, benefit_label)
