# Filter a vector by LASSO-selected variables

Returns elements of `x` that are in `lassokeep`.

## Usage

``` r
filter_by_lassokeep(x, lassokeep)
```

## Arguments

- x:

  Character vector.

- lassokeep:

  Character vector of selected variables.

## Value

Filtered character vector or NULL.

## Examples

``` r
filter_by_lassokeep(c("age", "size", "nodes"), lassokeep = c("age", "nodes"))
#> [1] "age"   "nodes"
filter_by_lassokeep(c("pgr", "er"), lassokeep = c("age", "nodes"))  # returns NULL
#> NULL
```
