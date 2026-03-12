# Count ID Occurrences in Bootstrap Sample

Counts the number of times an ID appears in a bootstrap sample.

## Usage

``` r
count_boot_id(x, dfb)
```

## Arguments

- x:

  ID value.

- dfb:

  Data frame of bootstrap sample.

## Value

Integer count of occurrences.

## Examples

``` r
df_boot <- data.frame(id = c(1, 2, 1, 3, 1), id_boot = 1:5)
count_boot_id(1, df_boot)  # returns 3
#> [1] 3
count_boot_id(4, df_boot)  # returns 0
#> [1] 0
```
