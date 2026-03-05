# Create Subgroup Indicator Columns from ForestSearch

Internal helper to create Qrecommend and Brecommend indicator columns.

## Usage

``` r
create_fs_subgroup_indicators(
  df,
  fs.est,
  col_names = c("Qrecommend", "Brecommend"),
  verbose = FALSE
)
```

## Arguments

- df:

  Data frame to modify.

- fs.est:

  A forestsearch object.

- col_names:

  Character vector of length 2. Names for the indicator columns: first
  for harm/questionable (treat.recommend == 0), second for
  benefit/recommend (treat.recommend == 1). Default: c("Qrecommend",
  "Brecommend")

- verbose:

  Logical. Print diagnostic messages.

## Value

Modified data frame with indicator columns.
