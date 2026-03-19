# Find Quantile for Target Subgroup Proportion

Determines the quantile cutpoint that achieves a target proportion of
observations in a subgroup. Useful for calibrating subgroup sizes.

## Usage

``` r
find_quantile_for_proportion(
  data,
  var_name,
  target_prop,
  direction = "less",
  tol = 1e-04
)
```

## Arguments

- data:

  A data.frame containing the variable of interest

- var_name:

  Character string specifying the variable name to analyze

- target_prop:

  Numeric value between 0 and 1 specifying the target proportion of
  observations to be included in the subgroup

- direction:

  Character string: "less" for values \<= cutpoint (default), "greater"
  for values \> cutpoint

- tol:

  Numeric tolerance for root finding algorithm. Default is 0.0001

## Value

A list containing:

- quantile:

  The quantile value (between 0 and 1) that achieves the target
  proportion

- cutpoint:

  The actual data value corresponding to this quantile

- actual_proportion:

  The achieved proportion (should equal target_prop within tolerance)

## Details

This function uses root finding (`uniroot`) to determine the quantile
that results in exactly the target proportion of observations being
classified into the subgroup. This is particularly useful when you want
to ensure a specific subgroup size regardless of the data distribution.

## See also

[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
