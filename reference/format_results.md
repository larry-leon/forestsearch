# Format results for subgroup summary

Formats results for subgroup summary table.

## Usage

``` r
format_results(
  subgroup_name,
  n,
  n_treat,
  d,
  m1,
  m0,
  drmst,
  hr,
  hr_a = NA,
  hr_po = NA,
  return_medians = TRUE
)
```

## Arguments

- subgroup_name:

  Character. Subgroup name.

- n:

  Character. Sample size.

- n_treat:

  Character. Treated count.

- d:

  Character. Event count.

- m1:

  Numeric. Median or RMST for treatment.

- m0:

  Numeric. Median or RMST for control.

- drmst:

  Numeric. RMST difference.

- hr:

  Character. Hazard ratio (formatted).

- hr_a:

  Character. Adjusted hazard ratio (optional).

- hr_po:

  Numeric. Potential outcome hazard ratio (optional).

- return_medians:

  Logical. Use medians or RMST.

## Value

Character vector of results.
