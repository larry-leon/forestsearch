# Region-consistency metrics on a simulated trial (ITT and optionally a subgroup)

Region-consistency metrics on a simulated trial (ITT and optionally a
subgroup)

## Usage

``` r
sim_region_metrics(sim, region_col, subgroup_index = NULL, pi_preserve = 0.5)
```

## Arguments

- sim:

  Data frame of one simulated trial with columns `y_sim`, `event_sim`
  and `treat_sim` (as from
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md))
  plus the region indicator.

- region_col:

  Character name of the 0/1 region indicator column in `sim`.

- subgroup_index:

  Optional logical vector of length `nrow(sim)` marking subgroup
  membership. When supplied, Cox HRs are also computed within the
  subgroup for the region and non-region strata. Default `NULL`.

- pi_preserve:

  Numeric threshold for the Method 1 consistency check: the check passes
  when log(HR region) / log(HR overall) is at least `pi_preserve`.
  Default `0.5`.

## Value

A list with numeric vectors `overall`, `non_region` and `region` (each
with elements `hr`, `lo`, `hi`, `n`, `d`; `NA` HR and limits when fewer
than 10 subjects or 5 events), `consistency_ratio`, logical
`method1_pass`, `method2_pass` and `region_hr_below_1`, and, when
`subgroup_index` is supplied, `region_subgroup` and
`non_region_subgroup`.
