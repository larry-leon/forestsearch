# Invert Effect Columns in Simulation Results for Benefit Search

Transforms simulation results from the switched treatment scale to the
original scale. For ratio-scale measures (HR, OR, IRR), inversion is
`1/x`; for identity-scale measures (MD, RD, IRD), inversion is `-x`.
Called automatically by table functions when
`subgroup_notation = "benefit"`.

## Usage

``` r
invert_hr_columns(res, effect_measure = NULL)
```

## Arguments

- res:

  Data frame of simulation results (from
  [`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)).

- effect_measure:

  Character or `NULL`. The effect measure used in the simulation. When
  `NULL` (default), ratio-scale inversion (`1/x`) is used for backward
  compatibility with survival analyses.

## Value

Data frame with all effect-related columns inverted. Non-effect columns
(classification metrics, subgroup sizes, etc.) are unchanged.
