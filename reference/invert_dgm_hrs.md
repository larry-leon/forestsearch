# Invert DGM Effect Values for Benefit Search

Transforms a DGM object's truth values from the switched treatment scale
to the original scale. For ratio-scale measures (HR, OR, IRR), fields
are reciprocated (`1/x`); for identity-scale measures (MD, RD, IRD),
fields are negated (`-x`). The effect measure is auto-detected from
`dgm$effect_measure` when not explicitly provided.

## Usage

``` r
invert_dgm_hrs(dgm, effect_measure = NULL)
```

## Arguments

- dgm:

  A DGM object (`gbsg_dgm`, `aft_dgm_flex`, or `glm_dgm`).

- effect_measure:

  Character or `NULL`. Overrides auto-detection from
  `dgm$effect_measure`. When `NULL` (default) and `dgm$effect_measure`
  is also absent, ratio-scale (`1/x`) is used for backward
  compatibility.

## Value

Modified DGM object with all effect fields inverted. Non-effect fields
(super-population data, model parameters, etc.) are unchanged.

## Details

Called automatically by table functions when
`subgroup_notation = "benefit"`.
