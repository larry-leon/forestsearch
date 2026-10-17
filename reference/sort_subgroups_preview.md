# Sort Subgroups by Focus (pre-consistency)

Sorts a data.table of candidate subgroups before consistency evaluation.
Mirrors
[`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
but operates on the pre-consistency column convention (`HR`, `n`, `K`)
and omits `Pcons` tiebreakers (consistency not yet available).

## Usage

``` r
sort_subgroups_preview(
  result_new,
  sg_focus,
  selection_rule = "neighborhood",
  effect_neighborhood = 0.1,
  effect_log_scale = FALSE
)
```

## Arguments

- result_new:

  A data.table of candidate subgroups with columns `HR`, `n`, `K`.

- sg_focus:

  Character. Sorting focus. One of `"hr"`, `"maxSG"`, `"hrMaxSG"`,
  `"minSG"`, `"hrMinSG"`.

- selection_rule:

  Character. See
  [`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md).
  Default `"neighborhood"`.

- effect_neighborhood:

  Numeric in `[0, 1)`. See
  [`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md).
  Default `0.10`.

- effect_log_scale:

  Logical. See
  [`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md).
  Default `FALSE`.

## Value

A sorted data.table.
