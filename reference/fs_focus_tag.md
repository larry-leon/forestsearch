# Stem tag for a (subgroup_method, sg_focus) pair

Returns the short label used in result-file stems and plot annotations
for the rule that a given `(subgroup_method, sg_focus)` combination
actually runs – not the spelling the caller passed. The two differ
whenever an alias or a behavioural synonym is in play, which is
precisely when a hand-written tag goes wrong.

## Usage

``` r
fs_focus_tag(subgroup_method, sg_focus)
```

## Arguments

- subgroup_method:

  Character scalar. One of `"consistency"`, `"dina"`, `"grf"`. Any value
  other than `"consistency"` takes the effect-argmax collapse.

- sg_focus:

  Character scalar. Any accepted `sg_focus` spelling, alias or
  canonical.

## Value

Character scalar; the stem tag.

## Details

The collapse is engine-specific:

- `"consistency"`:

  `eff`, `hr` and `maxcons` all tag as `"maxcons"` – they name the
  consistency argmax. `maxeff` and `maxeffCons` stay distinct, because
  the consistency floor genuinely separates them here.

- `"dina"`, `"grf"`:

  `eff`, `hr`, `maxcons`, `maxeff` and `maxeffCons` all tag as `"eff"`.
  Neither engine computes a Pcons, so the consistency qualifier has
  nothing to bind to and all five rank by `order(-eff)`.

The band foci collapse the same way on every engine: `hrMaxSG` tags as
`"effMaxSG"` and `hrMinSG` as `"effMinSG"`. `maxSG` and `minSG` are
canonical foci in their own right and pass through.

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for the per-engine semantics of each focus.

## Examples

``` r
fs_focus_tag("consistency", "eff")   # "maxcons"
#> [1] "maxcons"
fs_focus_tag("dina",        "eff")   # "eff"
#> [1] "eff"
fs_focus_tag("dina",        "maxeffCons")  # "eff" -- not "maxeffCons"
#> [1] "eff"
fs_focus_tag("grf",         "hrMaxSG")     # "effMaxSG"
#> [1] "effMaxSG"
```
