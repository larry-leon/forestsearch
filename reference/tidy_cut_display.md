# Tidy a cut / subgroup expression for display

Display-only cleanup for diagnostic printing of cut and
subgroup-definition strings. Collapses runs of whitespace to single
spaces (so padded GRF definitions like `"cd40 <= 320"` render cleanly)
and rounds each post-operator numeric threshold half-up to `digits`
places (matching the `collapse_cuts` representative rounding, e.g.
`"wtkg > 79.833600000000004"` -\> `"wtkg > 80"` at `digits = 0`).

## Usage

``` r
tidy_cut_display(s, digits = 0L)
```

## Arguments

- s:

  Character vector of expressions, e.g.
  `"cd40 <= 320 & wtkg > 79.8336"`.

- digits:

  Integer \>= 0. Rounding for embedded thresholds. Default 0 (nearest
  integer), matching the `collapse_cuts` default.

## Value

Character vector, tidied for display (same length as `s`).

## Details

Only numbers that follow a comparison operator (`<=`, `>=`, `==`, `<`,
`>`) are rounded, so integers embedded in variable names (`cd40`,
`prior_6mo`) are never altered. This function does NOT touch the
canonical, membership-exact labels used in the search; it is intended
only for the diagnostic strings passed to
[`cat()`](https://rdrr.io/r/base/cat.html).
