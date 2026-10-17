# True log-HR psi0(x) along x_pred with region histograms (paper Fig. 3 analogue)

True log-HR psi0(x) along x_pred with region histograms (paper Fig. 3
analogue)

## Usage

``` r
plot_region_modifier(dgm, file = NULL, main = NULL)
```

## Arguments

- dgm:

  An MRCT data-generating mechanism returned by
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md).
  Uses `df_super` (the effect modifier, the region indicator and
  `loghr_po`) and `hazard_ratios$AHR` for the reference line.

- file:

  Optional character path. When supplied, the plot is written to this
  PNG file (900 x 800 px) and the device is closed on exit; when `NULL`
  (default), it is drawn on the current device.

- main:

  Optional character title for the upper panel. When `NULL` (default),
  the title is `"True causal log-HR along <x_pred>"`.

## Value

`NULL`, invisibly. Called for its side effect of drawing a two-panel
plot.
