# Covariate involvement across the sweep

Reproduces the supplement's involvement figure: the share of detections
whose rule names the anchor, the true partner, and the proxy, against
the sweep dimension.

## Usage

``` r
plot_fs_identification_involvement(x, roles = NULL, palette = NULL)
```

## Arguments

- x:

  An `fs_identification` object.

- roles:

  Character vector. Which covariates to draw; defaults to the
  anchor/partner/proxy triple.

- palette:

  Named character vector of colours, or `NULL` for the colour-blind-safe
  default used in the supplement.

## Value

A `ggplot` object.
