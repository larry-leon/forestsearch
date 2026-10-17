# Resolve a de-biased SE from the IJ variance, with graceful fallback

Prefers the bias-corrected IJ variance; if it is non-positive (too few
draws), falls back to the raw IJ variance, then to the subgroup robust
SE.

## Usage

``` r
.fs_mr_se_from_ij(ij, se_fallback)
```

## Arguments

- ij:

  Output of
  [`.fs_mr_ij_var()`](https://larry-leon.github.io/forestsearch/reference/dot-fs_mr_ij_var.md).

- se_fallback:

  Robust subgroup SE (`sigma_D`) used as last resort.

## Value

List with `se`, `var`, and `source` (`"ij"`, `"ij_raw"`, or
`"wald_fallback"`).
