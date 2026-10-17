# Estimator Standard Deviation at a Given Trial Size

Converts a
[`fs_dgm_scale`](https://larry-leon.github.io/forestsearch/reference/fs_dgm_scale.md)
object into the sampling standard deviation of the within-region
difference-in-means estimator at trial size `n`.

## Usage

``` r
fs_scale_se(scale, n, region = NULL, jensen = FALSE)
```

## Arguments

- scale:

  An object of class `"fs_dgm_scale"`.

- n:

  Integer. Trial size.

- region:

  Character. Region name, matching a row of `scale$regions`. Default
  `NULL` returns every region.

- jensen:

  Logical. If `TRUE`, inflate by the exact factor accounting for the
  random region size \\n_g \sim \mathrm{Bin}(n, P(g))\\ and random arm
  split \\n_1 \mid n_g \sim \mathrm{Bin}(n_g, p)\\, returning the
  unconditional standard deviation. Default `FALSE`, which returns the
  idealised (fixed-count) value.

## Value

Named numeric vector of standard deviations.

## See also

[`fs_dgm_scale`](https://larry-leon.github.io/forestsearch/reference/fs_dgm_scale.md)
