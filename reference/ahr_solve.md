# Solve P(trigger) = target in one free parameter

Scans `f` over a grid of the free parameter, brackets every crossing of
`target`, and refines each with
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html): the pattern
of the mrct-consistency boundary and regional-entry scripts. `f` maps
the free parameter to the probability (typically
`ahr_design(...)$prob[["p_trigger"]]`); when the parameter enters the
build (the region prognostic log-HR `r`, the entry delay), `f` builds
its own population. Roots are only found where the scan brackets them:
two crossings between adjacent scan points are missed.

## Usage

``` r
ahr_solve(f, target, interval, ..., n_scan = 25L, scan = NULL, tol = 1e-08)
```

## Arguments

- f:

  Function of the free parameter (first argument) returning one number.

- target:

  The value to solve for.

- interval:

  Length-2 range of the free parameter.

- ...:

  Further arguments passed to `f`.

- n_scan:

  Number of equally spaced scan points on `interval` (default 25);
  ignored when `scan` is given.

- scan:

  Optional explicit scan points (sorted, within `interval`).

- tol:

  Tolerance of
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) in the free
  parameter.

## Value

A list: `roots` (sorted), `n_roots`, `values` (`f` at the roots), `scan`
(a data.frame of the scan points `x` and `value`), and `evaluations`
(calls to `f`).

## See also

[`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md)

## Examples

``` r
if (FALSE) { # \dontrun{
f <- function(a) ahr_design(ahr_population(dgm, 40, 11,
  entry_delay = list(var = "z_region", delay = a)), 0.55, 2, 1000)$prob[["p_trigger"]]
ahr_solve(f, 0.25, c(0, 24), scan = 0:24)
} # }
```
