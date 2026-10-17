# Random-benchmark subgroup specification

Describes the size-matched random benchmark subgroups regenerated in
every simulated trial. With the default `nested = TRUE`, a single index
draw of `max(sizes)` patients is taken per trial and the smaller
benchmarks are its prefixes (`random15` \\\subset\\ `random20`
\\\subset\\ `random40` \\\subset\\ `random60`), exactly as in the
extreme-subgroups vignettes. Membership columns are named
`paste0(prefix, sizes)` and take values 0/1.

## Usage

``` r
benchmark_spec(sizes = c(60L, 40L, 20L, 15L), nested = TRUE, prefix = "random")
```

## Arguments

- sizes:

  Integer vector of benchmark sizes. Order determines the
  column-creation order only; the draw uses `max(sizes)`.

- nested:

  Logical. `TRUE` (default): one draw, prefix nesting – the vignette
  scheme. `FALSE`: an independent draw per size, in the order given
  (consumes more RNG; not the vignette scheme).

- prefix:

  Column-name prefix, default `"random"`.

## Value

An object of class `"benchmark_spec"`.

## See also

[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md)

## Examples

``` r
benchmark_spec()
#> $sizes
#> [1] 60 40 20 15
#> 
#> $nested
#> [1] TRUE
#> 
#> $prefix
#> [1] "random"
#> 
#> attr(,"class")
#> [1] "benchmark_spec" "list"          
benchmark_spec(sizes = c(100, 50), prefix = "rnd")
#> $sizes
#> [1] 100  50
#> 
#> $nested
#> [1] TRUE
#> 
#> $prefix
#> [1] "rnd"
#> 
#> attr(,"class")
#> [1] "benchmark_spec" "list"          
```
