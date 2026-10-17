# Print Method for ForestSearch Bootstrap Results

Prints a concise but informative summary of an `fs_bootstrap` object
returned by
[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md).
Covers identification rate, agreement with the primary subgroup, top
identified subgroups, and bias-corrected effect estimates.

## Usage

``` r
# S3 method for class 'fs_bootstrap'
print(x, top_n = 3L, ...)
```

## Arguments

- x:

  An `fs_bootstrap` object.

- top_n:

  Integer. Maximum number of top identified subgroups to print. Default
  3.

- ...:

  Additional arguments (ignored; present for S3 consistency).

## Value

Invisibly returns `x`.

## Details

The print method is deliberately richer than the corresponding
`print.fs_kfold` / `print.fs_tenfold` methods because the bootstrap is
the primary inferential machinery of the package, while the CV routines
are complementary diagnostics. For the full per-factor breakdown,
consistency distribution, size distribution, and GRF cut tabulation,
call
[`summarize_bootstrap_subgroups`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_subgroups.md).

## See also

[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
to produce the object.
[`summarize_bootstrap_subgroups`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_subgroups.md)
for full tabulations.

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(
  df.analysis = mydata,
  confounders.name = c("age", "sex", "biomarker"),
  outcome.name = "time", event.name = "status", treat.name = "treat"
)
fs_bc <- forestsearch_bootstrap_dofuture(fs, nb_boots = 500)
print(fs_bc)          # rich default summary
print(fs_bc, top_n = 5)
} # }
```
