# Print Method for forestsearch Objects

Displays a concise summary of ForestSearch results including the
identified subgroup definition, consistency metrics, algorithm details,
and computation time.

## Usage

``` r
# S3 method for class 'forestsearch'
print(x, ...)
```

## Arguments

- x:

  A `forestsearch` object returned by
  [`forestsearch`](https://github.com/larry-leon/forestsearch/reference/forestsearch.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`summary.forestsearch`](https://github.com/larry-leon/forestsearch/reference/summary.forestsearch.md)
for detailed output,
[`plot.forestsearch`](https://github.com/larry-leon/forestsearch/reference/plot.forestsearch.md)
for visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(df.analysis = mydata, ...)
print(fs)
# or simply:
fs
} # }
```
