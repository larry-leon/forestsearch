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
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## Post-selection inference

When the fit carries multiplier-resampling results
(`mr_inference = TRUE`, attached as `x$mr_inference`) and the field
block ran (`ci_method = "field"`, the default), a further block reports:
the one-sided 95% lower bound for the identified subgroup (the field);
the one-sided 95% upper bound for its complement (field-s, the
studentized complement field); the Bonferroni pair for a joint
two-subgroup statement; the infinitesimal-jackknife (IJ) two-term
two-sided interval, as a secondary summary; and the re-selection
frequency \\\hat p(\widehat H)\\ as a descriptive diagnostic. Bounds are
read by their location relative to a clinically meaningful effect size,
not as tests at the null. With `ci_method = "ij"` only the two-sided
interval is shown; without MR results the output is unchanged. The
constructions and their simulation operating characteristics are
described in León and Anderson (2026),
[arXiv:2609.38361](https://arxiv.org/abs/2609.38361).

## See also

[`summary.forestsearch`](https://larry-leon.github.io/forestsearch/reference/summary.forestsearch.md)
for detailed output,
[`plot.forestsearch`](https://larry-leon.github.io/forestsearch/reference/plot.forestsearch.md)
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
