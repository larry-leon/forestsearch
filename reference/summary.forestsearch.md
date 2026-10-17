# Summary Method for forestsearch Objects

Provides a detailed summary of a ForestSearch analysis including input
parameters, variable selection results, consistency evaluation, and the
selected subgroup with key metrics.

## Usage

``` r
# S3 method for class 'forestsearch'
summary(object, ...)
```

## Arguments

- object:

  A `forestsearch` object returned by
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `object`.

## Post-selection inference

Carries the same block as
[`print.forestsearch()`](https://larry-leon.github.io/forestsearch/reference/print.forestsearch.md)
(see there) and adds the log-scale standard errors (field, field-s with
the naive complement SE beside it, and IJ two-term), the three largest
re-selection frequencies with their labels, and `p_hat_sum` over the
re-selection family. Absent MR results, the output is unchanged.

When the fit carries the membership-agreement diagnostics – that is,
under `mr_inference_args = list(field_recovery = TRUE)` – one further
line reports `sens_H`, the mean share of the identified subgroup that
the field's re-selections retain, beside \\\hat p\\. The pair separates
two very different low-\\\hat p\\ situations: \\\hat p\\ counts *exact*
re-selections of \\\hat H\\, so a low \\\hat p\\ with a high `sens_H`
says the draws re-selected near-twins of \\\hat H\\, while a low \\\hat
p\\ with a low `sens_H` says they re-selected unrelated regions. Like
\\\hat p\\ it is descriptive and no construction reads it, and it
answers a narrower question than a bootstrap or cross-validation
recovery rate: re-selection within the fixed candidate family under
perturbation, not re-discovery from resampled data. The line is absent –
and the rest of the output byte-identical – when the diagnostics were
not requested.
[`print.forestsearch()`](https://larry-leon.github.io/forestsearch/reference/print.forestsearch.md)
does not report it.

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(df.analysis = mydata, ...)
summary(fs)
} # }
```
