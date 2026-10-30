# Multiplier-resampling (MR) estimates table (bootstrap-table format)

Renders the selected subgroup's MR de-biased estimate
([`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
with `mr_inference = TRUE`) in the same `gt` layout as the
bias-corrected table produced by
[`summarize_bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md).
The descriptive columns (sample size, medians/RMST or rates) and the
naive effect do not depend on which de-biasing method was used, so they
are inherited from a full-bootstrap (FB) result used as a template; the
bias-corrected column is replaced with the MR values from
`fs$mr_inference` – the de-biased estimate for the selected subgroup
and, when MR was run with `include_complement = TRUE`, for the
complement as well (otherwise the complement cell is marked, since MR
de-biases the selected subgroup only).

## Usage

``` r
mr_estimates_table(fs, boot_results, est.scale = "hr", digits = 2)
```

## Arguments

- fs:

  A
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  result with a non-NULL `mr_inference` element.

- boot_results:

  A
  [`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
  result used as the template for the descriptive columns, the naive
  effect, the complement row, and the subgroup footnote. Must contain
  `FSsg_tab` and `SG_CIs`.

- est.scale:

  Character, `"hr"` or `"1/hr"`; passed to
  [`format_bootstrap_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_table.md).
  Should match the scale of `boot_results`.

- digits:

  Integer display precision (default `2`), applied with the same
  reformatter
  [`summarize_bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
  uses so the two tables align.

## Value

A `gt` table.

## Details

The selected-subgroup row is matched by its bias-corrected string
(`boot_results$SG_CIs$H_bc`) rather than by position, so the
substitution is robust to row ordering. The MR de-biased confidence
interval is a leading-order approximation of the full-bootstrap
infinitesimal-jackknife interval – it comes from the IJ variance
computed on the multiplier draws under `ci_method = "ij"` or `"field"`
(the forest-search path defaults to `"field"`, which adds the field
bounds; the DINA and GRF paths default to `"ij"`), and from the subgroup
robust SE only under `ci_method = "wald"`. Either way the bias-corrected
cell may differ slightly from the full bootstrap even when the naive
effect matches exactly.

## See also

[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for the `mr_inference` switch and the vocabulary section;
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
for the MR computation itself;
[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
for the full bootstrap that supplies the template;
[`fs_fdr_report()`](https://larry-leon.github.io/forestsearch/reference/fs_fdr_report.md)
for harm-confirmation threshold sweeps;
[`summarize_bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md),
[`format_bootstrap_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_table.md)
