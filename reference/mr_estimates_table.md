# Multiplier-resampling (MR) estimates table (bootstrap-table format)

Renders the selected subgroup's MR de-biased estimate
([`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
with `mr_inference = TRUE`) in the `gt` layout of the bias-corrected
table produced by
[`summarize_bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md).

## Usage

``` r
mr_estimates_table(fs, boot_results = NULL, est.scale = "hr", digits = 2)
```

## Arguments

- fs:

  A
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  result with a non-NULL `mr_inference` element.

- boot_results:

  `NULL` (default) to build the table from `fs` alone, or a
  [`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
  result used as the template for the descriptive columns, the naive
  effect, the complement row, and the subgroup footnote. A template must
  contain `FSsg_tab` and `SG_CIs`.

- est.scale:

  Character, `"hr"` or `"1/hr"`; passed to
  [`format_bootstrap_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_table.md).
  With a template it should match the scale of `boot_results`; without
  one, `"1/hr"` inverts every ratio (limits swap, and a one-sided lower
  bound becomes an upper bound) and is accepted for ratio measures only.

- digits:

  Integer display precision (default `2`), applied with the same
  reformatter
  [`summarize_bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
  uses so the two tables align.

## Value

A `gt` table. Without a template it carries the attribute `"mr_values"`:
a data frame, one row per subgroup (`row` `"H"`, `"Hc"`), of the numbers
the table displays, on the effect scale of `fs$mr_inference` (before any
`"1/hr"` inversion): `naive_est`, `naive_lower`, `naive_upper`,
`debiased_est`, `ij_lower`, `ij_upper`, and, when present, `field_bound`
and `bonferroni_bound` (`NA` where a complement quantity was not
computed; the complement's naive columns are then the descriptive
table's own robust Cox fit, not stored here).

## Details

**With a bootstrap template** (`boot_results` supplied): the descriptive
columns (sample size, medians/RMST or rates) and the naive effect do not
depend on which de-biasing method was used, so they are inherited from a
full-bootstrap (FB) result used as a template; the bias-corrected column
is replaced with the MR values from `fs$mr_inference` – the de-biased
estimate for the selected subgroup and, when MR was run with
`include_complement = TRUE`, for the complement as well (otherwise the
complement cell is marked, since MR de-biases the selected subgroup
only).

**Without a template** (`boot_results = NULL`, the default): the table
is built from `fs` alone. The descriptive columns are computed from
`fs$df.est` (the analysis data with the declared rule's
`treat.recommend`) by the builder the bootstrap uses
([`SG_tab_estimates()`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates.md);
for GLM outcomes its GLM counterpart), one row for the selected subgroup
and one for its complement in the analysis data. The naive effect is
`fs$mr_inference$naive` (and `$complement$naive`) with its robust
(sandwich) limits, and its header says so. The MR columns follow: the
de-biased estimate with its IJ interval (`debiased`); when
`fs$mr_inference$field` exists, the field one-sided 95% bound
(`field$lower_1s` for the selected subgroup, `field$complement$upper_1s`
for the complement); and when the joint block exists, the Bonferroni
pair (`field$joint$bonf_lower_H`, `field$joint$bonf_upper_Hc`). The
field and Bonferroni columns appear only when those quantities exist.
Complement cells carry MR values only when MR was run with
`include_complement = TRUE` (and the complement field with
`field_complement = TRUE`); otherwise they are marked with an em dash,
as with a template.

With a template, the selected-subgroup row is matched by its
bias-corrected string (`boot_results$SG_CIs$H_bc`) rather than by
position, so the substitution is robust to row ordering. The MR
de-biased confidence interval is a leading-order approximation of the
full-bootstrap infinitesimal-jackknife interval – it comes from the IJ
variance computed on the multiplier draws under `ci_method = "ij"` or
`"field"` (the forest-search path defaults to `"field"`, which adds the
field bounds; the DINA and GRF paths default to `"ij"`), and from the
subgroup robust SE only under `ci_method = "wald"`. Either way the
bias-corrected cell may differ slightly from the full bootstrap even
when the naive effect matches exactly.

For GLM outcomes the effect columns follow the fit's orientation: the
naive and MR values are those of `fs$mr_inference`, on the scale of the
fit's own estimates; the descriptive table's estimator (which gives the
complement's naive cell when MR did not cover the complement) is built
with the fit's stored `adverse_outcome` or, when none is stored, with
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
rule (`TRUE` for binary and count outcomes, `FALSE` otherwise); and a
template's effect column comes from
[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md),
which builds its estimator with the same rule. With the default
`est.scale = "hr"`, the descriptive columns `Rate(C)` or `Mean(C)`,
`Rate(T)` or `Mean(T)`, and their raw difference `Diff` (treated minus
control) are the raw outcome and do not follow it. So for an unadjusted
fit where a lower outcome is adverse (`adverse_outcome = FALSE`, the
default for continuous outcomes) the effect and `Diff` have opposite
signs: an MD is `-Diff`, an RD is `-Diff / 100`, and the log of an OR or
RR has the sign opposite to `Diff`. Count outcomes and the rate measures
(IRR, IRD) are not reoriented (see
[`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md)).

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
