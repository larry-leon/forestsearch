# Summarise ForestSearch Cross-Validation Diagnostics

Post-hoc diagnostic summary of
[`forestsearch_tenfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
output, aggregating the sim x fold grid captured in `fold_summary` into
publication-ready `gt` tables. Parallel to
[`summarize_bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
but focused on cross-validation stability and GRF candidate-discovery
behaviour.

## Usage

``` r
summarize_cv_results(
  cv_output,
  original_sg = NULL,
  original_grf_cuts = NULL,
  create_plots = FALSE,
  top_n = 15L
)
```

## Arguments

- cv_output:

  List of class `"fs_tenfold"`, as returned by
  [`forestsearch_tenfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md).
  Must contain a non-empty `fold_summary` data frame.

- original_sg:

  Character vector or `NULL`. Original full-data subgroup definition
  (e.g., `fs$sg.harm`). When supplied, `original_agreement` reports
  exact and partial-factor match rates against this definition.

- original_grf_cuts:

  Character vector or `NULL`. Original full-data `grf_cuts` (e.g.,
  `fs$grf_cuts`). When supplied alongside a `grf_cuts` column in
  `fold_summary`, `original_agreement` also reports GRF-cut match rate.

- create_plots:

  Logical. If `TRUE` and ggplot2 is installed, produce bar-chart
  diagnostics (default: `FALSE`).

- top_n:

  Integer. Maximum number of detail rows to retain in frequency tables
  (default: `15L`). If the true number of distinct values exceeds
  `top_n`, a single `"(N others)"` row aggregates the tail so that
  `n_folds` totals remain honest. Set to a large value (e.g., `Inf`) to
  disable truncation entirely. Plots use `min(top_n, 10L)`.

## Value

An object of class `"fs_cv_summary"`, a list with components:

- `identification_summary`, `grf_cut_summary`, `cut_vs_subgroup_xtab`,
  `no_subgroup_decomposition`, `pconsistency_distribution`,
  `fold_numeric_summary`, `original_agreement`, `metrics_table`:

  `gt` tables (or `data.frame` fallback if gt is unavailable). Some may
  be `NULL` depending on input columns and optional arguments.

- `plots`:

  `NULL` or a list of `ggplot` objects.

- `data`:

  List of raw `data.frame` tabulations underlying the tables
  (`identification`, `grf_cuts`, `cut_vs_subgroup`, `no_subgroup`,
  `pconsistency`, `fold_numeric_summary`, `original_agreement`).

- `n_sims`, `n_folds`, `total_pairs`:

  Grid dimensions.

- `has_grf_cuts`, `has_pconsistency`:

  Logical. Whether the corresponding `fold_summary` columns were
  available. The `fold_numeric_summary` slot adapts automatically to
  whichever numeric diagnostic columns are present and does not have a
  dedicated flag.

## Diagnostic components

The returned object includes (as `gt` tables by default):

- `identification_summary`: top identified subgroups across all sim x
  fold pairs, with counts and percentages.

- `grf_cut_summary`: top raw `grf_cuts` strings across training sets –
  how often GRF returned each output (and how often it returned none).
  `NULL` if `fold_summary$grf_cuts` is not available (e.g., objects from
  forestsearch prior to the GLM-extension branch).

- `cut_vs_subgroup_xtab`: cross-tabulation of unique `grf_cuts` strings
  against identified subgroups. Reveals whether GRF's discovery aligns
  with the pipeline's final selection.

- `no_subgroup_decomposition`: tabulation of folds that identified no
  subgroup, broken into "GRF returned no cut" vs "GRF returned a cut but
  consistency rejected it".

- `pconsistency_distribution`: binned distribution of the consistency
  probability (`Pcons`) achieved by the identified subgroup on each
  fold, restricted to folds that identified a subgroup; includes median
  and IQR. Reveals whether identifying folds are scraping
  `pconsistency.threshold` (e.g., Pcons just above 0.90) or are well
  above it. `NULL` if `fold_summary$pconsistency` is not available.
  Rejected- candidate Pcons values are not surfaced – only the selected
  candidate's Pcons is captured in `fold_summary`.

- `fold_numeric_summary`: tidy median / IQR / range table with one row
  per numeric diagnostic column present in `fold_summary` (`n_test`,
  `pconsistency`, `training_fs_hr`, `n_candidates_evaluated`). Absent
  columns simply don't produce rows, so the table adapts to
  pre-feature-branch objects without error. Training-fold subgroup
  effects are on the effect measure's natural scale (HR for survival;
  OR, RR, IRR, RD, IRD, or MD for GLM, exponentiated at capture for
  ratio measures so the column is always on natural scale regardless of
  outcome type); they are in-sample estimates (optimistically biased)
  and are surfaced for diagnostic comparison across folds only.

- `original_agreement`: if `original_sg` and/or `original_grf_cuts` are
  supplied, fraction of folds matching the original full-data analysis
  (exact subgroup, partial shared factor, GRF-cut match, both).

- `metrics_table`: formatted version of `cv_output$sens_summary` and
  `cv_output$find_summary` (delegates to
  [`cv_metrics_tables`](https://larry-leon.github.io/forestsearch/reference/cv_metrics_tables.md)).

- `plots`: if `create_plots = TRUE`, bar charts of top identified
  subgroups and top GRF cuts.

Raw `data.frame` versions of each tabulation are available at `$data$*`
for custom post-processing.

## Backward compatibility

The GRF-dependent slots (`grf_cut_summary`, `cut_vs_subgroup_xtab`,
`no_subgroup_decomposition`) require a `grf_cuts` column in
`fold_summary`; the `pconsistency_distribution` slot requires a
`pconsistency` column. Both were added on the `feature/glm-extension`
branch (`grf_cuts` first, `pconsistency` later). Older objects are
detected and the affected slots return `NULL` without error. Metadata
fields `has_grf_cuts` and `has_pconsistency` record which columns were
present.

## See also

[`forestsearch_tenfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
for running repeated K-fold CV;
[`cv_metrics_tables`](https://larry-leon.github.io/forestsearch/reference/cv_metrics_tables.md)
for the underlying metrics formatter;
[`summarize_bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
for the bootstrap analogue.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- survival::gbsg
df$grade3      <- as.integer(df$grade == "3")
df$time_months <- df$rfstime / 30.4375

fs <- forestsearch(
  df.analysis      = df,
  confounders.name = c("age", "meno", "size", "grade3",
                       "nodes", "pgr", "er"),
  outcome.name     = "time_months",
  event.name       = "status",
  treat.name       = "hormon"
)

tf <- forestsearch_tenfold(fs.est = fs, sims = 100, Kfolds = 10)

# Basic diagnostic summary
cv_diag <- summarize_cv_results(tf)
cv_diag$identification_summary
cv_diag$grf_cut_summary
cv_diag$no_subgroup_decomposition

# With original-analysis comparison and plots
cv_diag_full <- summarize_cv_results(
  tf,
  original_sg        = fs$sg.harm,
  original_grf_cuts  = fs$grf_cuts,
  create_plots       = TRUE
)
cv_diag_full$original_agreement
cv_diag_full$plots$identification
} # }
```
