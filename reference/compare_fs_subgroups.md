# Compare Subgroup Membership Across Two Analyses

Cross-tabulates subgroup assignments from two analyses and produces a
formatted `gt` summary table with concordance statistics. Accepts
`forestsearch` result objects, GRF result lists, or raw indicator
vectors. Non-detection (no subgroup found) is handled automatically: all
subjects are classified to the complement.

## Usage

``` r
compare_fs_subgroups(
  fs1,
  fs2,
  label1 = "Analysis 1",
  label2 = "Analysis 2",
  id.name = NULL,
  sg0_label = NULL,
  sg1_label = NULL,
  subgroup_notation = c("harm", "benefit"),
  title = "Subgroup Membership Comparison",
  subtitle = NULL,
  font_size = 12
)
```

## Arguments

- fs1:

  Either a `forestsearch` result object, a `grf.subg.harm` result (with
  `data$treat.recommend`), or an integer/logical vector of subgroup
  indicators (0 = in H/Q, 1 = in Hc/Qc).

- fs2:

  Same types as `fs1` for the second analysis.

- label1:

  Character. Display label for the first analysis.

- label2:

  Character. Display label for the second analysis.

- id.name:

  Character. Subject ID column for merge alignment. Only used when both
  inputs are result objects. If `NULL`, row-position alignment is
  assumed.

- sg0_label:

  Character or `NULL`. Label for H/Q (`treat.recommend == 0`). If
  `NULL`, auto-generated from `subgroup_notation`.

- sg1_label:

  Character or `NULL`. Label for Hc/Qc (`treat.recommend == 1`). If
  `NULL`, auto-generated.

- subgroup_notation:

  Character. `"harm"` uses H/Hc labels; `"benefit"` uses Q/Qc labels.
  Only used when `sg0_label`/`sg1_label` are `NULL`.

- title:

  Character. Table title.

- subtitle:

  Character or `NULL`. Auto-generated if NULL.

- font_size:

  Numeric. Font size for gt table.

## Value

A list with:

- table:

  A `gt` cross-tabulation table.

- crosstab:

  The raw cross-tabulation matrix.

- concordance:

  Agreement, kappa, counts.

- membership:

  Per-subject assignments from both analyses.

## Examples

``` r
if (FALSE) { # \dontrun{
# Harm search: Cox vs Poisson (both forestsearch objects)
result <- compare_fs_subgroups(fs_cox, fs_pois,
  label1 = "Cox PH", label2 = "Poisson",
  subgroup_notation = "harm")

# Benefit search: FS vs GRF (handles non-detection automatically)
result <- compare_fs_subgroups(fs_result, grf_result,
  label1 = "ForestSearch", label2 = "GRF",
  subgroup_notation = "benefit")

# Raw vectors
result <- compare_fs_subgroups(
  c(rep(0L, 30), rep(1L, 70)),
  c(rep(0L, 25), rep(1L, 75)),
  label1 = "Method A", label2 = "Method B")
} # }
```
