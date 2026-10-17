# Format Pareto Frontier of Candidate Subgroups

Renders the post-hoc Pareto frontier on (effect, N) – both maximized –
as a formatted gt table or returns it as a `data.table` for programmatic
use. Works uniformly across survival (Cox PH) and GLM outcome types: the
effect-column label and scale handling are derived from the
`forestsearch` object's `effect_measure`.

## Usage

``` r
pareto_frontier_table(
  fs,
  format = c("gt", "data.table"),
  digits = 2L,
  digits_effect = NULL,
  digits_pcons = NULL,
  digits_ci = NULL,
  include_dominated = FALSE,
  include_factor_columns = TRUE,
  highlight_selected = TRUE,
  ci_table = NULL
)
```

## Arguments

- fs:

  A `forestsearch` object returned by
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- format:

  Character. Either `"gt"` (default) for a formatted gt table or
  `"data.table"` for a plain `data.table` with effect on the natural
  scale and an `is_selected` flag.

- digits:

  Integer. Master decimal-places setting that drives all three
  column-specific defaults (`digits_effect`, `digits_pcons`,
  `digits_ci`) when those are left at their defaults. Default `2L`. Pass
  `digits = 3L` for a tighter-precision view; pass per-column overrides
  below for asymmetric formatting (e.g. `digits_effect = 3L` with
  `digits_ci = 2L`).

- digits_effect:

  Integer or `NULL`. Decimal places for the effect column. When `NULL`
  (default), inherits from `digits`. Pass an explicit integer to
  override.

- digits_pcons:

  Integer or `NULL`. Decimal places for `Pcons`. When `NULL` (default),
  inherits from `digits`. Pass an explicit integer to override.

- digits_ci:

  Integer or `NULL`. Decimal places used for the three CI columns
  (Naive, Split, FSBC) when `ci_table` is supplied. When `NULL`
  (default), inherits from `digits`. Pass an explicit integer to
  override.

- include_dominated:

  Logical. If `FALSE` (default), the table shows only the
  Pareto-non-dominated subgroups (the frontier). If `TRUE`, the table
  shows **all** passing candidates – frontier members AND dominated
  candidates – with additional columns indicating frontier membership
  and band eligibility (when the selection rule uses one). This is the
  "rich post-consistency view" that complements the in-flight summary
  printed by `forestsearch(show_candidate_summary = TRUE)`.

- include_factor_columns:

  Logical. If `TRUE` (default), include the `M.<factor>` columns
  identifying each subgroup's defining cuts. If `FALSE`, show only
  summary metrics.

- highlight_selected:

  Logical. If `TRUE` (default), the selected subgroup row is highlighted
  in the gt table. Ignored when `format = "data.table"`.

- ci_table:

  Optional data.table of frontier CIs from
  [`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md).
  When supplied, three CI columns are merged into the table; when `NULL`
  (default), no CI columns are shown. This function does not run any
  computation on `ci_table` – it simply displays the bounds already
  computed.

## Value

Depending on `format`:

- `"gt"`:

  A `gt_tbl` object. Returns `NULL` invisibly (with a message) if the
  frontier is unavailable or empty.

- `"data.table"`:

  A `data.table` of the frontier with effect on the natural scale, an
  `is_selected` logical column, and effect column renamed to the
  `effect_measure` label (e.g., `"HR"`, `"OR"`). Returns an empty
  `data.table` if the frontier is unavailable.

## Details

The frontier is a diagnostic: it lists candidate subgroups that are *not
dominated* on (effect, N) simultaneously. It is computed inside
[`sg_consistency_out`](https://larry-leon.github.io/forestsearch/reference/sg_consistency_out.md)
(see
[`compute_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/compute_pareto_frontier.md))
and attached to the result object as
`fs$grp.consistency$out_sg$pareto_frontier`. It is **not** used for
selection – the selected subgroup is chosen by `sg_focus` and may or may
not appear on the frontier.

**Scale handling.** For ratio measures stored on the log scale
internally (`"OR"`, `"RR"`, `"IRR"`), the effect column is exponentiated
for display. For `"HR"` (Cox PH, natural scale) and identity-scale
measures (`"RD"`, `"IRD"`, `"MD"`), values pass through unchanged.

**Selected-row identification.** The selected subgroup is identified by
matching the original-table row index `m` against the top row of
`fs$grp.consistency$out_sg$result`. This is robust to sorting and
outcome-type differences.

**For `sg_focus = "hrMinSG"`**, the selected subgroup may not appear on
the frontier – that focus deliberately prefers small subgroups, which
are typically N-dominated.

**Optional confidence intervals.** When `ci_table` is supplied (the
data.table returned by
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md)),
three CI columns are added to the table. All three CIs are computed by
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md);
this function only displays what is passed in. Pass the same `ci_table`
to
[`plot_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_frontier.md)
for consistent display across table and plot.

- `Naive 95% CI` – full-sample Wald CI from a Cox or GLM refit on each
  frontier member's data. Ignores subgroup-search selection;
  anti-conservative by construction.

- `Split ~ 95% CI` – subsample-derived approximation, computed from the
  empirical SD of averaged half-sample effects across the splits used by
  `compute_frontier_cis`.

- `FSBC ~ 95% CI` – bias-corrected interval following the bootstrap
  algorithm of Leon2024fs (eq 7) but treating the selected subgroup as
  fixed across half-jackknife replicates. The cell shows
  `est (lcl, ucl)` where `est` is the bias-corrected effect estimate
  \\2\hat\beta - \overline{\hat\beta^{(h)}}\\.

See
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md)
for the algebra.

## See also

[`compute_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/compute_pareto_frontier.md),
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md),
[`plot_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_frontier.md),
[`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md),
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Survival example
data(gbsg, package = "survival")
fs <- forestsearch(gbsg, ...)

# Basic table (no CIs)
pareto_frontier_table(fs)
pareto_frontier_table(fs, format = "data.table")

# With CIs (compute once, display anywhere)
ci_tab <- compute_frontier_cis(fs, n_splits = 1000, seed = 1)
pareto_frontier_table(fs, ci_table = ci_tab)
plot_pareto_frontier(fs,   ci_table = ci_tab)  # same intervals
} # }
```
