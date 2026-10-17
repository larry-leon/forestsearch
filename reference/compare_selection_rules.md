# Compare forestsearch Runs Across Selection-Rule Combinations

Runs
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
once per combination of `sg_focus[i]` + `selection_rule[i]` (tuple
semantics; the two vectors are paired element-by-element and must have
the same length). All other arguments are held fixed across runs and
passed through via `...`.

## Usage

``` r
compare_selection_rules(
  df.analysis,
  sg_focus,
  selection_rule,
  compute_cis = TRUE,
  n_splits = 1000L,
  ci_seed = 1L,
  plot_xlim = NULL,
  show_band = TRUE,
  combo_labels = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- df.analysis:

  Analysis-ready data frame; passed to
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- sg_focus:

  Character vector of `sg_focus` values to try.

- selection_rule:

  Character vector of `selection_rule` values to try. Must have the same
  length as `sg_focus`.

- compute_cis:

  Logical. If `TRUE` (default), runs
  [`compute_frontier_cis()`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md)
  for each fit so the plot includes CI bars. Set to `FALSE` to skip the
  (often expensive) CI step.

- n_splits:

  Integer. Passed to `compute_frontier_cis(n_splits = )`. Default
  `1000L`.

- ci_seed:

  Integer. Passed to `compute_frontier_cis(seed = )`. Default `1L`.

- plot_xlim:

  Numeric vector of length 2 (or `NULL`). Common x-axis range to apply
  to all plots so they're directly comparable. Default `NULL` = let each
  plot auto-range.

- show_band:

  Logical. Whether to draw the effect-band shading on each plot. Default
  `TRUE`.

- combo_labels:

  Character vector of length `length(sg_focus)`, used as the title for
  each plot and key for the returned lists. Default `NULL` = auto,
  e.g.\\ `"hrMaxSG / pareto"`.

- verbose:

  Logical. Print one line per combination during execution. Default
  `TRUE`.

- ...:

  Additional arguments passed through to every
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  call (data names, thresholds, parallel args, etc.). Do NOT pass
  `sg_focus` or `selection_rule` here – they're handled by the wrapper.

## Value

A list with class `"forestsearch_comparison"` and components:

- `combos`:

  `data.frame` with columns `sg_focus`, `selection_rule`, `label`.

- `fs`:

  Named list of `forestsearch` result objects, one per combination. NULL
  entries mark fits that errored.

- `ci_tab`:

  Named list of frontier-CI tables (or all NULL if
  `compute_cis = FALSE`).

- `plots`:

  Named list of `ggplot` objects, one per combination.

- `plot_grid`:

  A single patchwork object placing the plots side by side, or `NULL` if
  patchwork isn't installed.

- `plot_combined`:

  A single `ggplot` composing all configurations onto one Pareto plot,
  with each winner labeled `S1: <combo_label>`, `S2: <combo_label>`,
  etc. Returned only when ALL successful configurations share an
  identical passing set (see
  [`plot_pareto_combined`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_combined.md)
  for the equality criterion); `NULL` otherwise. See
  `plot_combined_subsets` below for the per-subset fallback when some –
  but not all – configurations match.

- `plot_combined_subsets`:

  Named list of `ggplot` objects, one per equivalence class of
  configurations sharing an identical passing set. Only groups of size
  \\\ge 2\\ are included (a singleton group has no peer to combine
  with). Names are derived from the shared `sg_focus` when constant
  within a group (e.g., `"effMaxSG"`); from the concatenated combo
  labels otherwise. When all valid combos share one passing set this
  list contains exactly one element, equal to `plot_combined`. The list
  is empty when no two combos match.

- `combined_skip_reason`:

  Character scalar or `NULL`. When `plot_combined` is `NULL` because the
  equality precondition failed, this captures the specific reason – size
  mismatch, definition-set mismatch, or value drift on a named column –
  as a single string. `NULL` when the combined plot was built or when
  fewer than two fits succeeded.

- `console`:

  Named list of character vectors – the captured stdout from each
  forestsearch run (full output).

- `diagnostics`:

  Named list of per-combo slices of the captured output: each element is
  a list with character-scalar fields `preview` (the
  `CANDIDATE EVALUATION PREVIEW` block), `summary` (the
  `CANDIDATE EVALUATION SUMMARY` block), and `full` (the entire captured
  stdout, same content as `console[[i]]` collapsed to a single string).
  `preview` and `summary` are `NA_character_` when the corresponding
  banner is absent (e.g.\\ `show_candidate_summary` was `FALSE`, or the
  run errored before reaching it). This slice is what a comparison
  document should display as the primary diagnostic; `full` (or
  `console[[i]]`) remains available for an expandable / on-demand view
  of the complete run output.

- `errors`:

  Named list of error messages (or NULL) for combos that failed.

## Details

For each combination the wrapper:

1.  Captures the in-flight console output from forestsearch (including
    the pre-consistency **Candidate Evaluation Preview** and
    post-consistency **Candidate Evaluation Summary** when
    `show_candidate_summary = TRUE` is passed).

2.  Optionally computes the frontier CIs via
    [`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md).

3.  Builds the Pareto-frontier plot via
    [`plot_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_frontier.md).

The captured stdout is returned per combination so the caller can
[`cat()`](https://rdrr.io/r/base/cat.html) it later for side-by-side
inspection in a Quarto chunk. Plot objects are returned individually and
(when patchwork is installed) composed into a single side-by-side
figure. The full `fs` object for each combination is returned so
downstream diagnostics
([`pareto_frontier_table`](https://larry-leon.github.io/forestsearch/reference/pareto_frontier_table.md),
[`explain_pareto_selection`](https://larry-leon.github.io/forestsearch/reference/explain_pareto_selection.md))
can be run without re-fitting.

## Tuple semantics

With `sg_focus = c("hrMaxSG", "hrMaxSG")` and
`selection_rule = c("pareto", "both")`, two runs are performed: one with
`(hrMaxSG, pareto)` and one with `(hrMaxSG, both)`. The length-mismatch
case raises an error.

## Console capture

All stdout from each
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
call is captured (not just the diagnostic tables) via
[`capture.output()`](https://rdrr.io/r/utils/capture.output.html). The
PREVIEW / SUMMARY blocks are clearly delineated by banner separators
within the captured text. Use `cat(out$console[[i]])` to replay the
output in a Quarto chunk; `results = "asis"` works but is not required
since the captured text already includes newlines.

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md),
[`plot_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_frontier.md),
[`plot_pareto_combined`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_combined.md),
[`pareto_frontier_table`](https://larry-leon.github.io/forestsearch/reference/pareto_frontier_table.md),
[`explain_pareto_selection`](https://larry-leon.github.io/forestsearch/reference/explain_pareto_selection.md),
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md).

## Examples

``` r
if (FALSE) { # \dontrun{
out <- compare_selection_rules(
  df.analysis      = actg_df,
  sg_focus         = c("hrMaxSG", "hrMaxSG"),
  selection_rule   = c("pareto",  "both"),
  # All other forestsearch args:
  confounders.name = confounders.name,
  outcome.name     = adverse_outcome,
  treat.name       = treat.name,
  id.name          = id.name,
  outcome_type     = "continuous",
  effect_measure   = "MD",
  adverse_outcome  = TRUE,
  hr.threshold     = 10,
  hr.consistency   = 5,
  pconsistency.threshold = 0.90,
  show_candidate_summary = TRUE,
  details          = TRUE
)

# Inspect captured console output for combo 1
cat(out$console[[1]])

# Side-by-side plot
print(out$plot_grid)

# Use fs object for downstream diagnostics
pareto_frontier_table(out$fs[[1]], ci_table = out$ci_tab[[1]])
explain_pareto_selection(out$fs[[2]], ci_table = out$ci_tab[[2]])
} # }
```
