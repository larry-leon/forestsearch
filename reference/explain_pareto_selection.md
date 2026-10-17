# Describe the Pareto-Frontier Selection in Words

Wraps
[`pareto_frontier_table`](https://larry-leon.github.io/forestsearch/reference/pareto_frontier_table.md)
with a verbal account of the selection: why the chosen subgroup
satisfies the configured `sg_focus` + `selection_rule` criterion, and,
for each non-dominated alternative, on which lexicographic axis it lost
to the winner.

## Usage

``` r
explain_pareto_selection(
  fs,
  ci_table = NULL,
  effect_label = NULL,
  digits = 3L,
  verbosity = c("medium", "compact", "detailed"),
  format = c("markdown", "character", "list")
)
```

## Arguments

- fs:

  A `forestsearch` result object.

- ci_table:

  Optional data.table from
  [`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md).
  When supplied, the Naive 95\\ Loser paragraphs do not reference CIs by
  design – the selection rule does not consider CIs.

- effect_label:

  Optional character override for the effect- measure label in prose,
  e.g. `"Hazard ratio"`. Defaults to a sensible string for the outcome
  type.

- digits:

  Integer. Decimal places for effect estimates and bound values in the
  prose. Default `3L`.

- verbosity:

  One of `"compact"`, `"medium"` (default), or `"detailed"`.

- format:

  One of `"markdown"` (default), `"character"`, or `"list"`.

## Value

Depending on `format`:

- `"markdown"`:

  A single string with embedded line breaks and bullet markers, ready
  for [`cat()`](https://rdrr.io/r/base/cat.html) in a Quarto chunk with
  `results = "asis"`.

- `"character"`:

  A character vector with one element per logical paragraph or bullet
  (no markdown markers).

- `"list"`:

  A structured list with named elements `criterion`, `selected`,
  `losers`, `closing`, and `meta`.

## Details

The wrapper is **descriptive**, not algorithmic: it reads the frontier
and the configuration off the returned `fs` object and renders prose. It
does not re-run the selection. A future change to the lexicographic key
in
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
would require updating this wrapper to remain accurate.

## Verbosity

- `"compact"`:

  Criterion sentence + selected line + a one-line bullet per loser.

- `"medium"` (default):

  Criterion sentence + selected paragraph (including the Naive CI if
  `ci_table` is given) + one short paragraph per loser, each naming the
  lex-key axis on which the loser was beaten.

- `"detailed"`:

  Medium + a closing paragraph that reproduces the formal lex key, the
  resolved indicator values, and a note about any ties broken
  downstream.

## Loser ordering

Losers are sorted by "closeness-to-winning", which is operationalised
as: **later loss-axis comes first**, then **smaller gap on that axis**.
Specifically, a loser that lost on `K` (the last tiebreaker) is reported
before one that lost on `N` (an early axis), and a loser that beat the
winner on every axis until the final one is reported before a loser that
fell off the band entirely. This puts the boundary cases first, which is
the most informative order for understanding how the rule made its
choice.

## See also

[`pareto_frontier_table`](https://larry-leon.github.io/forestsearch/reference/pareto_frontier_table.md),
[`plot_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_frontier.md),
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md).

## Examples

``` r
if (FALSE) { # \dontrun{
ci_tab <- compute_frontier_cis(fs, n_splits = 1000, seed = 1)

# Inline in a Quarto chunk
cat(explain_pareto_selection(fs, ci_table = ci_tab))

# Or get the structured object
explain_pareto_selection(fs, ci_table = ci_tab, format = "list")
} # }
```
