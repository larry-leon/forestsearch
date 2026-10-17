# GRF Subgroup Identification for Survival Data

Identifies subgroups with differential treatment effect using
generalized random forests (GRF) and policy trees. This function uses
causal survival forests to identify heterogeneous treatment effects and
policy trees to create interpretable subgroup definitions.

## Usage

``` r
grf.subg.harm.survival(
  data,
  confounders.name,
  outcome.name,
  event.name,
  id.name,
  treat.name,
  frac.tau = 1,
  n.min = 60,
  dmin.grf = 0,
  RCT = TRUE,
  details = FALSE,
  sg.criterion = "mDiff",
  maxdepth = 2,
  seedit = 8316951,
  return_selected_cuts_only = FALSE,
  tune_grf = FALSE,
  grf_selection = c("tree", "frontier"),
  frontier_rule = c("effMaxSG", "eff", "maxSG", "minSG", "effMinSG"),
  effect_neighborhood = 0.1,
  selection_rule = c("neighborhood", "pareto", "both")
)
```

## Arguments

- data:

  Data frame containing the analysis data.

- confounders.name:

  Character vector of confounder variable names.

- outcome.name:

  Character. Name of outcome variable (e.g., time-to-event).

- event.name:

  Character. Name of event indicator variable (0/1).

- id.name:

  Character. Name of ID variable.

- treat.name:

  Character. Name of treatment group variable (0/1).

- frac.tau:

  Numeric. Fraction of tau for GRF horizon (default: 1.0).

- n.min:

  Integer. Minimum subgroup size (default: 60).

- dmin.grf:

  Numeric. Minimum difference in subgroup mean (default: 0.0).

- RCT:

  Logical. Is the data from a randomized controlled trial? (default:
  TRUE)

- details:

  Logical. Print details during execution (default: FALSE).

- sg.criterion:

  Character. Subgroup selection criterion ("mDiff" or "Nsg").

- maxdepth:

  Integer. Maximum tree depth (1, 2, or 3; default: 2).

- seedit:

  Integer. Random seed (default: 8316951).

- return_selected_cuts_only:

  Logical. If TRUE, returns only cuts from the tree depth that
  identified the selected subgroup meeting `dmin.grf`. If FALSE
  (default), returns all cuts from all fitted trees (depths 1 through
  `maxdepth`).

- tune_grf:

  Logical. If `TRUE`, enables cross-validated hyperparameter tuning via
  `tune.parameters = "all"` in the causal survival forest. Default
  `FALSE`.

- grf_selection:

  Character, one of `"tree"` (default) or `"frontier"`. `"tree"` uses
  the GRF policy tree (standard behavior). `"frontier"` is an
  **experimental** alternative that, on the same doubly-robust scores,
  enumerates single-covariate thresholds (and depth-2 covariate-pair
  conjunctions when `maxdepth >= 2`), takes the Pareto frontier of
  (harm-effect, size), and selects one subgroup under `frontier_rule`.
  The selected subgroup is always a single conjunction. In small
  benchmarks the policy tree matched or beat the frontier on harm
  recovery, so `"frontier"` is provided for comparison and exploration,
  not as a recommended default.

- frontier_rule:

  Character, one of `"effMaxSG"` (default), `"eff"`, `"maxSG"`,
  `"minSG"`, or `"effMinSG"`; the rule applied to the frontier when
  `grf_selection = "frontier"`. Semantics match
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)'s
  `sg_focus`: `"eff"` takes the maximum harm-effect candidate;
  `"maxSG"`/`"minSG"` the largest/smallest eligible candidate (effect
  \>= `dmin.grf`); `"effMaxSG"`/`"effMinSG"` the largest/smallest
  candidate within `effect_neighborhood` (relative) of the max
  harm-effect.

- effect_neighborhood:

  Numeric in (0, 1); relative neighborhood for the `"effMaxSG"` /
  `"effMinSG"` rules. Default 0.10. Used only when
  `grf_selection = "frontier"`.

- selection_rule:

  Character, one of `"neighborhood"` (default), `"pareto"`, `"both"`;
  the inclusion band for the `eff*SG` frontier rules, matching
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  and
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md).
  Inert for the other rules.

## Value

A list with GRF results, including:

- data:

  Original data with added treatment recommendation flags

- grf.gsub:

  Selected subgroup information

- sg.harm.id:

  **Character vector of cut expressions** defining the identified
  subgroup (length = depth of the selected policy tree), or `NULL` if no
  subgroup was found. **Not** a per-subject membership indicator; see
  the *Field naming collision* section below.

- tree.cuts:

  Cut expressions - either all cuts from all trees (if
  `return_selected_cuts_only = FALSE`) or only cuts from the selected
  tree depth (if `return_selected_cuts_only = TRUE`)

- tree.names:

  Unique variable names used in cuts

- tree:

  Selected policy tree object

- tau.rmst:

  Time horizon used for RMST

- harm.any:

  All subgroups with positive treatment effect difference

- selected_depth:

  Depth of the tree that identified the subgroup (when found)

- return_selected_cuts_only:

  Logical indicating which cut extraction mode was used

Additional tree-specific cuts and objects (tree1, tree2, tree3) based on
maxdepth

## Details

The `return_selected_cuts_only` parameter controls which cuts are
returned:

- FALSE (default):

  Returns all cuts from all fitted trees (depths 1 to `maxdepth`). This
  provides the full set of candidate splits for downstream exploration
  and is the original behavior for backward compatibility.

- TRUE:

  Returns only cuts from the tree at the depth that identified the
  "winning" subgroup meeting the `dmin.grf` criterion. This is useful
  when you want a focused set of cuts associated with the selected
  subgroup, reducing noise from non-selected trees.

When `return_selected_cuts_only = TRUE` and no subgroup meets the
criteria, `tree.cuts` will be empty (character(0)).

## Field naming collision with forestsearch / subgroup.consistency

The field name `sg.harm.id` has **different semantics** on this object
versus on the `grp.consistency` list returned by
[`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
(and nested inside the
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
result):

|  |  |  |
|----|----|----|
| **Object** | **`sg.harm.id` contains** | **Length / type** |
| `grf.subg.harm.survival()` result (this function) | character vector of cut expressions | character, length = depth of selected tree |
| [`subgroup.consistency()`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md) result | per-subject 0/1 membership indicator | integer, length `nrow(data)` |

**Practical consequence.** On this object, pasting
`paste(obj$sg.harm.id, collapse = " & ")` correctly renders the
identified subgroup. On a `grp.consistency` list the same expression
silently concatenates a long 0/1 indicator and produces output like
`"0 & 0 & 1 & 0 & ..."`. When writing code that must handle both object
types, dispatch on class first, or use the top-level `sg.harm` field
that the
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
result exposes.

This naming collision is a documented CRAN-stable API for v0.1.x and
v0.2.x; it is expected to be resolved via a deprecation cycle in a
future minor release.

## Examples

``` r
if (FALSE) { # \dontrun{
# Return all cuts (default behavior)
result_all <- grf.subg.harm.survival(
  data = trial_data,
  confounders.name = c("age", "biomarker", "region"),
  outcome.name = "tte",
  event.name = "event",
  id.name = "id",
  treat.name = "treat",
  dmin.grf = 0.1,
  maxdepth = 2
)
result_all$tree.cuts
# Returns cuts from both depth 1 and depth 2 trees

# Return only cuts from the selected tree
result_selected <- grf.subg.harm.survival(
  data = trial_data,
  confounders.name = c("age", "biomarker", "region"),
  outcome.name = "tte",
  event.name = "event",
  id.name = "id",
  treat.name = "treat",
  dmin.grf = 0.1,
  maxdepth = 2,
  return_selected_cuts_only = TRUE
)
result_selected$tree.cuts
# Returns cuts only from the depth that identified the winning subgroup
} # }
```
