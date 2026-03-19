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
  return_selected_cuts_only = FALSE
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

## Value

A list with GRF results, including:

- data:

  Original data with added treatment recommendation flags

- grf.gsub:

  Selected subgroup information

- sg.harm.id:

  Expression defining the identified subgroup

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
