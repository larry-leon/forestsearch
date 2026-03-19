# Extract cuts from selected tree only

Extracts cut information only from the tree at the specified selected
depth. This provides a focused set of cuts from the tree that identified
the subgroup meeting the `dmin.grf` criterion, rather than cuts from all
trees.

## Usage

``` r
extract_selected_tree_cuts(trees, selected_depth, maxdepth)
```

## Arguments

- trees:

  List. Policy trees (indexed by depth)

- selected_depth:

  Integer. Depth of the selected tree (from best_subgroup\$depth)

- maxdepth:

  Integer. Maximum tree depth (for populating tree-specific slots)

## Value

List with cuts and names, structured similarly to extract_all_tree_cuts
but with only the selected tree's cuts in the `all` field

## Details

This function is used when `return_selected_cuts_only = TRUE` in
[`grf.subg.harm.survival()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md).
It returns:

- `tree1`, `tree2`, `tree3`: Individual tree cuts (still populated for
  reference)

- `names1`, `names2`, `names3`: Individual tree variable names

- `all`: Cuts from the SELECTED tree only (not union of all trees)

- `all_names`: Variable names from the SELECTED tree only

- `selected_depth`: The depth that was selected
