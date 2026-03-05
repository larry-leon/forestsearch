# Extract all cuts from fitted trees

Consolidates cut information from all fitted policy trees. This is the
default behavior that returns cuts from all trees regardless of which
tree identified the selected subgroup.

## Usage

``` r
extract_all_tree_cuts(trees, maxdepth)
```

## Arguments

- trees:

  List. Policy trees (indexed by depth)

- maxdepth:

  Integer. Maximum tree depth

## Value

List with cuts and names for each tree and combined
