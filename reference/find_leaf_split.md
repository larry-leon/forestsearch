# Find the split that leads to a specific leaf node (deprecated)

**Deprecated.**

Returns only the single split immediately above `leaf_node`, rendered as
`"var <= value"`. This is **not** a correct subgroup definition for a
leaf below depth 1: it drops the rest of the root-to-leaf path and
mis-renders right-turns (which are `>`, not `<=`). Use the path-based
`.grf_build_subgroup_definition()` instead, which
[`grf.subg.harm.survival()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
now calls internally. Retained (exported) only for backward
compatibility; it is no longer used inside the package and will be
removed in a future release.

## Usage

``` r
find_leaf_split(tree, leaf_node)
```

## Arguments

- tree:

  Policy tree object

- leaf_node:

  Integer. Leaf node identifier

## Value

Character string with a single split expression, or NULL.

## See also

[`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
for the standard entry point.

## Examples

``` r
if (FALSE) { # \dontrun{
# Deprecated; grf.subg.harm.survival() builds the full path internally.
} # }
```
