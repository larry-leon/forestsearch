# Cross-Validation Subgroup Match Summary

Summarizes the match between cross-validation subgroups and analysis
subgroups.

## Usage

``` r
CV_sgs(sg1, sg2, confs, sg_analysis)
```

## Arguments

- sg1:

  Character vector. Subgroup 1 labels for each fold.

- sg2:

  Character vector. Subgroup 2 labels for each fold.

- confs:

  Character vector. Confounder names.

- sg_analysis:

  Character vector. Subgroup analysis labels.

## Value

List with indicators for any match, exact match, one match, and
covariate-specific matches, one value per fold (`sg1[i]`, `sg2[i]` are
the factor labels of fold `i`'s subgroup, `NA` for an empty slot):

- any_found:

  1 if the fold found a subgroup (either slot not `NA`), else 0.

- exact_match:

  1 if the fold's subgroup is the full-data subgroup `sg_analysis`,
  else 0. For a two-factor `sg_analysis`, both factors appear in the
  fold's two slots, in either order. For a one-factor `sg_analysis`, the
  fold's subgroup is that factor alone: one slot holds it and the other
  is `NA`. A fold subgroup that combines it with a second factor is not
  exact; it counts in `one_match`.

- one_match:

  1 if at least one full-data factor appears in either slot, else 0.

- cov1_any, cov2_any:

  1 if the fold's subgroup uses the covariate of the first (second)
  full-data factor, with any cut (see
  [`find_covariate_any_match`](https://larry-leon.github.io/forestsearch/reference/find_covariate_any_match.md)),
  else 0.

- cov1_exact, cov2_exact:

  1 if the first (second) full-data factor appears in either slot, else
  0.

`cov2_any` and `cov2_exact` are `NA` unless `sg_analysis` has two
factors. A comparison with `NA` counts as 0. When `sg_analysis` has no
element, or more than two, `exact_match`, `one_match`, `cov1_any` and
`cov1_exact` are 0 for every fold.

## Examples

``` r
if (FALSE) { # \dontrun{
# CV_sgs() is called internally by forestsearch_KfoldOut().
# See forestsearch_Kfold() and forestsearch_KfoldOut() for the entry points.
} # }
```
