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
covariate-specific matches.
