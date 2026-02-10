# Assign data to subgroups based on selected node

Creates treatment recommendation flags based on identified subgroup

## Usage

``` r
assign_subgroup_membership(data, best_subgroup, trees, X)
```

## Arguments

- data:

  Data frame. Original data

- best_subgroup:

  Data frame row. Selected subgroup information

- trees:

  List. Policy trees

- X:

  Matrix. Covariate matrix

## Value

Data frame with added predict.node and treat.recommend columns
