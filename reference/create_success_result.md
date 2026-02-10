# Create result object for successful subgroup identification

Builds comprehensive result object when a subgroup is found

## Usage

``` r
create_success_result(
  data,
  best_subgroup,
  trees,
  tree_cuts,
  selected_tree,
  sg_harm_id,
  values,
  config
)
```

## Arguments

- data:

  Data frame. Original data with subgroup assignments

- best_subgroup:

  Data frame row. Selected subgroup information

- trees:

  List. All fitted policy trees

- tree_cuts:

  List. Cut information from trees

- selected_tree:

  Policy tree. The tree that identified the subgroup

- sg_harm_id:

  Character. Expression defining the subgroup

- values:

  Data frame. All node metrics

- config:

  List. GRF configuration

## Value

List with complete GRF results
