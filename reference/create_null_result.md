# Create result object when no subgroup is found

Builds result object for cases where no valid subgroup is identified

## Usage

``` r
create_null_result(data, values, trees, config)
```

## Arguments

- data:

  Data frame. Original data

- values:

  Data frame. Node metrics (may be empty)

- trees:

  List. Fitted policy trees

- config:

  List. GRF configuration

## Value

List with limited GRF results
