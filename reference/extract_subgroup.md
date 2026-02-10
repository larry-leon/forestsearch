# Extract Subgroup Information

Extracts subgroup definition and membership from results.

## Usage

``` r
extract_subgroup(df, top_result, index.Z, names.Z, confs_labels)
```

## Arguments

- df:

  Data.frame. Original analysis data.

- top_result:

  Data.table row. Top subgroup result.

- index.Z:

  Matrix. Factor indicators for all subgroups.

- names.Z:

  Character vector. Factor column names.

- confs_labels:

  Character vector. Human-readable labels.

## Value

List with sg.harm, sg.harm_label, df_flag, sg.harm.id.
