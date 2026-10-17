# Print detailed output for debugging

Displays detailed information about the GRF analysis

## Usage

``` r
print_grf_details(config, values, best_subgroup, sg_harm_id, tree_cuts = NULL)
```

## Arguments

- config:

  List. GRF configuration

- values:

  Data frame. Node metrics

- best_subgroup:

  Data frame row. Selected subgroup (or NULL)

- sg_harm_id:

  Character. Subgroup definition (or NULL)

- tree_cuts:

  List. Cut information

## Value

Called for its side effect of printing to the console. The value is
returned invisibly and is incidental: `tree_cuts$all` when a subgroup,
its definition and `tree_cuts` are all supplied, otherwise `NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
# print_grf_details() is called internally by grf.subg.harm.survival().
# See grf.subg.harm.survival() for the standard entry point.
} # }
```
