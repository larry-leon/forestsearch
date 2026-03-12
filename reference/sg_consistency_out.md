# Output Subgroup Consistency Results

Returns the top subgroup(s) and recommended treatment flags.

## Usage

``` r
sg_consistency_out(
  df,
  result_new,
  sg_focus,
  index.Z,
  names.Z,
  details = FALSE,
  plot.sg = FALSE,
  by.risk = 12,
  confs_labels
)
```

## Arguments

- df:

  Data.frame. Original analysis data.

- result_new:

  Data.table. Sorted subgroup results.

- sg_focus:

  Character. Sorting focus criterion.

- index.Z:

  Matrix. Subgroup factor indicators.

- names.Z:

  Character vector. Factor column names.

- details:

  Logical. Print details.

- plot.sg:

  Logical. Plot subgroup curves.

- by.risk:

  Numeric. Risk interval for plotting.

- confs_labels:

  Character vector. Human-readable labels.

## Value

List with results, subgroup definition, labels, flags, and group id.

## Examples

``` r
if (FALSE) { # \dontrun{
# sg_consistency_out is called internally by forestsearch().
# See forestsearch() for the standard entry point.
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
fs$grp.consistency$out_sg
} # }
```
