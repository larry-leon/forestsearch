# Create Factor Summary Tables from Bootstrap Results

Generates formatted GT tables summarizing factor frequencies from
bootstrap subgroup analysis. Creates two complementary tables: one
showing factor selection frequencies within each position (M.1, M.2,
etc.), and another showing overall factor frequencies across all
positions.

## Usage

``` r
create_factor_summary_tables(factor_freq, n_found, min_percent = 2)
```

## Arguments

- factor_freq:

  Data.frame or data.table. Factor frequency table from
  [`summarize_bootstrap_subgroups`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_subgroups.md),
  containing columns:

  - `Position`: Position identifier (e.g., "M.1", "M.2")

  - `Factor`: Factor definition string

  - `N`: Count of times factor appeared in this position

  - `Percent`: Percentage relative to times position was populated

- n_found:

  Integer. Number of successful bootstrap iterations (where a subgroup
  was identified). Used to calculate overall percentages.

- min_percent:

  Numeric. Minimum percentage threshold for including factors in the
  tables. Factors with selection frequencies below this threshold are
  excluded. Default is 2 (i.e., 2%).

## Value

A list with up to two GT table objects:

- `by_position`:

  GT table showing factor frequencies within each position. Percentages
  represent conditional probability of factor selection given that the
  position was populated. Within each position, percentages sum to
  approximately 100% (may not sum exactly to 100% after filtering).

- `overall`:

  GT table showing total factor frequencies across all positions.
  Includes additional columns indicating which positions each factor
  appeared in and how many unique positions used the factor. Percentages
  represent proportion of successful iterations where the factor
  appeared in any position.

If no factors meet the minimum threshold, the corresponding table
element will be NULL.

## Note

This function requires the gt package for table creation. The overall
table also requires dplyr for data aggregation. If dplyr is not
available, only the position-specific table will be created and the
overall element will be NULL.

Always check for NULL before using the returned tables:

    if (!is.null(factor_tables$by_position)) {
      print(factor_tables$by_position)
    }

If all factors have percentages below `min_percent`, both table elements
will be NULL.

## See also

- [`summarize_bootstrap_subgroups`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_subgroups.md)
  for generating the factor_freq input

- [`format_subgroup_summary_tables`](https://larry-leon.github.io/forestsearch/reference/format_subgroup_summary_tables.md)
  for creating all subgroup summary tables

- [`summarize_bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
  for complete bootstrap analysis workflow

- [`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
  for running bootstrap analysis
