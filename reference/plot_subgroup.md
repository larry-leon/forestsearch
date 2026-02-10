# Plot Subgroup Survival Curves

Plots weighted Kaplan-Meier survival curves for a specified subgroup and
its complement using the weightedsurv package.

## Usage

``` r
plot_subgroup(df.sub, df.subC, by.risk, confs_labels, this.1_label, top_result)
```

## Arguments

- df.sub:

  A data frame containing data for the subgroup of interest.

- df.subC:

  A data frame containing data for the complement subgroup.

- by.risk:

  Numeric. The risk interval for plotting (passed to
  [`weightedsurv::df_counting`](https://rdrr.io/pkg/weightedsurv/man/df_counting.html)).

- confs_labels:

  Named character vector. Covariate label mapping (not used directly in
  this function, but may be used for labeling).

- this.1_label:

  Character. Label for the subgroup being plotted.

- top_result:

  Data frame row. The top subgroup result row, expected to contain a
  `Pcons` column for consistency criteria.
