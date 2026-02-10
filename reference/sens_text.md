# Generate Cross-Validation Sensitivity Text

Creates formatted text summarizing cross-validation agreement metrics.

## Usage

``` r
sens_text(fs_kfold, est.scale = "hr")
```

## Arguments

- fs_kfold:

  K-fold cross-validation results from forestsearch_Kfold.

- est.scale:

  Character. "hr" or "1/hr".

## Value

Character string with formatted CV metrics.
