# Validate Dataset for MRCT Simulations

Checks that a dataset contains all required variables for MRCT
simulation functions and reports any issues. Required variables include
outcome (tte, event), treatment (treat), continuous covariates (age,
bm), and factor covariates (male, histology, prior_treat, regA).

## Usage

``` r
validate_mrct_data(df.case, verbose = TRUE)
```

## Arguments

- df.case:

  Data frame to validate

- verbose:

  Logical. Print detailed validation results. Default: TRUE

## Value

Logical. TRUE if all requirements met, FALSE otherwise (invisibly)

## Details

### Required Variables

The function checks for the following variables:

- **Outcome**: tte (time-to-event), event (0/1 indicator)

- **Treatment**: treat (0/1 indicator)

- **Continuous**: age, bm (biomarker)

- **Factor**: male (0/1), histology, prior_treat (0/1), regA (0/1)

The function also validates variable types and value ranges.

## See also

[`create_dgm_for_mrct`](https://larry-leon.github.io/forestsearch/reference/create_dgm_for_mrct.md)
for creating DGM from validated data
