# Read a shipped record of the MRCT consistency campaign

The records behind the MRCT consistency articles ship with the package,
in `extdata/mrct-consistency`, with their provenance in its
`PROVENANCE.md`. From a source tree the local copy is used; otherwise
the installed one.

## Usage

``` r
mrct_record(file, ...)
```

## Arguments

- file:

  The file's name, for instance `"rep5_trial.csv"`.

- ...:

  Passed to
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html).

## Value

A data frame, read with `stringsAsFactors = FALSE` and
`check.names = FALSE`.

## See also

[`mrct_case_study()`](https://larry-leon.github.io/forestsearch/reference/mrct_case_study.md),
[`mrct_trial_analysis()`](https://larry-leon.github.io/forestsearch/reference/mrct_trial_analysis.md).

## Examples

``` r
trial <- mrct_record("rep5_trial.csv")
dim(trial)
#> [1] 1000   18
mrct_record("mr_field_S2_rep5.csv")$sg_found
#> [1] "{z_cd40 <= 236} & !{z_preanti <= 196}"
```
