# Summarize Bootstrap Subgroup Analysis Results

Comprehensive summary of bootstrap subgroup identification results
including basic statistics, factor frequencies, consistency
distributions, and agreement with the original analysis subgroup.

## Usage

``` r
summarize_bootstrap_subgroups(results, nb_boots, original_sg = NULL, maxk = 2)
```

## Arguments

- results:

  Data.table or data.frame. Bootstrap results with subgroup
  characteristics including columns like Pcons, hr_sg, N_sg, K_sg, and
  M.1-M.k

- nb_boots:

  Integer. Total number of bootstrap iterations

- original_sg:

  Character vector. Original subgroup definition from main analysis
  (e.g., c("{age\>=50}", "{nodes\>=3}") for a 2-factor subgroup)

- maxk:

  Integer. Maximum number of factors allowed in subgroup definition

## Value

List with summary components:

- basic_stats:

  Data.table of summary statistics

- consistency_dist:

  Data.table of Pcons distribution by bins

- size_dist:

  Data.table of subgroup size distribution

- factor_freq:

  Data.table of factor frequencies by position

- agreement:

  Data.table of subgroup definition agreement counts

- factor_presence:

  Data.table of base factor presence counts

- factor_presence_specific:

  Data.table of specific factor definitions

- original_agreement:

  Data.table comparing to original analysis subgroup

- n_found:

  Integer. Number of successful iterations

- pct_found:

  Numeric. Percentage of successful iterations
