# Summarize Bootstrap Subgroup Analysis Results

Comprehensive tabulation of bootstrap subgroup identification results
including basic statistics, factor frequencies, consistency
distributions, GRF-cut frequencies, and agreement with the original
analysis subgroup. Returns raw `data.table` tabulations suitable for
custom post-processing, export, or plotting; for the full formatted
pipeline (gt tables, diagnostics, plots) call
[`summarize_bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
instead.

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

- grf_cut_freq:

  Data.table of GRF policy-tree cut frequencies across all bootstrap
  iterations (not just successful ones). Columns: `Rank`, `grf_cut`,
  `N`, `Percent`. `NULL` when `results` lacks the `grf_cuts_b` column
  (pre-Phase-C objects). Unlike subgroup-related tables, this tabulation
  is NOT restricted to iterations that identified a subgroup – it spans
  all bootstraps so that the analyst can distinguish "GRF produced no
  cut" from "GRF produced a cut that consistency rejected".

- original_agreement:

  Data.table comparing to original analysis subgroup

- n_found:

  Integer. Number of successful iterations

- pct_found:

  Numeric. Percentage of successful iterations

## See also

[`summarize_bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
for the full formatted summary pipeline (gt tables plus diagnostics and
plots).
[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
to run the bootstrap analysis that produces the `results` object
consumed here.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run bootstrap analysis
boot_out <- forestsearch_bootstrap_dofuture(fs.est, nb_boots = 500)

# Raw tabulations for custom post-processing
subg <- summarize_bootstrap_subgroups(
  results     = boot_out$results,
  nb_boots    = 500,
  original_sg = fs.est$sg.harm,
  maxk        = 2L
)

# Inspect the GRF cut frequencies across all bootstraps
subg$grf_cut_freq

# Export factor frequencies to CSV for an external report
data.table::fwrite(subg$factor_freq, "factor_freq.csv")
} # }
```
