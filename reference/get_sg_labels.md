# Get Subgroup Display Labels for Harm or Benefit Notation

Returns Unicode label components for H/Hc (harm search) or G/Gc (benefit
search) notation. Used internally by `build_classification_table`,
`build_estimation_table`, `interpret_estimation_table`, and
`format_oc_results` to switch between the two notation systems described
in Leon et al. (2024, Sections 2–4.2).

## Usage

``` r
get_sg_labels(notation = c("harm", "benefit"))
```

## Arguments

- notation:

  Character. `"harm"` (default) uses H/Hc notation for
  detrimental-subgroup searches. `"benefit"` uses G/Gc notation (G for
  "good"/"gain") for benefit-subgroup searches.

## Value

A named list with components:

- sg, sg_c:

  True subgroup and complement (e.g., "H", "Hc")

- sg_hat, sg_hat_c:

  Estimated subgroup and complement (Unicode hat)

- plain, plain_c:

  Plain-text labels for metric names

- word:

  Descriptive word: "harm" or "benefit"
