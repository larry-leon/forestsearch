# Covariate signatures of the realized rules

Collapses realized rules to the **set of covariates** they name,
discarding the cut values. A rules-by-frequency table lists every
distinct `sg_def` string and so is dominated by cut noise – a
500-replicate survival run produced 583 rows, in which
`{age <= 46} & {nodes <= 5}` and `{nodes <= 3} & {age <= 46}` are
separate entries despite being the same pairing. This reports how often
each covariate *combination* is selected, which is the structural
question.

## Usage

``` r
fs_rule_covariate_pairs(x, true_covariates = NULL, detected_only = TRUE)
```

## Arguments

- x:

  A results data frame carrying `sg_def`, a bundle with `results`, a
  path to an `.rds`, or a bare character vector of rules.

- true_covariates:

  Character vector. The true-region covariates, used to classify each
  signature as `"exact"` (the signature is exactly this set),
  `"partial"` (overlaps it) or `"none"`. `NULL` skips the
  classification.

- detected_only:

  Logical. Restrict to rows with `detected == 1` when a results frame is
  supplied.

## Value

A data frame ordered by frequency: `covariates`, `k`, `n`, `share`,
`n_rules`, and `match` when `true_covariates` is given.

## What a signature is

The covariate names in the rule, de-duplicated and sorted, joined by
`" & "`. Sorting is what makes the two orderings above collapse
together. Two consequences worth knowing:

- A rule naming one covariate twice — a band such as
  `!{nodes <= 3} & {nodes <= 5}` — has signature `"nodes"`, of size 1.
  The `k` column gives the number of distinct covariates, and `n_rules`
  the number of distinct rule strings that collapsed into the row, so
  bands are visible rather than hidden.

- Logically equivalent spellings of the same condition (`!{meno}` versus
  `{meno == 0}`) collapse together, which a raw rule table also splits.
