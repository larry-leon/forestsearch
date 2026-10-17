# Identification structure for one set of replicates

Classifies every detected replicate's realized rule against a planted
`anchor` / `partner` / `proxy` triple, and returns covariate
involvement, the structure composition, and classification accuracy.

## Usage

``` r
fs_identification_structure(
  results,
  anchor,
  partner,
  proxy,
  confounders = NULL,
  sg_def_col = "sg_def",
  covs_col = "covs"
)
```

## Arguments

- results:

  A per-replicate results data frame, a bundle carrying `results`, or a
  path to an `.rds` holding either.

- anchor:

  Character. The covariate the search is expected to recover.

- partner:

  Character. The true second covariate of the planted region.

- proxy:

  Character. A correlated covariate that may substitute for `partner`.

- confounders:

  Character vector. The analysis covariate pool, for the involvement
  vector. Defaults to the three named roles.

- sg_def_col, covs_col:

  Character. Columns holding the realized rule and its projected
  covariate names.

## Value

A list with `det_rate`, `n_detected`, `n_total`, `involvement`,
`structure` (a named share vector), and the four accuracy means.

## The classification

Each detected replicate falls in exactly one of five mutually exclusive
categories, tested in this order: the anchor is absent
(`anchor missed`); the rule names the anchor and nothing else
(`anchor only`); the rule pairs the anchor with the true `partner`; the
rule pairs it with the `proxy`; or the rule pairs it with something
else. The order matters – a rule naming both the partner and the proxy
counts as the true pairing.
