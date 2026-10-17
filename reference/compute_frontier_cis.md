# Compute Naive and Split-Derived 95% CIs for Pareto Frontier Members

For each candidate subgroup on
`fs$grp.consistency$out_sg$pareto_frontier`, computes three 95\\

## Usage

``` r
compute_frontier_cis(
  fs,
  n_splits = 1000L,
  ci_level = 0.95,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- fs:

  A `forestsearch` result object.

- n_splits:

  Integer. Number of 50/50 random splits per frontier member for the
  Split CI. Default `1000L`.

- ci_level:

  Numeric in `(0, 1)`. Default `0.95`.

- seed:

  Integer or `NULL`. Seed for reproducibility of the per-split halving.
  Default `NULL` (no explicit seed).

- verbose:

  Logical. Default `FALSE`. When `TRUE`, prints per-frontier-member
  progress.

## Value

A `data.table` keyed by frontier row index (`m`) with columns:

- `m`:

  Candidate id (matches `pareto_frontier$m`).

- `estimate`:

  Refit effect estimate on the natural scale.

- `naive_lcl`, `naive_ucl`:

  Naive Wald CI bounds.

- `split_sd`:

  Half-jackknife SE estimate: empirical SD of the \\2S\\ per-half effect
  estimates on the SE scale (log for ratio measures, natural for
  differences).

- `split_lcl`, `split_ucl`:

  Split-derived CI bounds.

- `fsbc_estimate`:

  Bias-corrected effect estimate on the natural scale (\\2\hat\beta -
  \overline{\hat\beta^{(h)}}\\ on the SE scale, then back-transformed).

- `fsbc_lcl`, `fsbc_ucl`:

  FSBC-mimic CI bounds.

- `n_valid_splits`:

  Number of splits that produced valid per-half estimates.

Rows where the refit fails (e.g., zero events in one arm) have `NA` in
all CI columns.

## Details

- **Naive CI**: full-sample Wald CI from a Cox or GLM model refit on the
  subgroup's data. This is the textbook CI that ignores the
  multiple-testing / subgroup-search context. Anti-conservative by
  construction.

- **Split CI**: a 95\\ estimate, with SE estimated by the empirical
  standard deviation of the \\2S\\ individual half-sample effect
  estimates produced by `n_splits` random 50/50 splits. This is a
  half-jackknife SE (Shao1996). Naming: the \\\sim\\ in column labels
  marks this as a resampling-derived approximation, not a model-based
  CI.

- **FSBC-mimic CI**: a bias-corrected interval inspired by the bootstrap
  bias-correction algorithm of Leon2024fs (eq 7), with the
  selection-on-bootstrap-data term \\\eta_b^\*(\hat H_b^\*)\\ zeroed out
  (the selected subgroup is treated as fixed across half-jackknife
  replicates). The bias-corrected estimate is
  \\\hat\beta^{\mathrm{FSBC}} = 2\hat\beta -
  \overline{\hat\beta^{(h)}}\\ where \\\overline{\hat\beta^{(h)}}\\ is
  the mean of the \\2S\\ half-sample estimates; the SE is the same
  half-jackknife SE as for the Split CI. See “Details — FSBC-mimic
  interpretation” below.

All three CIs are computed post-hoc on the returned `forestsearch`
object; this function does not modify any internal state.

## Details — Half-jackknife SE

For each of `n_splits` random 50/50 splits we fit the effect model
separately on each half, obtaining \\2S\\ per-half estimates
\\\\\hat\beta^{(s,h)} : s=1,\ldots,S; h=1,2\\\\. The half-jackknife SE
estimator is \$\$\widehat{\mathrm{SE}} =
\mathrm{sd}\bigl\\\hat\beta^{(s,h)}\bigr\\,\$\$ i.e., the plain
empirical SD of the \\2S\\ estimates. No \\\sqrt{2}\\ correction is
applied: simulation studies confirm that this estimator's CI achieves
nominal coverage in a no-selection setting (see the package's
pressure-test script). Earlier versions of this function used the
average \\(\hat\beta^{(s,1)} + \hat\beta^{(s,2)})/2\\ per split; that
approach is invalid because the two halves are complementary and their
average is nearly constant.

## Details — FSBC-mimic interpretation

The FSBC-mimic CI is a *quick approximation* to the bootstrap
bias-corrected estimator of Leon2024fs. Substituting the \\2S\\ per-half
estimates for bootstrap replicates, and treating the selected subgroup
\\\hat H\\ as fixed (i.e., omitting the \\\eta_b^\*(\hat H_b^\*) =
\hat\beta_b^\*(\hat H_b^\*) - \hat\beta(\hat H_b^\*)\\ term that
captures variability in FS selection on the bootstrap data), the
bias-corrected estimator from paper eq (7) collapses to
\$\$\hat\beta^{\mathrm{FSBC}}(\hat H) = \hat\beta(\hat H) -
\frac{1}{2S}\sum\_{s,h}\bigl\[\hat\beta^{(s,h)} - \hat\beta(\hat
H)\bigr\] = 2\hat\beta(\hat H) - \overline{\hat\beta^{(h)}}.\$\$ The SE
is the same half-jackknife SE as for the Split CI.

The FSBC-mimic CI is *not* a full FSBC CI: it omits the
selection-on-bootstrap-data variance. Use this column as a quick
bias-corrected diagnostic, not for inference. For full FSBC CIs use
[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md).

## See also

[`pareto_frontier_table`](https://larry-leon.github.io/forestsearch/reference/pareto_frontier_table.md),
[`frontier_member_flags`](https://larry-leon.github.io/forestsearch/reference/frontier_member_flags.md),
[`plot_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_frontier.md).
