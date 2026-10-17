# Per-Subject Membership Matrix for Pareto Frontier Members

For each candidate subgroup on
`fs$grp.consistency$out_sg$pareto_frontier`, returns a binary indicator
(one row per subject in `fs$df.est`, one column per frontier member)
marking which subjects belong to that subgroup.

## Usage

``` r
frontier_member_flags(fs)
```

## Arguments

- fs:

  A `forestsearch` result object.

## Value

A list with three elements:

- `flags`:

  Integer matrix of dimension `n_subjects` x `n_frontier_members`, with
  1 = subject in subgroup, 0 = not.

- `member_labels`:

  Character vector of length `n_frontier_members`, each entry the
  "`cut1 & cut2 & ...`" definition of that frontier member.

- `is_selected`:

  Logical vector of length `n_frontier_members`, marking which member is
  the selected subgroup.

## Details

This lets the user compute downstream quantities the frontier table does
not surface directly: per-frontier-member sample sizes by covariate
stratum, set comparisons between frontier members ("which subjects are
in the selected subgroup but not on the runner-up?"), and survival
summaries by frontier member.

## See also

[`pareto_frontier_table`](https://larry-leon.github.io/forestsearch/reference/pareto_frontier_table.md),
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md),
[`plot_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_frontier.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fmf <- frontier_member_flags(fs)
# How many subjects in each frontier member?
colSums(fmf$flags)
# Which subjects are unique to the selected subgroup?
sel <- which(fmf$is_selected)
rowSums(fmf$flags) == 1L & fmf$flags[, sel] == 1L
} # }
```
