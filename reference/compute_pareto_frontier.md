# Compute Pareto Frontier on (Effect, N)

Returns the subset of candidate subgroups that are not dominated on the
two-objective space (effect size, sample size), where both objectives
are maximized. A candidate \\i\\ is dominated iff another candidate
\\j\\ has \\hr_j \ge hr_i\\ and \\N_j \ge N_i\\, with at least one
inequality strict.

## Usage

``` r
compute_pareto_frontier(result_dt, effect_log_scale = FALSE)
```

## Arguments

- result_dt:

  A data.table with numeric columns `hr` and `N`. Typically the
  post-consistency result table.

- effect_log_scale:

  Logical. If `TRUE`, exponentiate `hr` before comparing (so the
  dominance test is on the natural scale for ratio measures). Default
  `FALSE`.

## Value

A data.table of non-dominated rows, sorted by `hr` descending. Returns
an empty 0-row data.table if input is empty.

## Details

The frontier is intended as a **post-hoc reporting** artifact, not a
selection criterion. The selected subgroup (under `sg_focus`) may or may
not appear on the frontier – in particular, `"hrMinSG"` may select an
N-dominated point by design (preferring small subgroups).
