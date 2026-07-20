---
title: "GBSG maxeff + Guo & He application — Claude Code brief"
bibliography: []
---

# Claude Code brief: GBSG maxeff analysis with Guo & He as a fourth correction

## Purpose

Add Guo & He (2021) Algorithm 3 as a fourth de-biasing correction (alongside
Naive, Full Bootstrap, Multiplier Resampling) on the GBSG maxeff-selected
subgroup. This is the "post-hoc identified subgroups" application of Guo & He,
using the full forestsearch machinery. Application only; no simulation, no
truth target. If it looks promising, Table 7 coverage (GH vs MR) is a possible
follow-on -- out of scope here.

## Files

Place in the repo:

```
quarto/GuoHe/fs_to_guohe.R                 # NEW: the FS -> GH adapter (sourced)
quarto/GuoHe/guohe_algorithm3.R            # existing (unchanged)
quarto/GuoHe/analysis_gbsg_cox_maxeff.qmd  # NEW: the cleaned single-track analysis
```

`fs_to_guohe.R` and `guohe_algorithm3.R` are SOURCED, not installed. Nothing in
the `forestsearch` package is modified. Do not modify `guohe_algorithm3.R`.

## What the adapter does

`fs_to_guohe(fit)` reads the maxeff candidate family from
`fit$grp.consistency$out_sg$result`, materializes each candidate as a 0/1
membership column by intersecting its active `names.Z` cut-indicator columns
(the same mechanism `extract_subgroup()` uses), asserts the selected
subgroup's exported column equals `fit$sg.harm.id` exactly, then runs
`guohe_algorithm3()` with `orient = -1`. It returns HR-scale naive / debiased /
lower-bound for the selected subgroup.

## Gates -- run in order, stop on any failure

### Gate 1 -- adapter sources and self-tests (seconds, no package needed)

```bash
Rscript -e 'source("quarto/GuoHe/guohe_algorithm3.R"); source("quarto/GuoHe/fs_to_guohe.R"); cat("sourced OK\n")'
```

### Gate 2 -- the maxeff fit returns an UNTRUNCATED family

The whole approach depends on `sg_focus = "maxeff"` exposing the full candidate
family. The code says it does (`subgroup.consistency()` forces
`stop_Kgroups = Inf` under maxeff), but confirm by execution. After the
`forestsearch-maxeff` chunk runs, check:

```r
fam <- fs$grp.consistency$out_sg$result
nrow(fam)                       # should be the full family (hundreds to ~1744),
                                # NOT capped at a small number like 50
```

If `nrow(fam)` is suspiciously small (e.g. exactly `max_subgroups_search`), the
pool was truncated and the de-biasing would be against a censored family --
STOP and report, do not proceed.

### Gate 3 -- membership integrity (the critical structural check)

`fs_to_guohe(..., verify = TRUE)` (the default) asserts the exported selected
column matches `fit$sg.harm.id`. If it errors with "membership disagrees",
that is a real mismatch between the adapter's reconstruction and forestsearch's
selection -- STOP and report the row count. Do NOT set `verify = FALSE` to get
past it; the check is the guarantee that GH is correcting the right subgroup.

### Gate 4 -- GH argmax coincides with the forestsearch selection

`fs_to_guohe()` warns if Guo & He's oriented argmax differs from forestsearch's
selected subgroup. Under maxeff they should coincide (both are the
effect-maximiser). If the warning fires, report it -- it means the reported GH
correction refers to a different subgroup than FB/MR, and the table would be
comparing different targets.

### Gate 5 -- render

```bash
GH_DIR=quarto/GuoHe quarto render quarto/GuoHe/analysis_gbsg_cox_maxeff.qmd
```

## Known-fragile fields to verify against the real summary schemas

I wired the comparison table from the multimethod document's field names, but
could not run them (no forestsearch in the authoring sandbox). Two spots need a
schema check against the ACTUAL objects on the analysis machine, and are the
most likely place a first render fails:

1. **FB summary fields.** The table reads `fb_sum$est_H`, `fb_sum$lo_H`,
   `fb_sum$hi_H`. Confirm `summarize_bootstrap_results()` returns those names
   for the selected (H) subgroup HR and CI. If the schema differs (e.g. the
   values live in `fb_sum$table` rows rather than top-level fields), adjust the
   `table` chunk to read the correct path. This is a display fix only -- it does
   not touch the estimates.

2. **MR gate fields.** The table reads `g$Hc_estimates$est/lo/hi` from
   `fs$debias_gate`. Confirm the gate's estimate object carries the
   selected-subgroup HR and interval under those names; adjust if not.

Fixing these two is in scope (they are display plumbing in the new qmd). Report
what the real schemas are so the fix is grounded in the actual objects, not
guessed. Everything else -- the fit, FB, MR, and GH computations -- is
unmodified package/adapter behaviour.

## Parallelism (already wired; stated so nothing is second-guessed)

- FB: `forestsearch_bootstrap_dofuture()` uses `%dofuture%`, inherits the plan,
  restores it on exit. Any plan works.
- MR: vectorized matrix algebra, no plan needed; benefits from threaded BLAS.
- GH: `guohe_algorithm3()`'s pair bootstrap parallelises via its own `parallel`
  argument. Because GH is SOURCED not installed, use forked workers
  (`plan(multicore)`, set automatically via `GH_PARALLEL` on Linux) so workers
  inherit the sourced function; multisession workers would not see it. On the
  Linux box multicore is available, so `GH_PARALLEL = TRUE`.

## After the render

Report the four-row comparison table (Naive / FB / MR / GH) and the GH family
size. If the GH debiased HR and lower bound are sensible relative to FB/MR
(same selected subgroup, correction toward the null), the application is
promising and Table 7 coverage becomes worth discussing. Do not commit yet --
send the rendered result back for review first.
```
