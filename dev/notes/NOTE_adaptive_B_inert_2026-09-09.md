# NOTE — `--adaptive-B` is inert in `guohe_reproduction_run.R` (2026-09-09)

Raised by: `quarto/GuoHe/REPORT_guohe_supp_stage0_2026-09-09.md` (`33d0add3`), re-verified in
`quarto/GuoHe/REPORT_guohe_supp_stage0_2026-09-09_v2.md`. Disposition set by Larry in
`dev/tasks/claude_cc_task_guohe_supplement_2026-09-09_v2.md` §A6: **record only, no repair.**

## The defect

- `b_adapt` is parsed from the command line at `quarto/GuoHe/guohe_reproduction_run.R:54`:

```r
b_adapt <- as.integer(opt("adaptive-B", as.character(b_boot)))
```

- It is printed to the console banner at `:67`:

```r
            if (adaptive) sprintf("   adaptive B = %d", b_adapt) else ""))
```

- It is written into the saved bundle's metadata at `:127`:

```r
      B = b_boot, adaptive_B = if (adaptive) b_adapt else NA_integer_,
```

- **It never reaches `gh_one_rep()`.** The replicate call at `:109` passes `B = b_boot`:

```r
        beta = s$beta, n = s$n, r_grid = GH_R_GRID, B = b_boot,
```

- `gh_one_rep()` has a single `B` formal (`quarto/GuoHe/guohe_reproduction_sim.R:131`) and
  forwards that one value to **both** the fixed-r Algorithm-3 loop and the adaptive call
  (`guohe_reproduction_sim.R:186`):

```r
      orient = -1, r_grid = r_grid, v = v, B = B, level = level,
```

- Inside `guohe_adaptive_r()` the coupling continues by design: one `B` serves the inner CV fits
  (`R/guohe_adaptive_r.R:255`) and the final refit (`:285`). That coupling is the function's own
  documented behavior and is **not** the defect here.

## Consequence

- The documented production commands, `quarto/GuoHe/guohe_reproduction_RUN.md:57` and `:60`:

```
Rscript guohe_reproduction_run.R --tables=35 --cores=120 --B=2000 --adaptive --adaptive-B=200
Rscript guohe_reproduction_run.R --tables=6  --cores=120 --B=2000 --adaptive --adaptive-B=200
```

  therefore recorded `adaptive_B = 200` in the committed Tables 3–6 bundles while the adaptive
  path **executed at `B = 2000`**.

- **This is a provenance-labelling defect, not a numerical error.** The Adaptive columns in those
  bundles are valid results at `B = 2000`; they are simply labelled 200. The published-comparison
  conclusions drawn from them are unaffected.

- `quarto/GuoHe/guohe_reproduction_RUN.md:98-99`:

```
| Adaptive at $B = 200$ | 184 | ~1.7 h |
| Adaptive at $B = 2000$ | 875 | ~8.1 h |
```

  is therefore a **projection ledger, not a record of what ran** — the runs sat on the 875 core-h
  line, not the 184 one. The same caution applies to the per-replicate figure at `:79`
  (`| $k=2$, $n=400$, $B=200$, adaptive marginal cost | 5.90 |`): it is labelled `B = 200` and may
  have been measured under either setting.

## Bearing on the current task

- `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` §5 pinned T2's `B = 200` as "the
  validated reproduction setting." That premise is false on this reading, and §A5 of the v2
  amendment supersedes it: **T2 runs at `B = 2000`** — the setting under which the function was
  actually validated, resolution parity with the fixed-r columns, and their method's strongest
  configuration.
- Gate 1b must therefore **measure** the adaptive cost rather than inherit the ledger's 5.90 s
  figure, per A5.

## Disposition

- **No repair in this task.** `guohe_reproduction_run.R` is not modified.
- **No re-run of committed work.** The Tables 3–6 bundles stand as they are.
- Whether to fix the flag, re-label the committed metadata, or leave both is **Larry's call**.
