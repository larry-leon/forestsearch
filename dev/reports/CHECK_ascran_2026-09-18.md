# CHECK — one-off `R CMD check --as-cran` certification run

**One-off, authorized by Larry on 2026-09-18 for this occasion only.** The standing `CLAUDE.md` check
policy is unchanged and was not edited: per task, the task's own acceptance tests and nothing else; the
full suite and any CRAN check only when Larry explicitly asks.

Queued behind Directive C, which completed fully green with no gate stopping it (`4447a6c6`).

## Surface and environment

| | |
|---|---|
| Surface | `rcmdcheck::rcmdcheck(path = ".", args = "--as-cran", error_on = "never")` — the certification surface, PDF manual and vignettes **included** |
| Tree | `4447a6c6` (Directive C complete), branch `feature/glm-extension`, nothing pushed |
| Pandoc | `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64`, also prepended to `PATH`, per the `CLAUDE.md` environment note |
| Wall clock | **10.3 min** (hard timeout 90 min, never approached) |
| R / platform | as recorded in the saved result object |

Phases that ran and passed, confirming this was the full certification surface and not a reduced one:
`checking examples ... OK` (28 s), `checking examples with --run-donttest ... OK` (211 s),
`checking package vignettes ... OK`, `checking re-building of vignette outputs ... OK` (50 s),
**`checking PDF version of manual ... OK`**.

## Counts

```
this run:   1 error  | 2 warnings | 2 notes
reference:  0 errors | 1 warning  | 2 notes
```

Reference set: `quarto/simulations/actg175/binary_020/REPORT_threshold_docs_2026-09-18.md` §6, recorded
pre- and post-change at `216f3405` / `458cef48` and identical in text between those two runs.

## Finding diff

### Carried over, unchanged (3)

| Kind | Finding | Status |
|---|---|---|
| WARNING | `checking code files for non-ASCII characters` — `R/fs_bias_coverage.R` | identical to reference |
| NOTE | `checking R code for possible problems` — 9 "no visible binding" in `fs_plot_bias_coverage` (`b`, `cell`, `cov`, `cov1`, `cov2`, `estimator`, `obs`, `r`, `ref`) | identical to reference |
| NOTE | `checking HTML version of manual` — validation skipped, no `tidy` on this machine (environmental) | identical to reference |

### (a) NEW findings (2) — both name files changed by these tasks

**NEW-1 — ERROR, `checking tests`** [281s/228s]. `FAIL 14 | WARN 22 | SKIP 66 | PASS 5097`.
**Flagged: names files added by the sync, A and C tasks.** All 14 failures sit in the three probe-based
test files this workstream added, and nothing else in the suite fails:

| File | Failures | Added by |
|---|---|---|
| `tests/testthat/test-threshold-sync.R` | 8 (`:128`, `:129`, `:144`, `:167`, `:168`, `:176`, `:177`, `:178`) | threshold sync (`a611928e`) |
| `tests/testthat/test-threshold-pair-directive-a.R` | 5 (`:98`, `:100`, `:106`, `:108`, `:139`) | Directive A (`4e1b6dd4`) |
| `tests/testthat/test-directive-c.R` | 1 (`:662`, an Error) | Directive C (`4447a6c6`) |

**Two distinct causes, both in the test harness, neither in the package.**

*Cause 1 — a relative path that only exists in the working tree (1 of the 14, and it is mine).*
`test-directive-c.R:662` reads `"../../dev/reports/baseline_directive_C_2026-09-18.csv"` for the Gate C
comparison. Under `devtools::load_all()` + `test_file()` the working directory is `tests/testthat/`, so it
resolves; under `R CMD check` the tests run from `<pkg>.Rcheck/tests/` and `dev/` is not in the built
tarball at all, so `file()` cannot open it: `Error in file(file, "rt"): cannot open the connection`.
Introduced by me in `4447a6c6`. The fix is one line — guard with
`skip_if_not(file.exists(...))` — **not applied**, because this run is report-only.

*Cause 2 — the statement-lifting probe does not survive installation (13 of the 14).*
All thirteen run through `.probe_resolver()` in `tests/testthat/helper-threshold-sync.R`, which rebuilds
`forestsearch()`'s threshold resolution by lifting statements out of `body(forestsearch)` and matching
their deparsed text with `fixed = TRUE`. Against the installed, byte-compiled package that matching does
not hold: **every one of the 175 cells errors** (`sum(!is.na(p$parent_error))` is 175, expected 18; every
resolution column comes back `NA`; `sum(b$violation)` is 0, expected 15). The consequence worth recording
is that **those three files contribute no coverage under `R CMD check`** — they pass only on the
`load_all()` surface they were written and gated on.

*Cross-check, run after the check:* all five task-owned files are green under `devtools::load_all()` at
this same tree — `test-threshold-sync.R` 16, `test-threshold-pair-directive-a.R` 86,
`test-binary-default-or-entry-points.R` 25, `test-binary-default-or.R` 66, `test-directive-c.R` 182;
**0 failures, 0 errors**. So the divergence is between the two surfaces, not a regression in resolved
behaviour. No failure names an `R/` file, and no failure asserts a wrong threshold, estimand, guard,
warning or display value.

**NEW-2 — WARNING, `checking Rd cross-references`.** **Flagged: names a file changed by Directive B.**

```
Missing link(s) in Rd file 'evaluate_combination_with_status.Rd':
  'subgroup_search'
Found the following Rd file(s) with Rd \link{} targets missing package anchors:
  evaluate_combination_with_status.Rd: subgroup_search
```

`man/evaluate_combination_with_status.Rd` is one of the three `.Rd` files Directive B regenerated
(`f391b816`, after the `7e4e4c87` `hr.threshold` scale rider touched `R/subgroup_search.R`'s roxygen). The
link target is `subgroup_search`, but the function's Rd is `subgroup.search.Rd` — the alias carries the
dot, the link does not. Not present in the reference set. Pre-existing in the source roxygen or introduced
by the rider is not established here; **nothing changed**, this is a report.

### Removed relative to reference

None. Every reference finding is still present.

## (b) Vignette verdict — the B session's open question, answered

**The vignettes built.** `* creating vignettes ... OK` at build time, then
`* checking package vignettes ... OK` and `* checking re-building of vignette outputs ... OK` [50s/53s].

This settles `REPORT_binary_default_or_2026-09-18.md` **F4** ("vignette build fails on an untouched tree"):
the cause was the **`pandoc` PATH**, not stale or broken vignette content. With
`RSTUDIO_PANDOC` exported, `vignettes/forestsearch.Rmd` builds and re-builds cleanly, and the PDF manual
builds too. F4's own text already suspected this; it is now confirmed on a tree four tasks later. The
`CLAUDE.md` environment note is correct as written and is what made this run succeed.

## Bottom line

No new finding implicates package behaviour. The certification surface reports the same one warning and
two notes it reported at the docs task, plus one Rd cross-reference warning and one test ERROR whose 14
failures are entirely test-harness artifacts of this workstream's own probe files — one of them a path bug
I introduced in `4447a6c6`. **No fix applied; this run was report-only, and the standing check policy in
`CLAUDE.md` is untouched.**
