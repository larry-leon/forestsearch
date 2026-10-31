## Submission

This is an update from forestsearch 0.4.0 (the current CRAN version) to
0.4.1. NEWS.md lists the changes under "forestsearch 0.4.1":

* Three arguments of `forestsearch()` are deprecated (`vi.grf.min`,
  `max_n_confounders`, `max.minutes`): the formals are kept at the same
  positions and defaults, a non-default value warns once and has no effect.
* Two per-candidate fits call the fitting engines directly
  (`survival::coxph.fit()`, `stats::glm.fit()`); results are bit-identical
  to the formula-interface fits and are covered by identity tests. The
  package now requires survival (>= 3.5-8).
* Bug fixes: the log axis of `plot_subgroup_results_forestplot()` (#7); the
  weighted fit for `ps_adjust_method = "iptw"` with a continuous outcome.
* Documentation: a new getting-started vignette reproducing a published
  analysis; the longer analyses are articles on the package website, so the
  package ships one vignette instead of two; the `?forestsearch` examples
  are corrected (one called `mr_estimates_table()` without the bootstrap
  result it requires).

## Test environments

* Local: Pop!_OS 24.04 (x86_64-pc-linux-gnu), R 4.6.1, `R CMD check --as-cran
  --run-donttest` on the built tarball (vignette built): Status OK
  (0 errors | 0 warnings | 0 notes)
* win-builder R-release: (Larry fills in)
* win-builder R-devel: (Larry fills in)

## R CMD check results

0 errors | 0 warnings | 0 notes (local, Status: OK).

Check time: the vignette rebuild takes about 75 s locally; the examples
about 40 s with `--run-donttest`; the test suite runs its long tests only
with `NOT_CRAN=true` (`skip_on_cran()`), so the tests take under a minute on
CRAN. CRAN's checks of 0.4.0 shipped two vignettes; 0.4.1 ships one, which
reproduces a published analysis at its settings.

## Reverse dependencies

There are no reverse dependencies (`tools::package_dependencies(reverse =
TRUE, which = "all")` against CRAN returns none).
