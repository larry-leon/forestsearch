## Submission

This is an update from forestsearch 0.1.0 (the current CRAN version) to
0.4.0. The intermediate version numbers recorded in NEWS.md (0.2.0 to 0.3.5)
were development versions only and were never submitted to CRAN.

The update changes the API. NEWS.md lists the changes under "Changes since the
last CRAN release (0.1.0)": the `showten_subgroups` argument was removed from
`forestsearch()` and `subgroup.consistency()`, several defaults changed, and
`save_leff_calibration()` now requires `output_dir`.

## Test environments

* Local: macOS 27.0.1, R 4.5.2 (aarch64-apple-darwin20), checked with
  `R CMD check --as-cran` on the built tarball (vignettes built), with
  `_R_CHECK_LIMIT_CORES_=TRUE` and CRAN weightedsurv 0.1.0
* win-builder, first round: R-release 4.6.1 and R-devel (r90650)
* win-builder, second round (after the fixes below): R-devel, pending
* mac-builder: unavailable (HTTP 502 on two attempts, 2026-10-09)

## R CMD check results

0 errors | 0 warnings | 0 notes (expected after the fixes)

* The first win-builder round (R-release and R-devel) gave one NOTE, from
  "checking HTML version of manual": an HTML validation problem in
  `inject_mrct_structure.Rd`, where a literal `<region>` was read as an HTML
  tag. This version fixes it.
* Local check status: `Status: 1 NOTE`. The NOTE is
  "checking HTML version of manual ... NOTE / Skipping checking HTML
  validation: 'tidy' doesn't look like recent enough HTML Tidy." It comes from
  the old HTML Tidy on the local check machine, not from the package.

Long-running tests are skipped on CRAN with `skip_on_cran()`; the full suite
runs locally.

## Reverse dependencies

There are no reverse dependencies (`tools::package_dependencies(reverse = TRUE,
which = "all")` against CRAN returns none).
