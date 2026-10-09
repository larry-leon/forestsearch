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
* win-builder (devel, release): TBD
* mac-builder: TBD

## R CMD check results

0 errors | 0 warnings | 1 note

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.

  The HTML Tidy installed on the local check machine is too old, so this NOTE
  comes from the environment and not from the package.

## Reverse dependencies

There are no reverse dependencies (`tools::package_dependencies(reverse = TRUE,
which = "all")` against CRAN returns none).
