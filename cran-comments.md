## Submission

This is an update from forestsearch 0.4.0 (the current CRAN version) to
0.5.0. NEWS.md lists the changes under "forestsearch 0.5.0".

This release removes `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()`,
which were deprecated in 0.4.0 (the first CRAN release to export them,
where they already warned that they were deprecated).

## Resubmission

This is a resubmission of 0.5.0. The first submission passed the incoming pre-tests on
Windows and Debian (Status: OK) and was flagged for the overall check time on Windows
(12 min > 10 min). In this version:

- the vignette is now a short introduction; the long analysis it ran is an article on the
  package website and is not part of the package;
- on CRAN a small core set of tests of the main function is run; the full test suite is run
  locally before each release.

The additional issue reported for 0.4.0 under MKL came from tests that required
bit-identical floating-point results. Those tests are no longer run on CRAN, and the
remaining CRAN-run comparison of repeated results uses a numerical tolerance.

## Test environments

- Local: Pop!_OS 24.04 LTS, x86_64-pc-linux-gnu, R 4.6.2 (2026-10-29);
  `R CMD check --as-cran` on the tarball built by `R CMD build`.
- win-builder: TO BE FILLED BY MAINTAINER

## R CMD check results

0 errors | 0 warnings | 1 note

The local check above ended with `Status: OK`.

On CRAN a small core set of tests of the main function `forestsearch()` is run. The full
test suite is run locally before each release.

## Reverse dependencies

There are no reverse dependencies: CRAN's package index (read on
2026-11-12) lists no package that depends on, imports, links to or suggests
forestsearch.
