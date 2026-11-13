## Submission

This is an update from forestsearch 0.4.0 (the current CRAN version) to
0.5.0. NEWS.md lists the changes under "forestsearch 0.5.0".

This release removes `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()`,
which were deprecated in 0.4.0 (the first CRAN release to export them,
where they already warned that they were deprecated).

## Test environments

- Local: Pop!_OS 24.04 LTS, x86_64-pc-linux-gnu, R 4.6.2 (2026-10-29);
  `R CMD check --as-cran` on the tarball built by `R CMD build`.
- win-builder (devel and release): TO BE FILLED BY MAINTAINER
- R-hub: TO BE FILLED BY MAINTAINER

## R CMD check results

0 errors | 0 warnings | 0 notes

The local check above ended with `Status: OK`.

## Reverse dependencies

There are no reverse dependencies: CRAN's package index (read on
2026-11-12) lists no package that depends on, imports, links to or suggests
forestsearch.
