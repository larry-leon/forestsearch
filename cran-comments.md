## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps: unable to verify current time
  This is a local network issue (time server unreachable during check)
  and is not related to the package.

## Test environments

* Local: macOS Tahoe 26.3.1, R 4.5.2 (aarch64-apple-darwin20)
* win-builder: R-devel
* win-builder: R-release

## Downstream dependencies

This is a new submission. There are no reverse dependencies.

## Resubmission

This is a resubmission addressing all reviewer feedback from the
initial submission, including:

* Removed trailing spaces in DESCRIPTION
* Added @return documentation to all exported functions
* Removed examples from unexported/internal functions
* Replaced all \\dontrun{} with \\donttest{} or removed
* Ensured output_dir defaults to tempdir()
* Added on.exit(par(oldpar)) for all par() changes in functions
* Restored par() in vignette code chunks
* Replaced options(warn = -1) with suppressWarnings()
* Reduced exported functions from 118 to 53
