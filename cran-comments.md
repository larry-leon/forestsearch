## R CMD check results

0 errors | 0 warnings | 1 note

* New submission
* Possibly misspelled words in DESCRIPTION:
  GRF, MRCT, et, al, combinatorial
  These are standard acronyms (Generalized Random Forests,
  Multi-Regional Clinical Trials) and citation fragments.

## Test environments

* Local: macOS Tahoe 26.3.1, R 4.5.2 (aarch64-apple-darwin20)
* win-builder: R-release (R 4.5.3, x86_64-w64-mingw32)

## Downstream dependencies

This is a new submission. There are no reverse dependencies.

## Resubmission

This is a resubmission addressing all reviewer feedback from the
initial submission, including:

* Removed trailing spaces in DESCRIPTION
* Added @return documentation to all exported functions
* Removed examples from unexported/internal functions
* Replaced all \dontrun{} with \donttest{} or removed
* Ensured output_dir defaults to tempdir()
* Added on.exit(par(oldpar)) for all par() changes in functions
* Restored par() in vignette code chunks
* Replaced options(warn = -1) with suppressWarnings()
* Reduced exported functions from 118 to 53
