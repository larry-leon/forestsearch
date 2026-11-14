# TASK — MRCT articles: general installation note (2026-11-14)

Session name: `fs-mrct-install-note`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**. Larry asked for
this on 2026-11-14. Read this whole file first.

## Why

Both MRCT articles open with a callout, "Installing forestsearch", written while the CRAN release
(0.4.0) lacked the functions they use. It tells every reader to install the development version
from GitHub. That will be out of date once 0.5.0 is on CRAN. The replacement below is general: it
states the version the article needs, installs from CRAN, and keeps a note on installing from
GitHub.

## Files

- `vignettes/articles/mrct-consistency.qmd`
- `vignettes/articles/mrct-consistency-simulations.qmd`

## Rules

- **Text only, in these two files.** Two edits per file, both given below. Nothing else changes
  in them, and no other file changes apart from this task file and the report.
- No render, no `R CMD check`, no test suite. No agents. One short report.
- Work on `master`. Commits are path-limited. **No push, pull, merge or tag.**

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate (the only stop):** branch `master`; `git status --porcelain` is empty; in **each** of
   the two files, each of these strings occurs exactly once:
   - `The CRAN release of forestsearch (0.4.0) does not yet contain the functions this article uses`
   - `This article needs the development version of forestsearch (see 'Installing forestsearch'); `

## Step 1 — edit 1: the callout

In each file, replace the whole callout headed `## Installing forestsearch` — from its opening
line `::: {.callout-note}` to its closing `:::` — with exactly the text between the two marker
lines below (the marker lines themselves are not part of it). The text is the same for both
files.

~~~~text
<<<BEGIN
::: {.callout-note}
## Installing forestsearch

This article needs forestsearch 0.5.0 or later. Run the steps in a fresh R session.

```r
# 1. a reliable CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 2. forestsearch from CRAN, which also installs every package it needs
install.packages("forestsearch")

# 3. two packages forestsearch only suggests, which this article uses (speff2trial for the
#    ACTG 175 data, latentcor for the copula expansion), and knitr for the tables
install.packages(c("speff2trial", "latentcor", "knitr"))

packageVersion("forestsearch")   # this article needs 0.5.0 or later
```

If a step asks whether to install packages from sources, answer no. The chunk below stops with a
message naming anything that is missing, before any computation starts.

**Installing from GitHub.** The GitHub repository holds the version under development, which can be
ahead of the CRAN release. Use it if you want changes that have not reached CRAN yet, or if the
chunk below reports that your installed version is too old and CRAN has nothing newer. Do steps 1 to
3 first, so that every package forestsearch needs is already installed, then:

```r
install.packages("remotes")
options(timeout = 600)   # R's default download limit of 60 seconds is too short for this repository
remotes::install_github("larry-leon/forestsearch", upgrade = "never")
packageVersion("forestsearch")
```

`upgrade = "never"` leaves the packages already installed as they are. The call downloads the whole
forestsearch repository, several gigabytes, and builds only the package from it, so it takes several
minutes. Restart R afterwards, so that the new version is the one loaded.
:::
END>>>
~~~~

## Step 2 — edit 2: the message of the `install-check` chunk

In each file, in the `install-check` chunk, replace

~~~~text
  stop("This article needs the development version of forestsearch (see 'Installing forestsearch'); ",
~~~~

with

~~~~text
  stop("This article needs forestsearch 0.5.0 or later (see 'Installing forestsearch'); ",
~~~~

Nothing else in the chunk changes.

## Step 3 — checks

1. `git diff --stat` lists only the two files. In each, the changed lines are the callout block
   and the one `stop(` line; quote the diff's hunk headers in the report.
2. In each file, the two gate strings of Step 0 no longer occur, and `forestsearch 0.5.0 or later`
   occurs three times (twice in the callout, once in the chunk).
3. The R code of each file still parses: extract it to scratch (`knitr::purl()` into `/tmp`) and
   `parse()` it. Nothing is rendered.
4. The CRAN tarball is unaffected: quote the `.Rbuildignore` line that excludes
   `vignettes/articles`.
5. **List only, change nothing:** every other file under `vignettes/` that mentions the
   development version, a GitHub install, or a stated CRAN release number in an installation
   note, with the file, the line and the sentence.

## Step 4 — report and close

1. One commit for the two files.
2. Write `dev/reports/REPORT_mrct_articles_install_note_2026-11-14.md`: first line DONE or NOT
   DONE; the Step 3 results; the Step 3.5 list; this session's commits. Commit it and copy it to
   `~/Downloads`.
3. **Final checks** (report the result): this session's commits touch only the two article files,
   this task file and the report; `git status --porcelain` is empty; nothing was pushed, merged
   or tagged.
4. Print the first line and the Step 3.5 list as the closing message.
