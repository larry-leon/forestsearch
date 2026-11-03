# TASK: MRCT article — tested install sequence, rendering note, platform tolerance (forestsearch, 2026-11-02)

**Repository:** forestsearch. **Branch:** create `article-install-r2` from `master` (at `8323c5d` or later; stop if `master` does not contain `8323c5d`). **Touches:** `vignettes/articles/mrct-consistency.qmd` and its tracked `.html`. Nothing under `R/`.

## 0. First action

Copy this file to `dev/tasks/TASK_fs-mrct-article-install-r2_2026-11-02.md` and commit it alone ("task: MRCT article install r2").

## Background (for the record)

On 2026-11-02 Larry ran a clean install on a Mac Studio (R 4.5.2, aarch64-apple-darwin20, Apple Accelerate BLAS) and rendered the article from a sparse clone. Findings:

- The install works in this order: mirror and timeout set; CRAN release; `speff2trial`, `latentcor`, `remotes`; `remotes::install_github(..., upgrade = "never")`. The default 60-second download limit made `install_github()` fail (the repository archive is several GB), and a slow CRAN mirror timed out on dependencies.
- The render stopped at `find-check`: the redraw of replicate 5 equals `rep5_trial.csv` in every integer and 0/1 column, but `y_sim` (425 values) and `z_wtkg` (357 values) differ by at most 4.7e-15 relative (largest absolute differences 1.2e-13 and 4.8e-13). With that comparison made to a tolerance, every other check in the article passed on the Mac.

## 1. Edits (verbatim)

Find each anchor by the quoted text (line wrapping may differ). If an anchor is not found exactly once, stop and say so in the notes.

### E1 — the install callout

Replace the whole callout that begins `::: {.callout-note}` / `## Installing forestsearch` (from that line through its closing `:::`, before the `install-check` chunk) with:

````markdown
::: {.callout-note}
## Installing forestsearch

The CRAN release of forestsearch (0.4.0) does not yet contain the functions this article uses, so the
last step below installs the development version from GitHub. Run the steps in a fresh R session.

```r
# 1. a reliable CRAN mirror, and a longer download limit than R's default of 60 seconds
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)

# 2. the CRAN release, which installs every package forestsearch needs
install.packages("forestsearch")

# 3. two packages forestsearch only suggests, which this article uses (speff2trial for the
#    ACTG 175 data, latentcor for the copula expansion), and remotes for the next step
install.packages(c("speff2trial", "latentcor", "remotes"))

# 4. the development version from GitHub; upgrade = "never" leaves the packages above as they are
remotes::install_github("larry-leon/forestsearch", upgrade = "never")
packageVersion("forestsearch")
```

Step 4 downloads the whole forestsearch repository, several gigabytes, and builds only the package
from it, so it takes several minutes; the longer download limit of step 1 is what lets it finish. If a
step asks whether to install packages from sources, answer no. The chunk below stops with a message
naming anything that is missing, before any computation starts.
:::
````

The `install-check` chunk stays as it is.

### E2 — rendering note

Insert directly after the closing fence of the `install-check` chunk (before `# Overview`):

````markdown
::: {.callout-tip collapse="true"}
## Rendering this article yourself

Rendering needs the Quarto command-line tools and about five minutes. To render the article outside any
working copy of the repository, fetch only the article and its data folder (a few megabytes) into a
scratch folder and render it there:

```bash
git clone --depth 1 --filter=blob:none --sparse https://github.com/larry-leon/forestsearch.git ~/scratch/fs-article
cd ~/scratch/fs-article && git sparse-checkout set vignettes/articles
cd vignettes/articles && quarto render mrct-consistency.qmd
```

The rendered `mrct-consistency.html` is written beside the `.qmd`. Delete the scratch folder when done.
:::
````

### E3 — `find-check`: a tolerance in place of exact equality

**E3a.** In the paragraph introducing step 3 of "Finding a case study", replace
`The redrawn trial must equal the shipped `rep5_trial.csv` in every shipped column, and`
with
`The redrawn trial must equal the shipped `rep5_trial.csv` in every shipped column, to floating-point precision (see the note below the check), and`

**E3b.** In the `find-check` chunk, replace
`all(vapply(names(trial), function(v) all(as.numeric(cs$trial[[v]]) == as.numeric(trial[[v]])), NA)),`
with
`all(vapply(names(trial), function(v) isTRUE(all.equal(as.numeric(cs$trial[[v]]), as.numeric(trial[[v]]), tolerance = 1e-12)), NA)),`

and add, as the first line of the chunk body after its comment line, the measured difference (used in E3c):

```r
redraw_rel <- max(vapply(names(trial), function(v) {
  a <- as.numeric(cs$trial[[v]]); b <- as.numeric(trial[[v]])
  max(abs(a - b) / pmax(abs(b), .Machine$double.eps))
}, 0))
```

(It needs `cs` and `trial`, both defined earlier; if the chunk's first statement must stay the `stopifnot()`, put these lines immediately before it.)

**E3c.** In the paragraph after the chunk, replace
`**Passed**: the redrawn trial equals `rep5_trial.csv`, and its`
with
`**Passed**: the redrawn trial equals `rep5_trial.csv` (largest relative difference `r format(redraw_rel, digits = 2)`), and its`

**E3d.** Insert after that paragraph:

```markdown
::: {.callout-note}
## Platform precision

`rep5_trial.csv` was written on the Linux machine that ran the campaign. On another platform the redraw
is the same trial, but continuous columns can differ in their last digits, because the linear-algebra
library and floating-point arithmetic differ. On macOS (arm64, R 4.5.2, Apple's Accelerate library),
two columns differ, `y_sim` and `z_wtkg`, by at most 5e-15 relative, and every integer and 0/1 column
is identical. The check therefore compares with a relative tolerance of 1e-12, not exactly.
:::
```

## 2. Gate and notes

- Render by the r3/r4 route (render, then `build_article()`, restore the local HTML). Stop on a render error or a failed check.
- On this machine `redraw_rel` should be 0 (the CSV was written here). Report its value; if it is not 0, say so, but do not stop.
- Commit the `.qmd` and the re-rendered `.html` together ("article: tested install sequence, rendering note, platform tolerance"), then `git checkout master && git merge --ff-only article-install-r2`. Do not push.
- Write `~/Downloads/NOTES_fs-mrct-article-install-r2_2026-11-02.md`: status, commits, new master SHA, `redraw_rel`, render time, any anchor not found, any error text.
