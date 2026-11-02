# TASK: MRCT article — installing forestsearch, and an early check (forestsearch, 2026-11-02)

**Repository:** forestsearch. **Branch:** create `article-install-notes` from `master` (after `track-article-html` is merged into `master`; stop if `master` does not contain that commit). **Touches:** `vignettes/articles/mrct-consistency.qmd` and its tracked `.html`. Nothing under `R/`.

## 0. First action

Copy this file to `dev/tasks/TASK_fs-mrct-article-install_2026-11-02.md` and commit it alone ("task: MRCT article install notes").

## 1. Verify first

From `NEWS.md`, confirm that the functions listed in §2.2 are under the development heading (added after the CRAN release 0.4.0). If any is in 0.4.0, drop it from the list; if the statement "the CRAN release does not yet contain the functions this article uses" is then false, stop and say so in the notes.

## 2. Insertions (verbatim)

### 2.1 After the opening paragraph (the one ending `is applied inside it.`), before `# Overview`

````markdown
::: {.callout-note}
## Installing forestsearch

The CRAN release of forestsearch (0.4.0) does not yet contain the functions this article uses. Until
the next CRAN release, install the development version from GitHub. The article also uses two
packages that forestsearch only suggests: `speff2trial`, for the ACTG 175 data, and `latentcor`, for
the copula expansion.

```r
# the CRAN release (does not yet include this article's functions)
install.packages("forestsearch")

# the development version from GitHub, with the suggested packages
install.packages("remotes")
remotes::install_github("larry-leon/forestsearch", dependencies = TRUE)

# or install the two suggested packages directly
install.packages(c("speff2trial", "latentcor"))
```

Rendering the article itself needs the Quarto command-line tools. The chunk below stops with a
message naming what is missing, before any computation starts.
:::
````

### 2.2 Immediately after that callout: an evaluated check chunk

```{r install-check}
#| code-fold: show
needed <- c("inject_mrct_structure", "scale_dgm_hazard", "ahr_population", "ahr_limit",
            "ahr_design", "ahr_solve", "mrct_region_sims", "mrct_select_case",
            "mrct_case_study", "mr_estimates_table")
missing_fns <- setdiff(needed, getNamespaceExports("forestsearch"))
if (length(missing_fns))
  stop("This article needs the development version of forestsearch (see 'Installing forestsearch'); ",
       "the installed version (", as.character(packageVersion("forestsearch")), ") lacks: ",
       paste(missing_fns, collapse = ", "), ".", call. = FALSE)
missing_pkgs <- c("speff2trial", "latentcor")[!vapply(c("speff2trial", "latentcor"), requireNamespace, NA, quietly = TRUE)]
if (length(missing_pkgs))
  stop("This article needs: ", paste(missing_pkgs, collapse = ", "),
       " (install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), ")).", call. = FALSE)
c(forestsearch = as.character(packageVersion("forestsearch")))
```

(Write it as a normal R chunk in the .qmd; the fence above is shown literally. Adjust the `needed` list per §1.)

The `install-check` chunk must come before the first chunk that uses any of these functions or `speff2trial`. If the existing `packages` chunk or a data chunk already uses one earlier, move the check to just after `library(forestsearch)` instead and say where in the notes.

**Exception for `latentcor`:** the collapsed "Reproducing this report's build exactly" path does not need it, but the reader's DGM build and the own-trial build do. Keep it in the check.

## 3. Gate and notes

- Render by the r3/r4 route (render, then `build_article()`, restore the local HTML). Stop on a render error or a failed check.
- Commit the `.qmd` and the re-rendered `.html` together ("article: install notes and early check"), then `git checkout master && git merge --ff-only article-install-notes`. Do not push.
- Write `~/Downloads/NOTES_fs-mrct-article-install_2026-11-02.md`: status, commits, the new master SHA, the §1 finding, where the check chunk sits, render time, any error text.
