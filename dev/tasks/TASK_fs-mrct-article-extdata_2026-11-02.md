# TASK: MRCT article — records shipped with the package, and §2.6 the S2 campaign in brief (forestsearch, 2026-11-02)

**Repository:** forestsearch. **Branch:** create `article-extdata` from `master` (stop if `master` does not contain `e99c636e`; on the Linux box, pull `master` first). **Touches:** the article's data folder (moved), `vignettes/articles/mrct-consistency.qmd` and its tracked `.html`, `dev/article/make_rep5_trial.R` (output path only), and any file that refers to the old data path (§1.3). **Nothing under `R/`** (if §1.3 finds a reference under `R/`, stop and report it).

**Package contents change:** the 29 article records become installed files (`inst/extdata/mrct-consistency/`, about 316 KB). No code changes behaviour.

## 0. First action

Copy this file to `dev/tasks/TASK_fs-mrct-article-extdata_2026-11-02.md` and commit it alone ("task: MRCT article records in extdata, §2.6").

## Why

The article must run for a reader who has only the installed package. Today its data folder `vignettes/articles/mrct-consistency/` is excluded from the build (`^vignettes/articles$` in `.Rbuildignore`), so an installed package has none of the 29 records, and the article stops in its first chunk (`stopifnot(file.exists(file.path(dd, "PROVENANCE.md")))`). Larry chose (2026-11-02) to ship the whole folder under `inst/extdata/`.

## 1. The records move to `inst/extdata/mrct-consistency/`

1.1 `git mv vignettes/articles/mrct-consistency inst/extdata/mrct-consistency` (all 29 files, `PROVENANCE.md` included).

1.2 In the article's `packages` chunk, replace the `dd <- ...` assignment and the `stopifnot(file.exists(...PROVENANCE.md...))` line after it with the following, and update the two comment lines above it as shown:

```r
# The article's records ship with the package (inst/extdata/mrct-consistency; provenance in its
# PROVENANCE.md). From a source tree the local copy is used; otherwise the installed one.
dd <- c(file.path("..", "..", "inst", "extdata", "mrct-consistency"),
        file.path("inst", "extdata", "mrct-consistency"),
        system.file("extdata", "mrct-consistency", package = "forestsearch"))
dd <- dd[nzchar(dd) & dir.exists(dd)][1]
stopifnot(!is.na(dd), file.exists(file.path(dd, "PROVENANCE.md")))
```

The `rd()` helper and every `rd()` and `file.path(dd, ...)` call stay as they are.

1.3 Grep the repository (outside `inst/extdata/` itself) for `vignettes/articles/mrct-consistency` and for `articles/mrct-consistency/`: tests, `dev/` scripts, `_pkgdown.yml`, `.Rbuildignore`, `.github/`, `PROVENANCE.md`, `NEWS.md`, the article text. Update each to the new location, using `system.file("extdata", "mrct-consistency", package = "forestsearch")` in tests (or the test's existing fixture convention). List every change in the notes. A reference under `R/`: stop.

1.4 `dev/article/make_rep5_trial.R`: its output path (`dd`) becomes `file.path("inst", "extdata", "mrct-consistency")`. Do not re-run it.

1.5 `PROVENANCE.md`: one line at the top stating that the folder is installed with the package at `system.file("extdata", "mrct-consistency", package = "forestsearch")`, moved from `vignettes/articles/mrct-consistency/` on 2026-11-02. Existing entries unchanged.

1.6 Article text: wherever the article says its records are "shipped with the article" or names the old folder, say "shipped with the package" / name `inst/extdata/mrct-consistency`. Report each sentence changed.

## 2. Install and rendering notes

2.1 In the "Installing forestsearch" callout, step 3's line becomes
`install.packages(c("speff2trial", "latentcor", "knitr", "remotes"))`
and its comment gains `knitr for the tables,` after the latentcor clause. Also add `"knitr"` to the `missing_pkgs` check in the `install-check` chunk (both occurrences of the package vector).

2.2 Replace the body of the "Rendering this article yourself" callout (keep its title and `collapse="true"`) with:

````markdown
Rendering needs the Quarto command-line tools and about five minutes. The article's records are
installed with the package, so only the `.qmd` file itself is needed. Download it into a scratch
folder and render it there:

```r
dir.create("~/scratch/mrct-article", recursive = TRUE, showWarnings = FALSE)
download.file("https://raw.githubusercontent.com/larry-leon/forestsearch/master/vignettes/articles/mrct-consistency.qmd",
              "~/scratch/mrct-article/mrct-consistency.qmd")
```

```bash
cd ~/scratch/mrct-article && quarto render mrct-consistency.qmd
```

The rendered `mrct-consistency.html` is written beside the `.qmd`. Delete the scratch folder when done.
````

## 3. New §2.6 — the S2 campaign in brief

Insert at the end of the case-study section (`#sec-case`), after the "Reading the bounds" callout and before `# Setting up simulations from your own trial`. Confirm in the render that it is numbered 2.6.

```markdown
## The S2 campaign in brief {#sec-case-campaign}

Replicate 5 is one trial of S2. The campaign ran 10,000 trials of the same design (@sec-sims); the table
gives how often each step of the pipeline above went the same way. The numbers are read from the
campaign's records, shipped with the package.
```

Then a chunk `case-campaign` reading, for S2:
- from `core_full_numeric.csv`, row `setting == "S2_rmlog4"`: `replicates`, `n_trigger`, `n_declared`, `n_declared_Ghat_lt1`;
- from `mr_field_core.csv`, row `lb_100`, column `S2` (Ĥ field lower bound > 1.0, % of declared triggered trials);
- from `entry_cells_core_numeric.csv`, quantity `Ghat_region`, column `S2` (Ĝ, % of region patients, over declared triggered trials).

Check with `stopifnot()` that the two records agree on S2: `n_trigger` equals `entry_cells_core_numeric.csv` `trigger_n`, and `100 * n_declared_Ghat_lt1 / n_trigger` equals its `joint` to 1e-9.

Show one table (`knitr::kable()`), three columns `step`, `S2`, `of`, rows in this order:

| step | S2 | of |
|---|---|---|
| trigger: ITT ub < 1.0 and region ITT > 1.0 | % (count) | all trials (replicates) |
| subgroup declared | % to two decimals (count) | triggered trials |
| Ĥ field lower bound (one-sided 95%) > 1.0 | % | declared triggered trials |
| declared, region(Ĝ) < 1.0 | % (count) | triggered trials |
| Ĝ, % of region patients | % | declared triggered trials |

Percentages as bare numbers to one decimal except where stated; counts with thousands separators. Caption: "S2 at n = 1,000, 10,000 trials: the pipeline of this section, as often as it went each way."

Then the reading (verbatim, with inline R for each number from the chunk):

```markdown
About one S2 trial in four reached the trigger (`r …`%). In the triggered trials a subgroup was
declared in `r …` of `r …`, region(Ĝ) was below 1.0 in about half (`r …`%), and Ĝ kept on average
`r …`% of the region's patients. The field lower bound for Ĥ was above 1.0 in `r …`% of the declared
triggered trials; in the others it lies below 1.0, and it is read by its location, as for replicate 5
in step 5. Replicate 5 is one of the triggered trials with region(Ĝ) < 1.0.
```

"on average": check `PROVENANCE.md` (or the record's own label) for how `Ghat_region` is computed. If it is a mean over declared triggered trials, keep "on average"; if it is another statistic, name that statistic instead; if the record does not say, drop "on average" and say so in the notes.

## 4. Gates and notes

- Render by the r3/r4 route (render, then `build_article()`, restore the local HTML). Stop on a render error or a failed check.
- Run the test files changed in §1.3, if any; stop on a failure.
- `R CMD build --no-build-vignettes --no-manual .`: confirm the tarball contains `forestsearch/inst/extdata/mrct-consistency/` with 29 files, and report the tarball size. (No full check in this task.)
- Commit on `article-extdata` (the move may be its own commit, before the article edits), then `git checkout master && git merge --ff-only article-extdata`. Do not push.
- Write `~/Downloads/NOTES_fs-mrct-article-extdata_2026-11-02.md`: status, commits, new master SHA, every path reference changed (§1.3, §1.6), the `Ghat_region` statistic finding, the §2.6 table as rendered, tarball size, render time, any error text.
