# REPORT — finish the articles, part B (2026-10-30)

Task: `dev/tasks/cc_task_articles_finish_partB_2026-10-30.md`. Branch `feature/articles-finish`, from `da17be3c`.
Record: `dev/reports/articles_finish_partB_2026-10-30/` (`logs/`, `results/`).
**Outcome: every gate passed. No article deviations.**

## §0 Provenance

- HEAD `da17be3c` on `feature/articles-finish`; no tracked changes. **UNTRACKED_0: empty.**
- The zip's md5 was `dc861f543191ecb2b756c1f730bff53e` (matches). `MANIFEST.md5`: all six files OK.
- Commit 1: `2150ba5c task: finish the articles, part B — task document and inputs`.
- R 4.6.1 (2026-06-24). forestsearch 0.4.1, pkgdown 2.2.1, quarto (R) 1.5.1, knitr 1.52,
  rmarkdown 2.32, roxygen2 8.1.1, devtools 2.5.2, xml2 1.6.0, gt 1.4.0, speff2trial 1.0.5,
  future 1.80.0. Quarto CLI `/usr/local/bin/quarto` 1.10.18 (on `PATH`). pandoc 3.1.3.
  urlchecker 2.0.0 installed; **spelling not installed**.
- `git diff --name-only c2967362 HEAD -- R/ DESCRIPTION NAMESPACE inst/`: empty.
  `Config/roxygen2/version: 8.1.1`. (`results/provenance.txt`)

## §1–3

- §1 pre-state: all seven md5s matched; the article did not exist. (`results/prestate_md5.txt`)
- §2: the article's md5 is `6782b084…` (matches). Applier md5 `1056eed3…`; `--check` and the apply both ended
  `RESULT: OK` with the six stated transitions (`logs/apply*.log`). Roxygen-only check: no
  non-roxygen lines; numstat `17 5 R/forestsearch_main.R`. (`results/gate2_roxygen.txt`)
- §3 `document()`: wrote `forestsearch.Rd` only; status ` M man/forestsearch.Rd`; numstat
  `17 5`; greps 1 / 1 / 0 / 1. (`results/gate3_rd.txt`)

## §4 `?forestsearch` Example 1

Exit 0, wall-clock **57 s** (cap 15 min). `example1_mr_table.html` is 15,029 bytes.
```
REGION {er <= 0} & {size <= 35}
MR corrected 1.444 (0.621, 3.360); complement 0.644 (0.402, 1.033)
BOOT nb_boots 20; H2 1.977
EXAMPLE1_OK
```

## §5–6 The article

| render | wall-clock | gate file | result |
|---|---|---|---|
| pkgdown `build_article(lazy = FALSE)` | 349 s | `results/gate_pkgdown.txt` | `RESULT: OK` |
| local `quarto render` | 328 s | `results/gate_local.txt` | `RESULT: OK` |

- `check_pkgdown()`: `No problems found`. §5.2: the page exists; the index links it ×2; the navbar label ×1.
- §5.3: the local HTML exists; `git check-ignore` → `.gitignore:39:vignettes/articles/*.html`.
- §5.5: `git status` is the same before and after the renders; the renders left no untracked paths.
- Gate (both renders): 7 tables. S16 and S17 match every published cell. The S16 Regions are
  `!{cd40 <= 380} & !{wtkg <= 86}` and `!{cd40 <= 507} & {gender}`. S3 matches, with FS Ĥ multiplier at
  `1.79 (0.63, 5.09)` (Gate C). Bootstrap column: 2.62 0.58 0.80 0.62 1.08 0.60. All 12 prose
  sentences appear once. `NaN` 0, `?@` 0, `Warning` 0.
- Key values (`results/keyvalues.txt`):
  ```
  Binary     !{cd40 <= 380} & !{wtkg <= 86}   N 72   fit 66 s (multiplier stage 11.5 s)
  Continuous !{cd40 <= 507} & {gender}        N 79   fit 102 s (multiplier stage 19.6 s)
  FS   !{cd40 <= 380} & !{wtkg <= 86}         N  72   fit  53.0 s (MR 11.3 s)  family 2794
  DINA {preanti >= 849.40000000000009} & {cd40 >= 338}            N 98  fit 2.4 s (MR 0.7 s)  family 147
  GRF  {wtkg > 84.369600000000005} & {cd40 > 368.20000000000005}  N 95  fit 16.7 s (MR 4.9 s) family 1289
  Cost:  FS 11.31 25 43.65 1.75 1746 | DINA 0.73 25 19.89 0.80 796 | GRF 4.91 10 26.08 2.61 2608
  ```
- Both renders took under 10 min, so there is no render-time OPEN ITEM.

**Article deviations: none.**

## §7 URLs and spelling

- `urlchecker::url_check()` reported 3 "invalid" URLs. All three are 301 Moved, and all are in
  `vignettes/articles/references.bib` (jthoracdis.com, R-project.org ×2 lines, bibsonomy.org), so none
  is in a gated file. **Gate: PASS.** (`logs/url_check.txt`, `logs/url_check_detail.txt`,
  `results/gate7_urls.txt`)
- spelling: not installed, so not run.

## §8 Release re-check

- `R CMD build` (with `RSTUDIO_PANDOC` exported): OK, 93 s. Tarball `forestsearch_0.4.1.tar.gz`,
  md5 `2db5afe8180a3e33981748c86e144474`.
- `R CMD check --as-cran --run-donttest`: **`Status: OK`**, 294 s.
  - Vignette re-build: [76s/74s]. `cran-comments.md` says "about 75 s".
  - Examples: [43s/44s]. `cran-comments.md` says "about 40 s".
  - Both are within half of the stated figures, so there is no OPEN ITEM.
- §8.3 `vignette_gate.R` on the shipped `inst/doc/forestsearch.html`: `RESULT: OK`.
- §8.4: no tarball path under `vignettes/articles/`, `dev/` or `docs/` (`logs/tarball_files.txt`).
- Execution note on §8.1: `git ls-files --cached` quotes two tracked non-ASCII paths
  (`"dev/glm/articles/An introduction to `grf` \342\200\242 grf.html"`,
  `"quarto/methodology/articles/Statistics in Medicine - 2023 - Luo - …\342\200\220….pdf"`). Because of the
  leading `"`, the `grep -vE '^(dev|quarto|…)/'` let them through, and `cp` failed on them.
  - Every other file was copied. I compared the copied tree with
    `git -c core.quotepath=off ls-files --cached` under the same filter: they are identical (627 files).
    The build and check then ran as written.
  - My `&&`-chained first attempt also wrote a 44-byte `build.log` into the repo root. That attempt never
    reached `/tmp/fs_partB`, so `R CMD build src` failed there. I deleted the stray `build.log`.
  - I extracted the shipped vignette to `/tmp/fs_partB/x/`, not `/tmp/fs_partB/`. It is the same file;
    the different directory avoided touching the check directory while the check was running.

## §9 Commits

```
2150ba5c task: finish the articles, part B — task document and inputs
17f0c718 docs: article — the ACTG 175 analyses (2026b S3, 2026a Table S3)
75285f8d docs: ?forestsearch examples — mr_estimates_table() with its bootstrap template; Example 2's source; a section reference
cd153dfe docs: vignette next steps, NEWS, README and cran-comments for the articles
(this report) dev: report — articles finish, part B (2026-10-30)
```
Post-conditions:
- `git diff --name-only 0bdc2036 HEAD | grep -vE '^(dev/|vignettes/articles/|\.gitignore$)'` prints
  `NEWS.md`, `R/forestsearch_main.R`, `README.md`, `_pkgdown.yml`, `cran-comments.md`, `man/forestsearch.Rd`,
  `vignettes/forestsearch.qmd`.
- Of these, the tarball carries `NEWS.md`, `R/forestsearch_main.R`, `README.md`, `man/forestsearch.Rd`,
  `vignettes/forestsearch.qmd`.
- `git ls-files 'vignettes/articles/*.html'` prints nothing.

## OPEN ITEMS (findings only)

1. **Floating-point region labels in the article.** The Part 3 printed output shows DINA's and GRF's
   regions as `{preanti >= 849.40000000000009}`, `{wtkg > 84.369600000000005}` and `{cd40 > 368.20000000000005}`.
   These are the package's own `sg.harm` labels; the gate does not check them, and the prose rounds them
   (849, 84.4 kg, 368). Fixing the display would need a change to either the label formatting in the
   DINA/GRF code or the article's assembly code. I made neither change, because this task forbids it.
2. **The §8.1 recipe and non-ASCII paths.** `git ls-files` quotes non-ASCII paths, so the exclusion grep
   misses them. Here the result was harmless (see the §8 note). For future tasks, use
   `git -c core.quotepath=off ls-files -z` or anchor `^"?`.
3. Three 301-moved URLs in `vignettes/articles/references.bib`. They are not in the tarball.
4. `spelling` is not installed, so the spell check has not been run.

## Verdict

```
Every gate of §0–§8 passed; the task ran top to bottom with no deviation under §5.4.
The ACTG 175 article renders in both routes (349 s / 328 s) and its gate reads RESULT: OK twice.
Table S16 and S17 reproduce every published cell; Table S3 matches with FS H at 1.79 (0.63, 5.09).
?forestsearch Example 1 now runs to the end (57 s), mr_estimates_table() with its bootstrap.
The Rd page changed by exactly the roxygen edit (17/5); no other generated file moved.
The shipped vignette carries the two new bullets and links; the old wording is gone.
R CMD check --as-cran --run-donttest on the tracked-files tarball: Status: OK.
Check timings (74 s vignette, 44 s examples) agree with cran-comments.md.
No articles, dev/ or docs/ paths enter the tarball; no HTML is tracked.
Ready for Larry's review and merge (§11); four findings above, none blocking.
```
