# REPORT — post-selection analysis to a pkgdown article (2026-10-28)

**Status: COMPLETE under v2** (all gates V1–V5, G10–G23 pass; five commits since `b2e43397`; nothing pushed). v1 stopped at G13 ("Pandoc not available"); the v1 stop record is kept below unchanged.

- Task: `dev/tasks/cc_task_postsel_to_article_2026-10-28.md`
- Start: 2026-10-28T13:31:58-07:00; stop: 2026-10-28T13:34 PT
- `future::availableCores()` = 128
- Record directory: `dev/reports/postsel_to_article_2026-10-28/` (`inventory.md`, `g10_purl_identity.txt`, `check_pkgdown.txt`, `check_pkgdown_diagnostic_with_RSTUDIO_PANDOC.txt`) — untracked.

## Failing gate and proposed fix

- **G13 FAIL.** The command as written, `timeout 120 Rscript -e 'pkgdown::check_pkgdown()'`, exited 1 with `Error in convert_markdown_to_html(): ! Pandoc not available` (full output in `check_pkgdown.txt`). Cause: pandoc is not on `PATH` in a non-interactive `Rscript` session on this machine (CLAUDE.md); the gate does not export `RSTUDIO_PANDOC`, unlike the Quarto runs under "Fixed values". Nothing in the edits is implicated.
- **Diagnostic, not a gate result:** the same command with `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64` exported exited 0 with `✔ No problems found.` (`check_pkgdown_diagnostic_with_RSTUDIO_PANDOC.txt`).
- **Proposed fix:** amend G13 to `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64 timeout 120 Rscript -e 'pkgdown::check_pkgdown()'` and resume from G13 on the current working tree (G1–G12, G14, G15 hold on it; the §7.1 commit then becomes the third commit as planned, so G20 is unaffected).

## Commits

| Hash | Subject |
|---|---|
| `0c07d261` | dev: task — post-selection vignette to pkgdown article (2026-10-28) |
| `119b58b4` | docs: move survival-post-selection.qmd to vignettes/articles/ (pkgdown article) |

Not made (stop): the §7.1 edits commit and the §7.2 report commit.

v2 (`dev/tasks/cc_task_postsel_to_article_2026-10-28_v2.md`):

| Hash | Subject |
|---|---|
| `50aa747a` | dev: task v2 — post-selection to pkgdown article, resume from G13 (2026-10-28) |
| `4768365e` | docs: post-selection analysis as a pkgdown article — note, cross-references, _pkgdown.yml, NEWS; primary vignette timing sentence |
| (this commit) | dev: report — post-selection to pkgdown article (2026-10-28) |

## Gates

| Gate | Result | Note |
|---|---|---|
| G1 | PASS | branch `feature/vignettes-0.4.1`; HEAD~1 = `b2e43397aa36…` |
| G2 | PASS | only the two `??` html entries |
| G3 | PASS | md5 `660b2f87bd26edd0a5e1fdd9801fd09e` |
| G4 | PASS | `git diff --quiet 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` exit 0 |
| G5 | PASS | `R100	vignettes/survival-post-selection.qmd	vignettes/articles/survival-post-selection.qmd` |
| G6 | PASS | all four *edit* hits are prose (blockquotes, one bullet) |
| G7 | PASS | `_pkgdown.yml`: only L45 navbar `href: articles/survival-post-selection.html` (kept) and L546 index entry |
| G8 | PASS | one bullet (L71–75) |
| G9 | PASS | 1 occurrence before the edit |
| G10 | PASS | purl byte-identical for all five edited files; no path-rewrite lines (none needed) |
| G11 | PASS | 0 / 0 / 0; `survival-post-selection` 4 = URL 4 (README.Rmd does not exist; README.md has no hits) |
| G12 | PASS | slug once, inside URL; only hunk `@@ -71,5 +71 @@`, development section (next heading now L73) |
| G13 | **FAIL** | "Pandoc not available" — see above |
| G14 | PASS | 0 / 1 on the source |
| G15 | PASS | `.Rbuildignore`, `NEWS.md`, `_pkgdown.yml`, the three *edit* articles, the moved article, `vignettes/forestsearch.qmd` |
| G16–G23 | NOT RUN | stop at G13 (v1) |
| V1 (v2) | PASS | HEAD~1 `119b58b4…`, HEAD~2 `0c07d261…` |
| V2 (v2) | PASS | `git diff --name-only` = the eight files; nothing cached |
| V3 (v2) | PASS | eight `M`, plus the four expected `??` entries |
| V4 (v2) | PASS | `git diff --quiet 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` exit 0 |
| V5 (v2) | PASS | restored line once in the development section, at L72, directly after the post-selection bullet; one-line insertion only |
| G10 (v2) | PASS | purl byte-identical for all five files (`g10_g15_v2.txt`) |
| G11 (v2) | PASS | 0 / 0 / 0; slug 4 = URL 4 |
| G12 (v2) | PASS | slug once, inside URL; only hunk `@@ -71,5 +71,2 @@` |
| G13 (v2) | PASS | `RSTUDIO_PANDOC` set; exit 0, `✔ No problems found.` (`check_pkgdown_v2.txt`) |
| G14 (v2) | PASS | 0 / 1 |
| G15 (v2) | PASS | the eight files |
| G16 (v2) | PASS | exit 0, wall 2:08; html 13:42:08 > start 13:39:24; no `_files/`; 0 `citeproc`/`not found` lines (`render_article.log`) |
| G17 (v2) | PASS | note 1; six HR/CI strings 4/4/3/3/3/3; `1.88` 2, `0.83` 1, `1.29` 1, `Efron` 1 (`g17_article_text_v2.txt`) |
| G18 (v2) | PASS | exit 0, wall 1:14; 1 / 0 / `post-selection article` 1 (`render_primary.log`, `g18_primary_text_v2.txt`) |
| G19 (v2) | PASS | no `HANDOFF_`, `dev/`, `vignettes/articles/` or `survival-post-selection` entries (`g19_tarball_v2.txt`, `g19_build_v2.log`) |
| G20 (v2) | PASS* | five commits `b2e43397..HEAD` |
| G21 (v2) | PASS | empty; no `HANDOFF_*.md` changed |
| G22 (v2) | PASS* | only the three `??` html entries |
| G23 (v2) | PASS | nothing at `vignettes/` top level begins `survival-post-selection` |

\* G20 and G22 depend on this report's own commit; they were checked immediately after it (results relayed in the close-out message, and identical to the state described here).

(G14 and G15 were checked in the same read-only pass as G11/G12, before G13 was run.)

## Before → after (verbatim, uncommitted)

**Article note** (`vignettes/articles/survival-post-selection.qmd`, inserted after the YAML close, L18; n_vig = 6 so the second sentence is appended):

> > **A website article.** In forestsearch 0.4.0 this was a CRAN vignette; from 0.4.1 it is an article on the package website, because the full analysis takes a few minutes to run, too long for CRAN's package checks. The text is otherwise unchanged; where it says "vignette", read "article".

**Cross-references (§5.3):**

- `vignettes/articles/forestsearch-2024.Rmd:29`
  - before: `> in `vignette("forestsearch")` and `vignette("survival-post-selection")`. The code below runs as written;`
  - after: `> in `vignette("forestsearch")` and the [post-selection article](https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html). The code below runs as written;`
- `vignettes/articles/methodology.Rmd:64`
  - before: `> current methodology is documented in `vignette("forestsearch")` and `vignette("survival-post-selection")`;`
  - after: `> current methodology is documented in `vignette("forestsearch")` and the [post-selection article](https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html);`
- `vignettes/articles/paper_simulations.Rmd:59`
  - before: `> `vignette("forestsearch")` and `vignette("survival-post-selection")`.`
  - after: `> `vignette("forestsearch")` and the [post-selection article](https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html).`
- `vignettes/forestsearch.qmd:611`
  - before: `- `vignette("survival-post-selection")` reads the bounds in more depth on the`
  - after: `- The [post-selection article](https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html) reads the bounds in more depth on the`

**Timing sentence** (`vignettes/forestsearch.qmd`): `about a minute and a half on one core` → `about a minute on one core`.

**NEWS bullet** (development section, L71–75 → L71):

- before:
  ```
  * `vignette("survival-post-selection")` now reproduces the GBSG analysis of León and Anderson (2026a):
    the three identifiers at the manuscript's settings, multiplier resampling in full and the full
    bootstrap at a small B for illustration, with `summarize_bootstrap_results()` and
    `mr_estimates_table()`. The methodology and simulation articles on the package website carry a
    note placing them as the 2024 workflow, with the deprecated arguments removed from their listings.
  ```
- after:
  ```
  * The post-selection analysis is now an article on the package website rather than a CRAN vignette, because the full analysis takes a few minutes to run: <https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html>. It reproduces the GBSG analysis of León and Anderson (2026a), illustrates the full bootstrap at a small number of replicates, and shows `mr_estimates_table()`.
  ```

**`_pkgdown.yml`** L546: `      - survival-post-selection` → `      - articles/survival-post-selection`.

**Path rewrites:** none (the article has no relative paths; inventory §2).

**`.Rbuildignore` lines added:** `^HANDOFF_.*\.md$` (`^vignettes/articles$` and `^dev$` already present).

**Moved, not staged:** `vignettes/survival-post-selection.html` (pre-move render, git-ignored by `.gitignore:15`) → `~/Downloads/postsel_to_article_2026-10-28/pre_move_render/`.

## Record-only references (for the release task)

- `vignettes/articles/survival-post-selection.qmd:113,127` — "companion vignette" (article text, unchanged by rule).
- `R/forestsearch_main.R:1435` — `@examples` `\dontrun{}` comment "…and the survival-post-selection ## vignette."; `man/forestsearch.Rd:1309` — generated copy.
- `HANDOFF_vignettes_release_0.4.1_2026-10-27.md`: L58, 74, 77, 81, 89, 94, 102.
- `HANDOFF_vignettes_release_0.4.1_2026-10-28.md`: L25, 35, 38, 39, 43, 61, 74, 99, 163.
- `quarto/`: `extreme_subgroups/gbsg_design_comparison.qmd:37,1251`; `extreme_subgroups/gbsg_survival_fixed.qmd:57,515,767`; `extreme_subgroups/gbsg_survival_resample.qmd:57,456,692`; `guides/extreme_subgroups.qmd:384,620`; `guides/forestsearch_diagnostic_scoping.qmd:52`; `presentations/my use of claude/extreme_subgroups.qmd:384,620`; `simulations/README.qmd:456`; `simulations/actg175/README.qmd:456`; `simulations/gbsg_020/REPORT_field_recovery_2026-09-09.md:446,484`; `simulations/gbsg_020/REPORT_print_vignette_2026-09-09.md:1,227`.
- `dev/`: many reports, tasks and logs (per-file counts in `inventory.md` §1).
- Tracked HTML with the slug: `vignettes/articles/forestsearch-2024.html`.

## Renders, G19

Not run (stop at G13). `availableCores()` = 128.

## Side notes (not acted on)

- The replaced NEWS bullet carried a second sentence — "The methodology and simulation articles on the package website carry a note placing them as the 2024 workflow, with the deprecated arguments removed from their listings." — which the replacement bullet, as specified, drops. If it should stay, it needs its own bullet.
- The moved article keeps its `vignette:` YAML block (`%\VignetteIndexEntry`, `%\VignetteEngine{quarto::html}`); harmless under `^vignettes/articles$`, but the release task may want it removed. Article L253 says "This vignette uses at most two, the limit a CRAN check imposes" — now covered by the note's "read 'article'", though the CRAN rationale no longer applies.
- The note is inserted after the YAML's closing `---` with a blank line before it; the original blank line that followed the header remains, so two blank lines precede the setup chunk (no rendering effect).

## Resumed under v2

- Run start 2026-10-28T13:39:24-07:00. `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64` set in every R/Rscript/Quarto/pandoc invocation.
- **Restored NEWS line** (L72, directly after the post-selection bullet): `* The methodology and simulation articles on the package website carry a note placing them as the 2024 workflow, with the deprecated arguments removed from their listings.` This closes the first v1 side note.
- **Render wall clocks:** article 2:08 (cap 6:00; measured 2:07 before); primary vignette 1:14 (cap 4:00; measured 1:14 before).
- **Article plain-text diff against the pre-move render** (`article_plaintext_diff_vs_pre_move.txt`, recorded, not gated): the added note, plus run-to-run timing values only (fit times, bootstrap time and projection, cost table timings). No estimate changed.
- **G19:** `forestsearch_0.4.0.9000.tar.gz`, 1,698,826 bytes. Top level: `build`, `DESCRIPTION`, `LICENSE`, `man`, `NAMESPACE`, `NEWS.md`, `R`, `README.md`, `tests`, `vignettes` — no flagged entries. `vignettes/`: `forestsearch.qmd` and `forestsearch.html`.
- Side note (not acted on): the `vignettes/forestsearch.html` in the tarball is the G18 render sitting in the working tree (git-ignored, not in `.Rbuildignore`); `--no-build-vignettes` copies it in. A release build that does build vignettes regenerates it, so this matters only for a build made with `--no-build-vignettes`.
- The other v1 side notes (the moved article's `vignette:` YAML block and "This vignette uses at most two" wording) stand for the release task.
