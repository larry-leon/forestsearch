# REPORT — Guo & He (S4) article: place, reference gate, render, pkgdown (2026-10-30)

Task: `dev/tasks/cc_task_guohe_article_2026-10-30_v2.md` (v2, authoritative; supersedes v1
`dev/tasks/cc_task_guohe_article_2026-10-30.md`). Record directory:
`dev/reports/guohe_article_2026-10-30/` (`s1_reference_gate.R`, `s2_text_gates.sh`, `results/`, `logs/`).

## 1. Provenance

| item | value |
|---|---|
| branch | `feature/articles-guohe` |
| BASE_SHA | `f2ddb5c0a4209451011f08857469526bbc670acd` (`git status --porcelain` empty) |
| commit 1 (v1 task doc alone) | `4ea458bb43930fbf38ab2e9bb94dd69084653107` |
| commit 1b (v2 task doc alone) | `4933e50f16f068b2c0c7d537e7b89d6199116474` |
| commit 2 (article + `_pkgdown.yml`) | `404abf72a7f9a049711debb2a73e0d49b9e51695` |
| commit 3 (record + this report) | see closing block |
| input | `~/Downloads/guohe-nested-continuum.qmd`, md5 `18091f48356588306e784c8580e388d4`, 939 lines; placed file identical (same md5) and committed unchanged |
| task doc filenames | `~/Downloads/cc_task_guohe_article_2026-10-30.md` → `dev/tasks/cc_task_guohe_article_2026-10-30.md`; `~/Downloads/cc_task_guohe_article_2026-10-30_v2.md` → `dev/tasks/cc_task_guohe_article_2026-10-30_v2.md` |
| forestsearch / availableCores | `0.4.1 128` |
| quarto | `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto`, 1.10.18 (the RStudio-bundled binary of earlier article tasks) |
| pandoc | `pandoc 3.1.3` (`/usr/bin/pandoc`); `RSTUDIO_PANDOC` set to the bundled tools dir for the render |
| render wall-clock | 5:24.18 (`/usr/bin/time -v`; 10:29:12 → 10:34:36 −07:00), max RSS 501 MB |
| plan the article resolved | 2 workers (`n_workers = min(2L, availableCores())`; caption "on 2 workers") |

**v1 → v2.** v1 stopped at its §5: its text gates counted the whole rendered page, and the
article's `code-fold: false` echoes every chunk, so gates 1, 2b, 2c and 5c came out one too
high (2, 2, 2, 1). Larry's disposition (option 1, 2026-10-30), recorded in v2: count the
rendered output only, no second render, article unchanged. v2 was copied and committed alone
(`4933e50f`) before §5–§7 were executed; the resume rule's first branch did not apply (§5–§7 had
not been completed under option 1). **Method to strip the echoed code:** every
`<pre …>…</pre>` element is deleted from the HTML (R, `gsub("(?s)<pre[ >].*?</pre>", "", …, perl = TRUE)`;
28 elements: 26 `sourceCode r` code blocks and 2 bare `<pre>` printed-output blocks — the
`design-one` summary vector and `sessionInfo()`, neither holding a gated string), then
`pandoc --wrap=none -f html -t plain` → `results/article_output.txt`. Script `s2_text_gates.sh`.

## 2. Reference gate (§3)

**Green: 504 comparisons, 504 passes, 0 failures. Corrections: none.**

| table | passes |
|---|---|
| ref_lower ↔ bound_t7 | 144/144 |
| ref_select ↔ ident_t7 | 108/108 |
| ref_location ↔ exp(beta2), exp(gamma_mean) | 12/12 |
| ref_upper ↔ comp_t7 | 72/72 |
| ref_joint ↔ joint_t7 | 60/60 |
| ref_16 ↔ pre16 | 96/96 |
| ref_gh_published ↔ GH52_TABLE_7 | 12/12 |

Payload meta beside the article's stated manuscript settings:

| meta | payload | article states | agree |
|---|---|---|---|
| nrep | 2000 | 2000 | yes |
| draws | 5000 | 5000 | yes |
| field_R_out | 1000 | 1000 | yes |
| field_R_in | 500 | 500 | yes |
| head_at_render | `d975a51eda0420775f6b6e4f55b424380e5657e8` | — | — |
| forestsearch_version | 0.3.5 | — | — |

Per-comparison record: `results/reference_gate.csv`; log `logs/s1_reference_gate.log`.

## 3. Render (§4)

Exit status 0, wall-clock 5:24.18, no `_files/` directory (embed-resources). No deviations.
The article's `stopifnot()` did not fire. Log: `logs/render.log`.

Transplanted from `results/article_output.txt`:

```
  --------------------------------------------------------------------------------------
    beta2 beta(30), estimated   anchor deviation / SE   beta(60)   basis
  ------- --------------------- ----------------------- ---------- ---------------------
      0.0 -0.0140               -1.48                   0.0000     exact zero
      0.3 0.2967                -0.36                   0.1484     isotonic projection
  --------------------------------------------------------------------------------------
  : The truth curve by cell (n_big = 200,000; 14.2 s).

  ------------------------------------------------------------------------------------
  cell             replicates   seconds   seconds per replicate (wall, across workers)
  -------------- ------------ --------- ----------------------------------------------
  beta_2 = 0.0             25     117.7                                            4.7
  beta_2 = 0.3             25     120.2                                            4.8
  ------------------------------------------------------------------------------------
  : Wall-clock of this article's run on 2 workers (Guo and He at B = 200; MR at 5000 draws,
    field 1000 x 500; truth curves 14.2 s in addition).

  ----------------------------------------------------------------------------------------
                                                                      hazard-ratio scale
  ------------------------------------------------------------------- --------------------
  selected cutpoint c-hat (true boundary 30)                          35.5
  selected subgroup |H-hat| (H is about 150)                          176 of 400
  hazard ratio in H-hat, the target exp{beta(H-hat)}                  1.28
  unadjusted: estimate, 95% lower bound                               1.58, 1.13
  Guo-He r = 1/12: bias-reduced estimate, 95% lower bound             1.41, 1.02
  MR (field): estimate, 95% lower bound                               1.47, 1.08
  MR (IJ two-term): estimate, 95% lower bound                         1.48, 0.80
  re-selection frequency p-hat(H-hat)                                 0.218
  complement, 95% upper bound: unadjusted / IJ two-term / field-s     1.22 / 1.70 / 1.31
  Bonferroni pair: lower on H-hat / upper on H-hat^c (0.975 a side)   1.01 / 1.35
  ----------------------------------------------------------------------------------------
  : Replicate 1, beta_2 = 0.3, at the manuscript's settings (141 candidates; 64 s).
```

## 4. Text gates (§5, v2: output only) — green

| gate | output-only | v1 whole-page | diff | expected | result |
|---|---|---|---|---|---|
| 1 walkthrough at manuscript settings | 1 | 2 | +1 | 1 | PASS |
| 2a own lower bounds | 1 | 1 | 0 | 1 | PASS |
| 2b own argmax selects | 1 | 2 | +1 | 1 | PASS |
| 2c own Bonferroni pair | 1 | 2 | +1 | 1 | PASS |
| 3 wall-clock on 2 workers | 1 | 1 | 0 | 1 | PASS |
| 4a `0.952 (0.941, 0.960)` | 2 | 2 | 0 | ≥ 2 | PASS |
| 4b `0.928 (0.916, 0.939)` | 1 | 1 | 0 | ≥ 1 | PASS |
| 4c `0.945 / 0.302 / +0.010` | 1 | 1 | 0 | 1 | PASS |
| 5a `NaN` | 0 | 0 | 0 | 0 | PASS |
| 5b `\bNA \(` | 0 | 0 | 0 | 0 | PASS |
| 5c whole-word `Inf` | 0 | 1 | +1 | 0 | PASS |
| 6 selection agreement | 1 | 1 | 0 | 1 | PASS |
| 7 numbered lines (source) | 9 inside `## Limits`…`## Running`, 0 outside | — | — | 0 outside | PASS |
| 8a lower-bound figure caption | 1 | 1 | 0 | ≥ 1 | PASS |
| 8b complement figure caption | 1 | 1 | 0 | ≥ 1 | PASS |

The whole-page counts exceed the output-only counts by exactly +1 on gates 1, 2b, 2c and 5c
and nowhere else, as v2 states. All in `results/keyvalues.txt` (v1 block, then v2 block).

## 5. `_pkgdown.yml` (§6) — green

```
@@ -43,6 +43,8 @@ navbar:
           href: articles/methodology.html
         - text: Post-Selection Inference
           href: articles/survival-post-selection.html
+        - text: Comparison with Guo & He (2021)
+          href: articles/guohe-nested-continuum.html
         - text: Treatment Effect Definitions
           href: articles/treatment_effect_definitions.html
@@ -544,6 +546,7 @@ articles:
     desc: Inference for trial estimands after data-driven subgroup discovery
     contents:
       - articles/survival-post-selection
+      - articles/guohe-nested-continuum
```

`git diff --stat`: `_pkgdown.yml | 3 +++` — 1 file changed, 3 insertions(+), 0 deletions.
`pkgdown::check_pkgdown()`: `✔ No problems found.` (`logs/check_pkgdown.txt`).

## 6. Article deviations; script notes

**Article deviations: none.** The article is committed byte-identical to the input.

Script notes:
- `s1_reference_gate.R` evaluates only the `LOWER_METHODS`/`UPPER_METHODS` assignments from chunk
  `summaries` and the body of chunk `reference-data` (minus `#|` lines) in a fresh environment,
  and evaluates only the `GH52_TABLE_7 <- …` assignment from `quarto/GuoHe/guohe_sec52_sim.R`
  (nothing else in that file is run). Printed precision is the task's per-column `d`; a
  both-`NA` pair (specificity, PPV at β₂ = 0) passes, a one-sided `NA` fails.
- `s2_text_gates.sh` first ran with a regex lacking the dot-all flag; its own guard
  (`stopifnot(!grepl("<pre", …))`) stopped it before `article_output.txt` was written, so that
  run's counts (all 0) were void. Fixed with `(?s)`; the recorded log is the corrected run.
- v1's §5 counts were taken by an inline command (whole page, task regexes verbatim), recorded
  in `keyvalues.txt`; `article.txt` was regenerated identically by `s2_text_gates.sh`.

## 7. OPEN ITEMS (not stops)

1. The article's NEWS line waits for the 0.4.1 release.
2. Tracking of the rendered `vignettes/articles/guohe-nested-continuum.html` — Larry's decision.
   Note: the sibling articles' HTML (`survival-post-selection.html`, `methodology.html`, …) is
   tracked, so leaving this one untracked is a departure from that pattern if it persists.
3. Surfaced by the render: in the truth table the β₂ = 0 row's estimated β(30) is −0.0140 with
   anchor deviation/SE −1.48 (the basis is the exact zero, so the estimate is diagnostic only);
   unremarkable, recorded for completeness.
4. v1's §5 gates were mis-specified for a page that echoes its code; any later article task
   with `code-fold: false` should count on output-only text from the start.

## 8. Verdict

1. Provenance gate green; branch and ancestry as stated; tree clean at BASE_SHA `f2ddb5c0`.
2. Input md5 matched; article placed and committed unchanged.
3. Reference gate green: 504/504 literals match the S4 payload and GH52_TABLE_7 at printed precision.
4. No literal corrected; payload meta agrees with the article's stated settings.
5. One render, exit 0, 5:24 on 2 workers (cap 40 min); no deviations; `stopifnot` held.
6. v1 stopped at §5 on four echo-inflated counts; v2 (output only) recorded and committed alone.
7. All v2 §5 gates green; v1−v2 differences exactly +1 on gates 1, 2b, 2c, 5c.
8. `_pkgdown.yml`: 3 insertions, `check_pkgdown()` clean.
9. No package code touched; the rendered HTML is untracked and unstaged.
10. No push, no merge. Ready for Larry's merge decision.

## Closing block

Commit 3 (record + report): `0d437f31570d67cf43c06c38d067277113b3c1a7`. This closing block is
added in a report-only follow-up commit (the Stage 0 precedent, `f2ddb5c0`), since it must
print the post-conditions after the last commit.

```
$ git diff --name-only c2967362 HEAD | grep -vE '^(dev/|vignettes/articles/guohe-nested-continuum\.qmd$|_pkgdown\.yml$)'
(nothing)
$ git status --porcelain | grep -v '^?? '
(nothing)
$ git status --porcelain | grep '^?? '
?? vignettes/articles/guohe-nested-continuum.html
```

Copied to `~/Downloads/guohe_article_2026-10-30/`: this report, `results/`, and the rendered HTML.
