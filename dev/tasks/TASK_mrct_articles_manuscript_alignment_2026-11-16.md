# TASK — MRCT articles: align with the manuscript and add the two region-interaction configurations (2026-11-16)

Session name: `fs-mrct-articles-align`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**.
Larry's decisions of 2026-11-16: the manuscript in the mrct-consistency repository is the
canonical source of ideas and terminology; the two MRCT articles are aligned with it, without
changing the articles much; the manuscript's two extended configurations ("S2, null" and
"S2, detrimental") are added; the wording is "false-positive increase in meeting consistency",
not "type-1 error". Read this whole file first.

## What chat delivers

Chat wrote the two revised articles. They are in `~/Downloads`:

| File in `~/Downloads` | Replaces | Lines | md5 |
|---|---|---|---|
| `mrct-consistency_aligned_2026-11-16.qmd` | `vignettes/articles/mrct-consistency.qmd` | 1,060 | `e1d0ac3edd31d5e01ba2f130ee4758dd` |
| `mrct-consistency-simulations_aligned_2026-11-16.qmd` | `vignettes/articles/mrct-consistency-simulations.qmd` | 1,617 | `d409ec28b872de31b0fe4f586592c4a4` |

They were made from the articles at `master` whose md5 are `71aa9b658466713efcd58c370427eeed`
and `2447dd482e9f0ccc5320bb5f6a0c9925`. Chat could not run R: the new code is untested, and
this task is where it is run.

## Rules

- **No `R/` change, no `inst/` change, no test change, no `_pkgdown.yml` change.** Only the two
  articles change (plus this task file and the report). forestsearch 0.5.0 already has what
  the new code uses (`REPORT_mrct_type1_scenarios_trace_2026-11-16.md`, Step 5).
- **Lean run.** No agents. One short report. No `R CMD check`, no test suite, no site build.
  The only compute is the render of the two articles.
- The render uses **the CRAN release**: forestsearch 0.5.0 installed from CRAN into a scratch
  library, `/tmp/mrctalign/lib`. Nothing is installed into the user or site library.
- Rendered `.html` never goes into `vignettes/`. Render output goes to `/tmp/mrctalign/`.
- Work on `master`. Commits are path-limited. **No push, pull, merge or tag.**
- In mrct-consistency: read only. No fetch, no checkout, no change.
- Comparisons of numbers use a tolerance. Never tighten a check to exact equality, and never
  loosen one to make it pass.
- Gates stop on failure; they never stop to ask.

## What CC may change in the two delivered files

1. **The new code** (chunks `ri-build`, `ri-table2`, `ri-table6`, the two added lines of the
   `campaign-code` chunk, and the `...` added to `build_cell()`): adapt it to the package's
   actual interface where it does not run as written. The intent is fixed: build S2 with
   `region_loghr_shift` and `region_harm_loghr_shift`; compute the analytic probabilities with
   `b_R` and `h_R`; check them against the manuscript's printed values.
2. **A stated value** that differs from the manuscript (Step 1.3): correct it to the
   manuscript's printed value.
3. **A break caused by chat's edits elsewhere** (a cross-reference that no longer resolves, a
   renamed label used by other code): the smallest fix.

Nothing else: no prose rewording, no restructuring, no change to an existing check. Every
change made under 1 to 3 is shown in the report as a diff.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate (stop):** branch `master`; no modified tracked files; the two articles at `HEAD`
   have the two base md5 above; the two files in `~/Downloads` have the md5 in the table.

## Step 1 — check what the delivered files assume (read-only)

1. **The build.** From the 0.5.0 source: does `mrct_build_dgm()` pass `region_loghr_shift` and
   `region_harm_loghr_shift` on to `inject_mrct_structure()`? If not, the new chunk must call
   `inject_mrct_structure()` (with `latent_R`) and `scale_dgm_hazard()` itself.
2. **The analytic call.** How do `ahr_limit()` and `ahr_design()` treat `b_R` and `h_R` for a
   mechanism built with the two shifts? Compare with how the manuscript's own code makes the
   call (mrct-consistency `docs/manuscript/sim/ext-predictions.R`, `ext-s2null.R`,
   `sims/region_interaction_cells_design.R`). The article's call must be the same computation.
3. **The stated values.** Compare with the manuscript at mrct-consistency `main`
   (`docs/manuscript/sim/sim-manuscript.qmd` and the records it reads:
   `docs/manuscript/sim/analytic-ext/`, `results/region_interaction_cells/summary/`,
   `results/mr_field_campaign/S2_E2/`, `S2_E1/`):
   - `ms2` in chunk `ri-table2`: Table 2, rows (1)-(7), four columns, three decimals;
   - `0.057` and `0.115` in the same chunk: Table 2 row (8) for the two configurations;
   - `t6` in chunk `ri-table6`: Table 6, all nine rows of the three columns;
   - the seeds `615400000` (S2, null) and `615500000` (S2, detrimental) in `campaign-code`;
   - in the application article: "its regional estimate sits in the lower quarter of S2's
     triggered trials" and "it lies above 78% of them" (manuscript Sections 4.1 and 4.2);
     "more than 60 patients and at least 10 events on each arm", "the mean, median and
     quartiles for the other continuous covariates" (Section 5.1).
   Which configuration is S2_E2 and which S2_E1 is to be confirmed from the records, not
   assumed. A value that differs is corrected to the manuscript's printed value (rule 2).

## Step 2 — install the CRAN release in the scratch library

`install.packages("forestsearch", lib = "/tmp/mrctalign/lib",
repos = "https://cloud.r-project.org")`. Confirm that the main process and one `multisession`
worker load forestsearch 0.5.0 from `/tmp/mrctalign/lib` (with `R_LIBS=/tmp/mrctalign/lib`).
If CRAN cannot be reached, build the scratch library from the tag instead
(`git archive v0.5.0`) and say so.

## Step 3 — put the files in place and render

1. Copy the two delivered files over the two articles. **Do not commit yet.**
2. Render each article alone with Quarto, with `R_LIBS=/tmp/mrctalign/lib`, in a copy under
   `/tmp/mrctalign/render/` so that no `.html` is written into `vignettes/`. Record the wall
   clock of each render.
3. Every `stopifnot()` in both articles must pass. If the render fails, fix within the limits
   above and render again.
4. **Gate (stop):** both renders exit 0. If that cannot be reached within the limits above,
   restore both articles (`git checkout -- vignettes/articles/mrct-consistency.qmd
   vignettes/articles/mrct-consistency-simulations.qmd`), write the report with first line
   `NOT INSTALLED: <the failing chunk and its error>`, and close. `master` never holds an
   article that does not render.

## Step 4 — post-conditions (all reported; the first three are gates)

1. In both rendered pages: no unresolved cross-reference (no literal `?@` and no literal
   `@sec-` or `@fig-` in the text).
2. The simulations page has the anchor `sec-region-interaction`; the application page links
   to `mrct-consistency-simulations.html#sec-region-interaction` (three links).
3. In the two sources: `oracle joint`, `P(trigger)`, `the report` and `report's` each occur
   0 times.
4. From the simulations render, quote: the output of chunk `ri-build` (the two regional
   limits); the rendered Table 2 and Table 6; the largest absolute difference between the
   computed rows (1)-(7) and the manuscript's printed values.
5. `git diff --stat v0.5.0 HEAD -- R/ inst/ DESCRIPTION NAMESPACE`: empty, so the released
   package is the one these articles describe.
6. If either render took more than twice the "about five minutes" its rendering note states,
   say so. Do not edit the note.

## Step 5 — commit, report and close

1. Commit the two articles (one commit, path-limited).
2. Copy the two rendered pages to `~/Downloads` as
   `mrct-consistency_aligned_2026-11-16.html` and
   `mrct-consistency-simulations_aligned_2026-11-16.html`, for Larry to read before he pushes.
3. Write `dev/reports/REPORT_mrct_articles_manuscript_alignment_2026-11-16.md`: first line
   `ARTICLES INSTALLED AND RENDER` or `NOT INSTALLED: ...`; Step 1's three findings; every
   change CC made to the delivered files, as a diff; render times; Step 4's results; this
   session's commits. Facts only.
4. Commit the report. Copy it to `~/Downloads`.
5. **Final checks** (report the result): this session's commits touch only the two articles,
   this task file and the report; no modified tracked files; no `.html` under `vignettes/`;
   nothing changed in mrct-consistency; nothing was pushed, pulled, merged or tagged.
6. Delete `/tmp/mrctalign/lib`. Keep `/tmp/mrctalign/render/`.

## For the record: what chat changed

**Added**

- Simulations article: a section "The region-interaction configurations"
  (`#sec-region-interaction`), before the campaign code. It builds S2, null and S2,
  detrimental, computes their regional limits and rows (1)-(7) of the manuscript's Table 2
  live, checks them against the manuscript's printed values, and shows rows (8)-(9) and the
  manuscript's Table 6 with the two configurations' campaign values as the manuscript states
  them (their records are not shipped with 0.5.0). The S2 column of Table 6 is recounted from
  the shipped records. Two lines in the campaign code; a pointer under the cells table.
- Application article: the two configurations defined after the extension; three pointers to
  the new section.

**Terminology and notation, both articles**

"oracle joint" to "oracle consistency"; "P(trigger)" to "P(non-consistency)"; "predicted" to
"analytic" and "frequencies" to "observed proportions" where campaign results are compared;
"population Cox limit" introduced as the name, "windowed Cox limit" kept as its gloss; "the
report" to "the manuscript", with section numbers; "conditional on the trigger" to "among
triggered trials"; "signal-generating, not confirmatory"; the manuscript's event notation in
place of `A`, `B`, `C*`; eta for the prognostic coefficients; upright R for the region; the
manuscript's proposition numbers in the simulations appendix; "a stress setting" for
`r = -log 4`; the identifier's steps under the manuscript's names (eligible, effect screen,
consistency screen, admissible family, selection and declaration).

**Larry's passages of 2026-11-09 that the manuscript's wording replaced** (application article)

| Passage | Was | Now |
|---|---|---|
| Overview, the question | "...the regional result is suggestive/compatible with a trial-wide subgroup. We therefore cast the MRCT regional consistency criteria in terms of the identified subgroup ... as if the subgroup was pre-specified. ... The strategy is signal generation, not confirmation." | The manuscript's Section 2 statement, keeping the point that the subgroup is not pre-specified while the procedure is |
| Overview, the list | "We consider two mechanisms for the trigger" | "Three explanations fit a triggered trial": noise, a region-specific effect, a trial-wide subgroup |
| (M2) | "(M2) is then gamma_R = gamma_H*R = 0" | "(M2) is then the testable null ..." |
| (M4) | "... per the observed configuration of ACTG 175 (we take patient's observed covariates as-is)" | Removed: the population is a copula expansion with a planted region, so "as-is" does not hold |
| The extension | "Writing the effect with treatment-by-region interactions extends (M2):" | "... makes (M2) a testable restriction rather than a built-in one:" |
| The pipeline, last paragraph | "The hold-out removes the subgroup selection ... While the analysis is post-hoc, we propose a pre-specified identifier algorithm ..." | The manuscript's Section 3.3 statement on the held-out region, and "the analysis is post hoc, but the pipeline is fixed" |

Kept as Larry wrote them: the sentence on why ACTG 175 is used, and "the regional consistency
failure as it presents in oncology multiregional clinical trials".

**Not done**

- The six-cell results table does not gain the manuscript's rows 6 and 7 (the identified
  complement's share of the region): they come from per-trial records the package does not
  ship.
- "The standard pathway" is kept: the manuscript's Table 1 uses it (its text also says
  "traditional pathway").
- The regulatory background and the model-based alternative are not brought in.
