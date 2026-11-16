ARTICLES INSTALLED AND RENDER

# REPORT — MRCT articles: aligned with the manuscript, the two region-interaction configurations (2026-11-16)

Session `fs-mrct-articles-align`. Host `pop-os`, `~/Documents/GitHub/forestsearch`, branch `master`.
Task: `dev/tasks/TASK_mrct_articles_manuscript_alignment_2026-11-16.md`.
Start SHA `a793c4040ae61e17e765e2c2d8a48ad8918bd174`.

**Summary.** The two delivered `.qmd` files are installed over the two MRCT articles **unchanged**
(CC made no change to them under rules 1-3). Both render with forestsearch 0.5.0 from CRAN in the
scratch library: simulations 237.5 s, application 194.6 s, exit 0, every `stopifnot()` passing.
Gates 4.1-4.3 hold (4.1 with the reading stated below). The live rows (1)-(7) of the manuscript's
Table 2 differ from its printed values by at most 4.82e-4.

## Step 0 — start

- Task file copied to `dev/tasks/` (`cmp` identical) and committed: `a4f936a1`.
- Gate: branch `master`; `git status --short` empty; at `HEAD` the two articles had md5
  `71aa9b658466713efcd58c370427eeed` (application) and `2447dd482e9f0ccc5320bb5f6a0c9925`
  (simulations); the two files in `~/Downloads` had `e1d0ac3edd31d5e01ba2f130ee4758dd` (1,060 lines)
  and `d409ec28b872de31b0fe4f586592c4a4` (1,617 lines). **Holds.**

## Step 1 — what the delivered files assume (read-only)

Sources: forestsearch tag `v0.5.0` (`git show v0.5.0:R/...`); mrct-consistency branch `main`
(`be7be00`) read with `git show main:...` only. The checkout there is on
`docs/workstream-index-2026-11-14`; `main` is an ancestor of it and `git diff --stat main HEAD --
docs/manuscript/sim sims results` is empty. No fetch, no checkout, no change.

### 1.1 The build

`mrct_build_dgm(data, ..., m = 1)` (`R/mrct_minor_helpers.R:120`) calls
`inject_mrct_structure(data, ...)` and then `scale_dgm_hazard(dgm, m = m)` when `m != 1`.
`inject_mrct_structure()` has the formals `region_loghr_shift = 0`, `region_harm_loghr_shift = 0`
and `latent_R = NULL` (`R/mrct_inject_dgm.R:371-388`), so both shifts and `latent_R` reach it
through `...`. **Finding: yes, `mrct_build_dgm()` passes both shifts on; the delivered `ri-build`
call (`build_cell(h = 2, region_loghr_shift =, region_harm_loghr_shift =)`, with `build_cell()`
forwarding `...` to `mrct_build_dgm()`) needs no direct `inject_mrct_structure()` /
`scale_dgm_hazard()` call.**

The shift (`R/mrct_inject_dgm.R:558-566`) is applied after the AFT fit to `theta_1`, `loghr_po`,
`lin_pred_1` and the treated rows of `lin_pred_obs` only; `theta_0`, the region, `flag_harm` and
`model_params` are not touched.

### 1.2 The analytic call

- `ahr_limit(pop, set, b, h, b_R, h_R)` and `ahr_design(pop, b, h, n, b_R, h_R)`: with `b` and
  `h` given, the per-patient log-HR is the four-cell value (`.ahr_psi()`, `R/ahr_analytic.R:283-303`:
  log b / log h outside the region, log b_R / log h_R inside it, by H*). `ahr_design()` always
  takes the cached route (`force_fast = TRUE`); `ahr_limit()` takes it when
  0.2 <= b, h, b_R, h_R <= 30 (here 0.55, 2, 1, 1 / 0.55, 2, 1, 2).
- The cached route (`.ahr_fast()`, `.ahr_cache_grid()`, `.ahr_cache_lt()`,
  `R/ahr_analytic.R:382-426`) reads `pop$theta0`, `pop$cell` (H* x region), `pop$group` (entry),
  `pop$mu`, `pop$tau` and the four multipliers. It never reads `pop$psi` (the DGM's `loghr_po`),
  and the grid depends on `range(theta0)` only.
- The manuscript's code: `ext-predictions.R` (S2_E1) and `ext-s2null.R` (S2_E2) call
  `ahr_design(pop_S2, b, h, n, sets = ..., b_R =, h_R =)` on **S2's base population**; both also
  check that the per-row route on a build with the two shifts agrees with the four-cell route
  (C4, to 1e-8). `sims/region_interaction_cells_design.R` builds the campaign cells with the two
  shifts (gamma_R = log(b_R/b), gamma_H*R = log(h_R/h) - log(b_R/b)) and computes truths with
  `ahr_limit(pop, set)` per row on the cell's own DGM.
- The delivered `ri-build` calls `ahr_design(pop_ri[[k]], b = 0.55, h = 2, n = 1000, b_R =, h_R =)`
  and `ahr_limit(pop_ri[[k]], "region", b = 0.55, h = 2, b_R =, h_R =)` on the **shifted build's**
  population. Since the shift leaves `theta_0`, region, H*, `mu` and `tau` as S2's (1.1) and the
  cached route does not read `psi`, this is the same arithmetic as the manuscript's call on S2's
  population. **Finding: the article's call is the same computation; no change.**
- The delivered `rows7()` forms row (3) as (1) - (4) and row (6) as (4) - (5); the manuscript's
  records compute them as their own orthants (`ext-s2null.R`: `p_A_notB`, `p_AB_notC`), and the
  manuscript's `dgm-predictions` chunk asserts (3)+(4)=(1) and (5)+(6)=(4) at full precision to
  <= 1e-12. Rows (1), (2), (4), (5) are `ahr_design()$prob`'s `p_itt_ub_lt1`, `p_region_itt_gt1`,
  `p_trigger`, `p_trigger_and_Gstar_lt1`: the quantities `ext-predictions.R` reads from
  `ahr_design()$prob` (`ext-s2null.R` rebuilds the same orthants from `ahr_design()`'s sets and
  exact correlations and reproduces the S2_E1 rows to 1e-12, its C2).

### 1.3 The stated values

Compared with the manuscript at `main`: the kept `docs/manuscript/sim/sim-manuscript.tex`
(rendered with the qmd's last commit, `041b0e7`) and the committed `sim-manuscript.pdf`
(`be7be00`, `pdftotext -layout`), which agree; Table 2 is `tab:predictions`, Table 6 is `tab:ri`.

| Delivered value | Manuscript | Records | Result |
|---|---|---|---|
| `ms2` (Table 2 rows (1)-(7), 4 columns, 28 values) | tex lines 545-551; PDF Table 2 | rows (1)-(6) are the qmd's `EXP6` gate strings | all 28 equal |
| `0.057`, `0.115` (row (8), S2, null / S2, detrimental) | tex line 553; PDF Table 2 | `results/region_interaction_cells/summary/core_numeric.csv`: 15.688% x 3,627 = 569 / 10,000; 23.504% x 4,880 = 1,147 / 10,000 | equal |
| row (9) check `c(0.725, 0.500, 0.353, 0.950)` | tex line 554 | | equal |
| `t6` (Table 6, 9 rows x 3 columns, 27 values) | tex lines 1622-1630; PDF Table 6 | `summary/core.md`, `summary/mr_field_core.md` | all 27 equal |
| seed `615400000` (S2, null) | | `build.csv`, `S2_E2_rmlog4/timing.csv`: S2_E2_rmlog4, b_R = 1, h_R = 1, seed 615400000 | equal |
| seed `615500000` (S2, detrimental) | | `build.csv`, `S2_E1_rmlog4/timing.csv`: S2_E1_rmlog4, b_R = 1, h_R = 2, seed 615500000 | equal |
| "its regional estimate sits in the lower quarter of S2's triggered trials" | §4.1 (PDF line 773) | | same wording |
| "it lies above 78% of them" | §4.2: "Among the 2,397 declared triggered trials of S2 it lies above the bounds of 78% of them" | | equal |
| "more than 60 patients and at least 10 events on each arm" | §5.1 | | same wording |
| "the mean, median and quartiles for the other continuous covariates" | §5.1: "at the mean, median and quartiles otherwise" | | equal |

Which configuration is which, from the records (not assumed): **S2_E2_rmlog4 = "S2, null"
(b_R = 1, h_R = 1, seed 615400000); S2_E1_rmlog4 = "S2, detrimental" (b_R = 1, h_R = 2, seed
615500000)** — `build.csv`, the two `timing.csv`, `sims/region_interaction_cells_design.R`, and the
manuscript's column heads (`P9_cols`: S2_E2 -> "S2, null", S2_E1 -> "S2, detrimental").
Section numbers cited by the article (4.1, 4.2, 5.1) are the manuscript's PDF headings.

**Finding: no stated value differs from the manuscript; nothing changed under rule 2.**

## Changes CC made to the delivered files

**None.** No change under rule 1 (the new code ran as written), rule 2 (no stated value differs,
Step 1.3) or rule 3 (no break: every cross-article anchor exists in its target, Step 4.2).

```
$ diff ~/Downloads/mrct-consistency_aligned_2026-11-16.qmd <(git show 2356d073:vignettes/articles/mrct-consistency.qmd)
$ diff ~/Downloads/mrct-consistency-simulations_aligned_2026-11-16.qmd <(git show 2356d073:vignettes/articles/mrct-consistency-simulations.qmd)
$
```

Both diffs are empty; the committed articles have md5 `e1d0ac3edd31d5e01ba2f130ee4758dd` and
`d409ec28b872de31b0fe4f586592c4a4`.

## Step 2 — the CRAN release in the scratch library

`R_LIBS=/tmp/mrctalign/lib Rscript` with `install.packages("forestsearch", lib = "/tmp/mrctalign/lib",
repos = "https://cloud.r-project.org")`: downloaded `forestsearch_0.5.0.tar.gz` (1.7 MB) from CRAN,
built from source, 21.3 s. Installed: Version 0.5.0, Repository CRAN, Date/Publication
2026-11-15 08:50:02 UTC. Only `forestsearch` was installed into `/tmp/mrctalign/lib`; every
dependency was already present in the user library and nothing was installed there or in the site
library.

Load check (`R_LIBS=/tmp/mrctalign/lib`, `future::plan(multisession, workers = 2)`):

| Process | pid | forestsearch | loaded from |
|---|---|---|---|
| main | 2245034 | 0.5.0 | `/tmp/mrctalign/lib` |
| multisession worker | 2245085 | 0.5.0 | `/tmp/mrctalign/lib` |

(A first try with `workers = 1` returned the main pid: future runs a one-worker `multisession`
plan in the main session, so it was repeated with two workers.) The user library also holds a
forestsearch 0.5.0 (present before this task, not touched); `/tmp/mrctalign/lib` is first on
`.libPaths()` in both processes.

## Step 3 — the render

The two delivered files were copied over the two articles, and from there to
`/tmp/mrctalign/render/`. Each article was rendered alone (`quarto render <file>.qmd --to html`,
Quarto 1.10.18, R 4.6.2, `R_LIBS=/tmp/mrctalign/lib`), simulations first. No `.html` was written
under `vignettes/`. The records were read from the installed package (the render folder is not a
source tree).

| Article | Wall clock | Exit | Article's own run-time line |
|---|---|---|---|
| `mrct-consistency-simulations.qmd` | 237.5 s (13:09:08-13:13:05) | 0 | "took 3.8 minutes to run" |
| `mrct-consistency.qmd` | 194.6 s (13:13:15-13:16:30) | 0 | "took 2.9 minutes to run" |

Both renders ran to the end, so every `stopifnot()` in both articles passed. Neither render log
has an error or warning line. The simulations page shows three "Ran out of iterations and did not
converge" messages in chunk output (the censoring model's fit, as in earlier site builds).

**Gate 3.4: holds.**

**A third render (observer only).** The article does not print the difference that Step 4.4 asks
for. The simulations article was therefore rendered once more, from an identical copy
(`cmp`), in `/tmp/mrctalign/render/capture/`, with `R_PROFILE_USER` pointing to a profile whose
only content is an exit finalizer that saves `t2`, `ms2`, `reg_ri`, `decl`, `t6`, `live_S2` and
`t_ri` from the knit's global environment to an `.rds` (Quarto runs `Rscript rmd.R` without init
flags and knits in `globalenv()`). No source changed. Exit 0, 235.0 s. Its page differs from the
first render's in 7 text lines, all timings (build seconds, section seconds, the `ri-build`
seconds 1.0 vs 1.1). The pages copied to `~/Downloads` are the first render's.

## Step 4 — post-conditions

### 4.1 No unresolved cross-reference (gate)

Visible text of each page (`<script>`/`<style>` removed), counted with and without the code
listings (`<pre>` blocks):

| Page | literal | in the page | outside code listings |
|---|---|---|---|
| application | `?@` | 0 | 0 |
| application | `@sec-` | 2 | **0** |
| application | `@fig-` | 0 | 0 |
| simulations | `?@` | 0 | 0 |
| simulations | `@sec-` | 1 | **0** |
| simulations | `@fig-` | 0 | 0 |

The raw HTML holds one `?@` per page, inside Quarto's embedded JavaScript (a character class of a
regular expression), not in the text. The three `@sec-` occurrences are in the code listings,
on lines that are unchanged from the base articles:

- application, chunk code `caption = "... counterparts of these rows are in @sec-case-campaign and in
  the simulations article."` — the caption renders as "... are in Section 2.6 and in the simulations
  article";
- application, chunk code `caption = "S2 at n = 1,000, ...; rows as in @sec-predictions, ..."` — renders
  "rows as in Section 1.3";
- simulations, R comment `# the simulation route against the case study of @sec-case-objects`.

The two captions are resolved references; the comment is code text, not a reference. Reading
applied: the gate concerns cross-references in the rendered text, and none is unresolved.
**Holds** (with that reading; the literal counts above include the code listings for the record).

### 4.2 The anchor and the three links (gate)

- Simulations page: one element with `id="sec-region-interaction"` (rendered as section 4.5,
  "The region-interaction configurations").
- Application page: three `<a href="mrct-consistency-simulations.html#sec-region-interaction">`
  (sources lines 175, 294, 850); no other cross-article `.html#` link on that page.
- The simulations source's links to the application article (`#sec-dgm-overview` x2,
  `#sec-window` x3) point at anchors present in the application source.

**Holds.**

### 4.3 Terms in the two sources (gate)

| Term | `mrct-consistency.qmd` | `mrct-consistency-simulations.qmd` |
|---|---|---|
| `oracle joint` | 0 | 0 |
| `P(trigger)` | 0 | 0 |
| `the report` | 0 | 0 |
| `report's` | 0 | 0 |

Counted with `grep -o -F`, case-sensitive and case-insensitive (same counts). **Holds.**

### 4.4 From the simulations render

Output of chunk `ri-build` (the two regional limits), as rendered:

```
       region margHR(40): S2, null region margHR(40): S2, detrimental 
                          1.000000                           1.193099 
                           seconds 
                          1.000000 
```

The text below it: "The region's population Cox limit, margHR(40), is 1.000 under S2, null, and
1.193 under S2, detrimental, against 0.868 in S2."

The rendered Table 2 ("The manuscript's Table 2"):

| | S2 | S2, null (b_R = 1.0, h_R = 1.0) | S2, detrimental (b_R = 1.0, h_R = 2.0) | N0 |
|---|---:|---:|---:|---:|
| (1) P(ITT ub < 1.0): the trial is positive | 0.844 | 0.804 | 0.732 | 0.999 |
| (2) P(region ITT > 1.0) | 0.326 | 0.500 | 0.727 | 0.055 |
| (3) P(consistency): positive, region ITT ≤ 1.0 | 0.602 | 0.443 | 0.238 | 0.944 |
| (4) P(non-consistency): the trigger | 0.241 | 0.361 | 0.494 | 0.054 |
| (5) P(oracle consistency) | 0.187 | 0.070 | 0.185 | 0.017 |
| (6) P(oracle non-consistency) | 0.054 | 0.291 | 0.309 | 0.038 |
| (7) P(consistency or oracle consistency) | 0.790 | 0.513 | 0.423 | 0.961 |
| (8) P̂(declared consistency) | 0.123 | 0.057 | 0.115 | 0.006 |
| (9) P̂(consistency or declared consistency) | 0.725 | 0.500 | 0.353 | 0.950 |

The rendered Table 6 ("The manuscript's Table 6"):

| | % of | S2 | S2, null (b_R = 1.0, h_R = 1.0) | S2, detrimental (b_R = 1.0, h_R = 2.0) |
|---|---|---:|---:|---:|
| trial power (ITT ub < 1.0) | all trials | 84.0 | 80.3 | 73.1 |
| non-consistency: trials reaching the trigger, % (count) | all trials | 24.0 (2,398) | 36.3 (3,627) | 48.8 (4,880) |
| subgroup declared, % (count) | triggered trials | 100.0 (2,397) | 100.0 (3,626) | 100.0 (4,879) |
| declared, region(Ĝ) < 1.0 | triggered trials | 51.2 | 15.7 | 23.5 |
| region(G*) < 1.0 | triggered trials | 77.7 | 19.4 | 37.3 |
| declared, region(Ĝ) ≥ 1.0 | triggered trials | 48.7 | 84.3 | 76.5 |
| MR lower bound on Ĥ > 0.90 | triggered trials | 18.7 | 18.9 | 18.5 |
| MR lower bound on Ĥ > 1.0 | triggered trials | 13.6 | 13.8 | 13.5 |
| MR upper bound on Ĥc < 1.0 | triggered trials | 100.0 | 100.0 | 100.0 |

Both tables equal the manuscript's printed Tables 2 and 6 in every cell (Step 1.3).

**Largest absolute difference, computed rows (1)-(7) minus the manuscript's printed values:
4.82e-4**, at row (3), S2 (computed 0.602482, printed 0.602). From the observer render's `t2`
and `ms2`; the article's own check is `max(abs(t2 - ms2)) < 6e-4`. Every one of the 28
differences is below 5e-4 in absolute value, so each computed value rounds to the printed one.
Per column, the largest: S2 4.82e-4, S2, null 2.77e-4 (row (6)), S2, detrimental 3.92e-4
(row (4)), N0 4.36e-4 (row (6)). The full-precision computed values:

| row | S2 | S2, null | S2, detrimental | N0 |
|---|---:|---:|---:|---:|
| (1) | 0.843721 | 0.804107 | 0.731688 | 0.998695 |
| (2) | 0.326120 | 0.500000 | 0.727339 | 0.054819 |
| (3) | 0.602482 | 0.442833 | 0.238080 | 0.944274 |
| (4) | 0.241239 | 0.361274 | 0.493608 | 0.054420 |
| (5) | 0.187410 | 0.069997 | 0.184777 | 0.016857 |
| (6) | 0.053829 | 0.291277 | 0.308831 | 0.037564 |
| (7) | 0.789892 | 0.512829 | 0.422857 | 0.961131 |

Also from that render: the S2 column of Table 6 against the shipped records, largest absolute
difference 0.0490 (declared, region(Ĝ) ≥ 1.0: 48.749 vs 48.7), within the chunk's 0.051.

### 4.5 The released package

`git diff --stat v0.5.0 HEAD -- R/ inst/ DESCRIPTION NAMESPACE` (HEAD `2356d073`, after the
articles commit): **empty**.

### 4.6 Render time against the rendering note

The note says "about five minutes". Simulations 237.5 s (4.0 min), application 194.6 s (3.2 min):
neither exceeds twice five minutes. The note was not edited.

## Commits (this session)

| SHA | Subject | Paths |
|---|---|---|
| `a4f936a1` | task: mrct-articles-manuscript-alignment 2026-11-16 (as received) | `dev/tasks/TASK_mrct_articles_manuscript_alignment_2026-11-16.md` |
| `2356d073` | docs: MRCT articles aligned with the manuscript; the two region-interaction configurations | `vignettes/articles/mrct-consistency.qmd`, `vignettes/articles/mrct-consistency-simulations.qmd` |

The report's own commits follow (see Final checks).

## Files outside the repository

- `~/Downloads/mrct-consistency_aligned_2026-11-16.html` (2,237,035 bytes) and
  `~/Downloads/mrct-consistency-simulations_aligned_2026-11-16.html` (1,730,296 bytes): the first
  renders, `cmp` identical to `/tmp/mrctalign/render/`.
- `/tmp/mrctalign/render/`: the two sources, the two first renders, and `capture/` (the observer
  render). Kept.
- `/tmp/mrctalign/logs/`: install, load-check and render logs, the capture profile, `sims_capture.rds`.
