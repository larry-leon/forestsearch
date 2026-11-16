NEVER IN THE ARTICLES: results exist in mrct-consistency `results/region_interaction_cells/` (commit 0108843, cells S2_E2_rmlog4 and S2_E1_rmlog4), reported in `docs/manuscript/sim/sim-manuscript.qmd` (§2 Table 2, §5.2 Table 6)

# REPORT — MRCT articles: where are the type-1 error scenarios of the extended model? (2026-11-16)

Task: `dev/tasks/TASK_mrct_type1_scenarios_trace_2026-11-16.md` (committed 20f204b1). Session
`fs-mrct-type1-trace`, host pop-os, forestsearch `master`. Read-only: no render, simulation,
install, check, test suite or agent. Facts only.

Search terms as in the task. Case-insensitive for all except `gamma_R`, `b_R`, `h_R`, `type I` and
`type-I`, where case carries the meaning (a case-insensitive `type I` matches "type is", "type int").
`b_R` and `h_R` were matched as whole tokens. Because `region_treat_loghr` turned out not to be the
argument that builds the extended model (Step 5), `region_loghr_shift` and `region_harm_loghr_shift`
were searched as well.

## Step 0 — start

- `master` before this session: `d059e054a70f548377c539b5fa41b263c92533d8`. Task commit:
  `20f204b19fc9154c7358d8080d4cd7cc8ae909ab` (path-limited, one file).
- `git status --short` (before and after the task commit): empty.

| File | Lines (work / HEAD) | md5 (work = HEAD) | Equals the table in the task |
|---|---|---|---|
| `vignettes/articles/mrct-consistency.qmd` | 1,001 / 1,001 | `71aa9b658466713efcd58c370427eeed` | yes |
| `vignettes/articles/mrct-consistency-simulations.qmd` | 1,428 / 1,428 | `2447dd482e9f0ccc5320bb5f6a0c9925` | yes |

## Step 1 — the history of the two articles

Refs searched (`--all`): `master`, `origin/master`, `origin/HEAD`, `origin/gh-pages`, tags v0.1.0,
v0.4.0, v0.5.0, `refs/stash`. `--reflog` adds no commit touching either article. `git fsck
--connectivity-only --unreachable` lists 62 unreachable commits (2023-03-29 to 2026-10-25); none
touches or holds either article or `inst/extdata/mrct-consistency/`.

**1.1 Commits.** Neither article had an earlier path: `mrct-consistency.qmd` starts as a new file in
43441a42; `--follow` on the simulations article reports dad519d6 as a rename from
`mrct-consistency.qmd` (+35 −479); with `--no-renames` it is a new file (+1,452 −0).

`vignettes/articles/mrct-consistency.qmd` (all on master; author Larry Leon throughout):

| SHA | Date | + | − | Subject (shortened) |
|---|---|---|---|---|
| 43441a42 | 2026-11-01 | 1025 | 0 | docs: article "MRCT regional consistency via subgroup identification" |
| 450b5492 | 2026-11-01 | 244 | 148 | docs: article r2 — the reader's path, "Finding a case study", Appendix equations |
| 0b1f5ba7 | 2026-11-01 | 1 | 1 | fix: MRCT series check — ASCII, … |
| 42f36b87 | 2026-11-02 | 532 | 161 | docs: article r3 — replicate 5 as the walkthrough; own-trial simulation setup |
| 172140b8 | 2026-11-02 | 57 | 50 | docs: article r4 — own-trial section built with inject_mrct_structure() observed-region mode |
| fa3b1397 | 2026-11-02 | 40 | 0 | article: install notes and early check |
| 8323c5d6 | 2026-11-02 | 18 | 5 | article: interpretive edits from review |
| af3f6329 | 2026-11-02 | 49 | 15 | article: tested install sequence, rendering note, platform tolerance |
| 98d5a801 | 2026-11-02 | 69 | 14 | article: records from the installed package; … §2.6 the S2 campaign in brief |
| a597d673 | 2026-11-04 | 186 | 59 | article: mrct-consistency appendix aligned with manuscript §3.5 |
| 6b81622b | 2026-11-05 | 99 | 0 | docs: Overview gains "The data-generating mechanism" and "Analytic predictions before the analysis" |
| 8a1ba513 | 2026-11-05 | 75 | 46 | docs: "The data-generating mechanism" carries the mechanism's mathematics (… the region-dependent extension with gamma_R and gamma_H*R …) |
| dad519d6 | 2026-11-05 | 22 | 1261 | docs: split into two independently executable articles |
| 9abec407 | 2026-11-06 | 61 | 57 | articles_wrappers A1: mrct-consistency on the new wrappers |
| 1d582bf7 | 2026-11-06 | 45 | 21 | five prose edits from the case-study review |
| 10c146b2 | 2026-11-06 | 16 | 9 | the (M1)-(M4) assumptions as a bulleted list |
| 7f23433e | 2026-11-06 | 71 | 35 | seven dense paragraphs as lists (R1-R7) |
| 9adc97a8 | 2026-11-07 | 13 | 5 | hazard-scale rationale and the two routes to consistency |
| db19388b | 2026-11-07 | 105 | 22 | align with the manuscript (A)-(E) |
| 488b351d | 2026-11-09 | 215 | 43 | Update mrct-consistency.qmd |
| 8a81ffa2 | 2026-11-09 | 7 | 8 | vignettes update manually |
| 9f91302b | 2026-11-14 | 24 | 13 | docs: MRCT articles, general installation note (Steps 1 and 2) |

`vignettes/articles/mrct-consistency-simulations.qmd`: dad519d6 2026-11-05 (+1452 −0, new file),
84d1cbd1 2026-11-06 (+100 −148, articles_wrappers A2: on the new wrappers), 71a7fc1e 2026-11-09
(+33 −20, update mrct consistency application), 9f91302b 2026-11-14 (+24 −13, installation note).

**1.2 Search terms** (`git log --all -S`, limited to `vignettes/ inst/extdata/ R/`; then every
article commit scanned with `-G`-style matching of removed lines):

| Term | Commits (-S) | In the articles' qmd |
|---|---|---|
| `region_treat_loghr` | 8 | added 43441a42, 450b5492; moved by dad519d6 (below) |
| `gamma_R` | 5 | added 8a1ba513, 10c146b2; never removed |
| `b_R` | 9 | added 8a1ba513, 10c146b2, db19388b; removed by 488b351d (below) |
| `h_R` | 8 | added 8a1ba513, 10c146b2; one line rewritten by 488b351d |
| `type-1` / `type 1` / `type I` / `type-I` | 4 / 4 / 1 / 4 | never; hits are `methodology.qmd`, `paper_simulations.Rmd` (and their HTML renders), `MRCT_slides_draft1b.pdf`, `R/grf_main.R`, `R/mrct_simulation.R`, `R/fs_oc_grid.R`, `R/simulation_tables.R` |
| `false-positive` | 6 | added db19388b; line rewritten by 488b351d |
| `false positive` | 13 | never |
| `no regional treatment` | 2 | added db19388b; line rewritten by 488b351d |
| `extension` | 33 | added 10c146b2, db19388b; never removed |
| `extended model` | 0 | never |
| `treatment-by-region` | 4 | added 8a1ba513, 10c146b2; line rewritten by 488b351d |
| `region_loghr_shift`, `region_harm_loghr_shift` (searched in `vignettes/ R/ tests/ inst/ dev/ man/`) | 2 each (0890b16f, 961a3de4) | never, in any version |

Commits that removed matching lines from an article:

- **dad519d6** (split into two articles), `mrct-consistency.qmd`: removed `region_treat_loghr = NULL,`
  twice, from the chunk `x3 = list(vars = "cd40", cuts = list(cd40 = 230), loghr = log(2 / 0.55)),  # h = 2 in H*`
  / `region_treat_loghr = NULL,` / `n_super = 20000, expand = "copula", seed = 7L,` and from the
  `build_own` chunk (`x3 = list(vars = "cd40", cuts = list(cd40 = cut), loghr = log(h / b)),` /
  `region_treat_loghr = NULL,`). The same lines were added to `mrct-consistency-simulations.qmd` in
  the same commit.
- **84d1cbd1** (articles_wrappers A2), simulations article: the `region_treat_loghr = NULL,` line of
  the `build_own` chunk re-indented (−1 +1).
- **488b351d** (2026-11-09, "Update mrct-consistency.qmd"), application article, in "We consider two
  explanations for the trigger" and "The data-generating mechanism":
  ```
  -  configurations with no regional treatment benefit ($b_R = 1.0$) quantify how often the strategy
  -  explains a trigger that a region-specific effect produced — the false-positive rate of the
  -  subgroup pathway (the accompanying manuscript's region-interaction cells).
  +  configurations with no regional treatment benefit (the true region effect is $> 1.0$) quantify how often the strategy
  +  falsely leads to meeting consistency in an identified subgroup — the false-positive rate of the
  +  subgroup pathway.
  -  effects inside the region become $(b_R, h_R)$, and (M2) is then the testable null
  -  $\gamma_R = \gamma_{H^*\!R} = 0$.
  +  effects inside the region become $(b_R, h_R)$, and (M2) is then $\gamma_R = \gamma_{H^*\!R} = 0$.
  -Writing the effect with treatment-by-region interactions makes (M2) a testable restriction rather
  -than a built-in one:
  +Writing the effect with treatment-by-region interactions extends (M2):
  ```
  The removed parenthetical "(the accompanying manuscript's region-interaction cells)" was added by
  db19388b (2026-11-07) from `dev/tasks/EDITS_fs-article-ms-alignment_2026-11-07.md:20` and is the
  only pointer either article ever had to region-interaction cells; it pointed to the manuscript.
  "Every cell of this article is the base model." was added by 8a1ba513 (2026-11-05) and is still in
  the article (`mrct-consistency.qmd:252`).

**1.3 Largest deletions.**

- `mrct-consistency.qmd`:
  - dad519d6 (−1,261): moved "Setting up simulations from your own trial", "How replicate 5 was generated: the report's simulation design" (with "The cells", "The analytic layer in brief", "Finding a case study", "Platform precision", "Because this trial is simulated"), "Simulation results" and "Appendix: the analytic tools" into the new simulations article. Every heading at 8a1ba513 is present in one of the two articles at dad519d6.
  - 42f36b87 (−161): replaced "Finding a case study" (top level) and "Replicate 5 of S2" (its subsections included "The identifier on the non-region patients…" and "The four sets") with "The case study: replicate 5 of S2", Steps 1–5. "Finding a case study" moved under "How replicate 5 was generated".
  - 450b5492 (−148): removed "A live demonstration", replaced by "The demonstration run" and "Finding a case study".
  - a597d673 (−59), 9abec407 (−57), 172140b8 (−50), 8a1ba513 (−46), 488b351d (−43): rewrote text and chunks in place. None removed a heading.
- `mrct-consistency-simulations.qmd`:
  - 84d1cbd1 (−148): rewrote the chunks on the new wrappers. No heading removed.
  - 71a7fc1e (−20): renamed "4.3a — Choosing the cells" … "4.3d — Reading the results" to the same titles without numbers.
  - 9f91302b (−13): rewrote the install note.

Heading texts that were ever in either article and are not at HEAD: "4.3a–4.3d …" (renumbered), "A
live demonstration", "Replicate 5 of S2", "The four sets", "The identifier on the non-region
patients, with post-selection inference" (all superseded as above). None concerns a regional
treatment effect or type-1 error.

**1.4 Answer: no.** At no commit on any ref, and in none of the 62 unreachable commits, did either
article run, tabulate or read records for a cell with a regional treatment effect. Across all 36
file states of the two articles at the 24 article commits: `region_treat_loghr` appears 0–2 times per version, always as
`region_treat_loghr = NULL`; `region_loghr_shift` / `region_harm_loghr_shift` never appear; `b_R` /
`h_R` appear only in prose and display math (from 8a1ba513), the only value being the prose
"($b_R = 1.0$)" of db19388b, removed by 488b351d; the only cell names in any version are N0, N1, S1,
S2, S3, S4, S2_D1 and S2_D2 (with the record suffixes `_rmlog4`).

## Step 2 — outside `master`

**2.1 Branches.** Local: `master` only. Remote: `origin/master` (= d059e054, merged), `origin/gh-pages`
(not merged; no common ancestor with `master`; 73 commits). gh-pages holds no `vignettes/` or
`inst/extdata/` path. It holds site renders `articles/mrct-consistency.{html,md}`,
`articles/mrct-consistency-simulations.{html,md}` (deploys 055401a4 2026-11-14, 622d322a and
c8c10716 2026-11-16) and 308/332-byte redirect pages `articles/articles/mrct-consistency*.html`.
Term matches in its renders: the application article's `gamma_R`, `b_R`, `h_R`, `false-positive`,
`no regional treatment`, `extension` and `treatment-by-region` prose (the same counts as the current
source) and the simulations article's two `region_treat_loghr` lines. Nothing else.

**2.2 Stash.** `stash@{0}` (9bd2fbe7, 2026-10-27, "On feature/sg-focus-transparency", quarto/simulations
logs) and `stash@{1}` (61d6965d, 2026-10-27, "On feature/glm-extension", rerender reports): neither
touches the articles or `inst/extdata/mrct-consistency/`.

**2.3 Copies on disk.** `find` over `~/Documents/GitHub` (13 repositories) and `~/Downloads`, names
`mrct-consistency*.qmd` / `*.html`: 23 files, all in forestsearch; none in `~/Downloads` or another
repository. "Blob" = first commit holding the identical file. Term counts as defined above.

| Path (forestsearch) | Size | Modified | Lines | Git | Blob | Matches |
|---|---|---|---|---|---|---|
| `vignettes/articles/mrct-consistency.qmd` | 59,765 | 2026-11-14 09:35 | 1,001 | tracked | 9f91302b | 18 (gamma_R 5, b_R 4, h_R 3, false-positive 1, no regional treatment 1, extension 2, treatment-by-region 2) |
| `vignettes/articles/mrct-consistency-simulations.qmd` | 86,074 | 2026-11-14 09:35 | 1,428 | tracked | 9f91302b | 2 (region_treat_loghr) |
| `dev/article_renders/pre_0.5.0/mrct-consistency.html` | 2,231,840 | 2026-11-09 13:30 | 4,438 | tracked | 8a81ffa2 | 18 (as the qmd) |
| `dev/article_renders/pre_0.5.0/mrct-consistency-simulations.html` | 1,683,496 | 2026-11-09 16:14 | 5,367 | tracked | 71a7fc1e | 2 (region_treat_loghr) |
| `docs/articles/mrct-consistency.html` | 599,895 | 2026-11-04 09:40 | 3,555 | ignored | none | 2 (region_treat_loghr) |
| `legacy_vignette/06Nov2026/mrct-consistency.qmd` / `.html` | 37,394 / 1,969,112 | 2026-11-06 16:46 | 661 / 4,075 | ignored | 9abec407 / 29ed666e | 9 / 9 |
| `legacy_vignette/06Nov2026_2/mrct-consistency.qmd` / `.html` | 39,495 / 1,971,435 | 2026-11-06 17:15 | 685 / 4,075 | ignored | 1d582bf7 | 9 / 9 |
| `legacy_vignette/06Nov2026_3/mrct-consistency.qmd` / `.html` | 39,765 / 1,971,947 | 2026-11-06 17:19 | 692 / 4,082 | ignored | 10c146b2 | 15 / 15 |
| `legacy_vignette/06Nov2026_4/mrct-consistency.qmd` / `.html` | 40,174 / 1,972,898 | 2026-11-06 19:40 | 728 / 4,129 | ignored | 7f23433e | 15 / 15 |
| `legacy_vignette/07Nov2026/mrct-consistency.qmd` / `.html` | 40,174 / 1,972,898 | 2026-11-06 20:22 / 19:40 | 728 / 4,129 | ignored | 7f23433e | 15 / 15 |
| `legacy_vignette/07Nov2026_2/mrct-consistency.qmd` / `.html` | 41,013 / 1,973,834 | 2026-11-07 09:53 | 736 / 4,129 | ignored | 9adc97a8 | 15 / 15 |
| `legacy_vignette/07Nov2026_3/mrct-consistency.qmd` / `.html` | 41,013 / 1,973,834 | 2026-11-07 12:34 | 736 / 4,129 | ignored | 9adc97a8 | 15 / 15 |
| `legacy_vignette/09Nov2026/mrct-consistency.qmd` | 47,665 | 2026-11-09 11:52 | 819 | ignored | **none** | 18 |
| `legacy_vignette/09Nov2026/mrct-consistency.html` | 1,990,623 | 2026-11-07 12:42 | 4,307 | ignored | db19388b | 19 |
| `legacy_vignette/09Nov2026/mrct-consistency-simulations.qmd` / `.html` | 84,342 / 1,683,496 | 2026-11-06 16:31 | 1,404 / 5,367 | ignored | 84d1cbd1 / 29ed666e | 2 / 2 |

The two copies not in history: `legacy_vignette/09Nov2026/mrct-consistency.qmd` is db19388b's text
with part of 488b351d's prose edits already made (its `diff` against db19388b covers lines 133–254:
the trigger paragraph, the two explanations, ACTG 175, (M2), (M4), hazard scale). It has no cell,
chunk or record beyond db19388b. `docs/articles/mrct-consistency.html` is a pre-split render
(2026-11-04) of the single article (Step 3). The 9-match copies predate 10c146b2 (the (M1)–(M4)
bullets), which added 6 matching occurrences.

## Step 3 — the earlier renders

Headings were taken from `<h1>`–`<h4>`, section numbers stripped and math markup normalised.

- `dev/article_renders/pre_0.5.0/mrct-consistency.html`: Overview: the question and the design;
  The trigger and the question; The data-generating mechanism; Analytic predictions before the
  analysis; The pipeline, for one trial; The case study: replicate 5 of S2; Step 1: the all-patient
  Cox fit; Step 2: the region's Cox fit, and the trigger; Step 3: the identifier on the non-region
  patients, with post-selection inference; Step 4: Ĥ applied to the region; the four sets; Step 5:
  the Pareto frontier and the bounds; The S2 campaign in brief; The one-call form; Appendix: the
  analysis window and the windowed Cox limit; Appendix: the manuscript's event notation.
- `dev/article_renders/pre_0.5.0/mrct-consistency-simulations.html`: The case study's objects,
  recomputed; Setting up simulations from your own trial; How replicate 5 was generated: the report's
  simulation design; The data-generating mechanism, built from ACTG 175; The cells; The analytic
  layer in brief; How r = −log 4 was chosen, and the later-enrolling regions; Finding a case study;
  Because this trial is simulated; Ĥ against H*, in the region and outside; Reproduction check
  against the committed record; The four sets beside their limits; Simulation results; The core
  operating characteristics at n = 1,000; The later-enrolling cells; Predicted against observed;
  Post-selection inference on Ĥ in the triggered trials; The campaign code and its wall-clock; The
  demonstration run; Appendix: the analytic tools.
- `docs/articles/mrct-consistency.html` (2026-11-04, pre-split): the union of the two lists above
  except "The data-generating mechanism", "Analytic predictions before the analysis", "The one-call
  form", the application article's two appendices and the simulations article's "The case study's
  objects, recomputed".
- The 12 legacy renders: every heading is in the current sources.
- The gh-pages renders (c8c10716): 15 and 20 headings, every one in the current sources (three
  differ only in MathML markup of Ĥ). Their `.html` and `.md` have 18 and 2 term matches, the
  same as the sources.

No heading in any render concerns the extended model or type-1 error. No table does either: the
simulations render's table captions run from "The cells. Every patient's hazard ratio is the same in
every region." to "The core rows of S2: the live demonstration beside the committed campaign".
The extended model appears in the renders only as the application article's prose.

**Headings present in a render but absent from the current sources: none.** The reverse also
occurs: the current application source has headings the pre_0.5.0 render lacks ("Installing
forestsearch", "Rendering this article yourself", "The region in this article", "The windowed Cox
limit", "The pipeline in one call", "How to use this article", "What the move below 1.0 shows, and
what it does not", "Reading the bounds").

## Step 4 — the shipped records

`inst/extdata/mrct-consistency/`: 29 tracked files, unchanged since 98d5a801/de222e41 (2026-11-02;
moved from `vignettes/articles/mrct-consistency/`, which held the same 29 names and nothing else):
`PROVENANCE.md`, `s2_example_definitions.md` and 27 CSVs. Distinct `cell` / `setting` values:

| CSV | Values |
|---|---|
| `build_full_summary.csv` (cell), `core_full_numeric.csv`, `timing_full.csv` (setting), `predictions_full.csv` (cell) | N0_rmlog4 N1_rmlog4 S2_rmlog4 S3_rmlog4 S4_rmlog4 |
| `core_n1000.csv` (setting) | N0_rmlog4 N1_rmlog4 S1_rmlog4 S2_rmlog4 S3_rmlog4 S4_rmlog4 |
| `core_S1_numeric.csv`, `timing_S1.csv` (setting); `predictions_S1.csv` (cell) | S1_rmlog4 S1_rmlog4_n700 |
| `predictions_n700.csv` (cell) | N0/N1/S2/S3/S4 `_rmlog4_n700` |
| `mr_field_S2_rep5.csv`, `s2_example_four_sets.csv`, `s2_example_trials.csv` (setting) | S2_rmlog4 |
| `timing_entry_cells.csv` (setting) | S2_D1 S2_D2 |
| `timing_mr_field_<cell>.csv` (setting) | one each of N0, N1, S1, S2, S3, S4 (`_rmlog4`), S2_D1, S2_D2 |

Files without such a column carry cells as columns: `mr_field_core.csv` (N0 … S4, S2_D1, S2_D2),
`entry_cells_core*.csv` (S2, S2_D1, S2_D2), `entry_roots.csv` (rows D1, D2 and two a = 0
references). No record concerns a cell other than the eight. No file in the directory matches any
search term. `PROVENANCE.md` says the records are "Records of the mrct-consistency repository, read
at commit `2341f96` (2026-11-01 …)". That commit is earlier than the region-interaction campaign
(0108843, 2026-11-05), and the file names no other cells.

## Step 5 — the package

**`region_treat_loghr`** occurs in `R/mrct_inject_dgm.R` only: header comment line 14 ("or
`region_treat_loghr` (Region-U-TE)"), roxygen 318–320, signature 378, checks 418–420 and 427–429,
use 495–498, record 587 (`treat_interaction_loghr = region_treat_loghr`). Help
(`man/inject_mrct_structure.Rd:68`): "NULL, or a log-HR for a direct treat x Region interaction (the
residual Region-U-TE pathway). Mutually exclusive with x3 (the package supports one interaction
term)." A value sets the subgroup variable to the region and puts the value on `treat_harm`
(lines 495–498), so it cannot be combined with the `x3` subgroup H*. Code errors: "x3 and
region_treat_loghr are mutually exclusive" and "region_treat_loghr and region_loghr_shift /
region_harm_loghr_shift are mutually exclusive". It does **not** build the `(b_R, h_R)` mechanism
the application article describes.

**The `(b_R, h_R)` mechanism** is built by `inject_mrct_structure(region_loghr_shift,
region_harm_loghr_shift)`, exported. Roxygen `R/mrct_inject_dgm.R:321–344`: γ_R and γ_{H*R}, with
"γ_R = log(b_R/b) and γ_{H*R} = log(h_R/h) − log(b_R/b)", applied to the treated potential
outcome after the AFT fit. The arguments are in the signature at 379–380, validated at 421–431,
applied at 559–560 and 607–610, and recorded in `$mrct$region_effect` at 599–600. Defaults 0. They
need finite numbers, and `region_harm_loghr_shift` needs `x3`. `ahr_limit()` and `ahr_design()`
(exported) take `b_R`, `h_R` (`R/ahr_analytic.R:503–510`, 531–539, 592–595, 611; "`NULL`
(default) means `b_R = b` and `h_R = h`, the two-cell model"). `ahr_truths_table()` and
`ahr_design_table()` pass them through `...` (`R/ahr_tables.R:28–29`). All of this was added by
0890b16f (2026-11-04, "feat: region-dependent treatment effect in inject_mrct_structure() and
four-cell ahr_limit()/ahr_design()", version 0.4.1.9001), under
`dev/tasks/TASK_fs-dgm-region-interaction_2026-11-04.md`, which says "No compute: this task runs no
simulation campaign. Downstream: TASK_mrct-sim-dgm-section-predictions (mrct-consistency repo)
depends on this task being green." `mrct_region_sims()` has no argument of its own for it and
consumes the built DGM.

**Tests.** `tests/testthat/test-mrct-region-interaction.R` (`skip_on_cran()`;
`mrct_skip_if_no_build()` in each test) exercises the shifts and `b_R`/`h_R`: T1 (line 38) default
shifts reproduce the pre-change build; T2 (65) `b_R = h_R = 1` gives ψ = 0 in the region; T3 (91)
the four-cell AHR identity; T4 (115) γ-form = (b_R, h_R)-form on a grid; a validation test (148). A
non-NULL `region_treat_loghr` occurs in tests only as an expected error (line 157). It is `NULL` in
`helper-mrct-analytic.R:46` and `test-mrct-minor-helpers.R:111`.

## Step 6 — the mrct-consistency repository (read-only)

Present. Checkout branch `docs/workstream-index-2026-11-14`, last commit `60a53bd`
(2026-11-15 00:05, "docs: index canonical shared reports with paired Quarto and PDFs at repository
root"); `git status --short` empty; `main` = be7be00 (2026-11-09).

**Terms** (`git grep` over 2,003 tracked files): 172 files match. Totals: b_R 159, false positive
154, h_R 149, extension 132, false-positive 47, region_treat_loghr 34, gamma_R 24, type-I 19,
treatment-by-region 17, type 1 15, type-1 2, extended model 1, no regional treatment 1, type I 1.
The ten files with most matches:

| File | Matches |
|---|---|
| `docs/manuscript/sim/ext-predictions.R` | 37 |
| `docs/manuscript/sim/sim-manuscript.qmd` | 34 |
| `docs/manuscript/sim/ext-s2null.R` | 31 |
| `sims/region_interaction_cells_design.R` | 30 |
| `reports/REPORT_mrct-companions-screen.md` | 22 |
| `docs/manuscript/sim/records.R` | 20 |
| `docs/literature_review/mrct_consistency_subgroup_review.html` | 16 |
| `docs/manuscript/sim/analytic-ext/s2null.csv` | 16 |
| `docs/manuscript/techreport/_book/search.json` | 16 |
| `docs/manuscript/sim/sim-manuscript.tex` | 15 |

**Simulations of cells with a regional treatment effect.** Built with
`inject_mrct_structure(region_loghr_shift = gamma_R, region_harm_loghr_shift = gamma_HR)`
(`sims/region_interaction_cells_design.R:64`). Base: S2 at r = −log 4 (b = 0.55, h = 2, H* = {CD4 ≤
230}), n = 1,000. Committed in **0108843** (2026-11-05 03:57, "mrct-sim-campaign-ext 2026-11-04:
the region-interaction cells S2_E2 (b_R = h_R = 1.0) and S2_E1 (b_R = 1.0, h_R = 2.0) at n = 1,000,
10,000 replicates each …"), under `dev/tasks/TASK_mrct-sim-campaign-ext_2026-11-04.md`, whose scope
line reads "mrct-consistency only … No forestsearch changes". 0108843 is on `main` and every other
local and remote branch.

| Path | Cell (b_R, h_R) | Replicates | n | Run |
|---|---|---|---|---|
| `results/region_interaction_cells/S2_E2_rmlog4/` (replicates, estimators, hhat_confusion, Hstar_nonregion, truth, timing, call_warnings) | S2_E2 (1.0, 1.0), seed 615400000 | 10,000 | 1,000 | 2026-11-04 22:43–23:05, 120 workers |
| `results/region_interaction_cells/S2_E1_rmlog4/` (same files) | S2_E1 (1.0, 2.0), seed 615500000 | 10,000 | 1,000 | 2026-11-04 23:05–23:27 |
| `results/region_interaction_cells/screen/` (same layout + `summary/screen_gate_G4.csv`) | both | 500 each | 1,000 | 2026-11-04 22:36–22:42 |
| `results/region_interaction_cells/summary/` (core, core_numeric, pred_obs, region_events, mr_field_core, mr_field_detail, timing), `build.csv`, `checks.csv`, `logs/` | both | — | — | records of the above |
| `results/mr_field_campaign/S2_E2/`, `S2_E1/` (mr_field.csv, timing.csv), `logs/run_region_interaction.log` | S2_E2 / S2_E1 | 3,627 / 4,880 triggered trials | 1,000 | 2026-11-04 23:27 – 2026-11-05 03:55 |

Scripts: `sims/250_region_interaction_cells_truth.R`, `251_region_interaction_cells.sh`,
`251_region_interaction_cells_call.R`, `252_…_records.R`, `253_…_tables.R`,
`254_region_interaction_mr_tables.R`, `region_interaction_cells_design.R`. The analytic (not
simulated) predictions of these configurations are in `docs/manuscript/sim/analytic-ext/` (written
by `ext-predictions.R`, `ext-s2null.R`; commits 31a6643, 29e3237, 1e3fe07, 2026-11-04). The
companion cells of `results/companions_10k/` and `results/companions_screen/` (N0, N1, S3, or3, n,
q variants) are base-model cells.

**Where they are reported.** SiM manuscript `docs/manuscript/sim/sim-manuscript.qmd`, read through
`docs/manuscript/sim/records.R` (`sim_ri_campaign()`, `sim_ri_table()`), as committed in
`sim-manuscript.pdf` (be7be00, 28 pp.):

- §2 "Setting, estimands and notation" (qmd from line 207): "**A region-dependent extension.**"
  (line 290), the sentence "Under the extended configurations the region outside H* is truly null
  ($b_R = 1.0$), so transportability-pathway consistency there is spurious: rows (5) and (8) read as
  the false-positive increase in meeting consistency via the oracle and the declared pathways"
  (line 362), and **Table 2** (PDF p. 6), "Analytic predictions (rows (1)–(7)) and campaign-based
  declared-consistency quantities (rows (8)–(9)) under the four configurations", with columns S2 |
  S2, null (b_R = 1.0, h_R = 1.0) | S2, detrimental (b_R = 1.0, h_R = 2.0) | N0. Row (8) of the two
  extended columns is campaign-observed (0.057, 0.115).
- §5.2 "Results" of §5 "Simulation study" (qmd 1751; table chunk 1877–1898, labels and caption set at
  1459–1468): **Table 6** (PDF p. 22), "Operating characteristics of S2 and the two
  region-interaction configurations at n = 1,000, 10,000 trials per cell". Non-consistency is 24.0%
  (2,398) in S2, 36.3% (3,627) in S2, null and 48.8% (4,880) in S2, detrimental. Declared,
  region(Ĝ) < 1.0 is 51.2 / 15.7 / 23.5, and Region(G*) < 1.0 is 77.7 / 19.4 / 37.3. Below the
  table (qmd 1900): "The oracle pathway's false-positive increment exceeds the identifier's in every
  null-type configuration".
- The same tables are in `sim-manuscript-complete.pdf`. The technical report
  (`docs/manuscript/techreport/*.qmd`) and the supplement (`sim-supplement.qmd`) do not name these
  cells; the techreport's term hits are `region_treat_loghr = NULL` (`dgm.qmd:49`), "extension" for
  the later-enrolling region (`mechanisms.qmd`, `discussion.qmd:351`), "a declared Ĥ is a false
  positive" for N0 (`campaign.qmd:34`), one "false positive" in `dgm.qmd` and "type 1-term"
  (`analytic-layer.qmd:996`).

## Step 7 — close

Commits of this session on forestsearch `master`: 20f204b1 (task file), 70ad388a (this report) and
the commit adding this section. The report is copied to `~/Downloads`.

Final checks (run after 70ad388a):

- `git diff --no-renames --name-status d059e054 HEAD`: `A dev/tasks/TASK_mrct_type1_scenarios_trace_2026-11-16.md`,
  `A dev/reports/REPORT_mrct_type1_scenarios_trace_2026-11-16.md`. Nothing else. `git status --short`
  is empty. Outside `.git/` and `quarto/`, no other file in forestsearch is newer than the session
  start (12:28).
- Other repositories: no file newer than 12:28 in any of the 12 under `~/Documents/GitHub/`.
  mrct-consistency HEAD is still `60a53bd`, with a clean status. `Cox_AHR` and `guidesubgroup` each
  show one `git status` entry, both older than this session.
- No push, pull, fetch, merge or tag. Tags: v0.1.0, v0.4.0 and v0.5.0, as at the start. The newest
  `origin/master` reflog entry is still `d059e054` (10:14:58), so neither session commit had reached
  the remote-tracking ref at the time of the check. GitHub Desktop runs on this host and has pushed
  `master` commits on its own before.
