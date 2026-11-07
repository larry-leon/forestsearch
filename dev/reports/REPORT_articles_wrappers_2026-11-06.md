# REPORT — the four articles on the new wrappers (2026-11-06)

Task: `dev/tasks/TASK_articles_wrappers_2026-11-06.md`, with Amendment 2.  
**Outcome: all four articles render on the new wrappers and every published number is reproduced, with 0 non-volatile differences in each article. One gate reading needs your review: 16 timing tokens moved outside the strict two-render volatile set (item 1 under "Read these first").**  
The work is on branch `fs-articles-wrappers`, unmerged and unpushed, for your review. The final article tree is `ae0af1bb`; the HTML, amendment and report commits sit on top of it. There are no `R/`, test or `man/` changes.

## Gates

| Gate | Status | Evidence |
|---|---|---|
| Gate 0: prerequisite | PASS | `f404556f` is an ancestor; no `R`, `tests`, `man`, `DESCRIPTION`, `NAMESPACE` or `inst` diff from it; 12 of 12 exports installed |
| Gate 0: control renders | PASS | 8 of 8 exit 0, so every internal check passes on the unedited tree |
| Gate 0: volatile sets, inventory | PASS | Sets recorded per article; 61 assertions inventoried (41 protected) |
| Gate R (a): renders | PASS | 4 of 4 exit 0, from `ae0af1bb` |
| Gate R (b): numeric tokens | PASS, see item 1 | 0 non-volatile differences in each article; 16 declared timing tokens moved outside the strict set |
| Gate R (c): one-call `stopifnot` | PASS | `mrct-consistency.qmd:650` ran; a failure would have stopped the render |
| D5: `stopifnot` counts | PASS | 14 to 14, 41 to 39, 5 to 5, 1 to 1, all explained by the listed removals and the one new check |
| Stage 3: confinement, as amended | PASS | 4 `.qmd`, 4 HTML, the task doc's amendment; no `.gitignore`, `R/`, `tests/`, `man/`, `NAMESPACE` or `DESCRIPTION` entry |
| Stage 3: suite and `check()` | Not run, as specified | All four articles are Rbuildignored (`.Rbuildignore:20`) |

## Read these first

1. **The volatile set is wider than D3's strict reading.** D3 takes the tokens that differ between the two control renders. A timing cell can match in both by chance (A1's "The call took 62 seconds", survival's `5.31`), and A1's "The article took" moves by design with the extra minute. So I declared the timer-fed positions from the article source and froze them in a file at 15:56, before the first edited render finished at 15:59. No position was added afterwards. The gate reports them apart from the strict set. **Under the strict set alone, Gate R would show 16 differences** (A1 2, A2 3, A3 5, A4 6). Each is listed verbatim below and each is a timing token. Two of them (A2 `#566`, `#567`) are a display effect: `c(replicates, workers, seconds)` printed `100.0` as `100` because seconds came out whole. Whether this satisfies D3 is your call.
2. **The one-call form is Section 2.7, not straight after step 5.** It sits after "The S2 campaign in brief" and before the closing timing line. Placed right after step 5 it would renumber Section 2.6, a numeric change the gate forbids. The task calls it a closing subsection, so I read this as within the task.
3. **Two visible changes that move no number.** KM panel 4's title reads "Region minus H-hat (G-hat)" where it read "Region minus Ĥ (Ĝ)", because the package's labels are ASCII. In A2 the three build warnings print as `Warning: Ran out of iterations and did not converge`, without the `survreg.fit(...)` call prefix, because `mrct_build_dgm()` muffles them and the chunk re-raises them from the attribute. The sentence "Its 3 warnings, shown above" still holds.
4. **Three A3 checks could not stay verbatim, against D2's wording,** because they name objects this task removes (`in_region()`, `g`, `cc`, `cal`). Each asserts the same thing on the wrapper's output, and none is a protected check against a record. They are `:240`, `:347` and `:477`, marked "re-expressed" in the A3 inventory.
5. **Thin helpers kept, for verbatim checks.** `design_tab()` stays as a map from `ahr_design_table()`'s columns to the committed records' names (`p_A` to `p_A_B_Cstar`), so the 1e-9 checks do not change. `build_own()` and `build_cell()` stay as one-call wrappers because later chunks call them. `demo-summary` calls `mrct_core_summary(structure(demo, estimators = est))` because the protected `find-select-check` chunk drops the attribute, and I left that chunk alone.
6. **Two small leftovers.** `get_dfpred` stays in A1's `needed` although A1 no longer calls it, since the task says extend; A2 still calls it directly. `s2_example_definitions.md` is not a CSV, so `mrct_record()` cannot read it and A2 reads it by `system.file()`: the installed copy, even from a source tree. `diff -rq` shows the installed records identical to `inst/extdata/mrct-consistency`.

## Amendment 2: rule or omission

**Finding: omission. No rule excludes the file, so `.gitignore` is unchanged.** `git check-ignore -v -n` reports no match (`::`, exit 1), with the file absent and with it in the tree, where `git status` lists it as `??`. The nearest rule, `vignettes/*.html` (`.gitignore:15`), does not reach `vignettes/articles/`; `vignettes/.gitignore` and `vignettes/articles/.gitignore` have no HTML rule, `.git/info/exclude` has no HTML line, and there is no global excludes file.

It is not quite "never added", though. The history:

| Commit | Date | What happened to `survival-post-selection.html` |
|---|---|---|
| `606a78ea` | 2026-10-28 | Added ("commit html vignettes") |
| `3e491054` | 2026-10-30 | Removed with four others, under a new rule `vignettes/articles/*.html` ("kept on disk, never tracked") |
| `259f079e` | 2026-11-02 | Rule replaced ("rendered article HTML is tracked"). Five HTMLs were added, none of them this file |

So the omission dates from `259f079e`: the rule went away, and the file was not re-added.

- **Item 2:** committed in `772a032d`, on its own. The file is byte-identical to its Gate R render (`cmp`) and its `.qmd` is unchanged since `ae0af1bb`, so nothing was re-rendered; the token result below is that render's.
- **Item 3:** see Stage 3.
- **The label:** this task doc had no Amendment 1. I recorded this one as received, numbered 2, in `3e3ab055`.
- **A stale commit message:** `29ed666e` says the file "is untracked and stays so". Amendment 2 supersedes that; I did not rewrite the commit.
- **Out of scope, untouched:** `extreme_subgroups`, `guohe-nested-continuum`, `methodology` and `paper_simulations` also have a tracked `.qmd` and no tracked HTML in `vignettes/articles/`.

## Commits

```
772a032d articles_wrappers: rendered HTML for survival-post-selection (Amendment 2)
3e3ab055 task: articles_wrappers Amendment 2 (as received)
29ed666e articles_wrappers: rendered HTML from the Gate R render
ae0af1bb articles_wrappers A4: survival-post-selection on the new wrappers
315497e3 articles_wrappers A3: actg175-glm-outcomes on the new wrappers
84d1cbd1 articles_wrappers A2: mrct-consistency-simulations on the new wrappers
9abec407 articles_wrappers A1: mrct-consistency on the new wrappers
7989b80f task: articles_wrappers 2026-11-06 (as received)
```

One commit per article, with no fix commits. This report's commit sits on top. Nothing is pushed.

## Stage 0 and Gate 0

```
f404556f is an ancestor of HEAD
git diff f404556f HEAD -- R tests man DESCRIPTION NAMESPACE inst | wc -l: 0
exported: 12 of 12 
installed version: 0.4.1.9002
```

| Article | Kind | Rendered output at baseline | `_freeze` |
|---|---|---|---|
| `mrct-consistency.qmd` | Rbuildignored website article | Tracked | None |
| `mrct-consistency-simulations.qmd` | Rbuildignored website article | Tracked | None |
| `actg175-glm-outcomes.qmd` | Rbuildignored website article | Tracked | None |
| `survival-post-selection.qmd` | Rbuildignored website article | Untracked, not ignored | None |

There is no `_quarto.yml` either.

| Render | mrct-consistency | simulations | actg175 | survival |
|---|---|---|---|---|
| Control A, exit 0, wall (s) | 139 | 250 | 364 | 145 |
| Control B, exit 0, wall (s) | 138 | 242 | 362 | 150 |
| Inner loop, exit 0, wall (s) | 172 | 233 | 344 | 139 |
| Gate R, exit 0, wall (s) | 193 | 244 | 356 | 145 |

That is 16 renders and about 33 minutes of wall time, because pairs ran side by side, against the task's estimate of 60 to 90.

## How the tokens were compared

- **What is read:** the text of `<main>` in the rendered HTML, without `<script>`, `<style>` and code blocks (`sourceCode`, `cell-code`), since code changes by design. Cell outputs, tables, captions and prose are kept.
- **A token:** a digit run with its decimal or thousands parts, a leading minus and an exponent, including digits inside names such as `z_cd40`.
- **Volatile set:** control A against control B, aligned with `difflib.SequenceMatcher(autojunk=False)`.
- **The gate:** final against control A. It passes only if every difference is an equal-length replacement at a strict-volatile or declared-timing position.
- **A1's new section:** the 73 tokens from "2.7 The one-call form" up to "The article took" are cut from the final render before comparing. The closing timing line stays in.
- **Not gating:** the visible text with digits blanked was reviewed for non-numeric changes, and embedded figures were hashed.
- **The tool:** a scratch script, not committed, as the diff is confined.

## Gate R

| Article | Tokens, baseline / final | Strict volatile (moved) | Declared timing | Moved outside strict | Non-volatile differences |
|---|---|---|---|---|---|
| A1 `mrct-consistency` | 346 / 346 | 2 (1) | 7 | 2 | **0** |
| A2 `mrct-consistency-simulations` | 1454 / 1454 | 5 (4) | 13 | 3 | **0** |
| A3 `actg175-glm-outcomes` | 669 / 669 | 19 (19) | 25 | 5 | **0** |
| A4 `survival-post-selection` | 530 / 530 | 19 (19) | 27 | 6 | **0** |

## A1: `mrct-consistency.qmd`

657 lines to 661.

**Chunks changed**

| Chunk | Now at | Lines before | Lines after | Change |
|---|---|---|---|---|
| `packages` | `:38` | 22 | 16 | `rd()`/`dd` block to `mrct_record()` |
| `install-check` | `:90` | 15 | 16 | `needed` extended |
| `pre-build` | `:226` | 27 | 27 | `rd(` to `mrct_record(` only |
| `case-read` | `:345` | 12 | 12 | `rd(` to `mrct_record(` only |
| `case-sets` | `:448` | 27 | 27 | `in_Hhat` to `fs_in_subgroup()`; identity check kept |
| `fig-case-forest` | `:480` | 20 | 8 | `plot()` on `mrct_four_sets()` |
| `case-km` | `:489` | 23 | 6 | `plot()` on `mrct_four_sets()` |
| `case-campaign` | `:583` | 36 | 36 | `rd(` to `mrct_record(` only |
| `one-call` | `:632` | none | 23 | new: the one-call form and its check |

**Prose edits, verbatim** (the `.qmd` outside chunks)

```diff
-region (1 = in the region) columns and the baseline covariates the identifier may use.
+region (1 = in the region) columns and the baseline covariates the identifier may use: they are the
+arguments of `mrct_trial_analysis()`, the one-call form of the steps below, which closes this section.
-*With your trial:* your time, event and treatment columns in the formula; if ITT ub $\ge$ 1.0 the
-pipeline stops here.
+*With your trial:* `mrct_trial_analysis()` fits this from your time, event and treatment columns; if
+ITT ub $\ge$ 1.0 the pipeline stops here, at stage `not_positive`.
-trigger. *With your trial:* your region column in `reg`; if region ITT $\le$ 1.0 the pipeline stops
-here.
+trigger. *With your trial:* `mrct_trial_analysis()` fits this from your region column; if region ITT
+$\le$ 1.0 the pipeline stops here, at stage `no_trigger`.
-$\hat H$ is read for every patient from the identifier's selected rule (`get_dfpred()`), checked against
+$\hat H$ is read for every patient from the identifier's selected rule (`fs_in_subgroup()`), checked against
-The package's `plot()` draws the forest plot and the Kaplan–Meier panels from an `mrct_case_study`
-object only (simulations article), so here they are drawn from these objects with the same code, kept short.
+`mrct_four_sets()` takes the trial's columns and this membership of $\hat H$, and its `plot()` draws the
+forest plot and the Kaplan–Meier panels, by the code `plot()` uses for an `mrct_case_study` object
+(simulations article).
-`get_dfpred()` applies your $\hat H$ to all your rows; drop the region($G^*$) row, as $H^*$ is unknown.
+`fs_in_subgroup()` applies your $\hat H$ to all your rows; without `truth.name`, `mrct_four_sets()` drops
+the region($G^*$) row, as $H^*$ is unknown.
+## The one-call form {#sec-one-call}
+
+`mrct_trial_analysis()` is steps 1 to 5 in one call: the two Cox fits and the trigger, the identifier on
+the non-region patients with multiplier resampling, $\hat H$ applied to every row, and the four sets. Here
+it runs on the same trial, with the identifier's arguments of step 3 and the record's seed.
+
+
+*With your trial:* your column names, your covariates in `fs_args$confounders.name` and your seed in
+`mr_inference_args`, without `truth.name`. `$stage` says how far the pipeline went (`not_positive`,
+`no_trigger`, `no_declaration` or `complete`), and none of these is an error; `plot()` on the result draws
+the two figures of step 4, and `mr_values()` and `mr_estimates_table()` read its `$fs`.
+
```

The blank lines in the new section are where the `one-call` chunk sits. "With your trial" notes on argument meanings (step 3, step 5) stay verbatim.

**Removed scaffolding assertions**

- `stopifnot(!is.na(dd), file.exists(file.path(dd, "PROVENANCE.md")))` (baseline `:50`) went with `dd`. Covered by Gate B5 of `REPORT_wrapper_batch_2026-11-06.md` (the absent-directory and absent-provenance errors) and `tests/testthat/test-mrct-minor-helpers.R:17`.

**Volatile set (control A / control B), verbatim**

```
tokens A 346 | tokens B 346 | volatile positions 2
  #75: 1.200000 / 1.100000  ::  s  |  -1.541238 -1.386294 3.000000  |  seconds  |  [[1.200000]]  | Seven design probabilit
  #170: 61.9 / 62.5  ::   candidate_family n_Hhat n_Hhat_complement  |  [[61.9]] 3003.0 73.0 733.0  | Code | 
```

**Declared timing positions outside the strict set, frozen before any edited render**

```
#171: 3003.0  ::  idate_family n_Hhat n_Hhat_complement  |  61.9 [[3003.0]] 73.0 733.0  | Code | [1] "{z
#172: 73.0  ::  amily n_Hhat n_Hhat_complement  |  61.9 3003.0 [[73.0]] 733.0  | Code | [1] "{z_cd40
#173: 733.0  ::   n_Hhat n_Hhat_complement  |  61.9 3003.0 73.0 [[733.0]]  | Code | [1] "{z_cd40 <= 23
#184: 62  ::  6 and prior ART over 196 days. The call took [[62]] seconds on one thread. W
#345: 2.2  ::  rials with region(Ĝ) < 1.0. | The article took [[2.2]] minutes to run.
```

**Gate R, verbatim**

```
baseline tokens 346 | final tokens 346 (cut '2.7 The one-call form' .. 'The article took') | strict volatile positions 2 (moved 1) | declared timing positions 7
timing-declared positions outside the strict set that moved: 2
    #184: 62 -> 60  ::  6 and prior ART over 196 days. The call took [[62]] seconds on one thread. W
    #345: 2.2 -> 3.1  ::  rials with region(Ĝ) < 1.0. | The article took [[2.2]] minutes to run.
non-volatile differences: 0
```

**Figures (not gating)**

```
figures baseline 3 | final 3
  fig 1: same (84b6e5c66249e67d / 84b6e5c66249e67d)
  fig 2: DIFFERENT (12a8c1999cff30ac / 65046017bd64e69d)
  fig 3: same (11d8717738c44266 / 11d8717738c44266)

fig 2, the Kaplan-Meier panels, by pixel:
size 2304x1728 (bytes per pixel 3) | 2304x1728 (bytes per pixel 3)
differing pixels 5274 of 3981312
row bands: [(895, 927)] | columns 1582 to 1975
```

The forest plot and the frontier figure are byte-identical. The KM figure differs only in one band, panel 4's title. The old KM chunk's second code fold (around `par(op)`) is gone as well.

**(c) The one-call check** is at `:650`: its H-hat, in_Hhat, the four sets at 1e-12 and `mr_values()` at 1e-12, timing excluded, against the step-by-step objects. The render exited 0 with stage `complete` printed in the section. The article's closing line moved from 2.2 to 3.1 minutes in this render, which ran beside the other three.

**Protected-check inventory** (`grep -c 'stopifnot('`: 14 before, 14 after). P = against a shipped record, committed value, published table or a reading's stated claim; I = internal identity; S = scaffolding.

| Baseline | Chunk | Kind | Asserts | Now | Status |
|---|---|---|---|---|---|
| `:50` | `packages` | S | `!is.na(dd), file.exists(file.path(dd, "PROVENANCE.md"))` | none | removed |
| `:270` | `pre-table` | I | `near(rows_pre[3] + rows_pre[4], rows_pre[1], 1e-12), near(rows_pre[5] + rows_pre[6], ...` | `:265` | verbatim |
| `:357` | `case-read` | P | `length(seed_rep5) == 1, !anyDuplicated(trial$id), all(c("id", "y_sim", "event_sim",  ...` | `:352` | verbatim |
| `:373` | `case-itt` | P | `itt_ub < 1` | `:369` | verbatim |
| `:390` | `case-region` | P | `trigger` | `:386` | verbatim |
| `:411` | `case-fs` | I | `RNGkind()[1] == "Mersenne-Twister"` | `:407` | verbatim |
| `:430` | `case-fs` | I | `!is.null(fs$sg.harm), mr$ci_method == "field", mr$settings$draws == 2000, mr$field$n ...` | `:426` | verbatim |
| `:455` | `case-sets` | I | `identical(in_Hhat, trial$id %in% c(fs$df.est$id[fs$df.est$treat.recommend == 0], fs$ ...` | `:451` | verbatim |
| `:469` | `case-sets` | I | `near(four$hr[1:2], c(fit_itt$conf.int[1, "exp(coef)"], hr_reg), 1e-12)` | `:465` | verbatim |
| `:560` | `frontier` | I | `nrow(ci_front) == nrow(pf), sum(pf$N == mr$n_selected & abs(pf$hr - mr$naive$est) <  ...` | `:529` | verbatim |
| `:572` | `mr-table` | I | `near(mv$naive_est, c(mr$naive$est, cc$naive$est), 1e-12), near(mv$field_bound, c(fl$ ...` | `:541` | verbatim |
| `:620` | `case-campaign` | P | `nrow(cc_core) == 1, sum(cc_mr$row == "lb_100") == 1, sum(cc_ent$quantity %in% c("tri ...` | `:589` | verbatim |
| `:631` | `case-campaign` | P | `camp$n_trig == cc_ev("trigger_n"), abs(camp$p_joint - cc_ev("joint")) <= 1e-9` | `:600` | verbatim |
| `:633` | `case-campaign` | P | `round(camp$p_trig / 25) == 1, abs(camp$p_joint - 50) < 5` | `:602` | verbatim |
| none | `one-call` | new | `identical(one_call$Hhat, Hhat_def), identical(one_call$in_Hhat, in_Hhat), identical( ...` | `:650` | new (D2) |

## A2: `mrct-consistency-simulations.qmd`

1452 lines to 1404.

**Chunks changed**

| Chunk | Now at | Lines before | Lines after | Change |
|---|---|---|---|---|
| `packages` | `:42` | 22 | 16 | `rd()`/`dd` block to `mrct_record()` |
| `install-check` | `:87` | 16 | 17 | `needed` extended |
| `bridge-case` | `:143` | 63 | 47 | `mrct_trial_analysis()`; downstream objects read from its result |
| `own-build` | `:211` | 37 | 33 | `mrct_build_dgm()`; warnings re-raised from the attribute |
| `own-analytic` | `:274` | 18 | 17 | `ahr_truths_table()`, `ahr_design_table()` |
| `own-summary` | `:343` | 11 | 4 | `mrct_core_summary()` (`eval = FALSE`) |
| `build-fn` | `:448` | 38 | 33 | `mrct_build_dgm()` |
| `build-check` | `:489` | 27 | 27 | `rd(` to `mrct_record(` only |
| `cells` | `:532` | 26 | 26 | `rd(` to `mrct_record(` only |
| `analytic-cells` | `:571` | 27 | 23 | `ahr_truths_table()` |
| `design-n1000` | `:604` | 20 | 18 | `ahr_design_table()` inside `design_tab()` |
| `design-n700` | `:625` | 10 | 10 | `design_tab()`'s new signature |
| `root-r` | `:646` | 13 | 13 | `rd(` to `mrct_record(` only |
| `find-select-check` | `:747` | 12 | 12 | `rd(` to `mrct_record(` only |
| `reproduction` | `:874` | 19 | 19 | `mrct_record()`; `defs` by `system.file()` |
| `mr-repro` | `:898` | 12 | 12 | `rd(` to `mrct_record(` only |
| `core-n1000` | `:947` | 32 | 32 | `rd(` to `mrct_record(` only |
| `core-entry` | `:999` | 18 | 18 | `rd(` to `mrct_record(` only |
| `mr-field` | `:1057` | 7 | 7 | `rd(` to `mrct_record(` only |
| `wall-clock` | `:1109` | 15 | 15 | `rd(` to `mrct_record(` only |
| `demo-summary` | `:1138` | 22 | 17 | `mrct_core_summary()` |

**Prose edits, verbatim** (the `.qmd` outside chunks)

```diff
-recomputes the same objects with the same calls, condensed; the application article is the
-walkthrough.
+recomputes the same objects by the one call of `mrct_trial_analysis()`, on the same arguments; the
+application article is the walkthrough.
-As a reading, not a check: the region's margHR(40) is `r f3(lim_own$region$hr)` here against
+As a reading, not a check: the region's margHR(40) is `r f3(lim_own$marghr[lim_own$set == "region"])` here against
-`ahr_limit()` on the sets you care about gives each set's four truths: read margHR at the analysis
+`ahr_truths_table()` on the sets you care about gives each set's four truths: read margHR at the analysis
-`ahr_design()` at your trial size gives P(ITT ub < 1.0), P(region ITT > 1.0), P(trigger) and the
+`ahr_design_table()` at your trial size gives P(ITT ub < 1.0), P(region ITT > 1.0), P(trigger) and the
+`ahr_truths_table()` and `ahr_design_table()` tabulate `ahr_limit()` and `ahr_design()` over the cells.
-             align = "llrrrr", caption = "The four truths of four sets in the six cells at r = −log 4 (ahr_limit()); equal to the committed predictions to 1e-9.")
+             align = "llrrrr", caption = "The four truths of four sets in the six cells at r = −log 4 (ahr_truths_table()); equal to the committed predictions to 1e-9.")
```

The last pair is the one caption edit in the four articles, inside `analytic-cells`. The `margHR(40)` line changes an inline `r` expression only; it renders the same.

**Removed scaffolding assertions**

- `stopifnot(!is.na(dd), file.exists(file.path(dd, "PROVENANCE.md")))` (baseline `:54`): as in A1, Gate B5 and `tests/testthat/test-mrct-minor-helpers.R:17`.
- `stopifnot(near(four$hr[1:2], c(fit_itt$conf.int[1, "exp(coef)"], hr_reg), 1e-12))` (baseline `:200`) went with `cox_set()`, `fit_itt` and `hr_reg`. This identity is not covered by the wrapper batch; it is covered by the task before it: `tests/testthat/test-mrct-trial-analysis.R:234` (`mrct_four_sets()` against direct `coxph()`) and Gate 2 of `REPORT_mrct_trial_analysis_2026-11-06.md` ("largest |HR, CI or SE difference|: 0"). A1 keeps the same assertion verbatim, and A1's one-call check compares `$sets` against the hand-fitted sets on every render.

**Volatile set (control A / control B), verbatim**

```
tokens A 1454 | tokens B 1454 | volatile positions 5
  #14: 63.7 / 55.6  ::  kthrough. | Code |  patients seconds_mr  |  1000.0 [[63.7]]  | [1] "{z_cd40 <= 236} & 
  #557: 20 / 19  ::  e analytic computations of this section took [[20]] seconds. | 3.4 Finding a c
  #568: 147.6 / 148.4  ::  ode | replicates workers seconds  |  100.0 100.0 [[147.6]]  | 2. Select the trial. mr
  #607: 6 / 7  ::  nce 0) are those of Section 1. The call took [[6]] seconds. | NotePlatform pr
  #1252: 4 / 3.9  ::  ting characteristics. The whole article took [[4]] minutes to run. | 5 Append
```

**Declared timing positions outside the strict set, frozen before any edited render**

```
#13: 1000.0  ::  the walkthrough. | Code |  patients seconds_mr  |  [[1000.0]] 63.7  | [1] "{z_cd40 <= 23
#38: 0.400000  ::  1.386294  |  build_warnings seconds  |  3.000000 [[0.400000]]  | Every patient’s log haz
#49: 0.4  ::  ation and 19.4% of the trial. The build took [[0.4]] seconds. Its 3 warnings,
#162: 3e-01  ::  lt here (build_cell()). | Code |  seconds rows  |  [[3e-01]] 2e+04  |  latent_R_max_abs
#163: 2e+04  ::  e (build_cell()). | Code |  seconds rows  |  3e-01 [[2e+04]]  |  latent_R_max_abs_diff 
#566: 100.0  ::  ion reads. | Code | replicates workers seconds  |  [[100.0]] 100.0 147.6  | 2. Select t
#567: 100.0  ::  ads. | Code | replicates workers seconds  |  100.0 [[100.0]] 147.6  | 2. Select the tri
#1206: 148  ::  n’s arguments and call seed, on 100 workers, [[148]] seconds), summarised int
```

**Gate R, verbatim**

```
baseline tokens 1454 | final tokens 1454 | strict volatile positions 5 (moved 4) | declared timing positions 13
timing-declared positions outside the strict set that moved: 3
    #566: 100.0 -> 100  ::  ion reads. | Code | replicates workers seconds  |  [[100.0]] 100.0 147.6  | 2. Select t
    #567: 100.0 -> 100  ::  ads. | Code | replicates workers seconds  |  100.0 [[100.0]] 147.6  | 2. Select the tri
    #1206: 148 -> 146  ::  n’s arguments and call seed, on 100 workers, [[148]] seconds), summarised int
non-volatile differences: 0
```

**Figures (not gating)**

```
figures baseline 1 | final 1
  fig 1: same (84b6e5c66249e67d / 84b6e5c66249e67d)
```

**Protected-check inventory** (`grep -c 'stopifnot('`: 41 before, 39 after). P = against a shipped record, committed value, published table or a reading's stated claim; I = internal identity; S = scaffolding.

| Baseline | Chunk | Kind | Asserts | Now | Status |
|---|---|---|---|---|---|
| `:54` | `packages` | S | `!is.na(dd), file.exists(file.path(dd, "PROVENANCE.md"))` | none | removed |
| `:154` | `bridge-case` | P | `length(seed_rep5) == 1, !anyDuplicated(trial$id), all(c("id", "y_sim", "event_sim",  ...` | `:149` | verbatim |
| `:158` | `bridge-case` | P | `itt_ub < 1` | `:167` | verbatim |
| `:163` | `bridge-case` | P | `trigger` | `:170` | verbatim |
| `:164` | `bridge-case` | I | `RNGkind()[1] == "Mersenne-Twister"` | `:151` | verbatim |
| `:183` | `bridge-case` | I | `!is.null(fs$sg.harm), mr$ci_method == "field", mr$settings$draws == 2000, mr$field$n ...` | `:174` | verbatim |
| `:186` | `bridge-case` | I | `identical(in_Hhat, trial$id %in% c(fs$df.est$id[fs$df.est$treat.recommend == 0], fs$ ...` | `:177` | verbatim |
| `:200` | `bridge-case` | S | `near(four$hr[1:2], c(fit_itt$conf.int[1, "exp(coef)"], hr_reg), 1e-12)` | none | removed |
| `:206` | `bridge-case` | I | `near(mv$naive_est, c(mr$naive$est, cc$naive$est), 1e-12), near(mv$field_bound, c(fl$ ...` | `:185` | verbatim |
| `:260` | `own-build` | I | `isTRUE(dgm_own$mrct$region$observed), dgm_own$mrct$region$z_name == "z_region", "z_c ...` | `:235` | verbatim |
| `:315` | `own-analytic` | P | `length(v) == 1` | `:289` | verbatim |
| `:547` | `build-check` | P | `nrow(br) == 1, setequal(setdiff(num_cols, names(ours)), c("prevalence_trials", "seed ...` | `:509` | verbatim |
| `:583` | `cells` | P | `!anyNA(cells$cell), near(bfull$h, cells$h[match(sub("_rmlog4$", "", bfull$cell), cel ...` | `:545` | verbatim |
| `:623` | `analytic-cells` | P | `length(v) == 1` | `:582` | verbatim |
| `:627` | `analytic-cells` | P | `near(unlist(truths[c("marghr_40", "marghr_inf", "ahr", "cde")]), tr_ref)` | `:586` | verbatim |
| `:631` | `analytic-cells` | P | `reg_S2$marghr_40 < 1, reg_S2$cde > 1` | `:590` | verbatim |
| `:657` | `design-n1000` | P | `all(abs(d1000[[q]] - vapply(d1000$ref, function(cc) pv(cc, q, ""), 0)) < 1e-9)` | `:614` | verbatim |
| `:670` | `design-n700` | P | `length(v) == 1` | `:627` | verbatim |
| `:673` | `design-n700` | P | `all(abs(d700[[q]] - vapply(d700$ref, pv7, 0, q = q)) < 1e-9)` | `:630` | verbatim |
| `:676` | `design-n700` | P | `all(d700$p_B[d700$cell %in% reg_lt1] > d1000$p_B[d1000$cell %in% reg_lt1])` | `:633` | verbatim |
| `:695` | `root-r` | P | `nrow(x) == 1` | `:652` | verbatim |
| `:698` | `root-r` | P | `sol_r$n_roots == 1, abs(sol_r$roots - ref0$r) < 1e-6` | `:655` | verbatim |
| `:717` | `root-entry` | P | `sol_a$n_roots == 1, abs(sol_a$roots - d1r$a) < 1e-6, sol_r6$n_roots == 1, abs(sol_r6 ...` | `:674` | verbatim |
| `:726` | `root-entry` | P | `abs(p[["p_itt_ub_lt1"]] - x$p_itt_ub_lt1) < 1e-6, abs(p[["p_region_itt_gt1"]] - x$p_ ...` | `:683` | verbatim |
| `:792` | `find-select-check` | P | `nrow(case) == 1, case$sim == 5` | `:749` | verbatim |
| `:794` | `find-select-check` | I | `nrow(demo) == nsims_demo, identical(demo$sim, seq_len(nsims_demo))` | `:751` | verbatim |
| `:798` | `find-select-check` | P | `near(demo$hr_itt[match(k5, demo$sim)], trials5$hr_itt[match(k5, trials5$sim)], 1e-10 ...` | `:755` | verbatim |
| `:833` | `find-check` | P | `RNGkind()[1] == "Mersenne-Twister", seed_S2 + case$sim == seed_rep5, cs$call_args$si ...` | `:790` | verbatim |
| `:891` | `hhat-vs-hstar` | I | `identical(cf_non, cs$confusion$nonregion), identical(cf_reg, cs$confusion$region)` | `:848` | verbatim |
| `:901` | `hhat-vs-hstar` | P | `hr4[["region_Gstar"]] < hr4[["region_Ghat"]], hr4[["region_Ghat"]] < 1, hr4[["itt_re ...` | `:858` | verbatim |
| `:923` | `reproduction` | P | `length(l) == 1` | `:880` | verbatim |
| `:927` | `reproduction` | P | `identical(rec$set, four$set), all(rec$sim == 5), identical(f3(four$hr), rec$hr_3dp), ...` | `:884` | verbatim |
| `:944` | `mr-repro` | P | `nrow(mf) == 1, mf$sim == 5, identical(mf$sg_found, Hhat_def), near(c(mv$naive_est[1] ...` | `:901` | verbatim |
| `:950` | `mr-repro` | P | `fl$lower_1s < 1, fl$complement$upper_1s < 1, mr$naive$est > mr$debiased$est` | `:907` | verbatim |
| `:976` | `rep5-limits` | P | `lims$itt_region$hr < 1, z[four$set == "itt_region"] > 0` | `:933` | verbatim |
| `:1005` | `core-n1000` | P | `nrow(x) == 1, x$replicates == 10000, lead(g("r1")) == p1(of, x$replicates), lead(g(" ...` | `:962` | verbatim |
| `:1016` | `core-n1000` | P | `!is.unsorted(xs$n_trigger), all(xs[c("S2", "S3", "S4"), "pct_declared"] > 99.9), all ...` | `:973` | verbatim |
| `:1046` | `core-entry` | P | `abs(gv("trigger_pct", k) - 100 * gv("trigger_n", k) / gv("n_trials", k)) < 1e-9, abs ...` | `:1003` | verbatim |
| `:1055` | `core-entry` | P | `gv("trigger_n", "S2") == x2$n_trigger, abs(gv("joint", "S2") - 100 * x2$n_declared_G ...` | `:1012` | verbatim |
| `:1105` | `mr-field` | P | `max(vapply(c("N0", "N1", "S1"), mv_, 0, row = "lb_090")) < 1` | `:1062` | verbatim |
| `:1161` | `wall-clock` | P | `identical(mft$cell, tm$cell)` | `:1118` | verbatim |

`demo-summary` keeps its exact row labels and the committed campaign column (`wide$S2`).

## A3: `actg175-glm-outcomes.qmd`

687 lines to 658.

**Chunks changed**

| Chunk | Now at | Lines before | Lines after | Change |
|---|---|---|---|---|
| `regions` | `:206` | 15 | 14 | `in_region()` to `fs_in_subgroup()` |
| `s16-helper` | `:236` | 45 | 44 | `s16_column()` on `mr_values()` |
| `declaration` | `:344` | 35 | 30 | `s17_column()` on `fs_declaration_table()` |
| `s3-selected` | `:474` | 13 | 13 | `fs_in_subgroup()` |
| `s3-columns` | `:532` | 23 | 9 | `mr_method_comparison()` |
| `cost` | `:588` | 17 | 11 | `mr_cost_table()` |

**Prose edits, verbatim** (the `.qmd` outside chunks)

```diff
-Every number of Table S16 but the whole-trial row is in `fit$mr_inference`:
-the unadjusted and corrected estimates and the selection bias in its top level,
-the bounds in `$field`, and the complement's quantities in `$complement` and
-`$field$complement`. The helper below reads them in the table's order, on the
-outcome's own scale.
+Every number of Table S16 but the whole-trial row is in `fit$mr_inference`, and
+`mr_values()` returns them as a row for the region and a row for its complement.
+The helper below reads them in the table's order, on the outcome's own scale.
-$\alpha = 0.10$. `fs_declaration_calibration()` reads the answer from the
-multiplier draws each fit already made; nothing is redrawn (see
+$\alpha = 0.10$. `fs_declaration_table()` tabulates `fs_declaration_calibration()`, which reads the
+answer from the multiplier draws each fit already made; nothing is redrawn (see
```

**Removed scaffolding assertions**

None.

**Volatile set (control A / control B), verbatim**

```
tokens A 669 | tokens B 669 | volatile positions 19
  #66: 75 / 74  ::  nary !{cd40 <= 380} & !{wtkg <= 86} N 72 fit [[75]] s (multiplier stage 13.0
  #67: 13.0 / 12.8  ::  {wtkg <= 86} N 72 fit 75 s (multiplier stage [[13.0]] s) | #> Continuous !{cd40 
  #71: 112 / 111  ::  ontinuous !{cd40 <= 507} & {gender} N 79 fit [[112]] s (multiplier stage 21.6
  #72: 21.6 / 21.8  ::   & {gender} N 79 fit 112 s (multiplier stage [[21.6]] s) | On the binary endpoin
  #304: 56.8 / 56.4  ::  > FS !{cd40 <= 380} & !{wtkg <= 86} N 72 fit [[56.8]] s (multiplier stage 12.0
  #311: 2.4 / 2.5  ::  849.40000000000009} & {cd40 >= 338} N 98 fit [[2.4]] s (multiplier stage 0.8 
  #318: 17.5 / 17.3  ::  0005} & {cd40 > 368.20000000000005} N 95 fit [[17.5]] s (multiplier stage 5.0 
  #334: 45.2 / 45.1  ::  ld be the cores you have. | #> FS DINA GRF  | #> [[45.2]] 20.7 28.6 | Naive, bootstr
  #336: 28.6 / 29.2  ::  cores you have. | #> FS DINA GRF  | #> 45.2 20.7 [[28.6]] | Naive, bootstrap and mul
  #485: 12.04 / 11.95  ::   per replicate | projected at B = 1,000 (s) | FS | [[12.04]] | 25 | 45.18 | 1.81 | 1807 | DINA | 
  #487: 45.18 / 45.15  ::  icate | projected at B = 1,000 (s) | FS | 12.04 | 25 | [[45.18]] | 1.81 | 1807 | DINA | 0.78 | 25 | 2
  #489: 1807 / 1806  ::  cted at B = 1,000 (s) | FS | 12.04 | 25 | 45.18 | 1.81 | [[1807]] | DINA | 0.78 | 25 | 20.72 | 0.83 | 
  #490: 0.78 / 0.81  ::  = 1,000 (s) | FS | 12.04 | 25 | 45.18 | 1.81 | 1807 | DINA | [[0.78]] | 25 | 20.72 | 0.83 | 829 | GRF | 5.
  #492: 20.72 / 20.70  ::  (s) | FS | 12.04 | 25 | 45.18 | 1.81 | 1807 | DINA | 0.78 | 25 | [[20.72]] | 0.83 | 829 | GRF | 5.03 | 10 | 28.
  #494: 829 / 828  ::  4 | 25 | 45.18 | 1.81 | 1807 | DINA | 0.78 | 25 | 20.72 | 0.83 | [[829]] | GRF | 5.03 | 10 | 28.61 | 2.86 | 2
  #495: 5.03 / 5.00  ::  18 | 1.81 | 1807 | DINA | 0.78 | 25 | 20.72 | 0.83 | 829 | GRF | [[5.03]] | 10 | 28.61 | 2.86 | 2861 | The m
  #497: 28.61 / 29.16  ::  1807 | DINA | 0.78 | 25 | 20.72 | 0.83 | 829 | GRF | 5.03 | 10 | [[28.61]] | 2.86 | 2861 | The multiplier
  #498: 2.86 / 2.92  ::  INA | 0.78 | 25 | 20.72 | 0.83 | 829 | GRF | 5.03 | 10 | 28.61 | [[2.86]] | 2861 | The multiplier corr
  #499: 2861 / 2916  ::  .78 | 25 | 20.72 | 0.83 | 829 | GRF | 5.03 | 10 | 28.61 | 2.86 | [[2861]] | The multiplier correctio
```

**Declared timing positions outside the strict set, frozen before any edited render**

```
#305: 12.0  ::  tkg <= 86} N 72 fit 56.8 s (multiplier stage [[12.0]] s) family 2794 | #> DINA {
#312: 0.8  ::  d40 >= 338} N 98 fit 2.4 s (multiplier stage [[0.8]] s) family 147 | #> GRF {wt
#319: 5.0  ::  000000005} N 95 fit 17.5 s (multiplier stage [[5.0]] s) family 1289 | #> FS and
#335: 20.7  ::   the cores you have. | #> FS DINA GRF  | #> 45.2 [[20.7]] 28.6 | Naive, bootstrap an
#488: 1.81  ::  projected at B = 1,000 (s) | FS | 12.04 | 25 | 45.18 | [[1.81]] | 1807 | DINA | 0.78 | 25 | 20.72 | 
#493: 0.83  ::   | 12.04 | 25 | 45.18 | 1.81 | 1807 | DINA | 0.78 | 25 | 20.72 | [[0.83]] | 829 | GRF | 5.03 | 10 | 28.61 | 2.
```

**Gate R, verbatim**

```
baseline tokens 669 | final tokens 669 | strict volatile positions 19 (moved 19) | declared timing positions 25
timing-declared positions outside the strict set that moved: 5
    #305: 12.0 -> 11.6  ::  tkg <= 86} N 72 fit 56.8 s (multiplier stage [[12.0]] s) family 2794 | #> DINA {
    #319: 5.0 -> 5.2  ::  000000005} N 95 fit 17.5 s (multiplier stage [[5.0]] s) family 1289 | #> FS and
    #335: 20.7 -> 19.6  ::   the cores you have. | #> FS DINA GRF  | #> 45.2 [[20.7]] 28.6 | Naive, bootstrap an
    #488: 1.81 -> 1.80  ::  projected at B = 1,000 (s) | FS | 12.04 | 25 | 45.18 | [[1.81]] | 1807 | DINA | 0.78 | 25 | 20.72 | 
    #493: 0.83 -> 0.79  ::   | 12.04 | 25 | 45.18 | 1.81 | 1807 | DINA | 0.78 | 25 | 20.72 | [[0.83]] | 829 | GRF | 5.03 | 10 | 28.61 | 2.
non-volatile differences: 0
```

**Protected-check inventory** (`grep -c 'stopifnot('`: 5 before, 5 after). P = against a shipped record, committed value, published table or a reading's stated claim; I = internal identity; S = scaffolding.

| Baseline | Chunk | Kind | Asserts | Now | Status |
|---|---|---|---|---|---|
| `:81` | `data` | P | `nrow(actg) == 1083, !anyNA(actg$cd420)` | `:81` | verbatim |
| `:211` | `regions` | P | `identical(in_bin,  df_bin$wtkg > 86 & df_bin$cd40 > 380), identical(in_cont, df_cont ...` | `:210` | verbatim |
| `:244` | `s16-helper` | I | `!anyNA(in_H), sum(in_H) == g$n_selected, sum(!in_H) == cc$n` | `:240` | re-expressed |
| `:352` | `declaration` | I | `identical(cal$c0$source, "capture")` | `:347` | re-expressed |
| `:486` | `s3-selected` | I | `identical(in_region(fit_s3_fs, df_bin), in_bin)` | `:477` | re-expressed |

The hard-coded published tables (`s16-published`, `s17-published`, `s3-published`) and the `*-values` chunks are byte-identical to baseline, as are the captions' framing and the `kable` column names.

## A4: `survival-post-selection.qmd`

504 lines to 482.

**Chunks changed**

| Chunk | Now at | Lines before | Lines after | Change |
|---|---|---|---|---|
| `three-columns` | `:307` | 26 | 10 | `mr_method_comparison()` |
| `cost` | `:366` | 17 | 11 | `mr_cost_table()` |

**Prose edits, verbatim** (the `.qmd` outside chunks)

None, as the task specifies.

**Removed scaffolding assertions**

None.

**Volatile set (control A / control B), verbatim**

```
tokens A 530 | tokens B 530 | volatile positions 19
  #58: 11.3 / 11.5  ::  e 0.4 s) family 84 | #> GRF {er <= 0} N 82 fit [[11.3]] s (multiplier stage 2.3 
  #81: 36.2 / 40.9  ::  ow to the cores you have. | #> FS DINA GRF  | #> [[36.2]] 26.3 24.9 | The package’s 
  #82: 26.3 / 26.8  ::   the cores you have. | #> FS DINA GRF  | #> 36.2 [[26.3]] 24.9 | The package’s own s
  #83: 24.9 / 25.1  ::  cores you have. | #> FS DINA GRF  | #> 36.2 26.3 [[24.9]] | The package’s own summar
  #91: 0.59 / 0.66  ::  ------- | #> Overall: | #> Total bootstrap time: [[0.59]] minutes (0.01 hours) | #> 
  #93: 0.02 / 0.03  ::  nutes (0.01 hours) | #> Average per iteration: [[0.02]] min (1.4 sec) | #> Project
  #94: 1.4 / 1.6  ::  1 hours) | #> Average per iteration: 0.02 min ([[1.4]] sec) | #> Projected for 10
  #96: 23.41 / 26.56  ::  2 min (1.4 sec) | #> Projected for 1000 boots: [[23.41]] min (0.39 hrs) | Treatment
  #97: 0.39 / 0.44  ::  sec) | #> Projected for 1000 boots: 23.41 min ([[0.39]] hrs) | Treatment Effect by
  #328: 36.16 / 40.91  ::  licate | projected at B = 1,000 (s) | FS | 5.31 | 25 | [[36.16]] | 1.45 | 1447 | DINA | 0.38 | 25 | 2
  #329: 1.45 / 1.64  ::   | projected at B = 1,000 (s) | FS | 5.31 | 25 | 36.16 | [[1.45]] | 1447 | DINA | 0.38 | 25 | 26.32 | 
  #330: 1447 / 1637  ::  ected at B = 1,000 (s) | FS | 5.31 | 25 | 36.16 | 1.45 | [[1447]] | DINA | 0.38 | 25 | 26.32 | 1.05 | 
  #333: 26.32 / 26.77  ::   (s) | FS | 5.31 | 25 | 36.16 | 1.45 | 1447 | DINA | 0.38 | 25 | [[26.32]] | 1.05 | 1053 | GRF | 2.32 | 10 | 24
  #334: 1.05 / 1.07  ::  S | 5.31 | 25 | 36.16 | 1.45 | 1447 | DINA | 0.38 | 25 | 26.32 | [[1.05]] | 1053 | GRF | 2.32 | 10 | 24.94 | 2
  #335: 1053 / 1071  ::  1 | 25 | 36.16 | 1.45 | 1447 | DINA | 0.38 | 25 | 26.32 | 1.05 | [[1053]] | GRF | 2.32 | 10 | 24.94 | 2.49 | 2
  #336: 2.32 / 2.31  ::  6 | 1.45 | 1447 | DINA | 0.38 | 25 | 26.32 | 1.05 | 1053 | GRF | [[2.32]] | 10 | 24.94 | 2.49 | 2494 | The m
  #338: 24.94 / 25.09  ::  447 | DINA | 0.38 | 25 | 26.32 | 1.05 | 1053 | GRF | 2.32 | 10 | [[24.94]] | 2.49 | 2494 | The multiplier
  #339: 2.49 / 2.51  ::  NA | 0.38 | 25 | 26.32 | 1.05 | 1053 | GRF | 2.32 | 10 | 24.94 | [[2.49]] | 2494 | The multiplier corr
  #340: 2494 / 2509  ::  38 | 25 | 26.32 | 1.05 | 1053 | GRF | 2.32 | 10 | 24.94 | 2.49 | [[2494]] | The multiplier correctio
```

**Declared timing positions outside the strict set, frozen before any edited render**

```
#46: 38.2  ::  lel. | #> FS {er <= 0} & {size <= 35} N 61 fit [[38.2]] s (multiplier stage 5.3 
#47: 5.3  ::  ize <= 35} N 61 fit 38.2 s (multiplier stage [[5.3]] s) family 1744 | #> DINA {
#53: 1.8  ::  #> DINA {grade3 >= 1} & {pgr <= 10} N 89 fit [[1.8]] s (multiplier stage 0.4 
#54: 0.4  ::  {pgr <= 10} N 89 fit 1.8 s (multiplier stage [[0.4]] s) family 84 | #> GRF {er 
#59: 2.3  ::   {er <= 0} N 82 fit 11.3 s (multiplier stage [[2.3]] s) family 858 | The multip
#92: 0.01  ::  rall: | #> Total bootstrap time: 0.59 minutes ([[0.01]] hours) | #> Average per it
#326: 5.31  ::   per replicate | projected at B = 1,000 (s) | FS | [[5.31]] | 25 | 36.16 | 1.45 | 1447 | DINA | 
#331: 0.38  ::   = 1,000 (s) | FS | 5.31 | 25 | 36.16 | 1.45 | 1447 | DINA | [[0.38]] | 25 | 26.32 | 1.05 | 1053 | GRF | 2
```

**Gate R, verbatim**

```
baseline tokens 530 | final tokens 530 | strict volatile positions 19 (moved 19) | declared timing positions 27
timing-declared positions outside the strict set that moved: 6
    #46: 38.2 -> 35.9  ::  lel. | #> FS {er <= 0} & {size <= 35} N 61 fit [[38.2]] s (multiplier stage 5.3 
    #47: 5.3 -> 5.2  ::  ize <= 35} N 61 fit 38.2 s (multiplier stage [[5.3]] s) family 1744 | #> DINA {
    #53: 1.8 -> 1.6  ::  #> DINA {grade3 >= 1} & {pgr <= 10} N 89 fit [[1.8]] s (multiplier stage 0.4 
    #59: 2.3 -> 2.4  ::   {er <= 0} N 82 fit 11.3 s (multiplier stage [[2.3]] s) family 858 | The multip
    #326: 5.31 -> 5.20  ::   per replicate | projected at B = 1,000 (s) | FS | [[5.31]] | 25 | 36.16 | 1.45 | 1447 | DINA | 
    #331: 0.38 -> 0.39  ::   = 1,000 (s) | FS | 5.31 | 25 | 36.16 | 1.45 | 1447 | DINA | [[0.38]] | 25 | 26.32 | 1.05 | 1053 | GRF | 2
non-volatile differences: 0
```

**Protected-check inventory** (`grep -c 'stopifnot('`: 1 before, 1 after). P = against a shipped record, committed value, published table or a reading's stated claim; I = internal identity; S = scaffolding.

| Baseline | Chunk | Kind | Asserts | Now | Status |
|---|---|---|---|---|---|
| `:94` | `data` | P | `nrow(df) == 686, sum(df$status) == 299` | `:94` | verbatim |

`table-s2` and `reading-values` are byte-identical to baseline, as is the narrative's `fit$mr_inference` reading.

## Inventory totals

61 assertions at baseline: 41 protected, 17 internal, 3 scaffolding. All 41 protected ones are retained verbatim. Of the 17 internal ones, 14 are verbatim and 3 re-expressed. The 3 scaffolding ones are removed with the code they tested. A1's one-call check is the only new assertion.

## Stage 3

`git diff --name-only 7989b80f 772a032d`, before this report's commit:

```
dev/tasks/TASK_articles_wrappers_2026-11-06.md
vignettes/articles/actg175-glm-outcomes.html
vignettes/articles/actg175-glm-outcomes.qmd
vignettes/articles/mrct-consistency-simulations.html
vignettes/articles/mrct-consistency-simulations.qmd
vignettes/articles/mrct-consistency.html
vignettes/articles/mrct-consistency.qmd
vignettes/articles/survival-post-selection.html
vignettes/articles/survival-post-selection.qmd

R/tests/man/NAMESPACE/DESCRIPTION entries: 0
.gitignore entries: 0
```

The four `.qmd` and four HTML are the set Amendment 2 permits. The task doc's entry is the amendment's own record, under the standing protocol; the amendment's list does not name it. All four HTML files are byte-identical to their Gate R renders (`cmp`).

`.Rbuildignore:20` is `^vignettes/articles$`, so all four articles are outside the built package. The test suite and `check()` were not run, as Stage 3 specifies.

## What this does not verify

- **`eval = FALSE` chunks did not run.** The new `mrct_core_summary(res_own)` line in `own-summary` is shown, not executed. The function itself ran in `demo-summary`.
- **Renders ran in scratch copies,** reading the installed records, for baseline and final alike. An in-tree render would read `../../inst/extdata`, which is identical by `diff -rq`, but no in-tree render was run. All renders ran on this one machine.
- **The token gate sees digits only.** Non-numeric text was reviewed by diff, not gated, and figures by hash and pixel diff.
- **No pkgdown or site build was run.**
