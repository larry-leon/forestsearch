# G1 — ACTG175 data anchors (verbatim)

**Spec:** `SPEC_grounded_actg175_cell.md` § G1 (received 2026-08-09 01:13).
**Source document:** `quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd`
(20854 bytes, mtime 2026-05-20 20:36, unmodified). Line citations are `:NN`.

G1 only. No DGM changes, no `R/` changes, no calibration. G2's picks wait.

**Items 1–3 are reported. Items 4–5 could not be produced; § 4 and § 5 say why,
and § 6 states what that blocks in G2.**

---

## 1. Analysis cohort N

Construction copied verbatim from the document (`:193`, `:205`):

```
N after arm subset (arms 1,3)      = 1083
N after dropping missing cd420     = 1083
  arm treat=1 (ZDV+ddI, arms==1)   = 522
  arm treat=0 (ddI,     arms==3)   = 561
  dropped for missing cd420        = 0
```

**N = 1083.** Note the `!is.na(cd420)` filter at `:205` removes nothing within
arms 1 and 3 — every subject in these two arms has a week-20 CD4. The filter is
not inert in general, only on this subset.

## 2. ITT on the adverse scale

`lm(y_decline ~ treat)` per `:244`, where `y_decline <- cd40 - cd420` (`:208`),
`treat = (arms == 1)` with **no switching** (`:197`). Verbatim:

```
          Estimate  Std. Error      t value        Pr(>|t|)
treat -27.59087836 7.888859547 -3.497448293 0.0004888206195
```

**Estimate −27.59087836, SE 7.888859547.**

Sign reading, since the whole cell depends on it: `y_decline` is adverse
(higher = worse), so a **negative** treat coefficient means the experimental arm
declines *less*. The ITT is a **benefit of 27.59 CD4 cells/mm³**, not a harm.

## 3. Residual SD of `y_decline`

```
residual SD from lm (sigma)        = 129.7229086190  (df = 1081)
marginal SD of y_decline (overall) = 130.3944911259
  SD within treat=1                = 144.2769190592  (n = 522)
  SD within treat=0                = 114.5336553143  (n = 561)
  mean y_decline treat=1           = -54.4482758621
  mean y_decline treat=0           = -26.8573975045
```

The arms are **not** homoscedastic: 144.28 vs 114.53, a variance ratio of 1.59.
The calibrator plants a mean shift on a homoscedastic Gaussian; this cell's real
data does not have that shape. Worth stating in the cell header rather than
discovering later.

## 4. Tidy selected-subgroups table, `maxeffCons` arm — NOT AVAILABLE

Two independent reasons, either sufficient:

**(a) The document has never been rendered in this tree.** There is no
`analysis_actg175_continuous_compare_all.html`, and the `_data/` directory it
writes to does not exist. The bundles the spec offers as the alternative source
— `selected_subgroups_continuous.rds` and `comparison_continuous.rds`
(`:599-602`, written to `output_dir` defined at `:147-150`) — are absent.
`quarto/applications/actg175/` contains `.rds` payloads for the *binary*
multimethod analyses only.

**(b) `maxeffCons` is not one of the document's arms.** The axes are, verbatim
(`:170-171`):

```r
sg_focus_vec       <- c("effMaxSG", "effMaxSG", "effMinSG", "effMinSG", "eff", "maxSG","minSG")
selection_rule_vec <- c("pareto",   "both", "pareto", "both", "neighborhood", "neighborhood", "neighborhood")
```

Seven combinations, none of them `maxeffCons`. So rendering the document as
committed would still not produce the requested row — it would need a new arm
added.

## 5. De-biased (MR/FB) estimate for a selected region — NOT AVAILABLE, BY DESIGN

The document excludes it explicitly (`:39-42`):

> Bootstrap, cross-validation, and the combined forest plot are
> **deliberately omitted** so the document stays focused on the
> *selection mechanism* and the comparison reads as a tight diagnostic
> of the rule axis alone.

There is no `forestsearch_bootstrap_dofuture()` call and no `mr_inference`
argument anywhere in it. The split-derived CIs it does compute
(`compute_cis = TRUE`, `:273`) are frontier-diagnostic CIs, not a bias-corrected
estimate, and are not a substitute.

## 6. What this blocks in G2 — needs Larry's call

G2 sets `target_effect` as (spec § G2):

> the de-biased value if G1 item 5 exists, else a stated fraction of the naive

**Neither operand exists.** Item 5 is absent by the document's design, and the
naive selected-region MD is absent because items 4(a) and 4(b) both bite. So the
magnitude decision has no anchor yet, and G2 cannot be written from this table
alone.

Three ways forward, cheapest first. This is a pick, not a recommendation I can
make for you:

1. **Anchor on the ITT instead of a selected region.** Available now: −27.59
   (benefit) with SE 7.89 and residual SD ~130. A planted harm on the
   `y_decline` scale stated as a multiple of the residual SD is defensible and
   needs no render. The committed −40 cell is 0.31 residual SD, for scale.
2. **Render compare_all as committed** (~7 combos, `use_lasso` and `use_grf`
   both `TRUE`, `fs.splits = 500L`) and take the naive selected-region MD from
   the `effMaxSG` arms as the closest available proxy to `maxeffCons`. Gives
   item 4 in proxy form; still no item 5.
3. **Add a `maxeffCons` arm to compare_all and render.** Gives item 4 exactly as
   specified. Still no item 5 unless a bootstrap arm is also added, which the
   document deliberately excludes.

---

## 7. Configuration anchors carried into G2

Recorded because G2's config block cites them. Verbatim.

```r
# :106-108  — the analysis's own floors, adopted by G2
fs_md_threshold    <- 10    # screening: MD >= 10 indicates harm
fs_md_consistency  <- 5     # consistency: MD >= 5

# :92-95, :288 — outcome framing
adverse_outcome    <- "y_decline"   # FS adverse-scale outcome
beneficial_outcome <- "y_improve"   # GRF beneficial-scale outcome
  adverse_outcome  = TRUE,          # at the call site

# :87-89 — covariates
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")

# :83
analysis_seed <- 8316951L
```

Two carry-over hazards, flagged not fixed:

**`y_improve` is declared but never built.** `:95` names it, but the
data-preparation chunk constructs only `y_decline` and sets
`outcome_vars <- c("y_decline")` (`:208-210`). Nothing in the document creates
or references `y_improve`. The dual-scale convention is a change to be
*written*, not one to be *adopted* from here.

**`use_lasso = TRUE` / `use_grf = TRUE` (`:119-120`) are MR-incompatible.** The
committed sweep pins both `FALSE` with an explicit note
(`mr_coverage_sweep_md_harm.qmd:128-132`): the MR alignment guard rejects the
call before fitting when either is `TRUE`, which produced the pilot's detected
0 of 50. G2 adopts this analysis's *thresholds*; it must not also adopt its
*front end*, or the MR arm dies silently.

**Treatment coding differs from the committed sim harness.** compare_all uses
`treat = (arms == 1)`, "NO switching" (`:197`); the committed twin uses
`treat = 1 - (arms == 1)`, ddI = 1. Combined with the outcome-column difference
(`cd40 - cd420` here vs `cd420 - cd40` there) and `adverse_outcome` `TRUE` vs
`FALSE`, the two reach the same harm orientation by two different routes. A
threshold or sign copied between them is only valid because those two
differences cancel.

---

## 8. Route 3 executed — the `maxeffCons` real-data anchors

Step 1 of the 2026-08-09 overnight set. `maxeffCons` added to the document's
axes; `consistency_method` pinned to `"resample"`; J-quantile grids specified on
all six continuous covariates. Run at `413c61c6`.

**Scoped render, not the full axes.** `compare_selection_rules()` accepts
length-1 axis vectors, so the `maxeffCons` arm was run alone rather than
re-running all eight. Wall time is per-arm and reported per row below.

**MR required the front ends off.** `.validate_mr_configuration()`
(`R/forestsearch_helpers.R:1909-1927`) `stop()`s under `consistency` when any of
`use_lasso`/`use_grf`/`use_dina` is TRUE. The document sets `use_lasso` and
`use_grf` TRUE, so the anchor runs with all three FALSE — the MR-compatible
configuration of the committed sweep. This is a real difference from the
document's other eight arms and is why these rows are labelled as their own
configuration, not as "the compare_all result".

### Both anchor rows, verbatim

| | **J = 10 on all six** | **default cuts** |
|---|---|---|
| definition | `!{cd40 <= 507} & {age <= 37}` | `!{cd40 <= 421} & {age <= 34}` |
| n (selected) | 66 | 137 |
| naive MD | **87.916666667** | **52.156436488** |
| naive 95% CI | 4.305358495, 171.527974839 | 1.805091251, 102.507781724 |
| Pcons | **0.95** | **0.93** |
| MR de-biased MD | **33.6147229870208** | **−3.59145136776183** |
| MR 95% CI | −72.0371931311213, 139.266639105163 | −62.5877541865567, 55.404851451033 |
| MR se_ij | 53.9050293533508 | 30.100707606951 |
| MR se_wald | 42.6596145804699 | 25.6899339140523 |
| ij_source / ij_draws | ij / 4858 | ij / 4176 |
| `harm_flag` | TRUE | **FALSE** |
| `selection_rate` | 0.9716 | 0.8352 |
| MR `n_family` | 4935 | 1450 |
| cut labels in family | 54 | 28 |
| candidate subgroups | 128, **truncated to 30** | 23 |
| candidates evaluated / passed | 1 / 1 | 2 / 1 |
| wall time | 57.5 s | 18.6 s |

Raw consistency rows:

```
##### J = 10 : n_cand_eval=1 n_cand_total=30 n_passed=1 #####
  Pcons          hr  N  E   g m K            M.1         M.2
1  0.95 87.91666667 66 66 121 1 2 !{cd40 <= 507} {age <= 37}

##### default : n_cand_eval=2 n_cand_total=23 n_passed=1 #####
  Pcons          hr   N   E  g m K            M.1         M.2
1  0.93 52.15643649 137 137 65 2 2 !{cd40 <= 421} {age <= 34}
```

### Four things to read before picking a magnitude

**(1) MR removes most of the effect, and under default cuts removes all of it.**
87.92 → 33.61 and 52.16 → **−3.59**. Both MR intervals cover zero. Under default
cuts the de-biased point estimate is *negative* and `harm_flag` is FALSE: after
correcting for selection there is no harm in the selected region at all. The
naive figures are the selection artefact the whole exercise is about.

**(2) The two cut specs are not small perturbations of each other.** Different
regions (`cd40 <= 507 / age <= 37` vs `cd40 <= 421 / age <= 34`), n 66 vs 137,
and de-biased estimates on opposite sides of zero. The cut grid is not a tuning
detail here; it selects a different scientific claim. Both rows are reported
because neither is obviously the right one.

**(3) The J = 10 family was truncated, so it did not do what it was for.**

```
Warning message:
max_subgroups_search = 30 truncated the candidate pool from 128 to 30.
Candidates ranked below 30 under the sg_focus = 'maxeffCons' preview ordering
were not evaluated.
```

The grids widened the family from 28 cut labels to 54, and the subgroup pool
from 23 to 128 — then `fs_max_subgroups <- 30` (`:116`, unchanged, as
instructed) discarded 98 of the 128 before evaluation. The stated intent, "so
the fixed family can explore the range of the continuous factors", is only
partly served: the range is *generated* but not *evaluated*. Raising
`max_subgroups_search` is a change nobody authorised, so it was not made.

**(4) `karnof` gets one cut at J = 10, not ten.** Its grid is
`karnof <= 90` alone, because the variable takes few distinct values (the
default spec produced two, `<= 95` and `<= 90`). J = 10 is a request, not a
guarantee; on near-discrete covariates the quantile grid collapses. So the
"ten thresholds each" description holds for age/wtkg/cd40/cd80, gives seven for
`preanti`, and one for `karnof`.

### What this unblocks

G1 item 5 now exists for the `maxeffCons` arm, so § 6's blocker is lifted: G2's
`target_effect` can be anchored on a de-biased value as the spec intends. Which
of the two — 33.61 or −3.59 — and on which cut spec, is the magnitude pick.
Note that the default-cuts row cannot anchor a *harm* magnitude at all, since
it is negative.

G2 remains not started.
