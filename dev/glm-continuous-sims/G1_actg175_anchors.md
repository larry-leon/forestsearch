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
