# REPORT — Directive C: the DINA guard, the frontier-key warning, honest displays

**Task:** `dev/tasks/TASK_directive_C_2026-09-18.md` (committed as received, `bef180ce`).
**Branch:** `feature/glm-extension`. **Not pushed.** Five commits, each revertable alone:

| Part | Commit | Class |
|---|---|---|
| 1 DINA identity-scale guard | `c5bd3b5d` | behaviour (new error on an opt-in combination) |
| 2 frontier-key warning | `b8243432` | behaviour (messaging) |
| 3 `Inf` display caps + trim warning | `31255426` | behaviour (display object only) |
| 4 frontier print retitle | `c7e95160` | messaging |
| 5 c2 / p\* echo annotation | `c8b475b4` | display |

The docs task's report (`REPORT_threshold_docs_2026-09-18.md`) is **not** in `dev/reports/`; everything
below was re-found by search from source, as the task requires.

## Part 1 — the route list

`m_diff` is derived from `hr.threshold` in exactly two places. Established by scanning every function body
in the namespace for the derivation (asserted in the tests, so a third site cannot appear unnoticed):

1. **`.forestsearch_dina_select()`** — `R/forestsearch_helpers.R:1525` — the `subgroup_method = "dina"`
   selection path (called from `R/forestsearch_main.R:2596`).
2. **`forestsearch()` SECTION 3B** — `R/forestsearch_main.R:3051` — the `use_dina = TRUE` +
   `dina_args$selected_only = TRUE` screening branch.

`dina_frontier()` takes `m_diff` from the caller and derives nothing (`R/dina_subgroup.R:1193`), so it is
not a third route. `dina_subgroup()`, `dina_bagged()`, `dina_subgroup_refit()` and
`dina_subgroup_bootstrap()` all receive `m_diff` as an argument.

The guard `.dina_assert_ratio_estimand()` (`R/forestsearch_helpers.R:1446`) is called at both sites
(`forestsearch_helpers.R:1524`, `forestsearch_main.R:3050`), **not** at the `subgroup_method` dispatch. At
site 1 the derivation was hoisted above the `dina()` fit block so the refusal costs no model run; at site 2
the derivation already precedes nothing but sits inside the screening `tryCatch`, so the refusal surfaces
as the existing `"DINA analysis failed: ..."` warning and the run continues with no DINA cuts — correct for
an optional screening path, and the consistency search (which handles RD properly) proceeds.

Condition: `effect_measure %in% c("RD", "IRD")` **and** `family != "gaussian"`. Gaussian is provably
untouchable — asserted over the full 4x7 (family, measure) grid: the guard fires on exactly 6 of 28 cells,
none of them gaussian.

### Exposure inventory 1 — identity-scale DINA callers

Every tracked file containing `subgroup_method = "dina"` or `use_dina = TRUE` (62 files) was checked for its
resolved `effect_measure`. Result: **zero committed callers reach DINA with an identity-scale estimand on a
non-gaussian family.** The measures found are `"OR"` (all `quarto/simulations/actg175/binary_methods/*`,
`quarto/dina/method_equivalence_checks.qmd`), `"MD"` on gaussian
(`quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd` — the mddina campaign),
and survival `HR` (the gbsg and `quarto/dina/*` files, which set no `effect_measure`). The per-outcome
defaults (`forestsearch_main.R:1612-1617`) are `binary -> "OR"`, `continuous -> "MD"`, `count -> "IRR"`, so
`"RD"`/`"IRD"` can only arrive by explicit request. **Gate passed; no pin-vs-proceed needed.**

## Part 2 — the frontier-key warning

The seven keys move into `.DINA_FRONTIER_KEYS` (`R/forestsearch_helpers.R:1055`) so `.resolve_dina_args()`'s
key list and the warning cannot drift. The warning fires once per fit in `.forestsearch_dina_select()`
(`R/forestsearch_helpers.R:1498-1507`), listing every offending key in one message. Inertness re-verified
from source: that function's only `dina_frontier()` call (`forestsearch_helpers.R:1553`) passes its own
`scope = "wide"` and `n_min = n.min` and never touches `da$frontier`.

## Part 3 — honest display caps

`dina_frontier()` defaults are now `Inf` / `Inf` (`R/dina_subgroup.R:1195-1196`). A finite cap that actually
removes rows warns via `.dina_warn_cap_trim()` (`R/dina_subgroup.R:1034`), naming kept vs. available, with
condition class `dina_frontier_cap_trim`.

### Exposure inventory 2 — `dina_frontier()` callers

All **20** committed in-repo call sites pass **both** caps explicitly, so **no committed rendered table
changes on re-render**:

`dev/identifier-alignment/sim_analyses/analysis_actg175_binary_multimethod.qmd:1163` ·
`.../analysis_gbsg_survival_multimethod.qmd:1449` · `dev/replication-check/legacy_v2_2A_reconstructed.qmd:1085` ·
`dev/replication-check/v2_2new_rendered_source_prerename.qmd:1086` ·
`dev/review/analysis_gbsg_survival_multimethod.qmd:1382` ·
`quarto/applications/actg175/analysis_actg175_binary_multimethod_{fixed_family:1154, frontend:1147, psi_v2_2:1132, psi_v3a:1132}.qmd` ·
`quarto/applications/actg175/_archive/20260730_..._{psi_v2_2A:1125, psi_v2_2A_mac:1125, psi_v2_2_mac:1125, psi_v2_2_mac_w2:1128}.qmd` ·
`quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd:1435` ·
`quarto/applications/gbsg/_archive/{2026-05-30_gbsg_survival_dina_cox:709, 2026-05-31_..._effMaxSG:876, 2026-05-31_..._maxSG:876, 2026-06-01_..._depth2_DINA:318, 2026-07-30_..._effMinSG:986, 2026-07-30_..._multimethod:1075}.qmd`.
**Gate passed.** Their `# default 3L` / `# default 10L` source comments are now stale; not edited (scope).

`.resolve_dina_args()` still supplies the finite `3L` / `10L` for `use_dina` screening
(`R/forestsearch_helpers.R:1160-1161`), so **the screening pool is byte-identical** — no selection logic
moved. Under the default `selected_only = TRUE` that frontier is computed and discarded, so a trim warning
there would name a table nobody sees; per your decision it is muffled by condition class at exactly that one
call (`R/forestsearch_main.R:3029-3038`) and stands under `selected_only = FALSE`, where the caps really do
shape the pool. No committed caller sets `selected_only = FALSE` in code.

## Part 4 — the retitle

`R/forestsearch_helpers.R:1569-1584`. `"DINA frontier candidates (per-covariate non-dominated):"` becomes
`"DINA frontier -- proposed single cuts (per-covariate non-dominated); display only, not the searched
family:"`, shown beside the unchanged `"Candidates searched / qualifying"` family counts; the empty case is
retitled to match. String change only. `forestsearch_main.R:3089`'s `"frontier candidates"` is the `use_dina`
**screening mode label**, where under `selected_only = FALSE` the frontier genuinely *is* the candidate pool
— correct as written, left alone.

## Part 5 — the echo-site list

Every place c2 or p\* is echoed, found from source before annotating any. One helper,
`.fs_c2_inert_note()` (`R/forestsearch_helpers.R:1039`), returns `""` on the consistency path.

| # | Site | File:line | Prints under |
|---|---|---|---|
| E1 | GLM config banner | `forestsearch_main.R:2227` | grf, consistency (already excludes dina) |
| E2 | survival config banner | `forestsearch_main.R:2280` | grf, consistency (already excludes dina) |
| E3 | Search Alignment Diagnostic | `interpret_search_config.R:181` | grf (already quiet under dina) |
| E4 | `summary.forestsearch()` | `forestsearch_methods.R:451` | dina, grf, consistency |
| E5 | bootstrap parameter banner | `bootstrap_dofuture_main.R:429` | grf, consistency (dina branch omits them) |
| E6 | `print_cv_params()` | `forestsearch_cross_validation.R:2062` | dina, grf, consistency |
| E7 | `fs_family_report()` consistency-screen row | `fs_family_report.R:326-333` | all; status now `inert` under dina/grf |

`interpret_search_config()` gains `subgroup_method = "consistency"` (a defaulted, backward-compatible
formal), so every existing caller prints exactly what it did. The motivating case —
`quarto/simulations/gbsg_020/summary_grfmr.qmd`, `hr.threshold = 0.90` against the inert default c2 = 1.0
under grf — now says so at E1/E3/E4. Not re-rendered.

`print.forestsearch()`'s `"Consistency: NN%"` (`forestsearch_methods.R:563`) is the *realized* split rate of
the selected subgroup, not a threshold echo; `fs_oc_predict`'s `c2` (`fs_oc_predict.R:364`) is that tool's
own setting. Neither annotated.

## Gates

- **Gate A — passed.** The guard fires on exactly 6 of the 28 (family, measure) cells: `{cox, binomial,
  poisson}` x `{RD, IRD}`. The two new warnings fire only where specified (frontier keys under
  `subgroup_method = "dina"`; a finite `dina_frontier()` cap that removes rows) and nowhere else: a clean
  consistency run, a clean dina run and the same `dina_args` under consistency/grf raise none of them.
- **Gate B — passed.** An `mddina`-shaped call (continuous, `MD`, `family = "gaussian"`,
  `subgroup_method = "dina"`) runs to completion with none of the three conditions raised; the guard is
  provably unable to fire on gaussian.
- **Gate C — passed.** All **175** cells of the threshold-resolution probe
  (`helper-threshold-sync.R`, one copy, unmodified) are **byte-identical** to the pre-edit baseline —
  147 consistency, 14 grf, 14 dina. `dev/reports/baseline_directive_C_2026-09-18.csv` vs
  `postC_directive_C_2026-09-18.csv`, `identical() == TRUE`.

## Acceptance tests and compute

`tests/testthat/test-directive-c.R` — **182 assertions, 0 failures, 5.6 s wall clock** (3-minute abort).
Two seeded micro-fits, both far inside the 5-minute cap: `dina()` gaussian n = 200, **0.27 s** (fitted once,
reused by Parts 3 and 4 through `dina_res`); `forestsearch()` + dina on continuous n = 120, **0.10 s**
(Gate B). Everything else is source inspection, direct helper calls, mocked-binding tripwires or string
assertions. No `R CMD check`, no full suite. No existing test file touches any changed surface (checked).

## Findings, no tasks attached

1. **The `details = TRUE` frontier print on the dina path is now long.** With `Inf` caps it shows the whole
   per-covariate non-dominated set — on continuous covariates that ran to ~200 rows in the micro-fit, where
   it used to show 10. Truthful, which is the point of Part 3, but the diagnostic is now verbose.
2. **The `use_dina` unoriented-floor residual is real and still present** (out of scope, unchanged). At
   `forestsearch_main.R:3054-3063` the screening `dina_subgroup()` call omits `tau_sign`, which
   `.forestsearch_dina_select()` computes via `.dina_tau_sign()` and passes
   (`forestsearch_helpers.R:1546`). Under `adverse_outcome = FALSE` on a binary or continuous outcome the
   two routes therefore orient the floor differently.
3. **`selection_rule` vocabularies disagree across the dina boundary.** `forestsearch()` accepts
   `"hr"`/`"maxSG"`/`"minSG"`/..., `dina_subgroup()` accepts `"neighborhood"`/`"pareto"`/`"both"`
   (`dina_subgroup.R`); a `subgroup_method = "dina"` call carrying a consistency-path `selection_rule`
   errors inside `dina_subgroup()`. Observed while writing the Part 4 test; nothing changed.
4. The 20 `dina_frontier()` call sites above carry `# default 3L` / `# default 10L` comments that Part 3
   made stale.
