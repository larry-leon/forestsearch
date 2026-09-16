# REPORT — Applied analyses under effMaxSG, ε = 0.20: Stage 0, read-only

Date: 2026-09-16. Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1, reference BLAS). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_actg175_applied_effmaxsg_stage0_2026-09-16.md` (committed as received, `3c6745ae`). Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix`. Read-only: no edit, render or install; the one computation is §3's two identification fits from a temporary purl outside the repo (removed at the end). Line numbers at HEAD `570bf8a5`.

## 1. Provenance — GATE PASS

```
pop-os
feature/glm-extension
570bf8a5
570bf8a5 ACTG175 continuous applied analysis: the field-s complement bound and pair ...
98643ba9 Add TASK_actg175_continuous_field_s_2026-09-16 as received
89df2cf3 Merge branch 'feature/glm-extension' of github.com:larry-leon/forestsearch into feature/glm-extension
[tracked modifications: none]
field-s task closed
[R / Rscript / quarto processes: none]
R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix
```
First commit: `3c6745ae Add TASK_actg175_applied_effmaxsg_stage0_2026-09-16 as received`.

## 2. S0.1 — Inventory of applied documents

Every tracked `.qmd` under `quarto/applications/` that calls `forestsearch()` (or `compare_selection_rules()`), generated from source by pattern (`sg_focus` / `effect_neighborhood` / `selection_rule` / `use_dina` / `use_grf` / `mr_inference` assignments in code lines); 48 further documents under `actg175/_archive/` (last commit `dd60b5fe`, 2026-08-13) and `gbsg/_archive/`, `gbsg/_broken/` (`9475181e`, 2026-08-13) are superseded working copies and are not listed row by row. "—" = not set in the document (the package default applies: `sg_focus = "hr"` → `maxcons` on the consistency engine, `effect_neighborhood = 0.10`, `selection_rule = "neighborhood"`).

| document | last commit | `forestsearch(` calls | `sg_focus` | `effect_neighborhood` | `selection_rule` | `use_dina` / `use_grf` | MR (`mr_inference`) |
|---|---|---|---|---|---|---|---|
| `actg175/analysis_actg175_binary_multimethod_fixed_family.qmd` | 43b051b6 2026-08-16 | 4 | `effMaxSG` (per-identifier variables) | 0.10 | `neighborhood` | DINA and GRF runs alongside FS | yes (`run_mr` knob) |
| `actg175/analysis_actg175_binary_multimethod_frontend.qmd` | cf4d6432 2026-08-15 | 4 | `effMaxSG` (per-identifier variables) | 0.10 | `neighborhood` (FS) / `both` | DINA and GRF runs alongside FS | yes (`run_mr` knob) |
| `actg175/analysis_actg175_binary_multimethod_psi_v2_2.qmd` | cf4d6432 2026-08-15 | 4 | `effMaxSG` (per-identifier variables) | 0.10 | `neighborhood` (FS) / `both` | DINA and GRF runs alongside FS | yes (`run_mr` knob) |
| `actg175/analysis_actg175_binary_multimethod_psi_v3a.qmd` | cf4d6432 2026-08-15 | 4 | `effMaxSG` (per-identifier variables) | 0.10 | `neighborhood` (FS) / `both` | DINA and GRF runs alongside FS | yes (`run_mr` knob) |
| `actg175/analysis_actg175_binary_sgfocus.qmd` | 889dc05a 2026-08-15 | 3 | `effMaxSG` / `maxeff` (a focus-sweep document) | 0.10 (swept) | `both` (swept) | FALSE / TRUE | yes (`run_mr` knob) |
| `actg175/analysis_actg175_continuous_compare_all.qmd` | cf4d6432 2026-08-15 | 4 (+ compare_selection_rules) | `maxeffCons` anchor + sweep `effMaxSG/effMinSG/eff/maxSG/minSG/maxeffCons` × `pareto/both/neighborhood` via `compare_selection_rules()` | no | `neighborhood` (anchor) + swept | FALSE / sweep sets `use_grf` | yes |
| `actg175/analysis_actg175_continuous_oc.qmd` | f87518e9 2026-09-08 | 1 | "maxeffCons" | no | "neighborhood" | FALSE / FALSE | yes |
| `actg175/analysis_actg175_continuous_oc_evaluation.qmd` | c73a6b06 2026-08-31 | 0 | "maxeffCons" | no | — | — / — | no |
| `actg175/template_actg175_continuous.qmd` | 61dd99df 2026-08-17 | 1 | `minSG` | no | `neighborhood` | — / FALSE | no |
| `count_data_demo.qmd` | b81e011d 2026-04-14 | 10 | no | — | no | — / TRUE (screening cuts) | no |
| `gbsg/analysis_gbsg_survival_effMaxSG.qmd` | 48c419eb 2026-08-15 | 5 | `effMaxSG` (focus sweep) | 0.10 (swept) | `neighborhood` (swept) | FALSE / FALSE | yes |
| `gbsg/analysis_gbsg_survival_frozen_family.qmd` | 37f38540 2026-09-06 | 4 | `maxeff` | no | — | FALSE / FALSE | yes |
| `gbsg/analysis_gbsg_survival_maxeff_mrconfirm.qmd` | 48c419eb 2026-08-15 | 3 | `maxeff` | no | — | FALSE / FALSE | yes (`run_mr` knob) |
| `gbsg/analysis_gbsg_survival_multimethod.qmd` | 6ad5a724 2026-09-01 | 7 | `effMaxSG` (per-identifier variables) | 0.10 | "neighborhood" | DINA and GRF runs alongside FS | yes (`run_mr` knob) |
| `gbsg/analysis_gbsg_survival_sgfocus.qmd` | 48c419eb 2026-08-15 | 5 | `effMaxSG` (focus sweep) | 0.10 (swept) | `both` (swept) | FALSE / FALSE | yes |
| `validation_glm_simulation_study.qmd` | b81e011d 2026-04-14 | 4 | no | — | no | — / TRUE (screening cuts) | no |
| `validation_hte_tests_crump.qmd` | ede6b683 2026-05-28 | 11 | no | — | no | — / — | no |

**Last recorded renders** (only one ACTG175 document has a render record):
- `actg175/analysis_actg175_continuous_oc.qmd`: `REPORT_actg175_continuous_intervals_2026-09-07.md:39` — Mac-Studio-3 (M4 Max, 36 GB), `-P n_workers:1`, "wall 20 min 49 s; the document's own clock … 20.2 min; evaluation loop 1123.4 s (18.7 min) over 20 jobs, 1 workers; peak summed R/quarto RSS 15,641 MB"; `:116` "≈ 11 GB per worker". **pop-os, 2026-09-16** (`REPORT_actg175_continuous_field_s_2026-09-16.md` §4): committed `n_workers = 14`, wall 3,214 s (53.6 min), evaluation loop 2,089.7 s, peak summed RSS 90,526 MB (≈ 6.5 GB per worker); CPU time not recorded.
- `actg175/analysis_actg175_continuous_oc_evaluation.qmd`: no render record; it reads its anchor and OC payload `S` from `dev/glm-continuous-sims/` (`:31`) and calls `forestsearch()` 0 times.
- `actg175/analysis_actg175_continuous_compare_all.qmd`: no render record; its payload manifest (`_payloads_2026-09-01/.../comparison_continuous.rds.MANIFEST.md`) records version 0.3.5 and a file mtime of 2026-09-01, no wall.
- `actg175/template_actg175_continuous.qmd` and the binary multimethod documents: no render records in the directory.

## 3. S0.2 — The anchor under both rules

`knitr::purl()` of `analysis_actg175_continuous_oc.qmd` into `/tmp/fs_s0_tqaQ` (outside the repo; removed); chunks run verbatim: `setup` (`:11`), `data-prep` (`:57`); then the `anchor` chunk's `forestsearch()` call (`:92–132`) with only these arguments changed: `mr_inference = FALSE` in both fits; (b) `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.20`, `selection_rule = "neighborhood"` (the `mdsgnb20` bundle `meta`: `sg_focus effMaxSG | effect_neighborhood 0.2 | selection_rule neighborhood`, computed here from the committed md40 n500 combined bundle). Data: "Analysis N = 1083; unadjusted ITT mean difference = -27.591 (SE 7.889)" (the document's own `data-prep` output). Orientation: the document works on `y_decline` (positive = harm).

| | (a) `maxeffCons`, as committed | (b) `effMaxSG`, ε = 0.20, `neighborhood` |
|---|---|---|
| selected subgroup | `{age <= 37} & !{cd40 <= 507}` | `!{cd40 <= 507} & {gender}` |
| n(Ĥ) | 66 | 79 |
| oriented MD estimate (the `hr` column) | 87.916667 | 71.095513 |
| consistency proportion (`Pcons`) | 0.95 | 0.91 |
| complement n | 1017 | 1004 |
| rows of `grp.consistency$out_sg$result` (consistency-qualifying candidates carried) | 1 (the winner only) | 8 |

- **GATE (a): PASS** — the extraction reproduces the committed anchor `{age <= 37} & !{cd40 <= 507}` with n = 66 (`setequal` on the clauses; identical label).
- **Family size.** With `details = FALSE` the returned object carries no `hr.subgroups` table, so the enumerated family is not in the fit; the MR gate's payload records `settings$n_family = 4935` for this call (the field-s record §2.1 / the committed payload), and the consistency-qualifying set the band is computed over has 8 members under (b).
- **(b) band:** max oriented MD among the consistency-qualifying candidates 87.916667 (the committed anchor itself); floor (1 − 0.20) × max = **70.333333**; **4 of 8** candidates in the band. The in-band candidates by size (all four; there is no fifth):
  1. `!{cd40 <= 507} & {gender}` — N 79, MD 71.0955, Pcons 0.91 (**selected**: largest in band)
  2. `{age <= 40} & !{cd40 <= 507}` — N 76, MD 78.7054, Pcons 0.93
  3. `{age <= 37} & !{cd40 <= 507}` — N 66, MD 87.9167, Pcons 0.95 (the `maxeffCons` anchor)
  4. `!{wtkg <= 73} & !{cd40 <= 507}` — N 61, MD 82.2688, Pcons 0.93
- **Overlap:** |Ĥ(a)| = 66, |Ĥ(b)| = 79, |Ĥ(a) ∩ Ĥ(b)| = **52**; **Ĥ(b) does not contain Ĥ(a)** (14 of the 66 fall outside; the two rules share the `cd40 > 507` clause and differ on `age ≤ 37` vs `gender`).

## 4. S0.3 — What depends on the anchor or the rule

**`analysis_actg175_continuous_oc.qmd`** (the headline document). The anchor feeds every later section:
- `:137–141` `H_def <- paste(fs_anchor$sg.harm, ...)`; `n_H <- sum(fs_anchor$df.est$treat.recommend == 0L)`; `win <- fs_anchor$grp.consistency$out_sg$result[1L, ]`; `T_obs <- as.numeric(win$hr)`; `p_cons <- as.numeric(win$Pcons)` — **the selected subgroup and quantities derived from it**.
- `:146–147` `stopifnot(setequal(fs_anchor$sg.harm, c("{age <= 37}", "!{cd40 <= 507}")), n_H == 66L, abs(T_obs - (87 + 11 / 12)) < 5e-7)` — **a literal assertion of the committed anchor**: under `effMaxSG` (Ĥ = `!{cd40 <= 507} & {gender}`, n = 79, T̂ = 71.10) the document stops here.
- `:148–150` `q_rungs[length(q_rungs)] <- T_obs; q_shared[...] <- T_obs; c1_ladder[abs(c1_ladder - (87 + 11 / 12)) < 1e-9] <- T_obs` and `:37–39` the literal `87 + 11 / 12` in the rung and ladder definitions — **derived from T̂** (§4–§12, every OC table and figure, the payload's `table`).
- `:154–156` the anchor prose (inline `H_def`, `n_H`, `T_obs`, `p_cons`) — §1.
- `:178–392` the intervals section (`g <- fs_anchor$mr_inference`, its `stopifnot` at `:197–198` on the label, `n_H`, `T_obs`; Tables 1–3; the reading; payload `iv`) — **the selected subgroup** and its MR gate (§2.1).
- `:431` `Q_primary <- list(age = 37, cd40 = list(type = "greater", value = 507))` — **the anchored truth Q, hard-coded to the committed anchor's clauses**; `:452` `dgm_top <- dgm_at(T_obs, Q_primary)`; `:456–482` the Q knob (`age_grid > 37`, `cd40_grid < 507`, the 2× / 3× prevalence variants) — §3 and everything downstream (DGM, family, OC grids, `Q_variants`).
- `:503–511` `fs_args` (`confounders.name`, `conf.cont_jcuts`, `cut_type`, `cont.cutoff`, `maxk`, `n.min`) — **carries no `sg_focus`, `effect_neighborhood` or `selection_rule`**; `:512`, `:573`, `:634` `fs_oc_family_enumerate(..., fs_args, ...)`; `:543` `fs_oc_predict`; `:579`, `:608` `fs_oc_grid` — the OC sections (§4–§12) depend on the **selection rule the OC functions model** (below), not on the document's `sg_focus`.
- `:1243–1313` the payload (`extras$anchor`, `intervals`, `table`) — derived.

**`analysis_actg175_continuous_oc_evaluation.qmd`**: no `forestsearch()` call; reads `S` from `dev/glm-continuous-sims/` (`:31`) and displays `S$meta$anchor$def`, `$n_H`, `$T_obs` (`:64–66`, `:106`, `:276`, `:304`); its prose hard-codes the committed anchor (`:85` "`age = 37`; `cd40 = list(type = "greater", value = 507)`", `:93`, `:385`); its payload records `sg_focus = "maxeffCons"` (`:422`). Depends on **the selected subgroup and T̂ through a `dev/` payload built outside the applied directory**; a switch would need that payload regenerated first.

**`analysis_actg175_continuous_compare_all.qmd`**: runs `compare_selection_rules()` over `sg_focus_vec = c("effMaxSG","effMaxSG","effMinSG","effMinSG","eff","maxSG","minSG","maxeffCons")` × `selection_rule_vec = c("pareto","both","pareto","both","neighborhood","neighborhood","neighborhood","neighborhood")` (`:202–203`) plus a `maxeffCons` anchor; a comparison document — the `effMaxSG` combos it runs are `pareto` and `both`, not the survival grid's `neighborhood`; **depends on the rule set it sweeps**, not on one anchor.

**The OC functions** (`getNamespaceExports`: `fs_oc_family_enumerate`, `fs_oc_grid`, `fs_oc_invert`, `fs_oc_predict`, `fs_family_report`). Formals: `fs_oc_family_enumerate(dgm, forestsearch_args, n, max_M = 2000L, verbose = FALSE)`; `fs_oc_grid(dgm = NULL, forestsearch_args = list(), n, c1, c2, family = NULL, consistency_method = c("resample", "split"), pconsistency = NULL, draws = 2e+05, block = 50000, seed = NULL, verbose = FALSE, ...)`; `fs_oc_invert(...)` and `fs_oc_predict(...)` likewise; `fs_family_report(x, data = NULL, outcome_type = NULL)`.
- **What they read from `forestsearch_args`:** `confounders.name` (`R/fs_oc_family.R:186`), the cut arguments, `maxk`, `minp`, `n.min` / `n.min.frac` (`:230–250`), `effect.threshold` / `consistency.threshold` / `pconsistency.threshold` (`R/fs_oc_grid.R:106–107`, `:285–286`, `:479`; `R/fs_oc_predict.R:173`, `:186–187`); and `sg_focus` **only** to test `identical(.arg("sg_focus"), "maxeff")`, which relaxes the floors (`R/fs_oc_family.R:236–238`). **Neither `effect_neighborhood` nor `selection_rule` is read anywhere in `fs_oc_family.R`, `fs_oc_grid.R` or `fs_oc_predict.R`** (grep); `fs_family_report()` reads `sg_focus`, `selection_rule` and `effect_neighborhood` for its report lines only (`R/fs_family_report.R:134–135`, `:164`, `:336–337`).
- **Their model of the selected subgroup is the `maxeffCons` pick and nothing else:** `R/fs_oc_predict.R:280–293` "S4: maxeffCons selection — `Bmask <- Bhat; Bmask[!pass] <- -Inf; winner <- max.col(Bmask, ties.method = "first")`" (argmax effect among the passing candidates); `R/fs_oc_grid.R:570–589` "per draw: `w = argmax_{g eligible} Bhat_g` … `w` is the maxeffCons winner whenever anything declares" (`.fs_oc_reduce`: `Bmask <- dr$Bhat; Bmask[!eligible] <- -Inf; w <- max.col(Bmask, ties.method = "first")`). No band, no size, no `effMaxSG` (largest in the band) is implemented.
- **Finding F1:** the OC sections cannot run under `effMaxSG` as the package stands; they would model the `maxeffCons` selection while the anchor were chosen by `effMaxSG`. Making them carry the rule is an `R/` change (out of scope; no proposal drafted).

## 5. S0.4 — Cost

| document | last render | machine / workers | wall | CPU time | memory per worker |
|---|---|---|---|---|---|
| `analysis_actg175_continuous_oc.qmd` | 2026-09-16 (field-s task) | pop-os / 14 | 3,214 s (53.6 min); OC loop 2,089.7 s | not recorded | ≈ 6.5 GB (peak summed 90,526 MB) |
| same | 2026-09-08 (record `:39`, `:116`) | Mac-Studio-3 / 1 | 20 min 49 s; loop 1,123.4 s | not recorded | ≈ 11 GB (one worker) |
| `analysis_actg175_continuous_oc_evaluation.qmd` | none recorded | — | — | — | — |
| `analysis_actg175_continuous_compare_all.qmd` | none recorded (payload mtime 2026-09-01) | — | — | — | — |

## 6. Facts for Larry's decision (facts only)

- **The anchor under each rule** (same data, same settings, MR off): `maxeffCons` → `{age <= 37} & !{cd40 <= 507}`, n = 66, oriented MD 87.92, Pcons 0.95; `effMaxSG` ε = 0.20 (`neighborhood`) → `!{cd40 <= 507} & {gender}`, n = 79, MD 71.10, Pcons 0.91 (band floor 70.33; 4 of 8 qualifying candidates in band; the `maxeffCons` anchor is in the band as the third-largest). The two subgroups share 52 patients; neither contains the other.
- **Documents and sections a switch changes:** `analysis_actg175_continuous_oc.qmd` — §1 (anchor), §2.1 (intervals, the MR gate on the new Ĥ), §3 (the anchored truth `Q_primary`, hard-coded to `age = 37`, `cd40 > 507`, and the Q knob), §4–§12 (every OC section keyed on `T_obs`), the payload; its literal assertion at `:146–147` stops the render under the new anchor until edited. `analysis_actg175_continuous_oc_evaluation.qmd` — its anchor comes from a `dev/` payload and its prose hard-codes the committed clauses. `analysis_actg175_continuous_compare_all.qmd` — a sweep; its `effMaxSG` cells use `pareto` / `both`, not `neighborhood`.
- **Whether the OC sections can run under `effMaxSG` without an `R/` change: no.** `fs_oc_family_enumerate()`, `fs_oc_grid()`, `fs_oc_predict()` and `fs_oc_invert()` take thresholds and family arguments only, read `sg_focus` solely for the `maxeff` floor relaxation, and model the selection as the argmax effect among passers (`maxeffCons`); `effect_neighborhood` and `selection_rule` are not read.
- **Cost of each re-render:** the headline document 53.6 min at 14 workers on this machine (≈ 6.5 GB per worker, 90.5 GB peak), 20.8 min at 1 worker on the Mac; the evaluation document has no record and depends on a `dev/` payload; the comparison document has no record.

## 7. Findings

- **F1.** The OC functions implement only the `maxeffCons` selection (`fs_oc_predict.R:291–293`, `fs_oc_grid.R:575–589`); no argument carries `effMaxSG`. An `R/` change would be needed for OC under `effMaxSG`; none is proposed.
- **F2.** The headline document asserts the committed anchor literally (`:146–147`) and hard-codes its clauses in `Q_primary` (`:431`), the Q-knob grid (`:473`) and the rung / ladder literals (`:37–39`); a rule switch is a document edit at each of those sites, not a knob.
- **F3.** The `effMaxSG` band on this data holds four candidates within 20% of the maximum effect; the largest (n = 79) is a different clause pair (`gender` for `age ≤ 37`) with a smaller effect (71.10 vs 87.92) and lower consistency (0.91 vs 0.95).
- **F4.** With `details = FALSE` the fitted object does not carry the enumerated family; the family size (4,935) is known only from the MR gate's payload.
- **F5.** `compare_all` already sweeps `effMaxSG`, but at `selection_rule` `pareto` / `both`, not the survival grid's `neighborhood`.
- **F6.** The inventory table is generated by pattern from code lines; where a document sweeps its focus through variables, the cell says so rather than listing every value.

## 8. Commits

```
3c6745ae Add TASK_actg175_applied_effmaxsg_stage0_2026-09-16 as received
<this record: the next commit>
```
