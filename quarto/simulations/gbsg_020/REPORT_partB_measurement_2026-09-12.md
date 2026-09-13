# REPORT — Part B initial measurement: what MR costs, and what the six criteria resolve to

- **Date:** 2026-09-12. **Executor:** Claude Code, unattended. **Branch:** `feature/glm-extension`, starting at `1a782e4a`; the task document was committed first as `1d1e6257`.
- **Task:** `dev/tasks/TASK_partB_measurement_2026-09-12.md`.
- **Outcome: STOPPED at Stage 0c.** The template has no `mr_inference` knob; it hard-codes `TRUE`. As the task requires, **Stage M was not run** and nothing was reached by other means: no MR-off render, no template edit, no `R/` change.
- **Gate P is therefore not computed.** No MR-off per-replicate cost exists to project from.
- Two sections below still stand. The Monte Carlo resolution needs no measurement. The cost figures from committed MR-on bundles are labelled **bounds**, not the measurement.
- Readout script: `scripts_dinamr/partB_stage0_readout.R`, which reads committed bundles only. Output: `scripts_dinamr/logs/partB_stage0_readout.log` and `scripts_dinamr/partB_stage0_readout.rds`.

---

## Stage 0a — the resolution map at HEAD

**Confirmed, no departure.** `R/fs_focus_tag.R:65–74` (consistency) and `:76–86` (dina / grf):

| Spelling | consistency | dina / grf |
|---|---|---|
| `eff`, `hr`, `maxcons` | `maxcons` (`:67–69`) | `eff` (`:77–79`) |
| `maxeff` | `maxeff` (default arm `:74`) | `eff` (`:80`) |
| `maxeffCons` | `maxeffCons` (default arm `:74`) | `eff` (`:81`) |
| `effMaxSG`, `hrMaxSG` | `effMaxSG` (`:70–71`) | `effMaxSG` (`:82–83`) |
| `effMinSG`, `hrMinSG` | `effMinSG` (`:72–73`) | `effMinSG` (`:84–85`) |
| `maxSG`, `minSG` | pass through (`:74`) | pass through (`:86`) |

### On DINA and GRF, `maxeffCons` and `maxeff` are the same run

The tag does **not** drive the run. The template passes the raw spelling to `forestsearch()` (template `:1009`, `sg_focus = sg_focus`). `focus_tag` only builds the output stem (`:390`, `:420`). The equivalence is real because the **engines** collapse the two spellings, as follows.

- **Neither spelling is an alias.** Both are canonical (`.FS_SG_FOCUS_ALIASES`, `R/forestsearch_helpers.R:905–908`). Both reach the engines unchanged after `.normalize_sg_focus()` (`R/forestsearch_main.R:1478`).
- **DINA, native ordering:** `hr`, `maxeff` and `maxeffCons` are all `order(-eff, idx)` (`R/dina_subgroup.R:515–517`).
- **DINA, as the template runs it** (`dina_select_statistic = "effect"`, template `:503–508`): selection goes through `.dina_reselect_on_effect()` (`R/forestsearch_helpers.R:1471–1474`).
  - Its ordering switch has no `maxeff` or `maxeffCons` arm, so both take the fallback `ok[order(-eff[ok], ok)]` (`:1263`).
  - That is the same key as its `hr` arm.
- **GRF:** `frontier_rule <- switch(sg_focus, ..., maxeff = "eff", maxeffCons = "eff")` (`R/forestsearch_main.R:2412–2415`).
  - `.grf_reselect_on_effect()` then calls `.grf_frontier_select(rule = "eff")` (`R/forestsearch_helpers.R:1637`; `R/grf_subgroup_labels.R:364`).
- **Admission does not depend on focus.** `.fs_admission_applies()` (`R/forestsearch_helpers.R:2268`) returns "effect floor, no consistency term" for `dina` / `grf` whatever the focus.
- **The one difference that is not gated on engine:** the `maxeff` override block (`R/forestsearch_main.R:1537–1571`).
  - It sets `pconsistency.threshold` to 0, `stop_threshold` to NULL, `use_twostage` to FALSE, `max_subgroups_search` to Inf and `minp` to 0.
  - **None of these reaches DINA or GRF selection.** The DINA section returns at `:2355` and the GRF section at `:2540`, both before the consistency search that consumes `minp` (`:3008`). None of the five appears inside either section. `pconsistency.threshold` feeds only a consistency floor, and `.fs_admission_applies()` excludes that floor on these engines.
  - Its only visible effect is a `warning()` (`:1552`, not gated on `quiet`). It fires under `maxeff` and not under `maxeffCons`. It is output only.
- **Under MR (not part of Part B)**, both spellings re-select as `"maxeff"` (`R/fs_mr_inference_methods.R:97`, `:103`).
- **Consequence for Part B.** Take the six criteria to be `maxeff`, `maxeffCons`, `effMaxSG`, `effMinSG`, `maxSG` and `minSG`.
  - Consistency then has 6 distinct runs, and DINA and GRF have 5 each: 6 + 5 + 5 = **288 cell-runs** (108 FS, 90 DINA, 90 GRF).
  - The count changes if `eff` / `hr` / `maxcons` is one of the six. It is a separate rule on consistency, but it collapses to `eff` on DINA and GRF along with `maxeff` and `maxeffCons`.
  - On DINA and GRF the two spellings also produce the **same output stem**. Run the pair once; a second run under the same campaign tag would overwrite the first.
- **Separate template blocker.** The focus guard `stopifnot(sg_focus %in% c("maxeffCons", "effMaxSG", "maxSG", "minSG"))` (template `:327`) admits only **four of the six**: `maxeff` and `effMinSG` are rejected. Part B needs this guard widened too. The change is add-only.

## Stage 0b — ε scope

**Confirmed on all three engines:** `effect_neighborhood` affects selection only under `effMaxSG` / `effMinSG`. The kickoff's citations (`forestsearch_main.R:602–606`, `dina_subgroup.R:214–222`, `grf_main.R:47–52`) are roxygen text. The operative code is:

- **Consistency.** `sort_subgroups()` (`R/subgroup_consistency_helpers.R:540`).
  - The `maxeffCons` key is `setorder(-hr, K)` with no band term (`:570–572`).
  - `.compute_inclusion_band()` is reached only inside `if (sg_focus %in% c("hrMaxSG", "hrMinSG"))` (`:587`). The preview sort does the same (`:659`, `:676`).
  - `forestsearch()` range-checks ε only for the band foci (`R/forestsearch_main.R:1519`).
- **DINA.**
  - Band only in the `hrMaxSG` / `hrMinSG` arms (`R/dina_subgroup.R:518–535`).
  - Band only in the same two arms of `.dina_reselect_on_effect()` (`R/forestsearch_helpers.R:1250–1261`).
- **GRF.**
  - `.grf_frontier_select()` computes the band only in its `eff*SG` arm (`R/grf_subgroup_labels.R:371`). The `eff`, `maxSG` and `minSG` arms never read `nbhd`.
  - The frontier path runs only when `grf_selection = "frontier"` (`R/grf_main.R:283–296`; the template literal at `:503` sets frontier).
- **MR re-selection.** `.inband()` is called only by `effMaxSG` / `effMinSG` (`R/fs_mr_inference.R:168–175`).

**GRF's default is 0.10.** It is stated at `R/grf_main.R:51` and set by the formal at `:170`. Under `forestsearch()`, GRF actually receives `forestsearch()`'s own `effect_neighborhood` (`R/forestsearch_main.R:2433`). That default is also 0.10 (`:1254`), as is the template's (`FS_S7_NBHD`, template `:335`). **So 0.20 must be set explicitly for GRF's band rules.** `grfmr.sh` does set it, and the GRF comparator's meta records `effect_neighborhood = 0.2`.

**The bundle metas confirm it: ε was inert on the nine 12.4% FS comparator cells.**
- All 27 `tier2` / `p12ext` bundles (18 batch, 9 combined) record `subgroup_method = consistency`, `sg_focus = maxeffCons` and `effect_neighborhood = 0.1`. The batch metas also record `selection_rule = neighborhood`.
- Under `maxeffCons`, ε enters neither the identifier nor MR's re-selection.
- **Those nine cells therefore differ from the 31% cells, and from DINA and GRF, in the rule, not in ε.**
- The only recorded quantity that reads ε there is the informational `band_n` column (template `:1090–1095`).

**Records correction. Reported, not edited.** The kickoff's premise needs a qualifier: the records do not say only "different ε". They name the rule; they are wrong where they present ε 0.10 as part of the difference.
- **Correct as they stand:** the comparator tables name the rule: `maxeffCons` ε 0.10 in `current_status.md:23–24`, `REPORT_fs_extraction_2026-09-11.md:18` and `REPORT_grfmr_completion_2026-09-12.md:169–170`.
- **Wrong where they make ε 0.10 operative:**
  - `REPORT_dinamr_blockC_2026-09-11.md:596` says "a different selection functional **at half the band width**". The half-width is not a property of that functional.
  - The captions of `summary_grfmr.qmd` gate `fs_matched` / `criterion_matched` on "sg_focus **and eps**".
  - Both `REPORT_grfmr_completion_2026-09-12.md:257` and `REPORT_grfmr_2026-09-11.md:418` quote "(`maxeffCons` ε 0.10 against `effMaxSG` ε 0.20)" as the criterion difference.
  - `fs_extraction.R`'s header carries eps as a criterion column.
- **Correct statement:** at 12.4% the FS comparator runs `maxeffCons`, the consistency-screened effect argmax, which has no band. ε is not a dimension of the confound. The recorded 0.10 is the inert default.

## Stage 0c — `mr_inference` is not reachable from the template

**No knob.** `mr_inference = TRUE` is a literal inside `base_args` (template `:1017`).
- None of the `FS_S7_*` knobs reads it: `NSIMS`, `FB`, `FB_PATH`, `JOIN_SKIP`, `NB_BOOTS`, `KNOISE`, `MODE`, `START`, `METHOD`, `FOCUS`, `NBHD`, `ER_JCUTS`, `HR`, `N`, `Z1Q`, `WINNER_ROWS`, `CAMPAIGN`, `QUICKRUN`, `SAVE_COMBINED`, `UNIFORM`, `FIELD_COMPLEMENT`, `IJ_RESIDUAL`, `FIELD_DECOMP`, `FIELD_SCALEC`, `FIELD_RECOV`, `RETURN_RESEL` and `WORKERS` (template `:249–674`).
- **Per the task: STOP.** The change is described below and not made.

**What an add-only, default-inert template change would touch:**
1. **Knob.** Add, say, `FS_S7_MR` (default `TRUE`), read beside the other MR controls (`:585–645`) and substituted at `:1017`.
2. **Stem.** Add a token only when MR is off (e.g. `_nomr`) in `rds_stem` (`:417–421`). MR-off bundles can then never share a stem or a `combine_glob` with MR-on ones. Default stems stay unchanged.
3. **Meta.** Add `mr_inference` to both saved metas (`:1516–1517`, `:1637–1638`) and to the combine poolability key (`:1588–1590`).
4. **Recorder.**
   - `detected`, `sg_def`, `covs`, and `sens` / `spec` / `ppv` / `npv` / `n_harm` (`:1246–1249`) are already recorded outside the MR block, so identification survives MR off.
   - **`n_sel`, `label`, `n_family`, `n_cons_qual`, `band_n` and `p_hat_*` are filled only inside `if (!is.null(g))` (`:1068–1243`), so they would be NA with MR off.** The identification check would compare `n_harm`, which is |Ĥ| from `sg.harm.id`. The alternative is an add-only new column; re-sourcing `n_sel` would not be add-only.
   - `fit_mr_secs` (`:1043`, commented "INCLUDING MR") would time the fit alone. The name can stay if the meta records the switch.
5. **Focus guard.** Widen `:327` to the six Part B criteria (Stage 0a).
6. **Combine mode.** The summary chunks read MR columns. Whether they tolerate all-NA MR columns is **untested** and needs a smoke render.

**A comparator problem Stage M would hit even with the knob.**
- The FS bundles at the Stage M coordinate are the designated `e1stud` and the same-criterion `p30sgnb20`. Both were built on **pop-os at 100 workers under R 4.6.1**, not on the Mac at 12 workers. The DINA and GRF comparators are Mac, 12 workers, R 4.5.2.
- On the two FS cells where both hosts ran the same seeds, pop-os per-replicate seconds are **3.9–4.4× the Mac's** (§3 of the readout, below).
- So an FS MR-off run on the Mac, set against `e1stud`, would mix the MR share with the host effect.
- **Recommendation:** when Stage M runs, pair the FS arm on the Mac, running MR-on and MR-off at 200 replicates each, rather than reading MR-on from `e1stud`. DINA and GRF can use their committed comparators as the task specifies.

---

## What the committed bundles do say: MR on only, as bounds

Definitions (the readout script's header):
- **Ceiling** = `fit_mr_secs`. MR off cannot cost more than MR on.
- **Bound** = `fit_mr_secs − fld_H_secs − fld_Hc_secs`, with NA counted as 0.
  - Both field timers are nested inside `fit_mr_secs` and disjoint (`walls.R:14–16`), so this removes only the field block. The de-biasing draws, the IJ and the re-selection record stay in.
  - **MR-off cost is below the bound, and 1 − bound/ceiling is a lower bound on the MR share.**

### 1. The three comparators at 31%, HR 1.50, n 500 (effMaxSG ε 0.20, 2,000 replicates, sim_id 1–2000)

| Engine | Host / workers | Ceiling q25 / median / q75 / p90 (s) | Bound q25 / median / q75 / p90 (s) | MR share, lower bound | KB per replicate (MR on) |
|---|---|---|---|---|---|
| consistency (`e1stud`) | pop-os / 100 | 64.74 / 72.88 / 81.22 / 87.73 | 29.06 / 33.88 / 38.78 / 42.68 | **≥ 0.535** | 0.753 |
| consistency (`p30sgnb20`) | pop-os / 100 | 64.06 / 72.46 / 80.82 / 87.86 | 28.44 / 33.54 / 38.44 / 42.42 | ≥ 0.537 | 0.630 |
| dina (`dinamr`) | Mac / 12 | 16.84 / 25.42 / 36.26 / 42.89 | 3.52 / 5.87 / 9.06 / 10.96 | **≥ 0.769** | 0.802 |
| grf (`grfmr`) | Mac / 12 | 15.23 / 15.93 / 16.55 / 17.04 | 3.75 / 3.83 / 3.93 / 4.11 | **≥ 0.760** | 0.796 |

- **Not available:** the MR-off quantiles, the measured MR share, the identification-unchanged check and the MR-off bundle size. All four need Stage M.
- **Detection:** 0.9995 / 0.9985 / 1.0000 (consistency / DINA / GRF).

### 2. FS cost across `sg_focus`: same cell (31%, HR 1.50, n 500), same host (pop-os / 100)

| Focus (campaign) | Ceiling median (s) | Bound median (s) | Field block median (s) | MR share, lower bound |
|---|---|---|---|---|
| `maxeffCons` (`p30`) | 53.98 | 33.33 | 21.12 | ≥ 0.383 |
| `effMaxSG` 0.10 (`p30sg`) | 71.66 | 33.28 | 38.47 | ≥ 0.536 |
| `effMaxSG` 0.20 (`p30sgnb20`) | 72.46 | 33.54 | 39.72 | ≥ 0.537 |
| `effMaxSG` 0.20 (`e1stud`) | 72.88 | 33.88 | 39.62 | ≥ 0.535 |
| `effMaxSG` 0.30 (`banddial`) | 73.31 | 33.68 | 40.49 | ≥ 0.541 |
| `maxSG` (`banddial`) | 54.79 | 33.27 | 22.06 | ≥ 0.393 |
| `minSG` (`banddial`) | 52.04 | 32.58 | 19.58 | ≥ 0.374 |

- **The bound, which contains the identification, is flat across sg_focus**: 32.6–33.9 s over seven settings.
- **The MR-on total is not flat.** The field block is roughly twice as expensive under the band rule as under `maxeffCons` / `maxSG` / `minSG`.
- **So the MR share depends on the focus, while the identification cost does not.** "Flat across sg_focus" is supported for Part B's MR-off cost and would be wrong for MR-on cost.
- **No committed bundle covers `maxeff` or `effMinSG`.**
  - `maxeff` is the least safe for the flat assumption: it disables truncation and early stopping and evaluates consistency for every candidate (`R/forestsearch_main.R:1537–1571`).

### 3. Host factor: FS `maxeffCons` at 12.4%, n 500, sim_id 1–1000, same seeds

| Cell | pop-os / 100 (`s7c`), ceiling / bound (s) | Mac / 12 (`tier2`), ceiling / bound (s) | Ratio, ceiling / bound |
|---|---|---|---|
| HR 1.00 | 35.46 / 19.77 | 9.18 / 3.50 | 3.86 / 5.64 |
| HR 1.75 | 45.83 / 26.72 | 10.38 / 4.04 | 4.41 / 6.62 |

- 14 of the 18 designated FS comparator cells are pop-os bundles.
- On the pop-os cells, FS per-replicate cost rises steeply with n (e.g. 45.6 → 96.0 → 164.9 s at 12.4% HR 1.50).
- The four Mac cells are nearly flat in n (10.3 → 11.7 → 12.1 s at HR 1.75).
- **FS seconds from the committed grid are therefore not Mac seconds.**

### 4. Bounds on the Part B sweep: NOT the Gate P projection

- **Method.** One criterion's 18-cell cost is taken from each engine's committed grid and multiplied by the number of distinct criteria: 6 for FS, 5 each for DINA and GRF. Wall time is worker-hours ÷ 12.
- **DINA and GRF** are all Mac / 12.
- **FS is 14 × pop-os / 100 plus 4 × Mac / 12.** Its hours inherit §3's inflation. The FS grid also mixes rules: `maxeffCons` at 12.4%, `effMaxSG` at 31%.

| Replicates per cell | FS ceiling / bound (wall h, 12 w) | DINA ceiling / bound | GRF ceiling / bound | All three: ceiling / bound |
|---|---|---|---|---|
| 2,000 | 435.5 / 256.6 | 67.9 / 16.5 | 67.3 / 20.5 | 570.7 / 293.6 |
| 1,000 | 217.7 / 128.3 | 33.9 / 8.2 | 33.7 / 10.2 | 285.3 / 146.7 |
| 500 | 108.9 / 64.1 | 17.0 / 4.1 | 16.8 / 5.1 | 142.7 / 73.3 |

- **Read the DINA and GRF columns as upper bounds on Mac wall time.** At 2,000 replicates their MR-off wall is **below** 16.5 h and 20.5 h respectively.
- **The FS column is an upper bound in pop-os seconds only.**
  - Applying §3's ratio of 3.9–6.6 would put FS on the Mac well below 100 h at 2,000 replicates.
  - That ratio comes from two cells at n 500. It is **not** a projection.
  - The FS sweep is the figure Stage M most needs to pin.

## Gate P — Monte Carlo resolution (needs no measurement)

Wilson 95% half-widths:

| Replicates per cell | Rate near 0.50 | Rate near 0.90 | Detected replicates at the thinnest cell* | Near 0.50 over those | Near 0.90 over those |
|---|---|---|---|---|---|
| 2,000 | ±0.0219 | ±0.0132 | 688 | ±0.0373 | ±0.0225 |
| 1,000 | ±0.0309 | ±0.0186 | 344 | ±0.0525 | ±0.0318 |
| 500 | ±0.0437 | ±0.0264 | 172 | ±0.0739 | ±0.0452 |

\*DINA, 12.4%, HR 1.00, n 1500, detection 0.344.

**The projected totals at 2,000 / 1,000 / 500 are not given.** The task defines them from the measured MR-off cost, and that measurement is blocked by Stage 0c.

## What is needed before Part B can be measured

1. **Approve the template change** listed in Stage 0c: the knob, the stem token, meta plus the poolability key, the focus guard widened to six criteria, and a combine smoke with MR off. It belongs in its own task.
2. **Decide the FS arm of Stage M:** a paired MR-on / MR-off run on the Mac (recommended), or MR-off on pop-os against `e1stud`.
3. **Confirm the six criteria.** 288 holds for {`maxeff`, `maxeffCons`, `effMaxSG`, `effMinSG`, `maxSG`, `minSG`}.

## Files

| Path | Status | Size |
|---|---|---|
| `dev/tasks/TASK_partB_measurement_2026-09-12.md` | committed (`1d1e6257`) | 6.4 KB |
| `quarto/simulations/gbsg_020/scripts_dinamr/partB_stage0_readout.R` | committed with this report | 10.2 KB |
| `quarto/simulations/gbsg_020/scripts_dinamr/logs/partB_stage0_readout.log` | committed with this report | 9.7 KB |
| `quarto/simulations/gbsg_020/scripts_dinamr/partB_stage0_readout.rds` | committed with this report | 2.8 KB |
| `quarto/simulations/gbsg_020/scripts_dinamr/README.md` | modified (one table row) | — |
| `quarto/simulations/gbsg_020/REPORT_partB_measurement_2026-09-12.md` | this file | — |

- **No measurement bundles and no renders.** Stage M did not run.
- Nothing is over 50 MB.
