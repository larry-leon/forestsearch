# REPORT — ACTG175 continuous (CD4 change): post-selection intervals under effMaxSG, ε = 0.20

Date: 2026-09-16. Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1, reference BLAS). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_actg175_continuous_intervals_2026-09-16.md` (committed as received, `e1c3c7a6`), on Larry's decision that the applied analysis runs `effMaxSG` at ε = 0.20. Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix`; no `R/` change. Source document `<oc>` = `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` at HEAD `700bdcdd`, read only. New document `<doc>` = `quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd`. The render ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` under `timeout 30m`, quarto 1.9.38.

## 1. Provenance — GATE PASS

```
pop-os
feature/glm-extension
700bdcdd
700bdcdd actg175/continuous closeout: generate current_status.md at f74271e4
f74271e4 actg175/continuous status_curated.md: open-work line -- DINA and GRF on this design: Stage 0 recorded, decisions pending
ee47ee61 DINA and GRF on the ACTG175 continuous (MD) design, Stage 0 (read-only): ...
[tracked modifications: none]
DINA/GRF Stage 0 in HEAD
no R/ change since the install
new path free
R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix
```
- Processes: the task's `ps | grep` pattern matched one 11-day-old `bash -c` wrapper (PID 2086198) from a `gbsg_020` `uburst` render, whose only child is `sleep` (its `free -m` memory-monitor loop was never killed; the renders it wrapped are long finished); by process name (`ps -eo comm`) no `R`, `Rscript`, `quarto` or `deno` process existed; load average 0.17. Gate read as passed; the orphan was left alone (finding F5).
- First commit: `e1c3c7a6 Add TASK_actg175_continuous_intervals_2026-09-16 as received`.

## 2. Read (at `700bdcdd`)

**2.1 Transplant sources from `<oc>`.**
- YAML header and `params`, `:1–9`:
  ```
  title: "ACTG175 Continuous Analysis: Applied Operating Characteristics, Self-Contained"
  format: html
  params:
    draws: 20000
    n_workers: 14
    c2_ratio: 0.8
    seed: 8316951
  ```
- Setup chunk `:11–49`: libraries (`forestsearch`, `speff2trial`, `gt`); `results_dir <- NULL` (`:24`), `dirout <- "analysis_actg175_continuous_oc_intervals"` (`:28`); the header note `:30–36` ("Everything below is computed in this document: the analysis, the anchored truths ..., the candidate family, and the operating characteristics. Nothing is read from disk. Re-render at other precisions via the YAML params (draws, n_workers, c2_ratio, seed) ..."); the OC rung/ladder lines `:37–40` (`q_rungs`, `q_shared`, `c1_ladder`, `c2_vec <- sort(unique(c(5, params$c2_ratio * c1_ladder)))`); helpers `fmt`, `fmt0`, `pct`, `pm`, `at`, `qlab` `:42–47`; `t_doc` `:48`.
- Data chunk `data-prep` `:57–78`: arms 1 and 3 of `speff2trial::ACTG175`, `y_decline <- cd40 - cd420`, six continuous and six binary confounders as factors, `stopifnot(N == 1083L)`, the unadjusted ITT `lm`.
- Anchor call `:92–135`, every argument (the field-s record quotes it in full at `:92–132` before the field-s edit; at `700bdcdd` it is the same call plus `field_scale_complement = "selected"` at `:133`):
  `df.analysis = actg_df, confounders.name = confounders.name, outcome.name = "y_decline", treat.name = "treat", id.name = "id", outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE, seedit = params$seed, sg_focus = "maxeffCons", selection_rule = "neighborhood", consistency_method = "resample", effect.threshold = 10, consistency.threshold = 5, pconsistency.threshold = 0.90, use_twostage = TRUE, conf.cont_jcuts = list(age = 10, preanti = 10, wtkg = 10, karnof = 10, cd40 = 10, cd80 = 10), cut_type = "default", maxk = 2L, n.min = 60L, d0.min = 10L, d1.min = 10L, fs.splits = 500L, use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE, is.RCT = TRUE, parallel_args = list(plan = "sequential"), details = FALSE, quiet = TRUE, mr_inference = TRUE, mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE, field_complement = TRUE, field_scale_complement = "selected", return_reselection = TRUE)`.
  Derived quantities `:137–143` (`H_def`, `n_H`, `win`, `T_obs`, `p_cons`, the `cat`); the literal anchor assertion `:146–147` `stopifnot(setequal(fs_anchor$sg.harm, c("{age <= 37}", "!{cd40 <= 507}")), n_H == 66L, abs(T_obs - (87 + 11 / 12)) < 5e-7)`; OC bookkeeping `:148–151` (rungs and ladder set to `T_obs`).
- Intervals chunks: prose `:160–176`; `intervals-objects` `:178–251` — guards `:180–184` and `:197–200`, `f/fc/jt/js/gc/rs` `:185`, `lab_of` `:190–194`, `thr_H <- c(0, 10, 20, 30, 40); thr_Hc <- c(0, 10, 20, 30)` `:201`, `p_hat_H` `:202`, `lo1` `:205`, `up1` with `field_s = fc$upper_1s_s` `:206–207`, the `iv` element `:210–250`; `intervals-table-H` (`tab1`) `:253–276`; `intervals-spread` `:278–285`; `intervals-table-Hc` (`tab2`) `:287–312`; `intervals-table-joint` (`tab3`) `:314–331`; `intervals-diagnostics` (`tab4`, prints `p̂(Ĥ)` in its first row) `:333–355`.
- Payload save lines: `export-payload` `:1254–1332`; `.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))` `:1265`; `intervals = iv` `:1303`; `saveRDS(payload, .payload_file)` `:1330`. Path: `quarto/applications/actg175/_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_payload.rds` — **tracked** (`git ls-files`), as is `analysis_actg175_continuous_oc.html`.
- Intervals interpretation: the callout "Reading the intervals (generated from the numbers above at render time)" `:357–408`, chunk `intervals-reading` `:360–405`, inline `cat(sprintf(...))` over `f`, `fc`, `jt`, `js`, `g`, `gc`: complement first (field-s `U`, unstudentized `U_u` named once), Ĥ, the IJ reference, the field-s Bonferroni pair, the price of selection, the regime (`p̂(Ĥ)`; the tie-regime branch `:397` cites "the simulation campaign mdf1"), the provenance line. `p̂(Ĥ)` is printed at `:335`/`:344` (tab4) and `:394–398`, so no extra line was needed.

**2.2 Rule.** `quarto/simulations/actg175/continuous/mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds`, `meta`: `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.2`, **`selection_rule = "neighborhood"`**.

**2.3 Reproduction target** (`REPORT_actg175_applied_effmaxsg_stage0_2026-09-16.md` §3, computed at `570bf8a5`): "(b) `effMaxSG`, ε = 0.20, `neighborhood`: selected subgroup `!{cd40 <= 507} & {gender}`; n(Ĥ) 79; oriented MD estimate (the `hr` column) 71.095513; consistency proportion (`Pcons`) 0.91; complement n 1004; rows of `grp.consistency$out_sg$result` 8." Band: "max oriented MD among the consistency-qualifying candidates 87.916667 (the committed anchor itself); floor (1 − 0.20) × max = **70.333333**; **4 of 8** candidates in the band": 1. `!{cd40 <= 507} & {gender}` N 79, MD 71.0955, Pcons 0.91 (selected); 2. `{age <= 40} & !{cd40 <= 507}` N 76, MD 78.7054, Pcons 0.93; 3. `{age <= 37} & !{cd40 <= 507}` N 66, MD 87.9167, Pcons 0.95; 4. `!{wtkg <= 73} & !{cd40 <= 507}` N 61, MD 82.2688, Pcons 0.93. Overlap with the `maxeffCons` anchor: 52 of 66; neither contains the other.

**2.4 Contrast citation** (`REPORT_actg175_continuous_field_s_2026-09-16.md` §4–§5, commit `570bf8a5`): Ĥ = `{age <= 37} & !{cd40 <= 507}`, `n_selected` 66 (§4). §5 table and comparators, oriented scale: "Ĥ field one-sided lower bound −43.7954"; "Hc one-sided 95% UPPER bound: unstudentized −18.5546, field-s −18.9041"; "Bonferroni pair: H lower −56.4409 / Hc upper −16.4175; Bonferroni joint prob. on aligned draws 0.9491; calibrated gamma 0.0250"; "Ĥᶜ naive upper −23.6581; Ĥᶜ IJ upper −8.7454". The record does not state the unadjusted and IJ one-sided lower bounds on Ĥ; from the committed payload at `570bf8a5` (`extras$intervals$H`): naive `lower_1s` 17.7478, IJ `lower_1s` −55.0512, field `lower_1s` −43.7954 (finding F2).

## 3. The transplant (`<doc>` at `da25abcc`)

The named changes against `<oc>`'s chunks (`:1–408` and the export chunk `:1244–1332`):
1. Header: the title; `params` keeps `seed` only (`draws`, `n_workers` dropped as OC-only; `c2_ratio` dropped too — it is read only by the OC rung/ladder lines, finding F3); the header note reworded.
2. Setup: `dirout` kept, `payload_suffix <- "_effmaxsg"` added; the OC rung/ladder/`c2_vec` lines `:37–40` dropped (they read `params$c2_ratio`); helpers and `t_doc` unchanged. Data chunk unchanged.
3. Anchor: `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.20`, `selection_rule = "neighborhood"`; MR arguments exactly as `<oc>` passes them; the `cat` prints `p.consistency` at `%.8f`; the assertion replaced by `stopifnot(setequal(fs_anchor$sg.harm, c("!{cd40 <= 507}", "{gender}")), n_H == 79L, abs(T_obs - 71.095513) < 1e-6, abs(p_cons - 0.91) < 1e-8)` (finding F1); the OC bookkeeping `:148–151` dropped.
4. New chunk `band`: the ε band from `fs_anchor$grp.consistency$out_sg$result` (columns `Pcons`, `hr`, `N`, `K`, `M.1`, `M.2`; `R/subgroup_consistency_helpers.R:1860–1862`), floor `(1 − ε) × max(hr)`, a `gt` table, and a guard that row 1 is in band and is `sg.harm`.
5. Intervals chunks unchanged in code (`lo1`, `up1` with `field_s`, `iv`, `tab1`–`tab4`, `intervals-spread`, `intervals-reading`) except: the three table titles and tab4's title gain "(effMaxSG, ε = 0.20)"; the reading's tie-regime string drops its `mdf1` clause (that campaign is `maxeffCons`'s, finding F4) and the provenance line names the rule; the prose `:160–176` drops its task-reference sentence.
6. New section "3. Interpretation" (inline R over the gate objects; the `maxeffCons` contrast quoted from §2.4) and the verbatim scope sentence.
7. Export chunk: `.payload_file` gains `payload_suffix`; the payload keeps `labels`, `meta`, `extras$anchor`, `extras$intervals = iv`, adds `extras$band`; the OC elements (`table = lad`, `type1`, `declared`, `calibration`, `q_variants`, ...) are absent. Nothing from `<oc>`'s OC sections (§3–§12) is transplanted.

The diff, `<oc>` `:1–408` + `:1244–1332` against `<doc>`:

```diff
--- analysis_actg175_continuous_oc.qmd (lines 1-408, 1244-1332 at 700bdcdd)
+++ analysis_actg175_continuous_intervals.qmd (da25abcc)
@@ -1,10 +1,7 @@
 ---
-title: "ACTG175 Continuous Analysis: Applied Operating Characteristics, Self-Contained"
+title: "ACTG175 continuous (CD4 change): post-selection intervals under effMaxSG, ε = 0.20"
 format: html
 params:
-  draws: 20000
-  n_workers: 14
-  c2_ratio: 0.8
   seed: 8316951
 ---
 
@@ -22,22 +19,19 @@
 #   dirout      = NULL -> the .qmd stem
 # ---------------------------------------------------------------------------
 results_dir <- NULL
-# TASK_actg175_continuous_intervals_2026-09-07: committed payloads are read-only, so
-# this render writes beside the committed OC payload rather than over it (the
-# tracked HTML is refreshed per the directory's convention).
-dirout      <- "analysis_actg175_continuous_oc_intervals"
-
-# Everything below is computed in this document: the analysis, the anchored
-# truths (the primary Q and its two supersets — the Q_variants knob, defined
-# visibly in section 3 once the population frame exists), the candidate
-# family, and the operating characteristics.  Nothing is read from disk.
-# Re-render at other precisions via the YAML params (draws, n_workers,
-# c2_ratio, seed); the archived companion
-# analysis_actg175_continuous_oc_evaluation.qmd holds the deep fixed-c2 run.
-q_rungs   <- c(0.01, 5, 10, 15, 20, 40, 60, 87 + 11 / 12)   # primary rungs; last = T_obs, gated below
-q_shared  <- c(0.01, 5, 10, 20, 40, 87 + 11 / 12)           # rungs run for each Q superset, incl. its own sub-threshold null
-c1_ladder <- sort(c(5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 87 + 11 / 12))
-c2_vec    <- sort(unique(c(5, params$c2_ratio * c1_ladder)))
+# TASK_actg175_continuous_intervals_2026-09-16: this document writes its payload
+# beside the maxeffCons intervals payload of analysis_actg175_continuous_oc.qmd
+# (same directory), with the file name suffixed _effmaxsg.
+dirout         <- "analysis_actg175_continuous_oc_intervals"
+payload_suffix <- "_effmaxsg"
+
+# Everything below is computed in this document: the data, one forestsearch()
+# call under effMaxSG (epsilon = 0.20, selection_rule = "neighborhood") with the
+# multiplier-resampling gate, the epsilon band the selection was made in, and the
+# interval constructions for H-hat and its complement.  Nothing is read from
+# disk.  The operating-characteristics evaluation of this design lives in
+# analysis_actg175_continuous_oc.qmd and remains under maxeffCons; nothing from
+# it is transplanted here.  Re-render at another seed via the YAML params.
 
 fmt  <- function(x, d = 3) formatC(x, format = "f", digits = d)
 fmt0 <- function(x) formatC(x, format = "d")
@@ -79,14 +73,16 @@
 
 ## 2. The analysis
 
-The search whose operating characteristics this document evaluates is the
-fixed-family `maxeffCons` configuration — the compare-all document's
-MR-anchor arm.  Multiplier-resampling inference and the lasso/GRF front ends
-are mutually exclusive under `subgroup_method = "consistency"` (the family
-would be estimated from the same data the search runs on), so the
-inferential anchor runs with the front ends off, and the operating
-characteristics evaluated here are the operating characteristics of exactly
-that procedure.
+The search is the fixed-family configuration of `analysis_actg175_continuous_oc.qmd`
+with its selection rule changed to **`effMaxSG` at ε = 0.20 under
+`selection_rule = "neighborhood"`**: among the consistency-qualifying candidates
+whose oriented mean difference is within 20% of the largest, the largest
+subgroup is selected.  This is the rule of the survival grid and of the
+`mdsgnb20` simulation campaign on this design.  Multiplier-resampling inference
+and the lasso/GRF front ends are mutually exclusive under
+`subgroup_method = "consistency"` (the family would be estimated from the same
+data the search runs on), so the inferential anchor runs with the front ends
+off.
 
 ```{r anchor}
 fs_anchor <- forestsearch(
@@ -99,7 +95,10 @@
   effect_measure   = "MD",
   adverse_outcome  = TRUE,
   seedit           = params$seed,
-  sg_focus         = "maxeffCons",
+  # TASK_actg175_continuous_intervals_2026-09-16: the effMaxSG rule at epsilon = 0.20,
+  # as the mdsgnb20 campaign bundles' meta record it (selection_rule "neighborhood").
+  sg_focus         = "effMaxSG",
+  effect_neighborhood = 0.20,
   selection_rule   = "neighborhood",
   consistency_method = "resample",
   effect.threshold       = 10,
@@ -139,29 +138,56 @@
 win    <- fs_anchor$grp.consistency$out_sg$result[1L, ]
 T_obs  <- as.numeric(win$hr)      # MD scale: the hr column holds the MD
 p_cons <- as.numeric(win$Pcons)
-cat(sprintf("H-hat = %s;  n(H-hat) = %d;  T_obs = %.6f;  p.consistency = %.2f\n",
+cat(sprintf("H-hat = %s;  n(H-hat) = %d;  T_obs = %.6f;  p.consistency = %.8f\n",
             H_def, n_H, T_obs, p_cons))
 
-# gates: the anchor must reproduce the compare-all / stage-0 result exactly
-stopifnot(setequal(fs_anchor$sg.harm, c("{age <= 37}", "!{cd40 <= 507}")),
-          n_H == 66L, abs(T_obs - (87 + 11 / 12)) < 5e-7)
-q_rungs[length(q_rungs)] <- T_obs
-q_shared[length(q_shared)] <- T_obs
-c1_ladder[abs(c1_ladder - (87 + 11 / 12)) < 1e-9] <- T_obs
-c2_vec <- sort(unique(c(5, params$c2_ratio * c1_ladder)))
+# gate: the anchor must reproduce the effMaxSG Stage 0 result
+# (REPORT_actg175_applied_effmaxsg_stage0_2026-09-16.md section 3, computed at 570bf8a5:
+#  !{cd40 <= 507} & {gender}, n = 79, oriented MD 71.095513, Pcons 0.91)
+stopifnot(setequal(fs_anchor$sg.harm, c("!{cd40 <= 507}", "{gender}")),
+          n_H == 79L, abs(T_obs - 71.095513) < 1e-6, abs(p_cons - 0.91) < 1e-8)
 ```
 
 The anchor: **Ĥ = `` `r H_def` ``**, n(Ĥ) = `r fmt0(n_H)` of `r fmt0(N)`,
-fitted mean difference **T̂~obs~ = `r fmt(T_obs, 6)`** CD4 units,
-consistency `r fmt(p_cons, 2)`.
+fitted mean difference **T̂~obs~ = `r fmt(T_obs, 6)`** CD4 units on the
+oriented `y_decline` scale (a CD4 change of `r fmt(-T_obs, 2)` on the raw
+scale), consistency `r fmt(p_cons, 2)`.
+
+### 2.0 The ε band the selection was made in
+
+`effMaxSG` picks the largest subgroup among the consistency-qualifying
+candidates whose oriented MD is at least (1 − ε) times the largest.  The
+qualifying set and the band are computed here from the fitted object's
+consistency table (rows are the candidates that passed the consistency
+threshold; the selected candidate is row 1).
+
+```{r band}
+eps  <- 0.20
+res  <- as.data.frame(fs_anchor$grp.consistency$out_sg$result)
+res$candidate <- ifelse(as.integer(res$K) == 2L, paste(res$M.1, res$M.2, sep = " & "), res$M.1)
+md_max   <- max(res$hr)
+md_floor <- (1 - eps) * md_max
+res$in_band  <- res$hr >= md_floor
+res$selected <- seq_len(nrow(res)) == 1L
+stopifnot(res$selected[1L], res$in_band[1L],
+          setequal(strsplit(res$candidate[1L], " & ", fixed = TRUE)[[1]], fs_anchor$sg.harm))
+n_qual <- nrow(res); n_band <- sum(res$in_band)
+bt <- res[order(-res$in_band, -res$N), c("candidate", "N", "hr", "Pcons", "in_band", "selected")]
+gt(data.frame(Candidate = bt$candidate, n = fmt0(bt$N), `MD (oriented)` = fmt(bt$hr, 4),
+              `CD4 change` = fmt(-bt$hr, 2), Pcons = fmt(bt$Pcons, 2),
+              `In band` = ifelse(bt$in_band, "yes", "no"), Selected = ifelse(bt$selected, "yes", ""),
+              check.names = FALSE)) |>
+  tab_header(title = sprintf("Consistency-qualifying candidates and the ε = %.2f band", eps),
+             subtitle = sprintf("max oriented MD %.6f; band floor (1 − ε) × max = %.6f; %d of %d qualifying candidates in band; effMaxSG selects the largest in band",
+                                md_max, md_floor, n_band, n_qual))
+```
 
 ### 2.1 Post-selection intervals for the found subgroup and its complement {#sec-intervals}
 
-The gate call above now runs multiplier resampling with `ci_method = "field"`,
-`field_complement = TRUE` and `return_reselection = TRUE`
-(`TASK_actg175_continuous_intervals_2026-09-07`; the mirror of the GBSG
-frozen-family intervals section), so the fitted object carries, for Ĥ and its
-complement Ĥᶜ, the interval constructions currently in use: naive, MR (IJ
+The call above runs multiplier resampling with `ci_method = "field"`,
+`field_complement = TRUE` and `return_reselection = TRUE` (the mirror of the
+GBSG frozen-family intervals section), so the fitted object carries, for Ĥ and
+its complement Ĥᶜ, the interval constructions currently in use: naive, MR (IJ
 two-term), MR (field) with its **one-sided bounds as the directional
 products**, and the Bonferroni / calibrated pair for a claim on both. The
 target throughout is the effect in the region the search found, with the
@@ -266,7 +292,7 @@
   `One-sided 95% bound (cd4_change): change ≤` = sprintf("%.2f", -lo1),
   check.names = FALSE)
 gt(tab1) |>
-  tab_header(title = "Table 1 — Ĥ: interval constructions for the found subgroup",
+  tab_header(title = "Table 1 — Ĥ: interval constructions for the found subgroup (effMaxSG, ε = 0.20)",
              subtitle = sprintf("Ĥ = %s, n = %d of %d; oriented = y_decline (positive = harm); cd4_change = sign flip",
                                 H_def, n_H, N)) |>
   tab_footnote(footnote = paste(
@@ -301,7 +327,7 @@
   `One-sided 95% bound (cd4_change): change ≥` = sprintf("%.2f", -up1),
   check.names = FALSE)
 gt(tab2) |>
-  tab_header(title = "Table 2 — Ĥᶜ: interval constructions for the complement (benefit claim)",
+  tab_header(title = "Table 2 — Ĥᶜ: interval constructions for the complement (benefit claim) (effMaxSG, ε = 0.20)",
              subtitle = sprintf("Ĥᶜ: n = %d; the one-sided UPPER bound is the exposed limit for a benefit claim (\"harm at most U\" = \"benefit at least −U\")",
                                 N - n_H)) |>
   tab_footnote(footnote = sprintf(paste(
@@ -325,7 +351,7 @@
                                          sprintf("%.3f", js$joint_prob)),
   check.names = FALSE)
 gt(tab3) |>
-  tab_header(title = "Table 3 — the joint pair (Ĥ lower bound, Ĥᶜ upper bound)",
+  tab_header(title = "Table 3 — the joint pair (Ĥ lower bound, Ĥᶜ upper bound) (effMaxSG, ε = 0.20)",
              subtitle = sprintf("unstudentized: γ = %.3f on the 0.025–0.050 grid, corr(Λ*, Λ*ᶜ) = %+.3f over %d aligned draws; field-s: γ = %.3f, corr = %+.3f over %d",
                                 jt$gamma, jt$corr, jt$n_joint_draws, js$gamma, js$corr, js$n_joint_draws))
 ```
@@ -350,7 +376,7 @@
             sprintf("%d; %.2f / %.2f = %.3f", fc$n_complement_fits, fc$se_field, gc$debiased$se_wald,
                     fc$se_field / gc$debiased$se_wald)),
   check.names = FALSE)
-gt(tab4) |> tab_header(title = "Re-selection and field diagnostics",
+gt(tab4) |> tab_header(title = "Re-selection and field diagnostics (effMaxSG, ε = 0.20)",
                        subtitle = "p̂ near 1 = settled selection; below ~0.5 = tie regime (the second-order field correction matters most)")
 ```
 
@@ -394,10 +420,10 @@
   "- **Regime:** p̂(Ĥ) = %.3f — the selected subgroup wins the re-selection on %.1f%% of multiplier draws, %s. ",
   "λ-SDᶜ / naive SE = %.3f on the complement.\n"),
   p_hat_H, 100 * p_hat_H,
-  if (p_hat_H >= 0.5) "a settled selection" else "the tie regime (family labels with near-identical membership compete; the simulation campaign mdf1 found the field's one-sided bounds at nominal in exactly this regime on this DGM family)",
+  if (p_hat_H >= 0.5) "a settled selection" else "the tie regime (family labels with near-identical membership compete)",
   fc$se_field / gc$debiased$se_wald))
 cat(sprintf(paste0(
-  "\nProvenance: forestsearch %s; ci_method = %s, %d multiplier draws (%s), field R_out/R_in = %d/%d, ij_residual = %s, seed = %d; ",
+  "\nProvenance: forestsearch %s; sg_focus = effMaxSG, effect_neighborhood = 0.20, selection_rule = neighborhood; ci_method = %s, %d multiplier draws (%s), field R_out/R_in = %d/%d, ij_residual = %s, seed = %d; ",
   "gate wall %.1f s (multiplier stage %.1f s, field %.1f s, complement field %.1f s).\n"),
   as.character(utils::packageVersion("forestsearch")), g$ci_method, g$settings$draws, g$settings$multiplier,
   f$R_out, f$R_in, g$ij_residual, params$seed,
@@ -406,12 +432,84 @@
 
 No claim is made beyond these numbers; the thresholds are the simulation design's reference points and may be replaced.
 :::
+
+## 3. Interpretation
+
+```{r interp-objects}
+#| echo: false
+band_in  <- bt[bt$in_band, ]
+band_txt <- paste(sprintf("%s (n %d, MD %.2f, consistency %.2f)", band_in$candidate, band_in$N, band_in$hr, band_in$Pcons),
+                  collapse = "; ")
+bonf_H_sup  <- say(thr_H[js$bonf_lower_H >= thr_H]);  bonf_H_not  <- say(thr_H[js$bonf_lower_H < thr_H])
+bonf_Hc_sup <- say(thr_Hc[js$bonf_upper_Hc <= thr_Hc]); bonf_Hc_not <- say(thr_Hc[js$bonf_upper_Hc > thr_Hc])
+regime <- if (p_hat_H >= 0.5) "a settled selection" else "the tie regime: family labels with near-identical membership compete for the pick, and the field's second-order correction is where the selection-adjusted bound is doing its work"
+```
+
+**The subgroup.** Under `effMaxSG` at ε = 0.20 the search selects
+**Ĥ = `` `r H_def` ``**, n = `r fmt0(n_H)` of `r fmt0(N)`, with fitted mean
+difference `r fmt(T_obs, 2)` on the oriented `y_decline` scale (a CD4 change of
+`r fmt(-T_obs, 2)` cells/mm³ on the raw scale) and consistency proportion
+`r fmt(p_cons, 2)`. It is the largest of the `r fmt0(n_band)` consistency-qualifying
+candidates whose oriented MD lies within 20% of the largest
+(`r fmt(md_max, 2)`; band floor `r fmt(md_floor, 2)`), out of `r fmt0(n_qual)`
+qualifying candidates in all. The in-band candidates it was chosen among:
+`r band_txt`.
+
+**Ĥ, by bound location.** The field's one-sided 95% lower bound on the oriented
+scale is `r fmt(L, 2)` (on the raw scale: CD4 change ≤ `r fmt(-L, 2)`). Against
+the ladder 0 / 10 / 20 / 30 / 40, "harm at least τ on Ĥ" is supported for
+τ = `r say(sup_H)` and not supported for τ = `r say(not_H)`. Beside it: the
+unadjusted lower bound `r fmt(lo1["naive"], 2)` (change ≤ `r fmt(-lo1["naive"], 2)`)
+and the IJ two-term lower bound `r fmt(lo1["ij"], 2)` (change ≤ `r fmt(-lo1["ij"], 2)`).
+The point estimates move from the naive `r fmt(g$naive$est, 2)` to `r fmt(g$debiased$est, 2)`
+(IJ, β̃) and `r fmt(f$est2, 2)` (field est₂).
+
+**Ĥᶜ, the benefit claim.** The field-s one-sided 95% upper bound on the oriented
+scale is `r fmt(U, 2)`, i.e. on the raw scale a CD4 change of at least
+`r fmt(-U, 2)` cells/mm³ on the complement. "Harm at most τ on Ĥᶜ" is supported
+for τ = `r say(sup_Hc)` and not supported for τ = `r say(not_Hc)`. The
+unstudentized complement field's upper bound is `r fmt(U_u, 2)` (change ≥
+`r fmt(-U_u, 2)`). Beside them: the unadjusted upper bound `r fmt(up1["naive"], 2)`
+and the IJ two-term upper bound `r fmt(up1["ij"], 2)`.
+
+**The pair.** The field-s Bonferroni pair (γ = 0.025 each) is Ĥ lower
+`r fmt(js$bonf_lower_H, 2)`, Ĥᶜ upper `r fmt(js$bonf_upper_Hc, 2)` on the oriented
+scale (raw: Ĥ change ≤ `r fmt(-js$bonf_lower_H, 2)`, Ĥᶜ change ≥
+`r fmt(-js$bonf_upper_Hc, 2)`), with joint probability `r fmt(js$bonf_joint_prob, 3)`
+on the `r fmt0(js$n_joint_draws)` aligned draws. Read jointly: "harm at least τ on
+Ĥ" is supported for τ = `r bonf_H_sup` and not for τ = `r bonf_H_not`; "harm at
+most τ on Ĥᶜ" is supported for τ = `r bonf_Hc_sup` and not for τ = `r bonf_Hc_not`.
+The calibrated field-s pair returns γ = `r fmt(js$gamma, 3)` and
+`r if (abs(js$gamma - 0.025) < 1e-9) "coincides with" else "is tighter than"` it.
+The unstudentized Bonferroni pair is (`r fmt(jt$bonf_lower_H, 2)`, `r fmt(jt$bonf_upper_Hc, 2)`),
+sharing the Ĥ side.
+
+**Re-selection.** p̂(Ĥ) = `r fmt(p_hat_H, 3)`: the selected subgroup wins the
+re-selection on `r pct(p_hat_H)` of the `r fmt0(g$settings$draws)` multiplier
+draws (selection rate `r fmt(g$selection_rate, 3)`, family of `r fmt0(g$n_family)`),
+which indicates `r regime`. The top re-selected labels are
+`r lab_of(names(top3)[1])` (`r fmt(top3[1], 3)`), `r lab_of(names(top3)[2])`
+(`r fmt(top3[2], 3)`) and `r lab_of(names(top3)[3])` (`r fmt(top3[3], 3)`).
+
+**Contrast with the maxeffCons analysis** (quoted from
+`REPORT_actg175_continuous_field_s_2026-09-16.md` §5 and the committed payload
+at commit `570bf8a5`; not recomputed here). Under `maxeffCons` the same data and
+settings select Ĥ = `{age <= 37} & !{cd40 <= 507}`, n = 66, with the field
+one-sided lower bound on Ĥ −43.7954 (unadjusted 17.7478, IJ −55.0512), the
+field-s one-sided upper bound on Ĥᶜ −18.9041 (unstudentized −18.5546), and the
+field-s Bonferroni pair (−56.4409, −16.4175) with joint probability 0.9491.
+The two selected subgroups share the `cd40 > 507` clause and 52 patients;
+neither contains the other. On this document's anchor the field lower bound on
+Ĥ is `r fmt(L, 2)` against the quoted −43.80, the field-s upper bound on Ĥᶜ is
+`r fmt(U, 2)` against the quoted −18.90, and the pair is
+(`r fmt(js$bonf_lower_H, 2)`, `r fmt(js$bonf_upper_Hc, 2)`) against the quoted
+(−56.44, −16.42).
+
+These are selection-adjusted bounds on one trial under one selection rule; the operating-characteristics evaluation in analysis_actg175_continuous_oc.qmd remains under maxeffCons.
+
 ```{r doc-clock}
 cat(sprintf("document compute wall-clock so far: %.1f min\n",
             (proc.time()[["elapsed"]] - t_doc) / 60))
-# TASK_actg175_continuous_intervals_2026-09-07: the OC loop needs about 11 GB per worker
-cat(sprintf("rendered with n_workers = %d (the OC loop needs about 11 GB per worker: one fs_oc_grid() job on %s draws x M = %d in one block)\n",
-            params$n_workers, format(params$draws, big.mark = ","), fam$M))
 ```
 
 # Reproducibility payload {#sec-export-payload}
@@ -427,13 +525,13 @@
 .base     <- if (is.null(results_dir)) file.path(.qmd_dir, "_payloads") else results_dir
 .out_dir  <- file.path(.base, .dirout)
 dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)
-.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))
+.payload_file <- file.path(.out_dir, paste0(.dirout, payload_suffix, "_payload.rds"))
 
 payload <- list(
-  table  = lad,
+  table  = NULL,
   labels = list(
     sg_harm  = H_def,
-    sg_focus = "maxeffCons",
+    sg_focus = "effMaxSG",
     focus    = NULL),
   meta = list(
     n_total    = N,
@@ -442,51 +540,18 @@
     event      = NULL,
     treat      = "treat",
     covariates = confounders.name,
-    c1         = 10,
-    c2         = params$c2_ratio * 10,
-    args_call  = list(consistency_method = "resample", pconsistency = 0.90,
-                      draws = params$draws, block = 5e4, seed = params$seed,
-                      M = fam$M, q_rungs = q_rungs, q_shared = q_shared,
-                      T_obs = T_obs, c1_ladder = c1_ladder, c2_vec = c2_vec,
-                      n_workers = params$n_workers, t_loop_secs = t_loop)),
+    args_call  = list(sg_focus = "effMaxSG", effect_neighborhood = eps,
+                      selection_rule = "neighborhood",
+                      consistency_method = "resample", pconsistency = 0.90,
+                      seed = params$seed, T_obs = T_obs)),
   extras = list(
-    type1 = list(diagonal = tI, c1_05_diag = c1_05_diag,
-                 c1_10_diag = c1_10_diag, at_analyst_diag = type1_at_10,
-                 surface = z, c1_05_fixed_c2_5 = c1_05_fix,
-                 homogeneous = list(diagonal = tIh, c1_05_hom = c1_05_hom,
-                                    c1_10_hom = c1_10_hom,
-                                    at_analyst_diag = hom_at_10,
-                                    at_comparator_c2_5 = hom_at_5,
-                                    surface = och),
-                 subthreshold_variants = tIv),
-    declared    = dec,
-    comparator_c2_5 = cmp,
-    calibration = lad[, c("q", "pT", "pTs", "pT5", "pT5s")],
-    anchor      = list(def = H_def, n_H = n_H, T_obs = T_obs,
-                       p_cons = p_cons, beta_treat = beta_treat),
-    # TASK_actg175_continuous_intervals_2026-09-07: the interval constructions
+    anchor      = list(def = H_def, n_H = n_H, T_obs = T_obs, p_cons = p_cons),
+    band        = list(eps = eps, md_max = md_max, md_floor = md_floor,
+                       n_qualifying = n_qual, n_in_band = n_band,
+                       candidates = bt),
+    # TASK_actg175_continuous_intervals_2026-09-07 schema: the interval constructions
     intervals   = iv,
-    purity      = list(nearest_lab = fam$lab[i_near],
-                       nearest_purity = purity_near,
-                       nearest_jaccard = jac[i_near]),
-    c2_policy   = list(headline = "c2 = c2_ratio * c1",
-                       c2_ratio = params$c2_ratio, comparator_c2 = 5,
-                       applies_to = "type-I and power (threshold policy)",
-                       calibration_read_at = c(analyst = params$c2_ratio * 10,
-                                               comparator = 5),
-                       reason = paste(
-                         "c2 sets candidate eligibility; the calibration",
-                         "interrogates the analysis as actually conducted,",
-                         "so it is read at the analyst's operating",
-                         "consistency floor, not on the policy diagonal")),
-    q_variants  = list(defs = Q_variants, prevalence = prev,
-                       table = knob, curves = cal,
-                       null_not_shared = list(dev_vs_primary = dev_vs_primary,
-                                              det_at_Tobs = det0)),
-    interpretation = list(
-      tail_table = tail_tab,
-      itt        = itt,
-      crossings  = knob[, c("variant", "P", "q05", "q05c", "q50", "q50c")])),
+    interpretation = list(itt = itt)),
   est_scale = "md",
   built_at  = Sys.time(),
   forestsearch_version = tryCatch(as.character(utils::packageVersion("forestsearch")),
```

## 4. Render — GATE PASS

`quarto render analysis_actg175_continuous_intervals.qmd` (params: `seed 8316951`), started 2026-09-16T21:21:24Z: **exit 0, wall 91 s, peak summed RSS (R + quarto + deno, 5 s sampling) 2,524 MB**. The document's own clock: 1.3 min; the MR gate 50.6 s (multiplier 19.8 s, field 23.5 s, complement field 7.4 s). Outputs: `analysis_actg175_continuous_intervals.html` (2,212,626 B) and `_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_effmaxsg_payload.rds` (2,395 B), both new.

**Anchor check** — the document's printed line:
```
H-hat = !{cd40 <= 507} & {gender};  n(H-hat) = 79;  T_obs = 71.095513;  p.consistency = 0.91000000
```
From the payload: `T_obs` 71.0955128205, |T_obs − 71.095513| = 1.8e-7 (< 1e-6); `p_cons` 0.9100000000, |Δ| = 0 (< 1e-8); the in-document assertion passed (the render ran through). The band chunk reproduces Stage 0 §3 exactly: max 87.916667, floor 70.333333, 4 of 8 in band, the same four candidates with the same n, MD and Pcons (finding F7).

**Payload checks:** 87 numeric fields of `iv`, all finite; `joint_s$bonf_lower_H` = `joint$bonf_lower_H` = −51.4511282860 (`identical()` TRUE); `settings$reselection = "effMaxSG"`, `n_family` 4935, `n_selected` 79, `field_scale_complement = "selected"`, `joint$n_joint_draws` = `joint_s$n_joint_draws` = 983.

**Facts:** `Hc$field_s$upper_1s` = **−17.6644** beside `Hc$field$upper_1s` = −17.3715 (oriented); **p̂(Ĥ) = 0.0638**.

## 5. The intervals (from the payload; oriented `y_decline`, positive = harm; raw `cd4_change` = sign flip)

**Ĥ = `!{cd40 <= 507} & {gender}`, n = 79 of 1083; T̂ = 71.10; consistency 0.91.**

| Method | est. (oriented) | SE | two-sided 95% (oriented) | one-sided 95% LOWER (oriented) | est. (cd4_change) | two-sided 95% (cd4_change) | change ≤ |
|---|---|---|---|---|---|---|---|
| Naive | 71.10 | 38.86 (robust) | (−5.07, 147.26) | 7.17 | −71.10 | (−147.26, 5.07) | −7.17 |
| MR (IJ two-term) | 21.79 (β̃) | 50.60 (IJ) | (−77.38, 120.96) | −61.44 | −21.79 | (−120.96, 77.38) | 61.44 |
| MR (field) | 16.74 (est₂); β̃ 21.79 | 32.31 (λ-SD) | (−51.45, 73.46) | −38.59 | −16.74 | (−73.46, 51.45) | 38.59 |

Two-term correction: selection term 48.07, fixed term 1.23; field Λ*-mean +5.05; field draws 983/1000 outer, 498.7/500 inner.

**Ĥᶜ, n = 1004.**

| Method | est. (oriented) | SE | two-sided 95% (oriented) | one-sided 95% UPPER (oriented) | est. (cd4_change) | two-sided 95% (cd4_change) | change ≥ |
|---|---|---|---|---|---|---|---|
| Naive | −36.15 | 7.78 (robust) | (−51.41, −20.90) | −23.35 | 36.15 | (20.90, 51.41) | 23.35 |
| MR (IJ two-term) | −32.74 (β̃ᶜ) | 15.06 (IJ) | (−62.25, −3.23) | −7.97 | 32.74 | (3.23, 62.25) | 7.97 |
| MR (field, unstudentized) | −31.10 (est₂) | 8.09 (λ-SDᶜ) | (−46.69, −15.38) | −17.37 | 31.10 | (15.38, 46.69) | 17.37 |
| MR (field-s) | −31.12 (est₂ˢ) | 7.89 (λ-SDᶜ, studentized) | (−46.32, −15.82) | −17.66 | 31.12 | (15.82, 46.32) | 17.66 |

Regime diagnostic: λ-SDᶜ / naive SE = 1.039 unstudentized, 1.014 studentized; IJ SE / naive SE = 1.935. Complement fits 1,149.

**The pair (983 aligned draws).**

| Pair | Ĥ lower (oriented) | Ĥᶜ upper (oriented) | Ĥ change ≤ | Ĥᶜ change ≥ | joint prob. |
|---|---|---|---|---|---|
| Bonferroni (γ = 0.025 each), unstudentized | −51.45 | −15.38 | 51.45 | 15.38 | 0.949 |
| Calibrated (γ = 0.025), unstudentized | −51.45 | −15.38 | 51.45 | 15.38 | 0.949 |
| **Bonferroni (γ = 0.025 each), field-s** | **−51.45** | **−15.82** | 51.45 | 15.82 | 0.949 |
| Calibrated (γ = 0.025), field-s | −51.45 | −15.82 | 51.45 | 15.82 | 0.949 |

corr(Λ*, Λ*ᶜ) = −0.038 unstudentized, −0.035 field-s.

**Re-selection:** p̂(Ĥ) = 0.064; selection rate 0.972; family 4,935; top re-selected `!{wtkg <= 81} & !{cd40 <= 413}` 0.077, `{age <= 37} & !{cd40 <= 507}` 0.077, `{age <= 40} & !{cd40 <= 507}` 0.066.

## 6. Reading (by bound location; thresholds are the simulation design's reference points)

- **Ĥ.** The field one-sided 95% lower bound is −38.59 oriented (CD4 change ≤ 38.59 raw): below every rung of the ladder 0 / 10 / 20 / 30 / 40, so "harm at least τ on Ĥ" is supported for no τ. The unadjusted lower bound 7.17 sits above 0 only; the IJ two-term bound −61.44 sits below the field's.
- **Ĥᶜ.** The field-s one-sided 95% upper bound is −17.66 oriented (CD4 change ≥ 17.66 raw): below 0, 10, 20 and 30, so "harm at most τ on Ĥᶜ" is supported at every threshold used, and a CD4 change of at least 17.66 cells/mm³ on the complement is supported. The unstudentized bound is −17.37; naive −23.35; IJ −7.97.
- **The pair.** The field-s Bonferroni pair (−51.45, −15.82) with joint probability 0.949: the Ĥ side supports no rung, the Ĥᶜ side supports 0–30; the calibrated split returns γ = 0.025 and coincides with it.
- **Regime.** p̂(Ĥ) = 0.064 — the tie regime; the top three re-selected labels are near 0.07–0.08 each, and the `maxeffCons` anchor is the second at 0.077.
- **Against the ladder, nothing reads differently from the `maxeffCons` analysis**: no rung on Ĥ under either rule, all four on Ĥᶜ under either.

## 7. The maxeffCons contrast (quoted, `570bf8a5`; not recomputed)

Under `maxeffCons`: Ĥ = `{age <= 37} & !{cd40 <= 507}`, n = 66; Ĥ one-sided lower bounds naive 17.7478, IJ −55.0512, field −43.7954; Ĥᶜ upper bounds field-s −18.9041, unstudentized −18.5546; field-s Bonferroni pair (−56.4409, −16.4175), joint probability 0.9491. Under `effMaxSG` (this document): Ĥ = `!{cd40 <= 507} & {gender}`, n = 79; field lower −38.59 (5.2 MD units closer to zero); field-s upper −17.66 (1.2 closer to zero); pair (−51.45, −15.82); p̂(Ĥ) 0.064 against 0.0872. The two subgroups share 52 patients.

## 8. Findings

- **F1.** The task's gate reads "MD within 1e-6 of 71.10"; Stage 0's value is 71.095513, and 71.10 is its rounding (|71.095513 − 71.10| = 4.5e-3 would fail a 1e-6 tolerance). The document asserts against 71.095513; observed |Δ| = 1.8e-7.
- **F2.** The field-s record does not state the unadjusted and IJ one-sided lower bounds on Ĥ (§5 gives the field bound and the Ĥᶜ comparators); those two were taken from the committed payload at the same commit `570bf8a5` and are cited as such in `<doc>`.
- **F3.** The setup chunk is not literally unchanged: `dirout`'s comment and the header note are reworded, `payload_suffix` is added, and the four OC rung/ladder/`c2_vec` lines are dropped with `params$c2_ratio` (the task names only `draws` and `n_workers` as dropped; `c2_ratio` is read nowhere else). The helper functions `at` and `qlab`, OC-only, were kept as definitions.
- **F4.** The reading callout's tie-regime sentence in `<oc>` cites the `mdf1` campaign (a `maxeffCons` campaign); that clause is removed in `<doc>` rather than re-attributed, since no `mdsgnb20` claim was verified for it.
- **F5.** An 11-day-old orphaned `bash` memory-monitor loop from a `gbsg_020` `uburst` render (PID 2086198, child `sleep`, writing `mem_uburst.txt` in an old scratchpad) survives on the host; no R or quarto process. Not killed (out of scope).
- **F6.** Under `effMaxSG` the top re-selected label on the multiplier draws is `!{wtkg <= 81} & !{cd40 <= 413}` (0.077), not the selected subgroup (0.064) and not any in-band candidate but one (`{age <= 40} & !{cd40 <= 507}`, third at 0.066); the tie regime is the same regime the `maxeffCons` gate reported (0.087).
- **F7.** The fitted object's consistency table (`out_sg$result`, 8 rows, label columns `M.1`/`M.2`) reproduces the Stage 0 band exactly; Stage 0's F4 (the enumerated family of 4,935 is not carried with `details = FALSE`) still holds — the band is over the 8 consistency-qualifying candidates only.
- **F8.** Cost: 91 s wall, 2.5 GB peak, against the OC document's 53.6 min and 90.5 GB — the intervals product is separable from the OC evaluation at negligible cost.

## 9. Commits

```
e1c3c7a6 Add TASK_actg175_continuous_intervals_2026-09-16 as received
da25abcc ACTG175 continuous: intervals document under effMaxSG, eps = 0.20 (TASK_actg175_continuous_intervals_2026-09-16) -- ...
<the closeout commit: the HTML, the payload, this record>
```
