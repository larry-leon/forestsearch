# REPORT — ACTG175 continuous applied analysis: the field-s complement bound and pair

Date: 2026-09-16. Machine: `pop-os` (AMD Ryzen Threadripper PRO 5995WX, 64 physical cores, 251 GB; R 4.6.1, reference BLAS/LAPACK 3.12.0). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_actg175_continuous_field_s_2026-09-16.md` (committed as received, `98643ba9`), on Larry's D5 disposition of the MD field re-run's Gate 0. Document: `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` (`<doc>`). Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix`; no `R/` change. The render ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` and the document's committed `params` (`draws: 20000`, `n_workers: 14`, `c2_ratio: 0.8`, `seed: 8316951`).

## 1. Provenance — GATE PASS

```
pop-os
feature/glm-extension
89df2cf3
[tracked modifications: none]
MD re-run closeout in HEAD
[R / Rscript / quarto processes: none]
R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix
```
- HEAD at kickoff `89df2cf3`, three merges past `b77ca649` (the MD re-run closeout): `git diff --stat b77ca649..89df2cf3` is two PDFs under `quarto/methodology/articles/` and `dev/tasks/TASK_md_field_rerun_stage0_2026-09-15.md` re-created upstream; nothing under `R/` or `quarto/applications/actg175/`.
- First commit: `98643ba9 Add TASK_actg175_continuous_field_s_2026-09-16 as received`.

## 2. Read (at `89df2cf3`)

**2.1 The anchor's focus.** `<doc>:92–132`:
```
fs_anchor <- forestsearch(
  df.analysis      = actg_df,
  confounders.name = confounders.name,
  outcome.name     = "y_decline",
  treat.name       = "treat",
  id.name          = "id",
  outcome_type     = "continuous",
  effect_measure   = "MD",
  adverse_outcome  = TRUE,
  seedit           = params$seed,
  sg_focus         = "maxeffCons",
  selection_rule   = "neighborhood",
  consistency_method = "resample",
  effect.threshold       = 10,
  consistency.threshold  = 5,
  pconsistency.threshold = 0.90,
  use_twostage     = TRUE,
  conf.cont_jcuts  = list(age = 10, preanti = 10, wtkg = 10,
                          karnof = 10, cd40 = 10, cd80 = 10),
  cut_type         = "default",
  maxk             = 2L,
  n.min            = 60L,
  d0.min           = 10L,
  d1.min           = 10L,
  fs.splits        = 500L,
  use_lasso        = FALSE,
  use_grf          = FALSE,
  use_dina         = FALSE,
  is.RCT           = TRUE,
  parallel_args    = list(plan = "sequential"),
  details          = FALSE,
  quiet            = TRUE,
  mr_inference     = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L,
                           include_complement = TRUE, field_complement = TRUE,
                           return_reselection = TRUE)
)
```
(comment lines `:124–127` omitted). The focus is `sg_focus = "maxeffCons"` with `selection_rule = "neighborhood"`; `effect_neighborhood` is not passed (the default 0.10, inert under this focus). MR arguments: `ci_method = "field"`, `draws = 5000L`, `include_complement = TRUE`, `field_complement = TRUE`, `return_reselection = TRUE`; `field_scale_complement` not passed (inherited `"selected"` since `fb62705c`; `"none"` when the committed payload was rendered on 2026-09-08).
- Mapping, `R/fs_mr_inference_methods.R:98–103`:
```
         # maxeffCons -> MR's "maxeff", which is
         # passers[which.max(beta[passers])]: argmax effect among PASSERS, and
         # `passers` is MR's consistency-qualifying set (driven by
         # p_star = pconsistency.threshold).  So the MR rule named "maxeff"
         # is the maxeffCons rule, NOT sg_focus = "maxeff" (which is ungated).
         maxeffCons = "maxeff",
```
- **`settings$reselection = maxeff` in the payload names MR's label for the identifier's `maxeffCons` focus** (the effect argmax over the consistency-qualifying passers), not a focus of `"maxeff"`; the identifier's focus is `maxeffCons`. (Stage 0 §9.3 of the MD re-run read this field as the applied gate re-selecting "under `maxeff`, not the simulation's `maxeffCons`"; that reading was wrong: the two are the same rule under MR's naming. Finding F1.)

**2.2 The identity reference.** `<doc>` writes one payload file: `:1243–1247` `.dirout <- if (is.null(dirout)) .qmd_stem else dirout` (`dirout <- "analysis_actg175_continuous_oc_intervals"`, `:28`), `.base <- file.path(.qmd_dir, "_payloads")`, `.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))`; `:1312` `saveRDS(payload, .payload_file)`. That file, `quarto/applications/actg175/_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_payload.rds` (13,875 B), is **tracked**, as is the rendered `analysis_actg175_continuous_oc.html` (2,237,401 B), which `quarto render` rewrites beside the source. No other file is written (no `write.csv`, `writeLines`, `ggsave` or `png()`; figures are embedded). The committed OC payload `_payloads/analysis_actg175_continuous_oc/analysis_actg175_continuous_oc_payload.rds` is tracked but not written by this document. Copies for §4: `/tmp/fs_ref_LG6E/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_payload.rds` and `/tmp/fs_ref_LG6E/analysis_actg175_continuous_oc.html`.

**2.3 The interpretation.** Section "Reading the intervals (generated from the numbers above at render time)", chunk `intervals-reading` (`:345–387`), **inline R** (`cat(sprintf(...))` over the gate objects `f`, `fc`, `jt`, `g`, `gc`; no literal number). Sentences reading the complement bound or the joint pair, as committed:
- `:346` `U <- fc$upper_1s; L <- f$lower_1s` and `:350–354` "**Complement first (the benefit claim).** The field's one-sided 95% upper bound on Ĥᶜ is %.2f on the oriented scale (a CD4 change of at least %.2f on the raw scale). Against the reading thresholds, "harm at most τ on Ĥᶜ" is supported for τ = %s and not supported for τ = %s. Comparators: MR (IJ) upper %.2f (change ≥ %.2f), naive upper %.2f (change ≥ %.2f)." — reads `U = fc$upper_1s`.
- `:359–364` "**MR (IJ two-term) as the conservative reference:** … upper bound on Ĥᶜ %.2f …; its bounds sit %s the field's on both blocks." — compares `up1["ij"] >= U`.
- `:365–369` "**A claim on both subgroups** uses the Bonferroni pair (Ĥ lower %.2f, Ĥᶜ upper %.2f; on the cd4_change scale: Ĥ change ≤ %.2f, Ĥᶜ change ≥ %.2f); the calibrated pair at γ = %.3f %s it (corr(Λ*, Λ*ᶜ) = %+.3f)." — reads `jt$bonf_lower_H`, `jt$bonf_upper_Hc`, `jt$gamma`, `jt$corr`.
- `:370–374` "**The price of selection:** … on Ĥᶜ the adjusted upper bounds span %.2f (field) to %.2f (IJ) against the naive %.2f." — reads `U`.
- `:375–380` "**Regime:** … λ-SDᶜ / naive SE = %.3f on the complement." — reads `fc$se_field` (a diagnostic, not a bound; left as it was, per E5's "no other prose changes").

## 3. Edits (`<doc>` only; the diff as applied, `git diff` at HEAD `98643ba9`)

```diff
diff --git a/quarto/applications/actg175/analysis_actg175_continuous_oc.qmd b/quarto/applications/actg175/analysis_actg175_continuous_oc.qmd
index 299727fc..9d2c505e 100644
--- a/quarto/applications/actg175/analysis_actg175_continuous_oc.qmd
+++ b/quarto/applications/actg175/analysis_actg175_continuous_oc.qmd
@@ -128,6 +128,9 @@ fs_anchor <- forestsearch(
   mr_inference     = TRUE,
   mr_inference_args = list(ci_method = "field", draws = 5000L,
                            include_complement = TRUE, field_complement = TRUE,
+                           # TASK_actg175_continuous_field_s_2026-09-16: the studentized
+                           # complement field (field-s), the package default, passed explicitly.
+                           field_scale_complement = "selected",
                            return_reselection = TRUE)
 )
 
@@ -179,7 +182,7 @@ stopifnot(is.list(g), identical(g$ci_method, "field"),
           is.list(g$field$complement), is.null(g$field$complement$note),
           is.list(g$field$joint), is.null(g$field$joint$note),
           is.list(g$reselection), is.list(g$complement), !is.null(g$complement$debiased))
-f  <- g$field; fc <- f$complement; jt <- f$joint; gc <- g$complement; rs <- g$reselection
+f  <- g$field; fc <- f$complement; jt <- f$joint; js <- f$joint_s; gc <- g$complement; rs <- g$reselection
 z95 <- qnorm(0.95)
 # Family labels are cut codes q<k>.<level>: k indexes fs_anchor$confounders.evaluated,
 # level 1 = the cut holds, level 0 = its negation -- decoded here for display.
@@ -201,7 +204,7 @@ top3    <- sort(rs$p_hat, decreasing = TRUE)[1:3]
 # one-sided bounds on the oriented scale (harm block: lower; complement: upper)
 lo1 <- c(naive = g$naive$est - z95 * g$debiased$se_wald, ij = g$debiased$lower_1s, field = f$lower_1s)
 up1 <- c(naive = gc$naive$est + z95 * gc$debiased$se_wald, ij = gc$debiased$est + z95 * gc$debiased$se_ij,
-         field = fc$upper_1s)
+         field = fc$upper_1s, field_s = fc$upper_1s_s)
 fmt2 <- function(l, h) sprintf("(%.2f, %.2f)", l, h)
 # the payload element (schema: TASK_actg175_continuous_intervals_2026-09-07)
 iv <- list(
@@ -210,7 +213,8 @@ iv <- list(
                     note = "simulation design's reference points (0/10/30/40 = no harm / consistency threshold / effect threshold / planted MD)"),
   settings = list(ci_method = g$ci_method, draws = g$settings$draws, multiplier = g$settings$multiplier,
                   reselection = g$settings$reselection, ij_residual = g$ij_residual, seed = params$seed,
-                  field_R_out = f$R_out, field_R_in = f$R_in, n_family = g$n_family, n_selected = g$n_selected),
+                  field_R_out = f$R_out, field_R_in = f$R_in, n_family = g$n_family, n_selected = g$n_selected,
+                  field_scale_complement = if (!is.null(fc$est2_s)) "selected" else "none"),
   H = list(naive = list(est = g$naive$est, lower = g$naive$lower, upper = g$naive$upper,
                         se_wald = g$debiased$se_wald, lower_1s = unname(lo1["naive"])),
            ij = list(est = g$debiased$est, lower = g$debiased$lower, upper = g$debiased$upper,
@@ -228,10 +232,17 @@ iv <- list(
                          lower_2s = fc$lower_2s, upper_2s = fc$upper_2s, lambda_sd = fc$se_field,
                          lambda_mean = fc$lambda_mean, q05 = fc$q05, q95 = fc$q95,
                          n_out_used = fc$n_out_used, n_complement_fits = fc$n_complement_fits,
-                         timing_seconds = fc$timing_seconds)),
+                         timing_seconds = fc$timing_seconds),
+            # field-s (TASK_actg175_continuous_field_s_2026-09-16): the studentized complement field, the same names as Hc$field
+            field_s = list(est2 = fc$est2_s, upper_1s = fc$upper_1s_s, lower_1s = fc$lower_1s_s,
+                           lower_2s = fc$lower_2s_s, upper_2s = fc$upper_2s_s, lambda_sd = fc$se_field_s,
+                           lambda_mean = fc$lambda_mean_s)),
   joint = list(gamma = jt$gamma, lower_H = jt$lower_H, upper_Hc = jt$upper_Hc, joint_prob = jt$joint_prob,
                bonf_lower_H = jt$bonf_lower_H, bonf_upper_Hc = jt$bonf_upper_Hc,
                bonf_joint_prob = jt$bonf_joint_prob, corr = jt$corr, n_joint_draws = jt$n_joint_draws),
+  joint_s = list(gamma = js$gamma, lower_H = js$lower_H, upper_Hc = js$upper_Hc, joint_prob = js$joint_prob,
+                 bonf_lower_H = js$bonf_lower_H, bonf_upper_Hc = js$bonf_upper_Hc,
+                 bonf_joint_prob = js$bonf_joint_prob, corr = js$corr, n_joint_draws = js$n_joint_draws),
   reselection = list(p_hat_H = p_hat_H, selection_rate = g$selection_rate,
                      top3_labels = vapply(names(top3), lab_of, character(1)), top3_codes = names(top3),
                      top3_p_hat = unname(top3)),
@@ -275,17 +286,18 @@ cat(sprintf(paste0(
 
 ```{r intervals-table-Hc}
 tab2 <- data.frame(
-  Method = c("Naive", "MR (IJ two-term)", "MR (field)"),
+  Method = c("Naive", "MR (IJ two-term)", "MR (field, unstudentized)", "MR (field-s)"),
   `Point est. (oriented)` = c(sprintf("%.2f", gc$naive$est), sprintf("%.2f (β̃ᶜ)", gc$debiased$est),
-                              sprintf("%.2f (est₂); β̃ᶜ %.2f", fc$est2, gc$debiased$est)),
+                              sprintf("%.2f (est₂); β̃ᶜ %.2f", fc$est2, gc$debiased$est),
+                              sprintf("%.2f (est₂ˢ); β̃ᶜ %.2f", fc$est2_s, gc$debiased$est)),
   `SE` = c(sprintf("%.2f (robust)", gc$debiased$se_wald), sprintf("%.2f (IJ)", gc$debiased$se_ij),
-           sprintf("%.2f (λ-SDᶜ)", fc$se_field)),
+           sprintf("%.2f (λ-SDᶜ)", fc$se_field), sprintf("%.2f (λ-SDᶜ, studentized)", fc$se_field_s)),
   `Two-sided 95% (oriented)` = c(fmt2(gc$naive$lower, gc$naive$upper), fmt2(gc$debiased$lower, gc$debiased$upper),
-                                 fmt2(fc$lower_2s, fc$upper_2s)),
+                                 fmt2(fc$lower_2s, fc$upper_2s), fmt2(fc$lower_2s_s, fc$upper_2s_s)),
   `One-sided 95% UPPER (oriented)` = sprintf("%.2f", up1),
-  `Point est. (cd4_change)` = sprintf("%.2f", -c(gc$naive$est, gc$debiased$est, fc$est2)),
+  `Point est. (cd4_change)` = sprintf("%.2f", -c(gc$naive$est, gc$debiased$est, fc$est2, fc$est2_s)),
   `Two-sided 95% (cd4_change)` = c(fmt2(-gc$naive$upper, -gc$naive$lower), fmt2(-gc$debiased$upper, -gc$debiased$lower),
-                                   fmt2(-fc$upper_2s, -fc$lower_2s)),
+                                   fmt2(-fc$upper_2s, -fc$lower_2s), fmt2(-fc$upper_2s_s, -fc$lower_2s_s)),
   `One-sided 95% bound (cd4_change): change ≥` = sprintf("%.2f", -up1),
   check.names = FALSE)
 gt(tab2) |>
@@ -293,26 +305,29 @@ gt(tab2) |>
              subtitle = sprintf("Ĥᶜ: n = %d; the one-sided UPPER bound is the exposed limit for a benefit claim (\"harm at most U\" = \"benefit at least −U\")",
                                 N - n_H)) |>
   tab_footnote(footnote = sprintf(paste(
-    "One-sided conventions: naive est + 1.645·SE(robust); MR (IJ) β̃ᶜ + 1.645·SE_IJ; MR (field) the gate's stored upper_1s = β̃ᶜ − q₀.₀₅(Λ*ᶜ).",
-    "Regime diagnostic: λ-SDᶜ / naive SE = %.3f (the IJ SE / naive SE = %.3f)."),
-    fc$se_field / gc$debiased$se_wald, gc$debiased$se_ij / gc$debiased$se_wald))
+    "One-sided conventions: naive est + 1.645·SE(robust); MR (IJ) β̃ᶜ + 1.645·SE_IJ; MR (field) the gate's stored upper_1s = β̃ᶜ − q₀.₀₅(Λ*ᶜ);",
+    "MR (field-s) the studentized complement field's upper_1s_s = β̃ᶜ − q₀.₀₅(Λ*ᶜₛ), each reading rescaled to the selected complement's influence-norm scale (field_scale_complement = \"selected\", the package default; the evaluated complement product).",
+    "Regime diagnostic: λ-SDᶜ / naive SE = %.3f unstudentized, %.3f studentized (the IJ SE / naive SE = %.3f)."),
+    fc$se_field / gc$debiased$se_wald, fc$se_field_s / gc$debiased$se_wald, gc$debiased$se_ij / gc$debiased$se_wald))
 ```
 
 ```{r intervals-table-joint}
 tab3 <- data.frame(
-  Pair = c("Separate one-sided 95% field bounds", "Bonferroni (γ = 0.025 each)",
-           sprintf("Calibrated from field$joint (γ = %.3f)", jt$gamma)),
-  `Ĥ lower (oriented)` = sprintf("%.2f", c(f$lower_1s, jt$bonf_lower_H, jt$lower_H)),
-  `Ĥᶜ upper (oriented)` = sprintf("%.2f", c(fc$upper_1s, jt$bonf_upper_Hc, jt$upper_Hc)),
-  `Ĥ: change ≤ (cd4_change)` = sprintf("%.2f", -c(f$lower_1s, jt$bonf_lower_H, jt$lower_H)),
-  `Ĥᶜ: change ≥ (cd4_change)` = sprintf("%.2f", -c(fc$upper_1s, jt$bonf_upper_Hc, jt$upper_Hc)),
+  Pair = c("Separate one-sided 95% field bounds (unstudentized)", "Bonferroni (γ = 0.025 each) (unstudentized)",
+           sprintf("Calibrated from field$joint (γ = %.3f) (unstudentized)", jt$gamma),
+           "Bonferroni (γ = 0.025 each), field-s", sprintf("Calibrated from field$joint_s (γ = %.3f), field-s", js$gamma)),
+  `Ĥ lower (oriented)` = sprintf("%.2f", c(f$lower_1s, jt$bonf_lower_H, jt$lower_H, js$bonf_lower_H, js$lower_H)),
+  `Ĥᶜ upper (oriented)` = sprintf("%.2f", c(fc$upper_1s, jt$bonf_upper_Hc, jt$upper_Hc, js$bonf_upper_Hc, js$upper_Hc)),
+  `Ĥ: change ≤ (cd4_change)` = sprintf("%.2f", -c(f$lower_1s, jt$bonf_lower_H, jt$lower_H, js$bonf_lower_H, js$lower_H)),
+  `Ĥᶜ: change ≥ (cd4_change)` = sprintf("%.2f", -c(fc$upper_1s, jt$bonf_upper_Hc, jt$upper_Hc, js$bonf_upper_Hc, js$upper_Hc)),
   `Joint prob. on the aligned draws` = c("< 0.95 by construction", sprintf("%.3f", jt$bonf_joint_prob),
-                                         sprintf("%.3f", jt$joint_prob)),
+                                         sprintf("%.3f", jt$joint_prob), sprintf("%.3f", js$bonf_joint_prob),
+                                         sprintf("%.3f", js$joint_prob)),
   check.names = FALSE)
 gt(tab3) |>
   tab_header(title = "Table 3 — the joint pair (Ĥ lower bound, Ĥᶜ upper bound)",
-             subtitle = sprintf("γ = %.3f on the 0.025–0.050 grid; corr(Λ*, Λ*ᶜ) = %+.3f over %d aligned draws",
-                                jt$gamma, jt$corr, jt$n_joint_draws))
+             subtitle = sprintf("unstudentized: γ = %.3f on the 0.025–0.050 grid, corr(Λ*, Λ*ᶜ) = %+.3f over %d aligned draws; field-s: γ = %.3f, corr = %+.3f over %d",
+                                jt$gamma, jt$corr, jt$n_joint_draws, js$gamma, js$corr, js$n_joint_draws))
 ```
 
 ```{r intervals-diagnostics}
@@ -343,15 +358,16 @@ gt(tab4) |> tab_header(title = "Re-selection and field diagnostics",
 ## Reading the intervals (generated from the numbers above at render time)
 
 ```{r intervals-reading, results = "asis"}
-U <- fc$upper_1s; L <- f$lower_1s
+U <- fc$upper_1s_s; U_u <- fc$upper_1s; L <- f$lower_1s   # field-s is the evaluated complement bound; U_u the unstudentized before/after
 sup_Hc <- thr_Hc[U <= thr_Hc]; not_Hc <- thr_Hc[U > thr_Hc]
 sup_H  <- thr_H[L >= thr_H];  not_H  <- thr_H[L < thr_H]
 say <- function(v) if (length(v)) paste(v, collapse = ", ") else "none"
 cat(sprintf(paste0(
-  "- **Complement first (the benefit claim).** The field's one-sided 95%% upper bound on Ĥᶜ is %.2f on the oriented scale ",
+  "- **Complement first (the benefit claim).** The field-s one-sided 95%% upper bound on Ĥᶜ is %.2f on the oriented scale ",
   "(a CD4 change of at least %.2f on the raw scale). Against the reading thresholds, \"harm at most τ on Ĥᶜ\" is supported for τ = %s ",
-  "and not supported for τ = %s. Comparators: MR (IJ) upper %.2f (change ≥ %.2f), naive upper %.2f (change ≥ %.2f).\n"),
-  U, -U, say(sup_Hc), say(not_Hc), up1["ij"], -up1["ij"], up1["naive"], -up1["naive"]))
+  "and not supported for τ = %s. The unstudentized field's upper bound, the paired before-and-after, is %.2f (change ≥ %.2f). ",
+  "Comparators: MR (IJ) upper %.2f (change ≥ %.2f), naive upper %.2f (change ≥ %.2f).\n"),
+  U, -U, say(sup_Hc), say(not_Hc), U_u, -U_u, up1["ij"], -up1["ij"], up1["naive"], -up1["naive"]))
 cat(sprintf(paste0(
   "- **Ĥ (the harm claim).** The field's one-sided 95%% lower bound is %.2f on the oriented scale (a CD4 change of at most %.2f). ",
   "\"Harm at least τ on Ĥ\" is supported for τ = %s and not supported for τ = %s. The naive lower bound is %.2f.\n"),
@@ -363,13 +379,15 @@ cat(sprintf(paste0(
   gc$debiased$se_ij / gc$debiased$se_wald,
   if (lo1["ij"] <= L && up1["ij"] >= U) "outside" else "not uniformly outside"))
 cat(sprintf(paste0(
-  "- **A claim on both subgroups** uses the Bonferroni pair (Ĥ lower %.2f, Ĥᶜ upper %.2f; on the cd4_change scale: Ĥ change ≤ %.2f, ",
-  "Ĥᶜ change ≥ %.2f); the calibrated pair at γ = %.3f %s it (corr(Λ*, Λ*ᶜ) = %+.3f).\n"),
-  jt$bonf_lower_H, jt$bonf_upper_Hc, -jt$bonf_lower_H, -jt$bonf_upper_Hc, jt$gamma,
-  if (abs(jt$gamma - 0.025) < 1e-9) "coincides with" else "is tighter than", jt$corr))
+  "- **A claim on both subgroups** uses the field-s Bonferroni pair (Ĥ lower %.2f, Ĥᶜ upper %.2f; on the cd4_change scale: Ĥ change ≤ %.2f, ",
+  "Ĥᶜ change ≥ %.2f): \"harm at least τ on Ĥ\" is supported for τ = %s and \"harm at most τ on Ĥᶜ\" for τ = %s; ",
+  "the calibrated field-s pair at γ = %.3f %s it (corr(Λ*, Λ*ᶜₛ) = %+.3f).\n"),
+  js$bonf_lower_H, js$bonf_upper_Hc, -js$bonf_lower_H, -js$bonf_upper_Hc,
+  say(thr_H[js$bonf_lower_H >= thr_H]), say(thr_Hc[js$bonf_upper_Hc <= thr_Hc]), js$gamma,
+  if (abs(js$gamma - 0.025) < 1e-9) "coincides with" else "is tighter than", js$corr))
 cat(sprintf(paste0(
   "- **The price of selection:** on Ĥ the adjusted one-sided lower bounds span %.2f (field) to %.2f (IJ) against the naive %.2f; ",
-  "on Ĥᶜ the adjusted upper bounds span %.2f (field) to %.2f (IJ) against the naive %.2f. The point estimates move from the naive %.2f ",
+  "on Ĥᶜ the adjusted upper bounds span %.2f (field-s) to %.2f (IJ) against the naive %.2f. The point estimates move from the naive %.2f ",
   "to %.2f (IJ) and %.2f (field est₂) on Ĥ.\n"),
   L, lo1["ij"], lo1["naive"], U, up1["ij"], up1["naive"], g$naive$est, g$debiased$est, f$est2))
 cat(sprintf(paste0(
```

## 4. Render — GATE

One render of `<doc>` (`quarto render analysis_actg175_continuous_oc.qmd`, committed `params`: `draws 20000`, `n_workers 14`, `c2_ratio 0.8`, `seed 8316951`; `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`; GNU `timeout 5400`; `scripts_mdsgnb20/mem_sampler.sh` at 5 s): **exit 0, wall 3,214 s (53.6 min), peak summed RSS 90,526 MB**; the document's own clock: evaluation loop 2,089.7 s over 20 jobs at 14 workers (`meta$args_call$t_loop_secs`); the MR gate 46.7 s (multiplier 20.5 s, field 19.4 s, complement field 6.8 s). Outputs rewritten in place: `analysis_actg175_continuous_oc.html` (2,894,253 B) and `_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_payload.rds` (13,648 B). The 53.6-min wall is above the task's 15–30-min estimate (finding F4) and under the 90-min timeout.

**Identity against the §2.2 copy** (`identity_fields.R`, flattened field by field; character/integer/logical `identical()`, numeric ≤ 1e-8 relative or ≤ 1e-10 absolute):

- Fields in the reference: **410**; compared: 410; missing in the new payload: 0; new fields: 17 (`settings$field_scale_complement`, the seven `Hc$field_s` fields, the nine `joint_s` fields).
- **Selected subgroup:** `{age <= 37} & !{cd40 <= 507}`, `n_selected` 66, identical; `anchor$def`, `n_H`, `T_obs` (87.9166666667) identical; `labels` identical; **no label-valued field differs** (no label tie to enumerate).
- **Numeric statistical fields: 352 compared, all within tolerance; largest relative difference 1.3e-9** (`extras$type1$diagonal$Enaive_bias`); the OC table's numeric columns agree to 9.95e-14; every interval field (`H`, `Hc`, `joint`, `reselection`) agrees to ≤ 1e-9 relative. The Mac copy was built under R 4.5.2 with Accelerate; this render under R 4.6.1 with reference BLAS — floating-point agreement, not bit identity, as the task expects.
- **Seven fields outside tolerance, all compute provenance, none statistical:** `meta$args_call$n_workers` (1 → 14: the committed payload came from the Mac render at `-P n_workers:1`, this one from the committed `params`), `meta$args_call$t_loop_secs` (1,123.4 → 2,089.7 s), and the five wall-clock fields `extras$intervals$H$field$timing_seconds`, `Hc$field$timing_seconds`, `timing_seconds$multiplier_stage`, `$field`, `$complement_field` (seconds). `built_at` differs by construction. **Classification:** these cannot reproduce on any re-render and carry no number the analysis reads; every prior identity gate in this repository excluded timings (`mdf1` Gate 2: "minus timings, messages and the FB columns"). They are enumerated here and excluded from the gate's tolerance rule; the task's rule read literally ("a numeric field outside tolerance" fails) would fail on wall-clock seconds, which is recorded as finding F2, and the gate is **passed on the statistical content**: same subgroup, every statistical numeric field within tolerance, no label tie.
- **New fields:** all seven `Hc$field_s` and all nine `joint_s` fields finite; `settings$field_scale_complement = "selected"`; `joint_s$n_joint_draws` = `joint$n_joint_draws` = 983, and `joint_s$bonf_lower_H` = `joint$bonf_lower_H` = −56.4408886659 exactly (the same harm draws).

**GATE 4: PASS** (on the classification above).

## 5. The field-s complement bound and pair beside the unstudentized values (oriented `y_decline`, positive = harm; raw `cd4_change` is the sign flip; from the new payload)

| quantity | unstudentized field | field-s |
|---|---|---|
| Hc est2 (oriented) | -31.9645 | -31.9737 |
| Hc lambda-SD | 8.0668 | 7.9578 |
| Hc one-sided 95% UPPER bound (oriented) | -18.5546 | -18.9041 |
| Hc two-sided 95% (oriented) | (-46.8086, -16.3473) | (-46.8281, -16.4175) |
| Bonferroni pair: H lower | -56.4409 | -56.4409 |
| Bonferroni pair: Hc upper | -16.3473 | -16.4175 |
| Bonferroni joint prob. on aligned draws | 0.9491 | 0.9491 |
| calibrated gamma | 0.0250 | 0.0250 |
| calibrated pair: H lower / Hc upper | -56.4409 / -16.3473 | -56.4409 / -16.4175 |
| corr(Lambda*, Lambda*c) | -0.0191 | -0.0169 |

Comparators (unchanged): Ĥ field one-sided lower bound −43.7954; Ĥᶜ naive upper −23.6581; Ĥᶜ IJ upper −8.7454. Thresholds: Ĥ lower 0 / 10 / 20 / 30 / 40; Ĥᶜ upper 0 / 10 / 20 / 30.

- **The field-s upper bound on Ĥᶜ is −18.90 on the oriented scale** (a CD4 change of at least 18.90 raw): "harm at most τ on Ĥᶜ" is supported for every threshold the document uses (τ = 0, 10, 20, 30). The unstudentized bound, the paired before-and-after, is −18.55: the studentization moves it by −0.35 MD (λ-SDᶜ 7.96 against 8.07; est₂ −31.97 against −31.96) and changes no threshold reading.
- **The field-s Bonferroni pair** is (Ĥ lower −56.44, Ĥᶜ upper −16.42): "harm at least τ on Ĥ" is supported for none of 0–40, "harm at most τ on Ĥᶜ" for all of 0–30; joint probability on the 983 aligned draws 0.949; the calibrated split returns γ = 0.025 and coincides with it. Its Ĥ side is identical to the unstudentized pair's (same harm draws); its Ĥᶜ side differs by −0.07.
- Between the comparators, the field-s upper bound sits between the naive −23.66 and IJ's −8.75, as the unstudentized bound did; nothing the document read by location changes under field-s on this data.

## 6. Findings

- **F1.** The payload's `settings$reselection = maxeff` is MR's label for the `maxeffCons` focus (`R/fs_mr_inference_methods.R:103`); the MD re-run's Stage 0 record (§9.3) read it as a different rule. Corrected here.
- **F2.** Identity gate classification: seven timing / worker-count fields differ from the Mac copy by construction and were excluded from the tolerance rule; the rule read literally would fail on wall-clock seconds. Every statistical field is within 1.3e-9 relative.
- **F3.** The committed payload's `n_workers` was 1 (the Mac render); this render used the committed 14. The OC loop took 2,090 s at 14 workers here against 1,123 s at 1 worker on the M4 Max: on this host the 14-way `mclapply` of `fs_oc_grid()` jobs (≈ 6 GB each, 90 GB peak) is slower than the Mac's serial loop.
- **F4.** Render wall 53.6 min, above the task's 15–30-min estimate, under the 90-min timeout.
- **F5.** Field-s equals the unstudentized complement field on this data to 0.35 MD in the upper bound (λ-SDᶜ ratio 0.986), as on the four `mdsgnb20` cells (≤ 0.003 in coverage).
- **F6.** The document's provenance line prints the gate's wall seconds, a displayed number that cannot reproduce on any render; every other displayed number reproduces (the tables read the payload fields that agree to ≤ 1e-9).

## 7. Commits

```
98643ba9 Add TASK_actg175_continuous_field_s_2026-09-16 as received
<the edit commit: <doc>, its HTML, the tracked payload, this record>
```
