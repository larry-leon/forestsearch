# REPORT — Reporting the certified survival products: the `print`/`summary` extension and the survival post-selection vignette

Date: 2026-09-09. Task: `dev/tasks/TASK_print_vignette_2026-09-09.md` (committed as received, 8b9da196). Executor: Claude Code (Linux), unattended. **Report-and-wait.** Classification: **adds code; existing output byte-identical when no MR results are present.** Compute: verification renders and one vignette build only.

Source of truth for every printed and written claim: `dev/notes/NOTE_survival_products_2026-09-09.md`.

## Stage 0 — Discovery

### 1. The MR results are reachable — no STOP

`forestsearch()` attaches the `fs_mr_inference()` return to its result. The gate runs at `R/forestsearch_main.R:3352` into a local `mr_out`:

```r
    mr_out <- tryCatch({
```

and the returned list carries it under the name **`mr_inference`** (`R/forestsearch_main.R:3479-3494`):

```r
  out <- list(
    grp.consistency = grp.consistency,
    ...
    # Multiplier-resampling result (NULL unless mr_inference = TRUE)
    mr_inference = mr_out,
    mr_harm_confirmed = if (!is.null(mr_out)) mr_out$harm_flag
                         else NA,
```

with `class(out) <- c("forestsearch", "list")` at `:3527`. So `x$mr_inference` is exactly the `fs_mr_inference()` return, and both methods can reach every certified product. The two other attachment sites (`:2337`, `:2522`, the DINA/GRF branches via `.fs_apply_mr()`) use the same name. **The object model is as the task assumed; no STOP.**

### 2. The MR return's structure, confirmed on a live object

Read from the GBSG fit used for Gate P (`R/fs_mr_inference.R:961-1002` for the construction):

```
names(g)                 selected_index selected_label measure log_scale ci_method naive
                         debiased selection_bias fixed_bias selection_rate mean_r mean_r_c
                         complement settings harm_flag n_family n_selected timing_seconds
                         reselection field ij_residual
names(g$debiased)        est lower upper lower_1s se se_ij se_wald var_ij ij_source
                         ij_draws se_ij_two_term se_ij_winner ... lower_1s_wf
names(g$field)           lambda_mean lambda_sd q05..q975 n_out_used n_in_used_mean est2
                         lower_1s lower_2s upper_2s se_field lower_se upper_se R_out R_in
                         seed_offset timing_seconds complement joint joint_s
names(g$field$complement) ... est2 upper_1s lower_1s lower_2s upper_2s se_field lower_se
                         upper_se ... lambda_mean_s se_field_s est2_s upper_1s_s lower_1s_s
                         lower_2s_s upper_2s_s lower_se_s upper_se_s
names(g$field$joint)     gamma joint_prob alpha lower_H upper_Hc bonf_gamma bonf_lower_H
                         bonf_upper_Hc bonf_joint_prob corr n_joint_draws grid_gamma
                         grid_joint_prob
names(g$field$joint_s)   (identical to joint)
names(g$reselection)     winner p_hat
```

**Discrepancy recorded.** The task names the joint members `bonf_loH` / `bonf_upHc`. The object's are **`bonf_lower_H` / `bonf_upper_Hc`** (`R/fs_mr_inference.R:1221-1226`):

```r
  list(gamma = gamma, joint_prob = jp, alpha = alpha,
       lower_H = to_eff(beta_deb - qh_g), upper_Hc = to_eff(bdc - qc_g),
       bonf_gamma = alpha / 2,
       bonf_lower_H = to_eff(beta_deb - qh_b), bonf_upper_Hc = to_eff(bdc - qc_b),
```

`bonf_loH` / `bonf_upHc` are the **simulation template's recorder column names** (`fld_joint_bonf_loH`, `fld_joint_s_bonf_upHc`), not the object's. The code uses the object's names.

`p_hat_sum` is likewise not an element: `reselection$p_hat` is the named frequency vector over the family, and the template derives `p_hat_sum <- sum(ph, na.rm = TRUE)` and `p_hat_H <- ph[[label]]` (template lines 1025-1033). The methods transplant that idiom exactly.

### 3. Vignette infrastructure — the task's assumption was wrong

The task asks to "confirm there is **no** `vignettes/` directory". **There is one**, and it is populated:

```
vignettes/forestsearch.Rmd            (the knitr-built package vignette)
vignettes/articles/*.Rmd              (7 files; ^vignettes/articles$ is .Rbuildignore'd)
```

DESCRIPTION at HEAD already declared:

```
VignetteBuilder: knitr
Suggests: DiagrammeR, doRNG, htmltools, MASS, sandwich, tidyr, forestploter,
    cubature, svglite, knitr, rmarkdown, katex, testthat (>= 3.0.0)
```

`quarto` was **not** a declared dependency, and the `quarto` R package was not installed on this host. A `.qmd` vignette needs the `quarto::html` vignette engine, which the `quarto` R package provides. **What building it adds to DESCRIPTION:** `VignetteBuilder: knitr, quarto` and `quarto` in `Suggests` — nothing else. The `quarto` R package installed here with **no new transitive dependencies** (processx, rmarkdown, jsonlite, yaml, cli, rlang, later, rstudioapi, xfun were all already present) and resolves the CLI at 1.9.38. Both the before- and after-check tallies below are on that basis.

### 4. Scale and orientation

`R/fs_mr_inference.R:578` and the line after it:

```r
  log_scale <- asm$log_scale
  to_eff    <- function(x) if (log_scale) exp(x) else x
```

For survival, `measure = "HR"` and `log_scale = TRUE`, so every bound returned through `to_eff()` is on the **HR scale**, with larger = more harm (`adverse_outcome = TRUE`). The harm block's certified product is a one-sided **lower** bound (`field$lower_1s = to_eff(beta_deb - qs[5])`); the complement's is a one-sided **upper** bound (`field$complement$upper_1s = to_eff(bdc - qs[1])`). The printed block labels each as such and prints the measure name from `g$measure`, so no bound is printed on the wrong scale or in the wrong orientation.

### 5. The two methods and their helpers

`print.forestsearch()` (`R/forestsearch_methods.R:86-153` at HEAD) and `summary.forestsearch()` (`:179-330`) were quoted in full during discovery; both are plain `cat()` sequences guarded by `is.null(x$sg.harm)`, ending in a timing block and `invisible()`. `.fs_get()` (`:19-34`) walks nested name paths and returns the first non-NULL; `.fs_sg_labels()` (`:46-56`) prefers `grp.consistency$out_sg$sg.harm_label` over `sg.harm`; `.fs_format_admission()` lives in `R/forestsearch_helpers.R:2369` and returns NULL for a NULL admission so `.print_param()` skips it. The extension reuses all three plus the `cat()` idiom; no formatting layer was invented.

## Part P — the reporting extension

Four internal helpers were added to `R/forestsearch_methods.R`, all `@noRd`: `.fs_mr_products()` (reads the MR object, returns NULL when MR did not run or the field block did not), `.fs_mr_caveats()` (the NOTE-sourced caveat strings), `.fs_cat_caveat()` (wrapping) and `.fs_print_mr_products()` (the shared block, `long = TRUE` for `summary`). Each method gained a four-line call site before its timing block. Nothing is recomputed: every number is read from `x$mr_inference`.

The complement's reported bound is **field-s** (`upper_1s_s`) whenever it is present, which is the default under `field_scale_complement = "selected"`; the unscaled `upper_1s` is the labelled fallback. The joint pair is taken from `field$joint_s` when field-s is what is being reported, so the pair and the marginal bound come from the same field.

`devtools::document()` regenerated `man/print.forestsearch.Rd` and `man/summary.forestsearch.Rd` only; NAMESPACE unchanged (both already exported). `tools::showNonASCIIfile()` reports the source ASCII clean.

### Gate P

## Pa -- absent-MR invariance

- object: forestsearch(mr_inference = FALSE) on GBSG; is.null(mr_inference) = TRUE
- print():   16 lines before, 16 after; identical() = TRUE
- summary(): 41 lines before, 41 after; identical() = TRUE
- differing lines: print none, summary none
- **Pa: PASS**

## Pb -- every printed number against its MR element

| printed quantity | MR element | element value | printed as | found in output |
|---|---|---|---|---|
| H one-sided 95% lower bound | `mr_inference$field$lower_1s` | 0.6097420234 | 0.610 | TRUE |
| Hc one-sided 95% upper (field-s) | `mr_inference$field$complement$upper_1s_s` | 0.8120160135 | 0.812 | TRUE |
| Joint Bonferroni H lower | `mr_inference$field$joint_s$bonf_lower_H` | 0.500556942 | 0.501 | TRUE |
| Joint Bonferroni Hc upper | `mr_inference$field$joint_s$bonf_upper_Hc` | 0.8477515058 | 0.848 | TRUE |
| Joint gamma each side | `mr_inference$field$joint_s$bonf_gamma` | 0.025 | 0.025 | TRUE |
| IJ two-sided lower | `mr_inference$debiased$lower` | 0.5504239795 | 0.550 | TRUE |
| IJ two-sided upper | `mr_inference$debiased$upper` | 2.900922774 | 2.901 | TRUE |
| se_field (H) | `mr_inference$field$se_field` | 0.4100666339 | 0.410 | TRUE |
| se_field_s (Hc) | `mr_inference$field$complement$se_field_s` | 0.1352705715 | 0.135 | TRUE |
| naive complement SE | `mr_inference$complement$debiased$se_wald` | 0.1343376846 | 0.134 | TRUE |
| se_ij (H, two-term) | `mr_inference$debiased$se_ij` | 0.424011697 | 0.424 | TRUE |
| p-hat(H) | `reselection$p_hat[[selected_label]]` | 0.006 | 0.006 | TRUE |
| p_hat_sum | `sum(reselection$p_hat)` | 1 | 1.000 | TRUE |
| top-1 mass | `sort(p_hat, dec)[1]` | 0.182 | 0.182 | TRUE |
| top-2 mass | `sort(p_hat, dec)[2]` | 0.156 | 0.156 | TRUE |
| top-3 mass | `sort(p_hat, dec)[3]` | 0.056 | 0.056 | TRUE |

- quantities checked: 16; all found at printed precision: TRUE
- top-3 labels present in summary(): TRUE
- **Pb: PASS**

### print(fs) with MR present

```
ForestSearch Results
====================

Selected Subgroup:
  Definition: {er <= 0} & {pgr <= 26} 
  sg_focus: hrMaxSG 
  N: 75 
  HR: 2.222 
  Pcons: 0.99 
  Algorithm: twostage 
  Candidate family: no-front-end 
  Admission set: effect floor 0; consistency floor 0 at p* = 0.8 
  Candidates evaluated: 120 
  Candidates passed: 16 

Post-selection inference (certified products):
  Harm subgroup H:        one-sided 95% lower bound on HR   0.610
  Complement Hc:          one-sided 95% upper bound on HR   0.812   [field-s]
  Joint (Bonferroni):     H lower 0.501, Hc upper 0.848  (gamma 0.025 each side) [field-s]
  Two-sided (IJ, secondary): H (0.550, 2.901)
  Re-selection frequency  p-hat(H) = 0.006

  Two-sided intervals are not certified: the IJ two-term interval's
  harm-block coverage falls to 0.913-0.917 at 12.4% prevalence with n >= 1000
  (0.971-0.981 at 31%).  Read two-sided statements at low prevalence and
  large n with that caveat. 

Computation time: 0.06 minutes
```

### summary(fs), post-selection block only

```
Post-selection inference (certified products):
  Harm subgroup H:        one-sided 95% lower bound on HR   0.610
  Complement Hc:          one-sided 95% upper bound on HR   0.812   [field-s]
  Joint (Bonferroni):     H lower 0.501, Hc upper 0.848  (gamma 0.025 each side) [field-s]
  Two-sided (IJ, secondary): H (0.550, 2.901)

  Standard errors (log scale):
    field (H)                 se_field    0.410
    field-s (Hc)              se_field_s  0.135   (naive complement SE 0.134)
    IJ two-term (H)           se_ij       0.424
  Re-selection frequency  p-hat(H) = 0.006
    top-3 re-selection mass:  q1.1 & q18.1 0.182 | q1.1 & q27.1 0.156 | q10.0 & q30.0 0.056
    p_hat_sum = 1.000 over a family of 1744 candidates

  Certified: the one-sided lower bound on H (field) and the one-sided
  upper bound on Hc (field-s), and the Bonferroni joint pair at gamma =
  0.025 each side.  Not certified: any two-sided interval.  p-hat(H) is a
  recorded diagnostic -- no construction reads it.  Read every bound by
  location against a clinically meaningful effect size, never as
  significance at the null.  Source: dev/notes/NOTE_survival_products_2026-09-09.md.

  Two-sided intervals are not certified: the IJ two-term interval's
  harm-block coverage falls to 0.913-0.917 at 12.4% prevalence with n >= 1000
  (0.971-0.981 at 31%).  Read two-sided statements at low prevalence and
  large n with that caveat. 

Computation time: 0.06 minutes
```

## Pc -- every printed caveat against its NOTE source line

| printed caveat (verbatim from the output) | NOTE source line |
|---|---|
| Two-sided intervals are not certified: ... harm-block coverage falls to 0.913-0.917 at 12.4% prevalence with n >= 1000 (0.971-0.981 at 31%). Read two-sided statements at low prevalence and large n with that caveat. | **Two-sided intervals are not certified.** The IJ two-term two-sided interval is retained and reported as the secondary, conservative option and is the only two-sided construction offered, but its harm-block coverage falls to 0.913–0.917 at 12.4% prevalence with n ≥ 1000 (0.971–0.981 at 31%). The fi |
| p-hat(H) >= 0.5 is the stable-pick regime, where the harm-block correction is under-corrected (+0.02 log units). The flag is directional, not calibrated. | **Analysis-time diagnostic.** p̂(Ĥ), the field's re-selection frequency, is recorded and reported; **no construction reads it**. Harm-block bias is a monotone increasing function of p̂, crossing zero near p̂ ≈ 0.5: **over-correction at low p̂** (bias −0.11 to −0.28 log units at 12.4%, n = 1500) and  |
| Certified: the one-sided lower bound on H (field) and the one-sided upper bound on Hc (field-s), and the Bonferroni joint pair at gamma = 0.025 each side. | **Joint two-subgroup claim:** Bonferroni, γ = 0.025 each side; 0.939–0.963 at every cell. |
| p-hat(H) is a recorded diagnostic -- no construction reads it. | **Analysis-time diagnostic.** p̂(Ĥ), the field's re-selection frequency, is recorded and reported; **no construction reads it**. Harm-block bias is a monotone increasing function of p̂, crossing zero near p̂ ≈ 0.5: **over-correction at low p̂** (bias −0.11 to −0.28 log units at 12.4%, n = 1500) and  |

- every printed claim traced to a NOTE line: TRUE
- p-hat threshold used: 0.5, from the NOTE's "crossing zero near p-hat ~ 0.5"; this analysis has p-hat(H) = 0.006, so the stable-pick line is correctly NOT printed
- **Pc: PASS**

## GATE P: PASS (Pa PASS, Pb PASS, Pc PASS)

## Part V — the vignette

`vignettes/survival-post-selection.qmd`, Quarto, 230 lines. Structure and fitting code transplanted from `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` — the same GBSG preparation (`time_months`, `grade3`, the seven confounders, `conf_force = c("er <= 0", "pgr <= 0")`, `conf.cont_jcuts = list(er = 10, pgr = 10)`) and the same `forestsearch()` call shape — cut to vignette size: `fs.splits` 1000 → 200, `mr_draws` → 500, `parallel_args = list(plan = "sequential")`, `plot.sg = FALSE`, and **no bootstrap, no cross-validation, no Guo & He, no leave-one-out**. The applications document's Pareto-frontier, provenance-reconciliation, CV and LOO blocks are absent.

Sections, in the task's order: the data and the fit; the identified subgroup; the certified products with `print()` shown; how to read a bound **by location** (with the naive / de-biased / field / field-s locations printed side by side against the 0.80 and 0.85 reading aids); p-hat as a diagnostic, both poles, with the 0.5 threshold stated as descriptive not calibrated; what is not certified (the two-sided interval at low prevalence with n >= 1000, plus the epsilon > 0.25 and `maxSG` caveats); and a pointer to the NOTE and the four campaign records.

**Measured build time: 34.55 seconds elapsed, 488 MB peak RSS** (`/usr/bin/time` on `quarto render`, the fit itself 30 s of that). Inside `R CMD check --as-cran` the vignette re-build step reports `[54s/57s] OK`.

DESCRIPTION plumbing, exactly what Stage 0 established as required and nothing else:

```
-VignetteBuilder: knitr
+VignetteBuilder: knitr, quarto
     katex,
+    quarto,
     testthat (>= 3.0.0)
```

The rendered HTML is a build artifact and is not committed; `R CMD build` produces it into `inst/doc`.

## `R CMD check --as-cran`, before and after

The certification surface (`rcmdcheck::rcmdcheck(args = "--as-cran")`, which builds the PDF manual), run on HEAD before any change and again with Part P and Part V in place. Pandoc is not on the Claude shell's PATH on this host; both runs used RStudio's bundled binary (`/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64`), without which `R CMD build` fails at "creating vignettes" before check begins.

| | ERRORS | WARNINGS | NOTES |
|---|---|---|---|
| **before** (HEAD 8b9da196) | 0 | 1 | 2 |
| **after** (this task) | 0 | 1 | 2 |

Item by item, both runs:

- WARNING `checking code files for non-ASCII characters` — `R/fs_bias_coverage.R`. Pre-existing; not touched by this task. `R/forestsearch_methods.R` is ASCII clean.
- NOTE `checking R code for possible problems` — nine `no visible binding for global variable` findings, all in `fs_plot_bias_coverage`. **Byte-identical before and after** once the timing string is stripped (`[41s/41s]` → `[44s/44s]`): the four new helpers add no new findings.
- NOTE `checking HTML version of manual` — "no command 'tidy' found", an artefact of this host, not of the package.

**Nothing new is attributable to this task.** The vignette built cleanly in both `R CMD build` (`* creating vignettes ... OK`) and check (`* checking re-building of vignette outputs ... [54s/57s] OK`), and `* checking installed files from 'inst/doc' ... OK`, so the `quarto::html` engine is correctly wired.

## Test suite

`devtools::test()` after the change: **FAIL 0 | WARN 32 | SKIP 3 | PASS 5051** — the same tally as the last recorded full run (`f4664aed`). No test was added or edited by this task.

## Reading

The GBSG worked example sits at the **low** p-hat pole: `p-hat(Hhat) = 0.006` over a family of 1,744 candidates, with the three largest re-selection frequencies 0.182 / 0.156 / 0.056 on *other* candidates. By the NOTE, that is the over-correction regime (bias −0.11 to −0.28 log units at 12.4% prevalence, n = 1500), so the harm-side lower bound printed here is the conservative one, and `print()` correctly withholds the stable-pick caveat.

Read by location: against a naive HR on Ĥ of 2.222 and a de-biased 1.263, the field one-sided lower bound is **0.610** — the adjusted evidence does not establish clinically meaningful harm in this subgroup, and most of the naive signal was selection. On the complement, the field-s one-sided upper bound is **0.812**, which sits below 0.85 but not below 0.80. The Bonferroni joint pair (0.501, 0.848) is what both statements cost when made together. The IJ two-sided interval (0.550, 2.901) is reported as the secondary, conservative option and is not certified.

## Scope and side issues

Kept: `R/` changes confined to `R/forestsearch_methods.R`; `man/` regenerated for the two methods; DESCRIPTION plumbing only; one new vignette. No change to `fs_mr_inference.R`, `forestsearch_main.R` or any construction. No campaign, no committed simulation bundle. The seven pre-existing untracked files are untouched.

Flagged, not done:

1. **`print()` carries only the two caveats the task enumerates** — the two-sided one and the high-p̂ stable-pick note. The NOTE also describes the **low**-p̂ pole (over-correction, bias −0.11 to −0.28), which is where this worked example actually sits and where a user is therefore quite likely to land. A third line, printed when p̂ < 0.5, would be equally NOTE-sourced. Not added, because the task enumerated two; the low pole is covered in `summary()`'s context and in the vignette's p-hat section. **Recommend deciding whether `print()` should carry it.**
2. **`quarto` is now a vignette-building dependency.** `VignetteBuilder: knitr, quarto` means any machine running `R CMD check` needs the `quarto` R package and the quarto CLI. That is the cost of a `.qmd` vignette and the task directed a `.qmd`; it is recorded here because it is a CRAN-facing change to the package's build requirements, not merely a Suggests addition.
3. **The task's Stage 0 assumptions were wrong in two places** — `vignettes/` exists, and the joint members are `bonf_lower_H` / `bonf_upper_Hc`. Both are recorded above; neither blocked the work.
4. **`.fs_apply_mr()`'s `ci_method` default** remains `"ij"` (`R/fs_mr_inference_methods.R:141`), so a DINA or GRF fit attaches an MR object with no field block and the new section prints nothing for it. Explicitly out of scope here; carried over from `REPORT_cimethod_flip_2026-09-09.md` side issue 1.
