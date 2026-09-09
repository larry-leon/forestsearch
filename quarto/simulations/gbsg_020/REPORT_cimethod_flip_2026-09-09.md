# REPORT — Part D2: the `ci_method` default becomes `"field"`, and the survival-products NOTE

Date: 2026-09-09. Task: `dev/tasks/TASK_cimethod_note_2026-09-09.md` (committed as received, ab9afb20). Executor: Claude Code (Linux). **Report-and-wait.** Classification: **changes behaviour**. Compute: verification renders only (3 × 5 replicates).

Predecessors: `REPORT_defaults_flip_2026-09-08.md` (Part D, which left `ci_method` as the open decision), `REPORT_cert20_2026-09-08.md`, `REPORT_tier2_2026-09-08.md`, `REPORT_fixedphat_ij2s_2026-09-09.md` (C-4), `dev/notes/NOTE_complement_product_2026-09-08.md`.

## Stage 0 — Discovery (quoted from HEAD before the change)

**1. The `ci_method` formal of `fs_mr_inference()`** (`R/fs_mr_inference.R:524`, inside the signature at 505–534):

```r
                           ci_method = c("ij", "wald", "field"),
```

resolved at line 539 by `multiplier <- match.arg(multiplier); ci_method <- match.arg(ci_method)`.

**2. Its roxygen `@param`** (`R/fs_mr_inference.R:301–316`):

```r
#' @param ci_method `"ij"` (default) bases the **de-biased** CI on the
#'   infinitesimal-jackknife variance (Leon et al. 2024, Eq. VInfJ_bc), computed
#'   from the same multiplier draws -- the leading-order analogue of the FB
#'   interval.  `"wald"` uses the subgroup robust SE (`sigma_D`).  The naive CI
#'   always uses the robust SE.  `"field"` computes everything the `"ij"` path
#'   computes -- the `debiased` element is identical -- and additionally runs
#'   the field-calibrated interval (method proposal,
#'   `dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md`), returned as a `field`
#'   element:
```

**3. The gate** (`R/fs_mr_inference.R:818`, opening the block whose header comment runs 797–816):

```r
  field <- NULL
  if (ci_method == "field") {
    t0f <- proc.time()
    if (!is.null(seed)) set.seed(as.integer(seed) + 900000L)
```

**4. The IJ SE is computed before and independent of that gate.** `R/fs_mr_inference.R:671` — unconditional, 147 lines above the gate:

```r
  ijH   <- .fs_mr_ij_var(Xi, r_H, ok_H)
  se_ij <- .fs_mr_se_from_ij(ijH, se_wald)
```

and the reported SE is chosen at `R/fs_mr_inference.R:686`, where only `"wald"` diverts it — `"field"` and `"ij"` take the same branch:

```r
  # "field" keeps the debiased element on the IJ interval (identical to "ij");
  # only "wald" switches it to the robust SE.
  se    <- if (ci_method == "wald") se_wald else se_ij_rep$se
```

The complement's SE follows the same convention at line 762: `sec_used <- if (ci_method %in% c("ij", "field")) se_ijc_rep$se else sec`. **This is the structural fact Gate D2c tests empirically: `"field"` is add-only with respect to the IJ path.**

**5. The `forestsearch()` fallback** (`R/forestsearch_main.R:3406`):

```r
        ci_method     = .g_mr(mr_inference_args$ci_method,   "ij"),
```

**6. The template does not expose a `ci_method` knob — it hard-codes the value.** `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:567`:

```r
mr_ci_method    <- "field"    # "ij" | "wald" | "field"
```

forwarded at line 627 in `mr_inference_args` and recorded in the bundle meta at line 1473. There is no `FS_S7_CI_METHOD` (every other MR knob — `FS_S7_FIELD_COMPLEMENT`, `FS_S7_FIELD_SCALEC`, `FS_S7_FIELD_DECOMP`, `FS_S7_IJ_RESIDUAL`, `FS_S7_RETURN_RESEL`, `FS_S7_UNIFORM` — is env-driven; this one is not). **Per the task, the template is left as it is:** campaigns set `ci_method` explicitly through this hard-coded line, and every committed bundle's `meta$ci_method` reads `"field"` because of it. No knob was added.

**7. `.fs_apply_mr()` (the DINA/GRF branches) falls back independently and does *not* inherit the new default.** `R/fs_mr_inference_methods.R:141`:

```r
      ci_method        = .g(mr_inference_args$ci_method,   "ij"),
```

This is a second, separate `"ij"` literal, in a third file. The task confines the change to two `R/` files, so **it is left untouched and reported here**: after this flip, `forestsearch()`'s consistency engine defaults to `"field"` while the DINA and GRF hooks still default to `"ij"`. That asymmetry is a side issue for a later decision, not fixed here.

## The change

Five files, exactly as scoped: two under `R/`, the two generated `man/` pages, `NEWS.md`, plus the NOTE (Part N2).

```
 NEWS.md                | 25 +++++++++++++++++++++----
 R/forestsearch_main.R  | 29 ++++++++++++++++++++++-------
 R/fs_mr_inference.R    | 30 +++++++++++++++++++-----------
 man/forestsearch.Rd    | 19 +++++++++++++------
 man/fs_mr_inference.Rd | 30 +++++++++++++++++++-----------
 5 files changed, 94 insertions(+), 39 deletions(-)
```

1. `R/fs_mr_inference.R:524` — `ci_method = c("field", "ij", "wald")`, so `match.arg` yields `"field"`. Roxygen `@param` rewritten: `"field"` is the recommended default and runs the field block (the certified one-sided products); the IJ elements are returned in every call regardless; `"ij"` restores the prior default and omits the block; the added per-fit Monte Carlo cost is named. The `@return` line that said "the IJ SE under the default `ci_method = "ij"`" now names `"field"` as the default and states that both take the IJ SE.
2. `R/forestsearch_main.R:3406` — `.g_mr(mr_inference_args$ci_method, "field")`, with an eight-line comment naming the new default, the reason (Part D's flips are inert without it), the fact that nothing is removed, and `"ij"` as the restore value. The `\item{\code{ci_method}}` roxygen block and the `mr_inference` return item updated to match.
3. Template: **left unchanged** (hard-coded, see Stage 0 item 6).
4. `NEWS.md`: one new bullet under the development header. The pre-existing Part D bullet said "the `ci_method` default is unchanged, so the complement field block still runs only under `ci_method = "field"`"; that clause is now false for the same unreleased version, so it was rewritten in the past tense with a pointer to the new bullet. This is the one edit beyond the four items the task enumerates, and it is recorded here rather than made silently.
5. `devtools::document()` regenerated `man/forestsearch.Rd` and `man/fs_mr_inference.Rd` only — `NAMESPACE` unchanged. `devtools::install(dependencies = FALSE)` succeeded.

**Installed-package read-back.**

```
formals(fs_mr_inference)$ci_method  ->  c("field", "ij", "wald")
match.arg resolves to               ->  field
body(forestsearch) contains         ->  ci_method = .g_mr(mr_inference_args$ci_method, "field")
                                        field_complement = .g_mr(mr_inference_args$field_complement, TRUE)
                                        field_scale_complement = .g_mr(mr_inference_args$field_scale_complement, "selected")
                                        return_reselection = .g_mr(mr_inference_args$return_reselection, TRUE)
```

`deparse()` of the installed bodies equals the source: `fs_mr_inference` body `TRUE`, `fs_mr_inference` formals `TRUE`, `forestsearch` body `TRUE`.

## Gate D2

Config: effMaxSG ε 0.20, HR 1.50, n = 500, `z1q` 0.60, `er_jcuts` 10, `seed_base = 8316951` (per-replicate seed `seed_base + sim_id`), `sim_id` 1–5, `FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none`, `mr_draws = 5000`, tags `cim_unset` / `cim_field` / `cim_ij`. Timing columns (`fb_secs`, `fit_mr_secs`, `fld_H_secs`, `fld_Hc_secs`, `fld_H_uniform_secs`) excluded from every comparison.

**Gate drivers.** Because the committed template hard-codes `mr_ci_method` and the task directs that it be left alone, the three arms were driven from three gate-local copies — `gate_d2_cim_unset.qmd`, `gate_d2_cim_field.qmd`, `gate_d2_cim_ij.qmd` — differing from the template in exactly two regions: the `mr_ci_method` assignment, and the `mr_inference_args` list, which passes `ci_method` only when a non-`NULL` `mr_ci_method_pass` is set. The `cim_unset` driver sets `mr_ci_method_pass <- NULL`, so **`ci_method` is never passed to `forestsearch()`** and the package default resolves it; its `mr_ci_method` (meta record only) is read back from the installed formal. Full diffs against the template are the two hunks quoted above and nothing else, verified by `diff`. The drivers and the three 5-replicate bundles are committed beside this record; the `cim_*` campaign tags sit outside every committed combine glob.

### Results

## Meta of the three gate renders

- cim_unset  ci_method=field  n_sims=5 sim 1-5 seed_base=8316951 hr=1.50 n=500 z1q=0.60 nbhd=0.20 focus=effMaxSG field_complement=TRUE field_decompose=TRUE field_scale_complement=selected ij_residual=two_term return_resel(meta absent)=n/a version=0.3.5
- cim_field  ci_method=field  n_sims=5 sim 1-5 seed_base=8316951 hr=1.50 n=500 z1q=0.60 nbhd=0.20 focus=effMaxSG field_complement=TRUE field_decompose=TRUE field_scale_complement=selected ij_residual=two_term return_resel(meta absent)=n/a version=0.3.5
- cim_ij     ci_method=ij     n_sims=5 sim 1-5 seed_base=8316951 hr=1.50 n=500 z1q=0.60 nbhd=0.20 focus=effMaxSG field_complement=TRUE field_decompose=TRUE field_scale_complement=selected ij_residual=two_term return_resel(meta absent)=n/a version=0.3.5

e1stud comparator meta: ci_method=field field_decompose=TRUE ij_residual=two_term field_scale_complement=selected nbhd=0.20

## D2a -- unset identical to explicitly "field"

- non-timing columns compared: 153 of 158 (timing excluded: fb_secs, fit_mr_secs, fld_H_secs, fld_Hc_secs, fld_H_uniform_secs)
- identical(): 153 / 153
- truth identical(): TRUE
- resolved meta$ci_method under unset: "field"
- **D2a: PASS**

## D2b -- the unset render reproduces the committed e1stud rows 1-5

- pre-existing non-timing columns compared: 153 (unset has 158 columns, e1stud 158; new-only: )
- identical(): 153 / 153
- truth identical(): TRUE
- **D2b: PASS**

## D2c -- IJ IS NOT LOST

(i)  IJ columns compared: 23 (mr_ok, mr_H_est, mr_H_lo, mr_H_hi, mr_H_se_ij, mr_Hc_est ...)
     identical("ij" render, unset render): 23 / 23
     every mr_H_/mr_Hc_ column finite on detected rows: unset TRUE, "ij" TRUE
(ii) field / field-s / joint numeric columns compared: 75; all NA in the "ij" render: 75 / 75
(iii) in the UNSET render, on all 5 detected rows: every mr_H_/mr_Hc_ finite = TRUE AND every key field column finite = TRUE (13 / 13)
- **D2c: PASS**

### D2c side-by-side: the IJ two-sided bounds, five replicates

| sim_id | detected | label | mr_H_est unset | mr_H_lo unset | mr_H_hi unset | mr_H_est "ij" | mr_H_lo "ij" | mr_H_hi "ij" | mr_Hc_lo unset | mr_Hc_hi unset | mr_Hc_lo "ij" | mr_Hc_hi "ij" | identical row |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | q8.1 & q27.0 | 0.9281 | 0.4815 | 1.7888 | 0.9281 | 0.4815 | 1.7888 | 0.4201 | 1.1765 | 0.4201 | 1.1765 | TRUE |
| 2 | 1 | q10.0 & q26.0 | 1.4413 | 0.5903 | 3.5193 | 1.4413 | 0.5903 | 3.5193 | 0.6184 | 1.5530 | 0.6184 | 1.5530 | TRUE |
| 3 | 1 | q24.0 & q27.0 | 0.8672 | 0.4318 | 1.7418 | 0.8672 | 0.4318 | 1.7418 | 0.4450 | 1.1924 | 0.4450 | 1.1924 | TRUE |
| 4 | 1 | q7.1 & q12.1 | 1.0255 | 0.5703 | 1.8441 | 1.0255 | 0.5703 | 1.8441 | 0.5045 | 1.4091 | 0.5045 | 1.4091 | TRUE |
| 5 | 1 | q21.0 & q28.0 | 1.4099 | 0.7046 | 2.8212 | 1.4099 | 0.7046 | 2.8212 | 0.4367 | 1.1835 | 0.4367 | 1.1835 | TRUE |

### The field products present only in the unset render

| sim_id | fld_H_est2 | fld_H_lo1s | fld_H_se | fld_Hc_up1s_s | fld_Hc_se_s | fld_joint_bonf_loH | fld_joint_s_bonf_upHc | same columns in the "ij" render |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.8331 | 0.4911 | 0.3117 | 0.9101 | 0.1529 | 0.4402 | 0.9685 | NA |
| 2 | 1.3069 | 0.5784 | 0.4065 | 1.2334 | 0.1220 | 0.4664 | 1.2589 | NA |
| 3 | 0.7647 | 0.4437 | 0.3266 | 0.9338 | 0.1402 | 0.3664 | 0.9616 | NA |
| 4 | 0.9323 | 0.5554 | 0.2891 | 1.0981 | 0.1542 | 0.5012 | 1.1506 | NA |
| 5 | 1.3626 | 0.8496 | 0.2789 | 0.9276 | 0.1373 | 0.7378 | 0.9569 | NA |

## GATE D2: PASS (D2a PASS, D2b PASS, D2c PASS)

## Reading

**D2a — unset ≡ explicitly `"field"`: PASS.** All 153 non-timing columns of 158 `identical()`, `truth` `identical()`, and the resolved `meta$ci_method` under *unset* reads `field`. The default resolves exactly as the explicit setting; there is no residual path difference.

**D2b — the unset render reproduces the committed `e1stud` rows 1–5: PASS.** All 153 pre-existing non-timing columns `identical()` and `truth` `identical()` against `fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_e1stud_combined_1_2000.rds`, which ran `ci_method = "field"` explicitly at package 0.3.5. The gate drivers add no columns (158 in both, no new-only names). So the flip reproduces the certified campaign bit-for-bit on its own configuration — the default now *is* what the campaigns were pinning.

**D2c — IJ IS NOT LOST: PASS on all three limbs.** (i) All 23 IJ and IJ-derived columns (`mr_ok`, `mr_H_est/lo/hi/se_ij`, the `_w` / `_wf` winner variants, `mr_Hc_*`, `ij_source`, `mr_harm_flag`) are `identical()` between the `"ij"` render and the unset (`"field"`) render, and every `^mr_H_` / `^mr_Hc_` column is finite on the detected rows of both. The field gate does not touch the IJ path — as Stage 0 item 4 predicted from the source, now confirmed on data. (ii) All 75 numeric field / field-s / joint columns are `NA` in the `"ij"` render: choosing `"ij"` still omits the block, so prior behaviour is genuinely restored, not merely re-labelled. (iii) In the unset render both families are present **simultaneously** on all 5 detected rows: every `mr_H_` / `mr_Hc_` column finite *and* all 13 key field columns (`fld_H_est2`, `fld_H_lo1s`, `fld_H_se`, `fld_Hc_est2`, `fld_Hc_up1s`, `fld_Hc_se`, the three `_s` companions, `fld_Hc_scale_ratio`, `fld_joint_gamma`, `fld_joint_bonf_loH`, `fld_joint_s_bonf_upHc`) finite.

The side-by-side table above is the concrete form of Larry's criterion: for each of the five replicates the IJ two-sided bounds on Ĥ and on Ĥᶜ are the same numbers under `"ij"` and under the new default, to the last recorded digit and to `identical()` on the underlying doubles — replicate 1's harm interval is 0.4815–1.7888 in both, replicate 5's 0.7046–2.8212 in both — while the unset render additionally carries the field products (replicate 5's field one-sided lower bound on β(Ĥ) is 0.8496 and its Bonferroni joint lower bound 0.7378, both absent under `"ij"`). **The default gains the certified one-sided products and loses nothing.**

Read against the reporting convention: on these five replicates the field one-sided lower bounds on β(Ĥ) run 0.44–0.85 against field point estimates (`fld_H_est2`) of 0.76–1.36, so no clinically meaningful harm is established on any of them; the field-s one-sided upper bounds on β(Ĥᶜ) run 0.91–1.23, none below a 0.80 or 0.85 benefit threshold. Five replicates of a verification render carry no operating characteristics — these are quoted only to show the products are populated and on the expected scale.

## Part N2 — the NOTE

`dev/notes/NOTE_survival_products_2026-09-09.md` written with the Part N2 text exactly as the task document specifies it (extracted programmatically from the committed task file's block quote, so the wording is the approved wording), under a four-line provenance header naming the source task, the superseded note and the evidence reports. A one-line pointer was prepended to `dev/notes/NOTE_complement_product_2026-09-08.md` marking it superseded and linking the new file; its body is unchanged.

The NOTE is committed together with the flip, as one commit — the flip is what makes the NOTE's first sentence ("Defaults are the recommendations") true, so separating them would leave one of the two wrong at every intermediate commit.

## Scope and side issues

Kept: one `R/` change in exactly two files, plus `man/`, `NEWS.md` and the NOTE. No campaign, no bundle beyond the three 5-replicate gate bundles, no change to the print/summary method, the vignette, or the field-recovery diagnostics. The seven pre-existing untracked files are untouched.

Flagged, not fixed:

1. **`.fs_apply_mr()` still defaults `ci_method` to `"ij"`** (`R/fs_mr_inference_methods.R:141`), so the DINA and GRF hooks do not inherit the new default while the consistency engine does. Out of scope here (third file); needs its own decision.
2. **The template has no `ci_method` knob.** Left as directed. Every arm of a future campaign that wants a non-`"field"` value must edit line 567 or use a driver copy, as this gate did.
3. **`NEWS.md`'s Part D bullet was amended**, not only appended to, because its "the `ci_method` default is unchanged" clause described the same unreleased version and is now false.
