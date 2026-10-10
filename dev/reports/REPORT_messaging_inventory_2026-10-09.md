# REPORT — Inventory of user-facing interpretive text (read-only)

- Date: 2026-10-09
- Task: `dev/tasks/TASK_messaging_inventory_2026-10-09.md` (commit `7e5b5a6d`)
- Branch: `feature/glm-extension`
- Type: read-only. Nothing outside this report was edited. No installs, no checks.

## Method

Scope searched: `R/` (code + roxygen), `vignettes/forestsearch.Rmd`,
`vignettes/survival-post-selection.qmd`, `README.md`, `inst/`.
`vignettes/articles/` excluded per the task.

Patterns:

- a. `dev/`, `REPORT_`, `NOTE_`, `TASK_`, `HANDOFF_` (plus `PROPOSAL_` and
  phrases pointing at unshipped records -- "batch documents", "these records" --
  reported as related hits)
- b. `certif*`, `source of truth`
- c. simulation OC figures quoted as text (coverage ranges, prevalence %,
  `n >=`, "log units", detection rates)
- d. `unstable-pick`, `stable-pick`, `exposed one`, `over-correct*`,
  `under-correct*`

"User-facing" was decided as follows. Plain `#` comments were dropped. Roxygen
hits were kept only when the block is not `@noRd`. Every kept roxygen hit was
then confirmed against the generated `man/*.Rd`, so each one listed here really
appears in an installed help page. String hits were traced to the printing
function.

Outcome-path facts used throughout (from `R/forestsearch_main.R` and
`R/fs_mr_inference_methods.R`):

- The forestsearch MR block (blocks R1–R5) prints only when `.fs_mr_products(x)`
  is non-NULL. That needs all of: `mr_inference = TRUE`, a subgroup identified,
  MR eligible, and the **field** block having run without a `note`, which means
  `ci_method = "field"`.
- MR eligibility (`forestsearch_main.R:3704-3708`) on the consistency route is:
  - **survival**: the Cox path (`estimator_fn` NULL). It runs under the default
    split consistency, and `ci_method` defaults to `"field"` (`:3807`).
  - **binary / continuous / count**: the GLM path. It requires
    `consistency_method = "resample"`, and `ci_method` defaults to `"field"`.
  - **DINA (`subgroup_method = "dina"`) and GRF (`subgroup_method = "grf"`)
    routes**: go through `.fs_apply_mr()`, which defaults `ci_method = "ij"`
    (`fs_mr_inference_methods.R:157`). So by default the field block does not
    run and **the block does not print** on these routes. It prints only with
    `mr_inference_args = list(ci_method = "field")`. This holds for any
    outcome.
- The caveat text (R3–R5) is **not outcome-specific**. It quotes survival/HR
  simulation figures and "log units" whenever the block prints, including on
  binary/continuous/count fits and on identity-scale measures (MD/RD). This is
  noted, not fixed.

---

## Part 1 — R/: text that reaches printed output

### R1. MR block header

- **Lines:** `R/forestsearch_methods.R:196`
- **Function:** `.fs_print_mr_products()`
- **Output path:**
  - `print.forestsearch()` → `.fs_print_mr_products(long = FALSE)` (`:375`)
  - `summary.forestsearch()` → `.fs_print_mr_products(long = TRUE)` (`:591`)
- **Categories:** b
- **Outcome paths:** see "Outcome-path facts" (survival default; GLM under
  `consistency_method = "resample"`; DINA/GRF only with `ci_method = "field"`).
- **Conditions:** `mr_inference = TRUE` and the field block ran.

```r
  cat("\nPost-selection inference (certified products):\n")
```

### R2. "Certified / Not certified" paragraph, with dev/ source path (summary only)

- **Lines:** `R/forestsearch_methods.R:230-235`
- **Function:** `.fs_print_mr_products()` (`long = TRUE` branch)
- **Output path:** `summary.forestsearch()` → `.fs_print_mr_products(long = TRUE)`.
  Not printed by `print()`.
- **Categories:** a, b
- **Outcome paths:** as R1.
- **Conditions:** as R1, and `summary()` only.

```r
    cat("\n  Certified: the one-sided lower bound on H (field) and the one-sided\n")
    cat("  upper bound on Hc (field-s), and the Bonferroni joint pair at gamma =\n")
    cat("  0.025 each side.  Not certified: any two-sided interval.  p-hat(H) is a\n")
    cat("  recorded diagnostic -- no construction reads it.  Read every bound by\n")
    cat("  location against a clinically meaningful effect size, never as\n")
    cat("  significance at the null.  Source: dev/notes/NOTE_survival_products_2026-09-09.md.\n")
```

### R3. Two-sided caveat with simulation coverage figures

- **Lines:** `R/forestsearch_methods.R:146-151`
- **Function:** `.fs_mr_caveats()`
- **Output path:**
  `print.forestsearch()` / `summary.forestsearch()` → `.fs_print_mr_products()`
  → `.fs_mr_caveats()` → `.fs_cat_caveat()` (`:237-238`)
- **Categories:** b, c
- **Outcome paths:** as R1. Survival-derived figures are printed on every outcome.
- **Conditions:** as R1, and `is.finite(p$ij_lo)`, which holds whenever the IJ
  interval was computed. In practice that is always the case when the block
  prints.

```r
  if (!is.null(p$ij_lo) && is.finite(p$ij_lo))
    out <- c(out, paste(
      "Two-sided intervals are not certified: the IJ two-term interval's",
      "harm-block coverage falls to 0.913-0.917 at 12.4% prevalence with",
      "n >= 1000 (0.971-0.981 at 31%).  Read two-sided statements at low",
      "prevalence and large n with that caveat."))
```

### R4. Stable-pick caveat (p-hat >= 0.5)

- **Lines:** `R/forestsearch_methods.R:155-159`
- **Function:** `.fs_mr_caveats()`
- **Output path:** as R3.
- **Categories:** c, d
- **Outcome paths:** as R1.
- **Conditions:** as R1, and `p_hat_H >= 0.5`.

```r
  if (is.finite(p$p_hat_H) && p$p_hat_H >= 0.5)
    out <- c(out, paste(
      "p-hat(H) >= 0.5 is the stable-pick regime, where the harm-block",
      "correction is under-corrected (+0.02 log units).  The flag is",
      "directional, not calibrated."))
```

### R5. Unstable-pick caveat (p-hat < 0.20), with dev/ path

- **Lines:** `R/forestsearch_methods.R:166-172`
- **Function:** `.fs_mr_caveats()`
- **Output path:** as R3.
- **Categories:** a, c, d
- **Outcome paths:** as R1.
- **Conditions:** as R1, and `p_hat_H < 0.20`. Nothing prints on [0.20, 0.5).

```r
  if (is.finite(p$p_hat_H) && p$p_hat_H < 0.20)
    out <- c(out, paste(
      "p-hat(H) < 0.20 is the unstable-pick pole, where the harm-block",
      "correction runs the other way -- over-correction, -0.11 to -0.28 log",
      "units at 12.4% prevalence with n = 1500 -- so the one-sided lower bound",
      "on H is conservative here and the two-sided interval is the exposed",
      "one.  See dev/notes/NOTE_survival_products_2026-09-09.md."))
```

### R6. Guo–He Algorithm 3 print: [CERTIFIED] / [EXTENSION] / [HEURISTIC] labels

- **Lines:** `R/guohe_algorithm3.R:515-520` (the match is `:515`; the sibling
  labels are included as one block)
- **Function:** `print.guohe_a3()` (exported S3 method)
- **Output path:** `print()` / autoprint of a `guohe_algorithm3()` result
- **Categories:** b
- **Outcome paths:** `guohe_algorithm3(outcome = c("survival", "binary",
  "continuous"))`. All three print it. Count is not supported by this function.
  This is independent of `forestsearch()`.
- **Conditions:** always, when a `guohe_a3` object is printed.

```r
  cat(sprintf("  %s bound (one-sided, %.3f) : %.4f   [CERTIFIED]\n",
              x$bound_side, 1 - x$level, x$bound_one_sided))
  cat(sprintf("  two-sided CI (%.2f, quantile inversion): (%.4f, %.4f)   [EXTENSION]\n",
              1 - x$level, x$ci_two_sided[1], x$ci_two_sided[2]))
  cat(sprintf("  Wald CI (%.2f, Zhao-style)             : (%.4f, %.4f)   [HEURISTIC]\n",
              1 - x$level, x$ci_wald[1], x$ci_wald[2]))
```

No other `cat`/`message`/`warning`/`stop`/`sprintf` string in `R/` matched a–d.
`.mr_announce()`/`.mr_skip()` and `print.fs_mr_oc()` were checked and are clean.

---

## Part 2 — R/: roxygen text that becomes man/ pages

For each block, the `man/` file and line where the text appears are given in
brackets.

### D1. `print.forestsearch()` — section "Post-selection inference"

- **Lines:** `R/forestsearch_methods.R:252-289`
  - matches at 255, 257, 264, 276, 278, 280
  - 282-286 repeat "the NOTE"
- **Man page:** `man/print.forestsearch.Rd:27-52`
- **Categories:** a, b, d
- **Outcome paths:** documentation; applies to all outcomes.

```
#' @section Post-selection inference:
#' When the fit carries multiplier-resampling results (`mr_inference = TRUE`,
#' which attaches the [fs_mr_inference()] return as `x$mr_inference`) and the
#' field block ran, a further block reports the **certified** survival
#' post-selection products of
#' `dev/notes/NOTE_survival_products_2026-09-09.md`, in that document's order
#' of standing: the one-sided 95% **lower** bound on the harm subgroup
#' \eqn{\beta(\widehat H)} (the field) and the one-sided 95% **upper** bound
#' on the complement \eqn{\beta(\widehat H^c)} (**field-s**, the studentized
#' complement field) first; the Bonferroni **joint** pair at
#' \eqn{\gamma = 0.025} each side second; the IJ two-term **two-sided**
#' interval third, explicitly labelled secondary because no two-sided interval
#' is certified; and the re-selection frequency \eqn{\hat p(\widehat H)}
#' last, as a recorded diagnostic that no construction reads.
#'
#' The field and field-s bounds exist only under `ci_method = "field"`, which
#' is the default; `ci_method = "ij"` omits the field block and with it this
#' entire section apart from the two-sided line.  When MR results are absent
#' the printed output is unchanged from a build without this section.
#'
#' Caveat lines are printed only when they apply and are compressions of the
#' NOTE, not of the analysis: a two-sided caveat whenever a two-sided interval
#' is shown, and one of two directional \eqn{\hat p} notes at the poles of
#' the harm-block bias.  At \eqn{\hat p(\widehat H) \ge 0.5} -- the
#' stable-pick regime -- the correction under-corrects, so the bounds are
#' somewhat optimistic.  At \eqn{\hat p(\widehat H) < 0.20} the correction
#' runs the other way and over-corrects, which makes the one-sided lower bound
#' on \eqn{\widehat H} the conservative product and leaves the two-sided
#' interval as the exposed one.
#'
#' **0.20 and 0.5 are descriptive band edges taken from the NOTE, not
#' calibrated thresholds.** 0.5 is where the NOTE puts the bias crossing zero;
#' 0.20 is the top of the band it characterises as low.  They say which pole an
#' analysis sits at, and the flag is directional, not a decision rule.
#' \eqn{\hat p} in \eqn{[0.20, 0.5)} prints no note at all, because the NOTE
#' supports no claim there.  Bounds are reported by location and must be read
#' against a clinically meaningful effect size, never as significance at the
#' null.
```

### D2. `summary.forestsearch()` — section "Post-selection inference"

- **Lines:** `R/forestsearch_methods.R:396-406` (matches at 397, 404, 405)
- **Man page:** `man/summary.forestsearch.Rd:25-34`
- **Categories:** a, b

```
#' @section Post-selection inference:
#' Carries the same certified-products block as [print.forestsearch()] -- see
#' its documentation for the ordering, the `ci_method = "field"` requirement,
#' the by-location reading convention and the descriptive (not calibrated)
#' \eqn{\hat p} threshold -- and adds, on the log scale, the field SE for the
#' harm block, the field-s SE for the complement with the naive complement SE
#' beside it, and the IJ two-term SE; the three largest re-selection
#' frequencies with their labels; `p_hat_sum` over the re-selection family;
#' and one paragraph naming what is certified and what is not, sourced from
#' `dev/notes/NOTE_survival_products_2026-09-09.md`.  Absent MR results, the
#' output is unchanged from a build without this section.
```

### D3. `fs_mr_inference()` — `@param return_reselection`

- **Lines:** `R/fs_mr_inference.R:290-300` (match at 298)
- **Man page:** `man/fs_mr_inference.Rd:132`
- **Categories:** d

```
#' @param return_reselection Logical. When `TRUE`, the return list gains a
#'   `reselection` element exposing the per-draw re-selection already computed
#'   for the bias term: `winner` (integer vector of length `draws`, the index
#'   of the re-selected candidate in the kept family on each draw, `NA` where
#'   no draw winner exists) and `p_hat` (named numeric over the kept family,
#'   `tabulate(winner) / draws` -- so `sum(p_hat) == selection_rate`, and the
#'   frequencies sum to 1 exactly when every draw produced a winner). Default
#'   `TRUE`: the re-selection diagnostics are part of the recommended output
#'   (they carry the stable-pick frequency the coverage caveat is read by).
#'   Nothing in the arithmetic depends on this switch, so `FALSE` reproduces
#'   the previous, smaller return object exactly.
```

### D4. `fs_mr_inference()` — `@param ci_method`

- **Lines:** `R/fs_mr_inference.R:301-324` (matches at 305 and 316)
- **Man page:** `man/fs_mr_inference.Rd:103`, `:114`
- **Categories:** a, b

```
#' @param ci_method `"field"` (default, the recommendation) computes everything
#'   the `"ij"` path computes -- the `debiased` element, its IJ SE and its IJ
#'   intervals are identical, and are returned in every call regardless of this
#'   argument -- and **additionally** runs the field-calibrated interval, whose
#'   one-sided products are the certified ones (the lower bound on the harm
#'   subgroup, the studentized complement's upper bound, and the Bonferroni
#'   joint pair).  Nothing is removed by choosing it; it adds a per-fit Monte
#'   Carlo cost (`field_R_out` x `field_R_in`, defaults 1000 / 500, plus the
#'   complement's own block under `field_complement = TRUE`).  `"ij"` restores
#'   the previous default: the **de-biased** CI from the infinitesimal-jackknife
#'   variance (Leon et al. 2024, Eq. VInfJ_bc), computed from the same
#'   multiplier draws -- the leading-order analogue of the FB interval -- with
#'   the field block omitted.  `"wald"` uses the subgroup robust SE
#'   (`sigma_D`) for the de-biased CI and also omits the field block.  The
#'   naive CI always uses the robust SE.  The field block (method proposal,
#'   `dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md`) is returned as a `field`
#'   element: Gaussian-multiplier perturbations of the shrunk candidate-effect
#'   field (`w = beta_hat` with the winner's entry replaced by the de-biased
#'   estimate) are pushed through the configured re-selection map, with an
#'   inner Monte Carlo estimating the selection drift, and the resulting
#'   `Lambda*` distribution is inverted around the de-biased estimate
#'   (basic-bootstrap form).  Drawn under `seed + 900000L` when `seed` is
#'   given, after -- never inside -- the main multiplier stream, so `"ij"` and
#'   `"wald"` output is unaffected byte-for-byte.
```

### D5. `fs_mr_inference()` — `@param field_uniform`

- **Lines:** `R/fs_mr_inference.R:328-347` (match at 332)
- **Man page:** `man/fs_mr_inference.Rd:144`
- **Categories:** a

```
#' @param field_uniform Logical (default `FALSE`); only consulted under
#'   `ci_method = "field"`.  When `TRUE`, after the field block completes the
#'   gate additionally runs the uniform (kappa) calibration
#'   (`fs_mr_field_uniform()`, method proposal
#'   `dev/tasks/TASK_mr_field_uniform_2026-09-05.md`) and attaches its result
#'   as `field$uniform`: the smallest widening factor `kappa` in `[1, 2]`
#'   such that the widened two-sided interval attains 95% coverage uniformly
#'   over the **winner-profile protection family** -- hypothetical true
#'   fields equal to the shrunk field `w` with the winner's entry set to the
#'   runner-up level plus `delta` winner-SEs, `delta` on a 0-4 grid --
#'   computed from the trial's own influence structure (`Sigma_hat` via
#'   `B_eff`), restricted to the **mass-carrying candidate set** (smallest
#'   `M <= 12` candidates holding >= 99% of the re-selection mass, winner
#'   always included).  Guarantee as documented there: the one-sided field
#'   bound is uniformly valid; the two-sided interval widened by the
#'   reported `kappa` is uniformly valid over the winner-profile family; the
#'   plain quantile interval is approximate.  The sweep draws under
#'   `seed + 910000L` after -- never inside -- the field's own stream, so
#'   every other output, the field block included, is byte-identical whether
#'   or not it runs.
```

("Guarantee as documented there" points the reader at the unshipped `dev/tasks` file.)

### D6. `fs_mr_inference()` — `@param field_complement`

- **Lines:** `R/fs_mr_inference.R:353-375` (matches at 358, 366, 367)
- **Man page:** `man/fs_mr_inference.Rd:172`, `:180-181`
- **Categories:** a, d

```
#' @param field_complement Logical (default `TRUE`, the recommended
#'   construction); consulted only under
#'   `ci_method = "field"` with `include_complement = TRUE`.  When `TRUE`,
#'   after the harm field block completes, the same field procedure is
#'   evaluated for the complement subgroup (method proposal
#'   `dev/tasks/TASK_mr_field_complement_2026-09-06.md`) and attached as
#'   `field$complement`.  **Orientation:** selection acts on the harm-subgroup
#'   effects; the complement is its by-product, and its naive estimate is
#'   pushed toward *benefit*.  For a benefit claim ("the effect in the
#'   complement is at most U") the exposed limit is therefore the **upper**
#'   bound on the complement's working-scale effect, and the primary product
#'   here is the one-sided 95% upper bound `upper_1s`.  The guarantee is the
#'   mirror of the harm side's: an upper bound is vulnerable to
#'   under-correction of the beneficial optimism (U too low over-claims) and
#'   immune to over-correction.  Selection is the harm field's own -- every
#'   outer and inner draw's winner is re-read, never re-selected -- and the
#'   complement's perturbations use the harm field's own multipliers projected
#'   through the complement's influence, so its noise carries its correct
#'   correlation with every candidate's noise; complements of draw winners the
#'   multiplier stage did not fit are fit lazily.  The block draws nothing, so
#'   every other output -- the harm `field` block and its `uniform` sub-block
#'   included -- is byte-identical whether or not it runs, so `FALSE`
#'   reproduces the pre-complement output exactly.
```

### D7. `fs_mr_inference()` — `@param field_decompose` (related: `PROPOSAL_`)

- **Lines:** `R/fs_mr_inference.R:376-383` (match at 381)
- **Man page:** `man/fs_mr_inference.Rd:196`
- **Categories:** a (related prefix `PROPOSAL_`)

```
#' @param field_decompose Logical (default `FALSE`); consulted only when the
#'   complement field block runs.  When `TRUE`, `field$complement` gains a
#'   `decomp_fields` list of add-only scale diagnostics (the selected
#'   complement's influence-norm scale, the draw-weighted mean and CV of the
#'   outer winners' scales, their ratio `scale_ratio_c`, and the variance /
#'   covariance of the two pieces of `Lambda*c`; PROPOSAL_complement_field_scale
#'   2026-09-08 v2, Stage 1); it reads existing objects only, so every
#'   pre-existing output is unchanged whether or not it runs.
```

### D8. `fs_mr_inference()` — `@param field_scale_complement` (related: `PROPOSAL_`)

- **Lines:** `R/fs_mr_inference.R:384-397` (match at 391)
- **Man page:** `man/fs_mr_inference.Rd:206`
- **Categories:** a (related prefix `PROPOSAL_`)

```
#' @param field_scale_complement `"selected"` (default, the recommended
#'   construction: the studentized complement field, "field-s") or `"none"`;
#'   consulted only when the complement field block runs.  Under `"selected"` the
#'   complement's field readings are studentized to the selected complement's
#'   own scale before differencing -- each outer reading `zeta^c_{r, G_r}` is
#'   multiplied by `s_sel / s_{G_r}` and each inner reading by `s_sel /
#'   s_{G(v_r + zeta'_j)}`, with `s_g = sqrt(sum(Bc[, g]^2))` the candidate's
#'   complement influence-norm scale (the percentile-t root; PROPOSAL_complement_
#'   field_scale_2026-09-08_v2 s5, variant R1) -- and `field$complement` gains
#'   the `_s` companions of its bounds (`est2_s`, `upper_1s_s`, `lower_1s_s`,
#'   `lower_2s_s`, `upper_2s_s`, `lower_se_s`, `upper_se_s`, `se_field_s`,
#'   `lambda_mean_s`) with `field$joint_s` beside `field$joint`.  Add-beside:
#'   no existing field changes under either setting; the default reproduces
#'   prior output exactly.
```

### D9. `fs_mr_inference()` — `@param ij_residual`

- **Lines:** `R/fs_mr_inference.R:398-417` (match at 401)
- **Man page:** `man/fs_mr_inference.Rd:217`
- **Categories:** a

```
#' @param ij_residual Which IJ residual populates the **reported** de-biased
#'   SE and interval (`debiased$se_ij`/`lower`/`upper`/`lower_1s`, and the
#'   complement's), method proposal
#'   `dev/tasks/TASK_complement_refinements_2026-09-06.md`.  `"two_term"`
#'   (default) is the paper's residual
#'   \eqn{r_b = (bias_{sel} + bias_{fix}) - D_{H^*_b}(b) - D_{\hat H}(b)},
#'   whose limit variance is \eqn{(p + e_{\hat H})' \Sigma (p + e_{\hat H})} --
#'   \eqn{4\sigma^2} when one candidate dominates, a doubled SE exactly where
#'   there is least to correct (the complement is effectively always in that
#'   regime).  `"winner"` uses the winner-only (total-derivative) residual
#'   \eqn{r^w_b = bias_{sel} - D_{H^*_b}(b)}, exact when one candidate
#'   dominates but under-covering at ties on the harm side (the delta method
#'   does not hold in the local regime); `"winner_floor"` floors it at the
#'   naive robust variance, \eqn{\max(V^w, \hat\sigma^2_D)}, which repairs the
#'   tie regime at zero cost.  All three SEs are computed from the same draws
#'   and returned side by side in every call (`se_ij_two_term`,
#'   `se_ij_winner`, `se_ij_winner_floor` with their bounds `*_w` / `*_wf`);
#'   the argument only selects the reported one, and the point estimates
#'   never change.  Add-only: the default reproduces prior output exactly
#'   (the new elements aside).
```

### D10. `fs_mr_inference()` — `@param field_recovery`

- **Lines:** `R/fs_mr_inference.R:418-435` (match at 421)
- **Man page:** `man/fs_mr_inference.Rd:238`
- **Categories:** a

```
#' @param field_recovery Logical (default `FALSE`); consulted only under
#'   `ci_method = "field"`.  When `TRUE`, `field` gains a `recovery` list of
#'   **descriptive membership-agreement diagnostics**
#'   (`dev/tasks/TASK_field_recovery_2026-09-09.md`, R-1): over the outer
#'   draws that produced a winner, each re-selected candidate's patient set
#'   `kept[[G_r]]` is compared with the observed pick's `kept[[sel]]`, and the
#'   four classification metrics of the CV membership cross-tab
#'   (`forestsearch_Kfold()`'s `sens_H` / `sens_Hc` / `ppv_H` / `ppv_Hc`) are
#'   averaged over draws, with `npv_Hc` reported as the standing-convention
#'   name of the same quantity as `ppv_Hc`.  These are computed **from draws
#'   already made**: no new fits, no new randomness, no RNG consumption, and
#'   **no construction reads them** -- every bound, bias and SE is
#'   byte-identical whether or not the block runs.  They answer a narrower
#'   question than the full bootstrap (FB) or cross-validation (CV) recovery
#'   diagnostics: re-selection *within the fixed kept family* under
#'   multiplier perturbation, not re-discovery of a subgroup from scratch on
#'   resampled data.  The two are related but not interchangeable, and a
#'   report that quotes one should say which.
```

### D11. `forestsearch()` — `mr_inference_args$ci_method` item

- **Lines:** `R/forestsearch_main.R:989-998` (match at 990)
- **Man page:** `man/forestsearch.Rd:934`
- **Categories:** b

```
#'     \item{\code{ci_method}}{\code{"field"} (default) additionally runs the
#'       field-calibrated block, whose one-sided products are the certified
#'       ones; the de-biased element, its infinitesimal-jackknife SE and its
#'       IJ intervals are returned in every call regardless, so nothing is
#'       removed -- only a per-fit Monte Carlo cost is added.  \code{"ij"}
#'       restores the previous default: the de-biased interval from the
#'       infinitesimal-jackknife variance, the leading-order analogue of the
#'       full-bootstrap interval, with the field block omitted.
#'       \code{"wald"} uses the subgroup robust SE and also omits the field
#'       block.  The naive interval always uses the robust SE.}
```

### D12. `fs_sim_bias_coverage()` — description paragraph

- **Lines:** `R/fs_bias_coverage.R:26-36` (match at 32)
- **Man page:** `man/fs_sim_bias_coverage.Rd:78`
- **Categories:** a

```
#' Conventions follow the committed template's diagnostics section: the
#' one-sided lower bound is the stored `fld_<block>_lo1s` for the field row
#' and `exp(log est - qnorm(level) * SE)` for naive and MR (IJ); two-sided
#' coverage uses each estimator's stored bounds; the field estimator is the
#' shrunk-field `est2` with the two-sided Lambda-quantile interval.  Under
#' `side = "upper"` (the complement's benefit-claim orientation,
#' TASK_mr_field_complement_2026-09-06) the one-sided coverage is
#' `truth <= upper bound`, with the field's stored `fld_<block>_up1s` and
#' `exp(log est + qnorm(level) * SE)` for naive and MR (IJ), and the Gaussian
#' reference is `Phi(z_level * r + b)` -- a positive retained bias helps an
#' upper bound where it hurts a lower one.
```

("the committed template's diagnostics section" also refers to an unshipped record.)

### D13. `fs_sim_bias_coverage()` — `@param estimators`

- **Lines:** `R/fs_bias_coverage.R:44-48` (match at 47)
- **Man page:** `man/fs_sim_bias_coverage.Rd:29`
- **Categories:** a

```
#' @param estimators Subset of `c("naive", "mr", "fld", "mr_w", "mr_wf")`
#'   (default the first three).  `"mr_w"` / `"mr_wf"` are MR (IJ) with the
#'   winner-only / winner-floor SE (`mr_<block>_se_w` / `_se_wf` and their
#'   `_lo_w/_hi_w` / `_lo_wf/_hi_wf` bounds; TASK_complement_refinements_2026-09-06),
#'   available only on bundles that recorded them.
```

### D14. `fs_sim_bias_coverage()` — `@param scale`

- **Lines:** `R/fs_bias_coverage.R:58-67` (match at 64)
- **Man page:** `man/fs_sim_bias_coverage.Rd:50`
- **Categories:** a

```
#' @param scale `"log"` (default) or `"identity"`.  The working scale of the
#'   estimate columns: `"log"` for ratio measures (HR/OR/RR; bias, SD and the
#'   normal-based one-sided bound `exp(log est -/+ z * SE)` are computed on
#'   the log scale, and non-positive estimates are dropped), `"identity"`
#'   for difference measures (MD/RD; bias, SD and the bound `est -/+ z * SE`
#'   on the natural scale, no positivity guard).  Add-only
#'   (TASK_continuous_field_mac_2026-09-07): the default reproduces the
#'   previous output exactly.  Under `"identity"` the `bias_log` column
#'   holds the identity-scale bias (the column name is kept so stacked
#'   tables and [fs_plot_bias_coverage()] read either scale unchanged).
```

### D15. `fs_oc_predict()` — `"resample"` gate item

- **Lines:** `R/fs_oc_predict.R:56-75` (match at 71)
- **Man page:** `man/fs_oc_predict.Rd:133`
- **Categories:** a

```
#'   \item{\code{"resample"} -- the package's production screen}{The
#'     consistency stage of \code{\link{forestsearch}} (default
#'     \code{consistency_method = "resample"}) represents the random 50/50
#'     split's half effects as \eqn{\hat\beta \pm D} and computes the
#'     consistency rate in closed form as
#'     \eqn{2\Phi((\hat\beta - c_2)/\sigma_D) - 1}, a candidate passing when
#'     that rate is at least \code{pconsistency}
#'     (\code{R/consistency_resample.R}; \code{R/subgroup_consistency_helpers.R}
#'     drops a candidate when \code{p.consistency < pconsistency.threshold}).
#'     Inverting: \eqn{\hat\beta \ge c_2 + z_p \sigma_D} with
#'     \eqn{z_p = \Phi^{-1}((1+p)/2)}.  \eqn{\sigma_D} is
#'     \code{sqrt(sum(dfbeta[i, treat]^2))}, the sandwich standard error of the
#'     subgroup treatment coefficient; in the draw model above the analogous
#'     quantity is \code{se_g}, and on a simulated MD40 trial the two agree
#'     within a few percent with no prevalence trend
#'     (\code{dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R}), so the
#'     wrapper identifies \eqn{\sigma_D = se_g}.  The gate is then a single
#'     threshold on the full-sample draw,
#'     \code{(Bhat >= c1) & (Bhat - c2 >= z_p * se_g)}, and only
#'     \code{Bhat ~ N(beta_g, Sg)} is drawn -- one matrix, not two.}
```

### D16. `fs_mr_oc_summary()` — sections citing "batch documents" and a quoted detection/coverage figure

- **Lines:** `R/fs_mr_oc_summary.R:18-23`, `:25-31`, `:48-59`
  - "analytic documents" at 22
  - "batch documents" at 30 and 53
  - figures at 58-59
- **Man page:** `man/fs_mr_oc_summary.Rd:85-86` for the figures
- **Categories:** a (related), c

```
#' Summarises a simulation payload produced by the multiplier-resampling (MR)
#' harness into the estimation, coverage and classification quantities used by
#' the analytic documents, returning plain numbers rather than a formatted
#' table.
...
#' \code{orient = if (adverse_outcome) 1 else -1}, applied here exactly as the
#' batch documents apply it.  Getting this backwards silently negates every
#' bias.
...
#' @section Two classification conventions:
#' Both are returned, because the project contains both and they disagree
#' whenever detection is not near certain:
#' \describe{
#'   \item{\code{conditional}}{Averaged over DETECTED replicates only. This is
#'     what the batch documents report.}
#'   \item{\code{unconditional}}{Averaged over ALL replicates, following
#'     \code{build_classification_table()} and Leon et al. (2024, Table 1),
#'     with non-detection scored as sens = 0, ppv = 0, spec = 1, npv = 1.}
#' }
#' With detection at 999/1000 the two agree to three decimals; in a null cell
#' at 0.930 they do not. Quote which one you mean.
```

### D17. `guohe_algorithm3()` `@return` and `print.guohe_a3()` description

- **Lines:**
  - `R/guohe_algorithm3.R:191-200` (match at 194)
  - `R/guohe_algorithm3.R:497-501` (match at 501)
- **Man pages:** `man/guohe_algorithm3.Rd:89`, `man/print.guohe_a3.Rd:20`
- **Categories:** b
- **Note:** this is the companion documentation to R6.

```
#' @return An object of class `"guohe_a3"`: a list with the selected candidate
#'   (`selected`), the naive and de-biased effects on the reporting scale
#'   (`naive`, `debiased`), the certified one-sided bound (`bound_one_sided`,
#'   with its conservative side in `bound_side`), the two-sided
#'   quantile-inversion interval (`ci_two_sided`), the heuristic Wald interval
#'   (`ci_wald`), and bookkeeping fields (...)
...
#' Print method for guohe_a3 objects
#'
#' Prints the family size, the selected candidate, the naive and de-biased
#' effects, and the three interval objects returned by [guohe_algorithm3()],
#' labelled by their inferential status (certified / extension / heuristic).
```

### D18. `fs_betaHhat_one()` — description

- **Lines:** `R/betaHhat_truth.R:208-213` (match at 212)
- **Man page:** `man/fs_betaHhat_one.Rd:53`
- **Categories:** b. This is the paper's Theorem 2 in its mathematical sense; it
  is listed for completeness.

```
#' Conditional estimand `beta(Hhat)` for one realized rule
#'
#' Computes the population value of the standard within-subgroup analysis on a
#' realized selection region and its complement, on a fixed evaluation frame.
#' This is the estimand Theorem 2 certifies coverage for, evaluated at the
#' rule a detector actually returned -- not the marginal true-subgroup target.
```

### D19. `.consistency_glm_pieces()` — internal man page (`@keywords internal`, not `@noRd`)

- **Lines:** `R/consistency_resample.R:217-222` (match at 222)
- **Man page:** `man/dot-consistency_glm_pieces.Rd:70` (installed, hidden from the index)
- **Categories:** a

```
#' Earlier versions of this documentation claimed that propensity-adjusted
#' effects "return `NULL` (the caller falls back to splitting)". That was never
#' true and was not implementable as written. Whether to extend the influence
#' path, add the argument purely so the claim can be honoured, or refuse the
#' combination outright is unresolved; see
#' `dev/identifier-alignment/f4_ps_effect_gap_findings.qmd`.
```

### D20. `.fs_mr_reselection_from_focus()` — internal man page

- **Lines:** `R/fs_mr_inference_methods.R:71-79` (match at 79)
- **Man page:** `man/dot-fs_mr_reselection_from_focus.Rd:22` (internal)
- **Categories:** b. This is the engineering sense of "single source of truth";
  it is listed for completeness.

```
#' Map `sg_focus` to MR's re-selection rule, faithfully per engine.
#'
#' MR must re-select under the same rule the search used.  The only rule
#' that differs by engine is the `"hr"`/`"eff"` focus: the consistency search
#' ranks by consistency rate, whose MR analog is `"maxcons"`; DINA/GRF rank by
#' effect, whose analog is `"maxeff"`.  Every other focus maps identically
#' (`maxSG`/`minSG` pass through; `hrMaxSG`/`hrMinSG` -> the size-within-effect-
#' neighborhood rules `effMaxSG`/`effMinSG`).  This makes `sg_focus` the single
#' source of truth, so callers never need to specify `reselection` by hand.
```

### Excluded R/ hits (not user-facing)

- `@noRd` blocks:
  - `forestsearch_methods.R:60-81` (`.fs_mr_products`), `:126-139` (`.fs_mr_caveats`), `:186-192` (`.fs_print_mr_products`)
  - `consistency_resample.R:125` (`.dfbeta_glm`)
  - `fs_mr_inference.R:1151`, `:1239`, `:1430`
  - `forestsearch_helpers.R:895`, `:2195`
  - `bootstrap_calculations_helpers.R:8`
- Plain comments:
  - `forestsearch_methods.R:116-118`, `:142-165`, `:370-373`
  - `guohe_algorithm3.R:436`
  - `fs_mr_inference.R:895-896` and similar
- `inst/patches/2026-05-21_extract_grf_args.md:18` ("single source of truth"):
  `^inst/patches$` is in `.Rbuildignore`, so it is not shipped.
- `README.md`: no hits.
- Example data in `calibrate_L_eff` (`sim_fpr = c(0.185, ...)`): these are
  illustrative inputs, not OC claims, so they were not counted.

---

## Part 3 — `vignettes/survival-post-selection.qmd`

This file is not in `.Rbuildignore`, so it ships as a package vignette.

### V0. Rendered output of R1–R5

- **Lines:** chunks `print` (`:107-109`) and `summary` (`:134-136`)
- The vignette's rendered HTML contains R1, R3 and R5 from `print(fs)`. Per
  `:203`, this fit has p-hat = 0.006, so R5 prints rather than R4. It also
  contains R1, R2, R3 and R5 from `summary(fs)`.
- So the `dev/notes/...` path appears twice and the coverage figures appear
  twice in the built vignette's output, in addition to the prose below.

### V1. Title and index entry

- **Lines:** `:2`, `:13`
- **Categories:** b

```
title: "Survival post-selection inference: the certified products"
  %\VignetteIndexEntry{Survival post-selection inference: the certified products}
```

### V2. Introduction, paragraph 1

- **Lines:** `:26-30` (match at 30)
- **Section:** (preamble)
- **Categories:** b

```
`forestsearch()` identifies a subgroup by searching; anything estimated on that
subgroup afterwards is estimated on a subgroup the data chose. This vignette
runs the pipeline end to end on a public dataset and shows what the package
reports about that selection, and -- as importantly -- what it declines to
certify.
```

### V3. Introduction, paragraph 2

- **Lines:** `:32-36` (matches at 33, 34)
- **Section:** (preamble)
- **Categories:** a, b

```
Every claim about coverage below comes from
`dev/notes/NOTE_survival_products_2026-09-09.md`, which is the source of truth
for the certified products; the simulation evidence behind it is named at the
end. Nothing here is new methodology, and no number in this document is
computed by it: the reporting block reads elements the fit already carries.
```

### V4. Section "The certified products": heading and three coverage paragraphs

- **Lines:** `:103`, `:113-115`, `:117-121`, `:123-125`
- **Section:** "The certified products"
- **Categories:** b, c

```
## The certified products
...
**The harm subgroup's one-sided lower bound** on `HR(Hhat)` is the *field*
construction. Its coverage is 0.944--0.974 across ten harm cells at 12.4% and
31% prevalence, n = 500 / 1000 / 1500, and it is flat in n.

**The complement's one-sided upper bound** on `HR(Hhat^c)` is *field-s*, the
studentized complement field. Its small-sample shortfall closes with n:
0.912 -> 0.942 -> 0.947 at 31% prevalence with HR 1.50, and 0.919 -> 0.942 ->
0.946 at HR 1.75; it is at or above 0.94 at every cell of the 12.4% regime.
field-s is at or above the unscaled field in all thirteen certification cells.

**The joint two-subgroup claim** is the Bonferroni pair at gamma = 0.025 each
side, covering 0.939--0.963 at every cell. Use it when both statements are
made together.
```

### V5. Lead-in to `summary()`

- **Lines:** `:131-132`
- **Section:** "The certified products"
- **Categories:** b

```
`summary()` adds the standard errors, the re-selection mass and the
certified/not-certified paragraph:
```

### V6. Reading aids "in these records"

- **Lines:** `:148-152` (match at 150)
- **Section:** "How to read a bound: by location"
- **Categories:** a (related: refers to unshipped records)

```
The complement-side question is a benefit claim: *is the one-sided upper bound
below a threshold that would matter?* Name the threshold you are reading
against -- 0.80 and 0.85 are the usual reading aids in these records -- and say
whether the bound clears it. A threshold quoted this way is a reading aid, not
a decision rule.
```

### V7. p-hat bias bullets

- **Lines:** `:184-195` (matches at 188-192)
- **Section:** "p-hat as a diagnostic"
- **Categories:** c, d

```
Harm-block bias is a monotone increasing function of p-hat, crossing zero near
p-hat ~ 0.5:

- **Low p-hat, below 0.20** -- the pick is unstable and the correction runs the
  other way, *over*-correcting (bias -0.11 to -0.28 log units at 12.4%
  prevalence, n = 1500). The one-sided lower bound on Hhat is the conservative
  product here; the two-sided interval is the exposed one. `print()` says so.
- **High p-hat, 0.5 and above** -- the stable-pick regime, where the correction
  *under*-corrects (+0.02 log units) and the bounds are somewhat optimistic.
  `print()` says so.
- **In between, 0.20 to 0.5** -- no line is printed, because the NOTE supports
  no claim there.
```

### V8. Band edges "from the NOTE"

- **Lines:** `:197-201`
- **Section:** "p-hat as a diagnostic"
- **Categories:** a (refers to the NOTE)

```
Those two edges come from the NOTE: 0.5 is where it puts the bias crossing
zero, and 0.20 is the top of the band it characterises as low. They are
**descriptive band edges, not calibrated thresholds** -- they tell you which
pole an analysis sits at, and the flag is directional, not a decision rule.
The complement's coverage, by contrast, is a stable function of p-hat.
```

### V9. This fit's pole

- **Lines:** `:203-205` (match at 204)
- **Section:** "p-hat as a diagnostic"
- **Categories:** d

```
This fit sits at the low pole (p-hat = 0.006), so the printed output carries
the unstable-pick line and not the stable-pick one, and the harm-side lower
bound should be read as the conservative product.
```

### V10. Section "What is not certified": two-sided paragraph

- **Lines:** `:249-256` (matches at 249, 251, 253-255)
- **Section:** "What is not certified"
- **Categories:** b, c

```
## What is not certified

**No two-sided interval is certified.** The IJ two-term two-sided interval is
retained and reported as the secondary, conservative option, and it is the only
two-sided construction offered -- but its harm-block coverage falls to
0.913--0.917 at 12.4% prevalence once n is 1000 or more (against 0.971--0.981
at 31%). The field's own two-sided interval is lower still, 0.878--0.930. Read
two-sided statements at low prevalence and large n with that caveat.
```

### V11. "Also on record"

- **Lines:** `:258-260`
- **Section:** "What is not certified"
- **Categories:** c (0.876 of the error SD); "on record" refers to unshipped records

```
Also on record: `effect_neighborhood > 0.25` is not adoptable; under a pure-size
pick (`sg_focus = "maxSG"`) the naive SE mis-calibrates, at 0.876 of the error
SD, and field-s inherits that.
```

### V12. One-sided vs two-sided rationale

- **Lines:** `:262-265`
- **Section:** "What is not certified"
- **Categories:** b

```
The one-sided products certify and the two-sided one does not for a reason
worth stating: each one-sided bound is exposed to one pole of the p-hat bias
only, and in the conservative direction. A two-sided interval is exposed to
both.
```

### V13. Section "Where the evidence lives"

- **Lines:** `:267-276` (matches at 269-275)
- **Section:** "Where the evidence lives"
- **Categories:** a, b, c

```
## Where the evidence lives

- `dev/notes/NOTE_survival_products_2026-09-09.md` -- the certified products
  and their caveats; the source of truth for every claim above.
- `REPORT_cert20_2026-09-08.md` -- the 31%-prevalence sample-size profile.
- `REPORT_tier2_2026-09-08.md` -- the 12.4% dominated regime.
- `REPORT_fixedphat_ij2s_2026-09-09.md` -- the fixed-p-hat diagnostics, the
  harm-block two-sided decay and detection conditioning.
- `REPORT_cimethod_flip_2026-09-09.md` -- why `ci_method = "field"` is the
  default.
```

The task named the arXiv/JRSS-B methods paper as the intended citation. It is
not cited anywhere in this vignette.

---

## Part 4 — `vignettes/forestsearch.Rmd`: known issues

No a–d hits in this file. The four known issues, verbatim:

1. **GBSG treatment described as "versus chemotherapy"**
   - **Lines:** `vignettes/forestsearch.Rmd:115`, `:121`
   - **Section:** "Data: German Breast Cancer Study Group Trial" → "Study Background"

   ```
   The GBSG trial evaluated hormonal treatment (tamoxifen) versus chemotherapy 
   - **Treatment**: Hormonal therapy (tamoxifen) vs. chemotherapy
   ```

2. **"included in forestsearch package" comment**
   - **Lines:** `vignettes/forestsearch.Rmd:126`
   - **Section:** "Data Preparation", chunk `data-setup`

   ```
   # Load GBSG data (included in forestsearch package)
   ```

3. **Key Findings line using `exp()` of the `hr` column**
   - **Lines:** `vignettes/forestsearch.Rmd:864`
   - **Section:** "Summary and Interpretation" → "Key Findings", chunk `summary-findings`

   ```
   cat("  Unadjusted HR:", sprintf("%.2f", exp(fs$grp.consistency$out_sg$result$hr[1])), "\n")
   ```

4. **"95%95%"**
   - **Lines:** `vignettes/forestsearch.Rmd:792`
   - **Section:** "KM Difference plots: ITT and subgroups"

   ```
   The solid black line denotes the ITT Kaplan-Meier treatment difference estimates along with 95%95% CIs (the grey shaded region). K-M differences corresponding to subgroups are displayed.
   ```

---

## Count summary

| Surface | Distinct blocks | Blocks touching each category |
|---|---|---|
| R/ printed output (R1–R6) | 6 | a: 2 (R2, R5); b: 4 (R1, R2, R3, R6); c: 3 (R3, R4, R5); d: 2 (R4, R5) |
| R/ roxygen → man/ (D1–D20) | 20 | a: 15 (D1, D2, D4–D10, D12–D16, D19); b: 7 (D1, D2, D4, D11, D17, D18, D20); c: 1 (D16); d: 3 (D1, D3, D6) |
| survival-post-selection.qmd (V0–V13) | 14 (incl. V0 rendered output) | a: 4 (V3, V6, V8, V13); b: 8 (V1–V5, V10, V12, V13); c: 5 (V4, V7, V10, V11, V13); d: 2 (V7, V9) |
| forestsearch.Rmd known issues | 4 issues (5 lines) | — |
| README.md, inst/ (shipped) | 0 | — |

Totals:

- **40 distinct user-facing blocks** (6 printed + 20 man + 14 vignette) plus the
  4 known Rmd issues.
- References to unshipped files:
  - `dev/notes/NOTE_survival_products_2026-09-09.md` appears in printed output
    (R2, R5), in 2 man pages (D1, D2), and in the survival vignette (V3, V13).
  - Other `dev/`/`TASK_`/`PROPOSAL_` paths appear in 5 more man pages
    (`fs_mr_inference`, `fs_sim_bias_coverage`, `fs_oc_predict`,
    `dot-consistency_glm_pieces`, plus `forestsearch` via "certified").
  - Four `REPORT_*` files are named in V13.

## Side observations (not in scope; not acted on)

1. The R3–R5 caveats quote survival/HR simulation figures and "log units", and
   they print unchanged on binary/continuous/count MR fits. That includes
   identity-scale MD/RD, where "log units" is wrong.
2. On the DINA/GRF routes `.fs_apply_mr()` defaults `ci_method = "ij"`, while
   the consistency route defaults to `"field"`. So by default the
   certified-products block appears on one route and not the others. This is
   documented as deliberate at `fs_mr_inference_methods.R:137-139`.
3. `survival-post-selection.qmd` hard-codes values from this particular fit in
   prose:
   - p-hat = 0.006 (`:203`, `:211`, `:232`, `:234`)
   - "1,744 candidates" (`:209`)
   - "more than half" / "about a quarter" (`:233-234`)

   If the fit output drifts, the prose will disagree with the rendered chunks.
4. `forestsearch.Rmd:858` hard-codes the ITT result
   `cat("  HR (95% CI): 0.69 (0.54, 0.89)\n\n")` rather than computing it.
