# REPORT — `dinamr` Stage 0: discovery, and a hard STOP at item 3

**Date:** 2026-09-10. **Executor:** Claude Code (Linux, pop-os). **Task:** `dev/tasks/TASK_dina_mr_2026-09-10.md`.
**Outcome: STOP at Stage 0 item 3.** Stage 1 was not run. Stage 2 and Stage 3 were not run. No `R/` file was changed.
**Decision requested: O-1, which is Larry's.**

---

## 0. Framing — what this campaign would and would not have established

This campaign **does not certify DINA and must not be described as doing so.** The manuscript's
fixed-family condition (Section 2.1) requires a candidate family that depends on no outcome, no fitted
effect and no learned surface. DINA's candidates are read off a cross-fit surface that a genuine
bootstrap would regenerate on every resample. After the alignment repair — re-ranking qualifiers on
the inferential coefficient beta-hat(g) — the de-biased estimate and its interval are first-order
exact **conditional on the proposed family**, with the omitted family-generation component standing as
the gap to the unconditional target. The estimand of every number the campaign would produce is
therefore the **conditional-on-proposed-family estimand**, and every coverage figure would have to be
labelled as coverage of *that*, beside each table, not of the unconditional target.

**What would have been new.** The existing DINA evidence (Supplementary Table S5) predates the current
product set: it carries no field-s complement bound, no Bonferroni joint pair, no p-hat, and no
recovery diagnostics. None of those has ever been measured on a model-generated family. That gap is
real and this report does not close it.

**Premise correction, for the record.** DINA is **not** `maxeff`-only. Per the per-engine `sg_focus`
resolution table in `R/forestsearch_main.R` (quoted verbatim at section 2 below), for DINA and GRF
`maxeff` and `maxeffCons` are **synonyms** — neither engine has a consistency floor, neither computes
`Pcons`, and no DINA or GRF sort key contains one — while `effMaxSG`, `effMinSG`, `maxSG` and `minSG`
remain available and remain distinct rules. On DINA the set {`eff`, `hr`, `maxcons`, `maxeff`,
`maxeffCons`} is one rule under five spellings: `order(-eff)`.

One refinement to that correction, measured rather than quoted (section 2.3): the synonymy is about the
*sort key*, not about the *admission set*. On DINA, `maxeff` carries the **effect floor** and no
consistency floor; on the consistency engine, `maxeff` carries **neither** floor. So `maxeff` does not
mean the same admission set across the two engines, and the FS reference runs (which used `maxeffCons`
and `effMaxSG`) differ from a DINA `maxeff` run in admission as well as in identifier and family
construction.

---

## 1. Item 1 — Does a committed document drive `subgroup_method = "dina"` under the M1 DGM?

**Yes, literally — and no, in substance.** The literal STOP in item 1 ("if none exists, STOP, do not
author one") is therefore **not** triggered; the STOP in this report comes from item 3. But the
documents that exist cannot carry this campaign, and that is reported here rather than at item 3
because it changes what a transplant would have to be.

### 1.1 What exists

Seven committed drivers under `quarto/simulations/gbsg_020/` set `subgroup_method <- "dina"` on the M1
DGM, all at HR 1.00 (null), n = 500, `knoise0`:

| Document | `sg_focus` | lines |
|---|---|---|
| `sim_dina_maxeff_mr_m1_h10_knoise0_n500_batch_1_500.qmd` | `maxeff` | 1500 |
| `sim_dina_maxeffCons_mr_m1_h10_knoise0_n500_batch_1_500.qmd` | `maxeffCons` | 1502 |
| `sim_dina_eff_mr_m1_h10_knoise0_n500_batch_1_500.qmd` | `eff` | 1500 |
| `sim_dina_eff_fb_mr_m1_h10_knoise0_n500_batch_1_200.qmd` | `eff` | 1504 |
| `sim_dina_eff_fb_mr_m1_h10_knoise0_n500_batch_201_250.qmd` | `eff` | 1504 |
| `sim_dina_eff_fb_mr_m1_h10_knoise0_n500_batch_251_500.qmd` | `eff` | 1504 |
| `sim_dina_eff_fb_mr_m1_h10_knoise0_n500_combine_1_500.qmd` | `eff` | 1504 |

Their knobs, quoted from `sim_dina_maxeff_mr_m1_h10_knoise0_n500_batch_1_500.qmd`:

```r
subgroup_method <- "dina"  # "consistency" | "dina" | "grf"
sg_focus           <- "maxeff"
target_hr_harm  <- 1.0    # calibrate k_inter to this Cox HR in the harm subgroup
n_sample        <- 500L   # M1 sample size
n_sims       <- 500L
nb_boots     <- 0L
k_random_noise <- 0
mr_draws   <- 5000L
quickrun     <- TRUE
dgm_model       <- "alt"
analysis_time   <- 84
cens_adjust     <- log(1.5)
n_super         <- 100000L
dina_args       <- list()
dina_select_statistic <- "effect"
hr_threshold       <- 0.90
hr_consistency     <- 0.80
pconsistency       <- 0.90
seed_base  <- 8316951L
confounders_base <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")
```

### 1.2 Why they cannot carry this campaign

**They predate the entire field apparatus.** Measured, not inferred:

| | DINA drivers | field template |
|---|---|---|
| occurrences of `FS_S7_` | **0** | **42** |
| recorder columns in the saved bundle | **57** | **158** |
| `fld_*` columns | 0 | 80 |
| `p_hat_*` columns | 0 | 6 |
| `n_family` column | absent | present |

The committed DINA bundles under `results/` carry the pre-field 57-column recorder
(`forestsearch_version` 0.2.0 / 0.2.1). There is no `ci_method` knob, no field block, no p-hat, no
joint pair, no `n_family`. Running the `dinamr` cells from these documents would produce none of the
products the task is about.

For completeness, the committed null-cell DINA evidence (`dina_maxeff_mr_m1_h10_knoise0_n500_quickrun_res_1_500.rds`,
500 replicates, HR 1.00, n 500): **350/500 = 70.0% detection**, 302 distinct selected subgroups among
the 350 detections, 49 distinct selected sizes. The identifier is behaving normally and is not
degenerate. (An earlier pass of this analysis reported a single distinct subgroup; that was an indexing
error on my side — `detected` is a 0/1 **integer**, so `r[r$detected, ]` indexes by position, not
logically. The corrected figures are the ones above.)

---

## 2. Item 2 — The per-engine `sg_focus` table, and the DINA branch's family

### 2.1 The resolution table, verbatim (`R/forestsearch_main.R:559-586`)

```
#'   \strong{Per-engine resolution.}  Every accepted spelling, and the rule
#'   that actually runs for it on each engine:
#'
#'   \tabular{lll}{
#'     \strong{sg_focus} \tab \strong{consistency} \tab \strong{dina / grf} \cr
#'     \code{eff}, \code{hr}, \code{maxcons} \tab consistency argmax
#'       (\code{Pcons}-primary, effect breaks ties) \tab effect argmax,
#'       \code{order(-eff)} \cr
#'     \code{maxeff} \tab effect argmax, no floors \tab effect argmax,
#'       \code{order(-eff)} \cr
#'     \code{maxeffCons} \tab effect argmax with the consistency floor \tab
#'       effect argmax, \code{order(-eff)} \cr
#'     \code{effMaxSG}, \code{hrMaxSG} \tab largest N within the effect band,
#'       per \code{selection_rule} \tab as documented in
#'       \code{\link{dina_subgroup}} \cr
#'     \code{effMinSG}, \code{hrMinSG} \tab smallest N within the effect band,
#'       per \code{selection_rule} \tab as documented in
#'       \code{\link{dina_subgroup}} \cr
#'     \code{maxSG} \tab largest N among qualifiers \tab largest N \cr
#'     \code{minSG} \tab smallest N among qualifiers \tab smallest N \cr
#'   }
#'
#'   The collapse sets, stated explicitly.  On \code{"consistency"},
#'   \{\code{eff}, \code{hr}, \code{maxcons}\} is one rule, and \code{maxeff}
#'   and \code{maxeffCons} are two further \emph{distinct} rules.  On
#'   \code{"dina"} and \code{"grf"}, \{\code{eff}, \code{hr}, \code{maxcons},
#'   \code{maxeff}, \code{maxeffCons}\} is \strong{one} rule -- five spellings,
#'   one \code{order(-eff)}.
```

`fs_focus_tag()` confirms the collapse at the stem level (measured):

```
  dina/maxeff     -> eff        consistency/maxeff     -> maxeff
  dina/maxeffCons -> eff        consistency/maxeffCons -> maxeffCons
  dina/eff        -> eff        consistency/eff        -> maxcons
  dina/hr         -> eff        consistency/hr         -> maxcons
  dina/maxcons    -> eff        consistency/maxcons    -> maxcons
  dina/effMaxSG   -> effMaxSG   consistency/effMaxSG   -> effMaxSG
  dina/maxSG      -> maxSG      consistency/maxSG      -> maxSG
  dina/minSG      -> minSG      consistency/minSG      -> minSG
```

A `dinamr` run at `sg_focus = "maxeff"` would therefore write under the stem token `dina_eff_...`,
**not** `dina_maxeff_...`.

### 2.2 The DINA branch's `.fs_mr_family_from_table()` call, verbatim (`R/forestsearch_main.R:2327-2349`)

```r
      # MR's family is DINA's qualifying candidates -- the whole of it.  This
      # is what Algorithm Step 2 specifies, and what the consistency engine
      # already does.  A former native-statistic band (.fs_mr_restrict_native)
      # narrowed it so MR's re-selection neighbourhood would mirror the one
      # the full bootstrap explores; MR is not required to mimic FB, and the
      # band's premise -- DINA ranking on tau-hat while MR ranks on perturbed
      # beta-hat -- no longer holds now that select_statistic = "effect" is
      # the default.
      .mr_fam  <- .fs_mr_family_from_table(.mr_df, dsel$candidates,
                                           op_right = ">=", n_min = n.min)
      out$mr_inference <- .fs_apply_mr(
        df = .mr_df, candidates = .mr_fam,
        selected_members = which(dsel$grp.consistency$sg.harm.id == 1L),
        spec = .mr_spec,
        # Admission resolved once, not rebuilt here.  DINA and GRF have no
        # consistency screen, so no consistency term may enter their admission
        # set -- previously `c_consistency = 0` with `p_star` still set meant
        # t_g carried a z * sigma_D term these engines never applied.
        admission = admission_resolved,
        effect_neighborhood = effect_neighborhood,
        reselection_default = .fs_mr_reselection_from_focus(sg_focus, engine = "effect"),
        selection_rule_default = selection_rule,
        mr_inference_args = mr_inference_args, seedit = seedit)
```

**What family MR receives, and that it is frozen before MR runs.** `dsel` is the completed DINA
selection: the cross-fit fit, the qualifier table and the selected subgroup are all final before this
block executes. `.fs_mr_family_from_table()` (`R/fs_mr_inference_methods.R:42-69`) then walks
`dsel$candidates` row by row, rebuilds each one- or two-cut conjunction as an explicit membership
vector on `dsel$df.est`, and keeps every candidate with at least `n_min` members. It re-derives
nothing, re-ranks nothing and fits nothing. So the family MR de-biases over is exactly DINA's own
qualifying set, fixed at the moment of the call — which is precisely why the resulting estimand is
conditional on it, and why a bootstrap that regenerated the surface would not target the same thing.

### 2.3 Admission, measured per engine (`.fs_resolve_admission`)

```
consistency  maxeff      effect_floor=NULL      consistency=NULL
consistency  maxeffCons  effect_floor=-0.1054   consistency=c_cons=-0.2231 p_star=0.90
consistency  effMaxSG    effect_floor=-0.1054   consistency=c_cons=-0.2231 p_star=0.90
consistency  maxSG       effect_floor=-0.1054   consistency=c_cons=-0.2231 p_star=0.90
dina         maxeff      effect_floor=-0.1054   consistency=NULL
dina         maxeffCons  effect_floor=-0.1054   consistency=NULL
dina         effMaxSG    effect_floor=-0.1054   consistency=NULL
dina         maxSG       effect_floor=-0.1054   consistency=NULL
```

(`hr.threshold = 0.90`, `hr.consistency = 0.80`, `pconsistency.threshold = 0.90`, on the log scale.)
DINA keeps the harm floor under every focus and never carries a consistency term — as the branch
comment above states. Consistency-engine `maxeff` uniquely carries neither floor.

---

## 3. Item 3 — `.fs_apply_mr()`'s argument list, and the STOP

### 3.1 The wrapper, verbatim (`R/fs_mr_inference_methods.R:119-148`)

```r
.fs_apply_mr <- function(df, candidates, selected_members, spec,
                                  admission,
                                  effect_neighborhood, reselection_default,
                                  selection_rule_default = "neighborhood",
                                  mr_inference_args = list(), seedit = NULL) {
  .g <- function(a, b) if (is.null(a)) b else a
  if (is.null(mr_inference_args)) mr_inference_args <- list()
  tryCatch(
    fs_mr_inference(
      df               = df,
      candidates       = candidates,
      spec             = spec,
      selected_members = selected_members,
      admission        = admission,
      t_confirm        = mr_inference_args$t_confirm,          # NULL -> near-null
      confirm_rule     = .g(mr_inference_args$confirm_rule, "point"),
      reselection      = .g(mr_inference_args$reselection, reselection_default),
      effect_neighborhood = effect_neighborhood,
      selection_rule   = .g(mr_inference_args$selection_rule, selection_rule_default),
      draws            = .g(mr_inference_args$draws,       2000L),
      multiplier       = .g(mr_inference_args$multiplier,  "poisson"),
      include_complement = .g(mr_inference_args$include_complement, TRUE),
      ci_method        = .g(mr_inference_args$ci_method,   "ij"),
      seed             = .g(mr_inference_args$seed,        seedit)),
    error = function(e) { ... NULL })
}
```

It forwards **fifteen** arguments. `fs_mr_inference()` accepts **twenty-five**. The ten it never
forwards are:

`return_reselection`, `field_R_out`, `field_R_in`, `field_uniform`, `field_M_cap`,
`field_complement`, `field_decompose`, `field_scale_complement`, `ij_residual`, `field_recovery`.

### 3.2 The five arguments the task asks about, stated plainly

| argument | forwarded by `.fs_apply_mr()`? | value the DINA branch actually gets | campaign asks for |
|---|---|---|---|
| `ci_method` | **YES** | whatever `mr_inference_args$ci_method` says (wrapper default `"ij"`; the consistency branch's default is `"field"`) | `"field"` — **reachable** |
| `field_complement` | **NO** | `fs_mr_inference()` default `TRUE` | `TRUE` — value coincides, **not controllable** |
| `field_scale_complement` | **NO** | default `"selected"` | `"selected"` — value coincides, **not controllable** |
| `field_decompose` | **NO** | default **`FALSE`** | `TRUE` — **UNREACHABLE** |
| `field_recovery` | **NO** | default **`FALSE`** | `TRUE` — **UNREACHABLE** |

Two further Stage-2 knobs land the same way: `ij_residual` is not forwarded but its default
(`"two_term"`) coincides with `FS_S7_IJ_RESIDUAL=two_term`; `return_reselection` is not forwarded but
its default (`TRUE`) means p-hat is produced anyway.

**The distinction that matters: coinciding is not controlling.** Three of the campaign's field knobs
(`FIELD_COMPLEMENT`, `FIELD_SCALEC`, `IJ_RESIDUAL`) would appear in the run's meta record as set, and
would in fact be inert — the branch would produce the right values for the wrong reason, and flipping
any of them to test the alternative would silently do nothing. Two (`FIELD_DECOMP`, `FIELD_RECOV`)
cannot produce their products at all.

### 3.3 Confirmed by measurement, not only by reading

A single pilot fit per replicate, DINA / M1 / HR 1.75 / n 500, passing the campaign's
`mr_inference_args` **verbatim** (`ci_method="field"`, `field_complement=TRUE`, `field_decompose=TRUE`,
`field_scale_complement="selected"`, `field_recovery=TRUE`, `ij_residual="two_term"`,
`return_reselection=TRUE`, `draws=2000`), then the *same* arguments on the consistency engine:

```
===== DINA sim_id=1  (4.4 s) =====         ===== CONSISTENCY sim_id=1  (16.3 s) =====
  field (harm)                 PRESENT       field (harm)                 PRESENT
  field$complement             PRESENT       field$complement             PRESENT
  field-s (est2_s)             PRESENT       field-s (est2_s)             PRESENT
  field$complement$decomp_fields ABSENT      field$complement$decomp_fields PRESENT
  field$recovery                 ABSENT      field$recovery                 PRESENT
  field$joint                  PRESENT       field$joint                  PRESENT
  field$joint_s                PRESENT       field$joint_s                PRESENT
  reselection$p_hat            PRESENT       reselection$p_hat            PRESENT
  grp.consistency$out_sg$result  ABSENT      grp.consistency$out_sg$result  data.frame 965x9
```

Identical arguments, identical DGM, identical seed. The only difference is the branch. Over the
20-replicate probe of section 4, `decomp_fields` and `recovery` were absent in **16 of 16** replicates
that produced an MR object, and the harm field, complement field, field-s, joint and p-hat were present
in **16 of 16**.

### 3.4 Is there a route that avoids an `R/` change?

**No.** The only channel from a driver into `fs_mr_inference()` on the DINA branch is
`mr_inference_args`, and `.fs_apply_mr()` drops the field arguments out of it. There is no environment
variable route: the only `Sys.getenv()` call anywhere in `R/` is the `_R_CHECK_LIMIT_CORES_` guard in
`forestsearch_main.R:59`.

Nor can a driver call `fs_mr_inference()` itself, post hoc, on DINA's family: the returned object does
not expose the family. The DINA branch returns `dina_res` — which is the fitted cross-fit model
(`coefficients`, `vcov`, `family`, `n`, `d`, `cross_fitting`, `n_folds`, `call`) — and sets
`dina_cuts = NULL`; the qualifier table `dsel$candidates` that `.fs_mr_family_from_table()` consumes is
local to the branch and is never returned. `grp.consistency$out_sg$result` is absent on this path too.
A driver-side reconstruction would have to re-derive DINA's qualifier table through unexported
internals and hope it matched — which is a re-implementation of package internals inside a document,
not a workaround.

### 3.5 STOP

> **The field products cannot reach the DINA branch without an `R/` change. That is decision O-1, it is
> Larry's, and it is not taken here.**

Per the task's protocol ("**No `R/` change.** If the campaign appears to need one, STOP and report"),
Stage 1 was not run, Stage 2 was not run, Stage 3 was not run, and no `R/` file was touched.

For the record, the shape of the change O-1 would authorize is small and additive — forwarding the ten
missing arguments through `.fs_apply_mr()` with `fs_mr_inference()`'s own defaults, so that every
existing DINA and GRF call stays byte-identical — but proposing, sizing or writing it is out of scope
here, and this report makes no recommendation.

---

## 4. Item 4 — The typical proposed-family size per replicate

Measured, since no committed run records it (`n_family` does not exist in the 57-column DINA recorder).
DINA / M1 / HR 1.75 / n 500 / `knoise0` / `sg_focus = "maxeff"` / `dina_select_statistic = "effect"`,
seeds `8316951 + sim_id`, sim_id 1–20, template DGM defaults (`z1_quantile` 0.25, super-population harm
prevalence **0.1242**, `k_inter` 1.32025).

**Detection: 16/20 (80%).** All 16 detections produced an MR object.

**Proposed-family size (`n_family`), over the 16 replicates with an MR object:**

| statistic | value |
|---|---|
| mean | 336.5 |
| min | 15 |
| 10th pct | 25.0 |
| 25th pct | 45.3 |
| median | **91** |
| 75th pct | 364.5 |
| 90th pct | 971.0 |
| max | 1745 |

**The family is not degenerate.** The task flagged that a family of one would make selection
deterministic and the correction vacuous; that is not what happens. The smallest family observed was
15 candidates, the median 91. p-hat mass is correspondingly spread: over the three fits inspected in
detail the top re-selection probability was 0.654 (family 35), 0.364 (family 126) and 0.670 (family
1745) — the selected subgroup is usually the modal re-selection but rarely a certainty.

The distribution is **strongly right-skewed** — two orders of magnitude between the 10th percentile and
the maximum — and per-replicate cost tracks it closely (0.13 s for a non-detection, 3.7 s at family 15,
27.8 s at family 1745; median 5.1 s, mean 7.6 s over the 20). Any future timing projection for this
engine must be built from the family-size distribution rather than from a mean, and must not be
projected from the FS walls: DINA fits a model per replicate and its field cost scales with a family
size that varies by 100x across replicates. No Gate 1 projection is offered here, because Stage 1 was
not reached.

---

## 5. The STOP report the task asks for

### (a) Closest transplant source

`quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` — as the task expected.
It is already **engine-generic**: `subgroup_method` is a setup-chunk knob, and the call assembles a
per-engine block that already includes DINA:

```r
  method_args <- switch(
    subgroup_method,
    consistency = list(
      consistency_method = consistency_method,
      use_lasso = use_lasso, use_grf = use_grf, use_twostage = use_twostage,
      use_dina = use_dina,
      conf_force = fs_conf_force, conf.cont_jcuts = fs_conf.cont_jcuts,
      fs.splits = fs_splits, maxk = maxk, d0.min = d0_min, d1.min = d1_min),
    dina = list(dina_select_statistic = dina_select_statistic,
                dina_args             = dina_args),
    grf  = list(grf_selection = grf_selection, grf_depth = grf_depth,
                dmin.grf = dmin.grf, grf_select_statistic = grf_select_statistic),
    stop("unknown subgroup_method: ", subgroup_method)
  )
```

The recorder is likewise engine-blind: every field read is `%||% NA`, so a DINA run records the same
158 columns and simply leaves absent blocks NA.

### (b) Exact list of changes a transplant would require

1. `subgroup_method <- "consistency"` → `"dina"` (line 300).
2. `sg_focus` (line 307, `.env_chr("FS_S7_FOCUS", "maxeffCons")`) → `"maxeff"` for this campaign.
   Note the tag consequence: `fs_focus_tag("dina", "maxeff")` is `"eff"`, so the output stem token
   becomes `dina_eff_...`. This does **not** collide with the committed pre-field bundles
   (`dina_eff_fb_mr_m1_...`), because the field template's stem carries the extra `_field_` token —
   `dina_eff_fb_mr_field_m1_...`. Worth verifying on the first render rather than assuming.
3. The DINA selector knobs already present and already correct: `dina_select_statistic <- "effect"`,
   `dina_args <- list()` (lines 496–497). No change needed.
4. Summary-side globs and comparator labels in the Stage 3 transplant (`summary_cert20.qmd` →
   `summary_dinamr.qmd`), per the task: globs and labels only.
5. **Nothing else in the driver.** The transplant is genuinely a knob change — which is exactly why
   the blocker is not in the document at all, but in `R/`.

### (c) Which `FS_S7_*` knobs go inert on a DINA transplant

| knob | status on DINA |
|---|---|
| `FS_S7_FIELD_DECOMP` | **inert — product unreachable.** Not forwarded; `field_decompose` defaults `FALSE`. |
| `FS_S7_FIELD_RECOV` | **inert — product unreachable.** Not forwarded; `field_recovery` defaults `FALSE`. |
| `FS_S7_FIELD_COMPLEMENT` | inert, but the default (`TRUE`) coincides with the campaign value. |
| `FS_S7_FIELD_SCALEC` | inert, but the default (`"selected"`) coincides with the campaign value. |
| `FS_S7_IJ_RESIDUAL` | inert, but the default (`"two_term"`) coincides with the campaign value. |
| `FS_S7_UNIFORM` | inert (not forwarded); default `FALSE`, which the campaign does not override. |
| `FS_S7_RETURN_RESEL` | inert (not forwarded); default `TRUE`, so p-hat is produced regardless. |
| `FS_S7_ER_JCUTS` | inert — feeds `fs_conf.cont_jcuts`, which lives in the **consistency-only** `method_args` block and is never passed to DINA. The task's "J = 10 where applicable" does not apply. |
| `FS_S7_NBHD` | inert under `maxeff` on either engine (band half-width; no band term in `order(-eff)`). |
| `FS_S7_HR`, `FS_S7_N`, `FS_S7_KNOISE`, `FS_S7_Z1Q`, `FS_S7_NSIMS`, `FS_S7_START`, `FS_S7_WORKERS`, `FS_S7_MODE`, `FS_S7_FB`, `FS_S7_CAMPAIGN`, `FS_S7_QUICKRUN`, `FS_S7_SAVE_COMBINED`, `FS_S7_JOIN_SKIP`, `FS_S7_WINNER_ROWS`, `FS_S7_FB_PATH` | **live** — DGM, batching, campaign tagging and reporting, all engine-independent. |

`ci_method` is not an `FS_S7_` knob but is the one field-related setting that **is** live on DINA.

### (d) Which recorder columns lose meaning without a consistency screen

Measured on the pilot, not inferred:

- **`n_cons_qual`** and **`band_n`** — both read `fs.est$grp.consistency$out_sg$result`, which is
  **absent** on the DINA branch (present on consistency as a 965x9 frame with columns
  `Pcons, hr, N, E, g, m, K, M.1, M.2`). Both would be NA in every DINA row. They are the
  consistency identifier's own qualifier table and its observed-effect band count; DINA has neither.
- **`p_star`** is not a recorder column. It is the admission set's consistency term, and on DINA it is
  deliberately `NULL` (section 2.3) — the branch comment records why: "DINA and GRF have no
  consistency screen, so no consistency term may enter their admission set -- previously
  `c_consistency = 0` with `p_star` still set meant `t_g` carried a `z * sigma_D` term these engines
  never applied."
- The `Pcons`-derived reporting in any transplanted summary would have to be dropped rather than
  rendered as NA rows, which is a `summary_dinamr.qmd` concern, not a recorder one.

### (e) Do the field / field-s / joint / p-hat / recovery columns reach the DINA branch at all?

Measured, 16/16 replicates with an MR object:

| product | recorder columns | reaches DINA? |
|---|---|---|
| harm field | `fld_H_*` | **YES** (requires `ci_method="field"`, which is forwarded) |
| complement field | `fld_Hc_*` | **YES** (by `fs_mr_inference()` default, not by the knob) |
| field-s | `fld_Hc_*_s` | **YES** (by default `"selected"`, not by the knob) |
| Bonferroni joint pair | `fld_joint_*`, `fld_joint_s_*` | **YES** |
| p-hat | `p_hat_*`, `n_family` | **YES** (by `return_reselection` default `TRUE`) |
| complement scale decomposition | `fld_Hc_scale_*` | **NO — always NA** |
| recovery diagnostics | `fld_recov_*` | **NO — always NA** |

So the campaign would have delivered five of its seven product families, with three of them switched on
by accident rather than by configuration, and two absent. Reporting those two as "absent" rather than
silently skipping them is exactly what the task's Gate 2 requires — which is why running the cells and
discovering it at Gate 2 would have been the more expensive way to learn this.

---

## 6. Side issues, flagged not fixed

1. **Prevalence mismatch against the Stage 3 comparator.** The template's default `FS_S7_Z1Q` is
   `0.25`, which on this DGM gives super-population harm prevalence **0.1242** (measured) — the
   `tier2` / `p12ext` prevalence, not `cert20`'s. `cert20` ran at `z1q60` (31%), visible in its own
   stems (`fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_z1q60_nb20_cert20_...`). The task specifies no
   `FS_S7_Z1Q` for `dinamr` but directs the Stage 3 transplant from `summary_cert20.qmd`. Unset, the
   campaign would run at 12.4% and be reported in a frame transplanted from a 31% study. This needs a
   decision before any future Stage 2, and it is not taken here.
2. **`cert20`'s bundles carry no `fld_recov_*` columns** (0 of 158) — it predates `field_recovery`
   (2026-09-09). Any "FS reference beside DINA" table for the recovery diagnostics has no FS column to
   place beside it from `cert20`.
3. **The `.fs_apply_mr()` gap is not DINA-specific** — GRF reaches `fs_mr_inference()` through the same
   wrapper (`R/forestsearch_main.R:2534`) and is affected identically. GRF is explicitly out of scope
   for this task; noted only because O-1 would resolve both at once.

---

## 7. Provenance

- Installed package **0.3.5**, matching `DESCRIPTION` at HEAD; `R/` and `DESCRIPTION` clean at commit
  `4b64d3b5`. `fs_mr_inference()` formals confirmed to include `field_recovery`, so the installed copy
  carries the 2026-09-09 work. `load_all()` was not used.
- Pilot fits were run as standalone `Rscript` against the installed package, writing nothing into the
  repository. The seven pre-existing untracked files were not touched.
- Measurements in sections 1.2, 3.3 and 4 are from those pilot fits and from committed `.rds` bundles
  under `quarto/simulations/gbsg_020/results/`; every code block in sections 2 and 3.1 is quoted
  verbatim from `R/` at this commit.

**No acceptance criteria were pre-registered, none are proposed, and this report makes no
recommendation. Report and wait.**
