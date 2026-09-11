# REPORT — O-1: DINA and GRF made operational for the current product set

**Date:** 2026-09-10. **Executor:** Claude Code (Linux, pop-os). **Task:** `dev/tasks/TASK_o1_forwarding_2026-09-10.md`.
**Resolves:** the STOP recorded in `REPORT_dinamr_stage0_2026-09-10.md` §3.5.
**Outcome: Gate F PASSES on both engines** (Fa, Fc, Fd literal; Fb and Fe substituted with the reason stated and the substitute defended).
**One `R/` change, confined to `R/fs_mr_inference_methods.R`; adds code, byte-identical defaults. No template change. No default change on any branch.**
**Report and wait: no recommendation, no campaign, no acceptance criteria.**

---

## Part F — The forwarding fix

### F1. Before and after, quoted

**Before** — `.fs_apply_mr()` forwarded **15** of `fs_mr_inference()`'s **25** arguments:

```
df, candidates, spec, selected_members, admission, t_confirm, confirm_rule,
reselection, effect_neighborhood, selection_rule, draws, multiplier,
include_complement, ci_method, seed
```

**After** — it forwards all **25**. The ten added:

```
return_reselection, field_R_out, field_R_in, field_uniform, field_M_cap,
field_complement, field_decompose, field_scale_complement, ij_residual,
field_recovery
```

**The formals the defaults are read from**, quoted verbatim from `R/fs_mr_inference.R:549-571`:

```r
fs_mr_inference <- function(df, candidates, spec, selected_members,
                           admission,
                           t_confirm = NULL, confirm_rule = c("point", "ci"),
                           reselection = c("maxcons", "maxeff", "maxSG",
                                           "minSG", "effMaxSG", "effMinSG"),
                           effect_neighborhood = 0.10,
                           selection_rule = c("neighborhood", "pareto", "both"),
                           draws = 2000L,
                           multiplier = c("poisson", "gaussian", "rademacher"),
                           include_complement = FALSE,
                           ci_method = c("field", "ij", "wald"),
                           seed = NULL,
                           return_reselection = TRUE,
                           field_R_out = 1000L,
                           field_R_in = 500L,
                           field_uniform = FALSE,
                           field_M_cap = NULL,
                           field_complement = TRUE,
                           field_decompose = FALSE,
                           field_scale_complement = c("selected", "none"),
                           ij_residual = c("two_term", "winner", "winner_floor"),
                           field_recovery = FALSE) {
```

**The defaults are not restated in the wrapper.** They are read from those formals at call time:

```r
  .mr_fml <- formals(fs_mr_inference)
  .d <- function(nm) eval(.mr_fml[[nm]], envir = environment(fs_mr_inference))
```

and each forwarded argument uses the existing `.g()` idiom, e.g.

```r
      field_decompose  = .g(mr_inference_args$field_decompose,
                            .d("field_decompose")),
```

`match.arg()`-style defaults are forwarded as their whole candidate vector, which
`fs_mr_inference()` resolves exactly as it resolves a missing argument. `ci_method` is
untouched and keeps this wrapper's own `"ij"` default — deliberately *not* the consistency
branch's `"field"`; changing it is a separate decision and is not taken here.

### F2. Documentation

Roxygen `@section Forwarded argument set:` added to `.fs_apply_mr()` stating that the wrapper
now forwards the full set, that the DINA and GRF branches are therefore controllable through
`mr_inference_args`, and that defaults are unchanged. `devtools::document()` regenerated
`man/dot-fs_apply_mr.Rd` (the only Rd affected; `NAMESPACE` unchanged — the change adds no
export and no import). One `NEWS.md` bullet under the development header.

> **Pre-existing tooling defect, flagged not fixed.** `devtools::document()` aborts with
> `fs_bias_coverage.R:17: @description failed to evaluate inline markdown code` —
> `` `r = se_mean / sd_emp` `` is parsed by roxygen as inline R. **Falsified as unrelated:**
> reverting `R/fs_mr_inference_methods.R` to HEAD and re-running `document()` reproduces the
> identical error, and the block predates this task (`git log` shows its last touch at
> `2f118042`). It does not prevent the affected Rd from being written, and `R CMD check` does
> not re-run roxygen. It is in the same file as the pre-existing check WARNING (§Tallies).

### F3. Install and read-back

`devtools::install(dependencies = FALSE)`. Verified on the **installed** package:

```
deparse(installed .fs_apply_mr) == deparse(source)            TRUE
call node arguments                                           25
fs_mr_inference() formals                                     25
forwarded set == accepted set                                 TRUE   (missing: none)
```

Resolved defaults, read back from the installed wrapper's call node:

| argument | wrapper expression | resolves to |
|---|---|---|
| `return_reselection` | `.g(mr_inference_args$return_reselection, .d("return_reselection"))` | `TRUE` |
| `field_R_out` | `.g(mr_inference_args$field_R_out, .d("field_R_out"))` | `1000L` |
| `field_R_in` | `.g(mr_inference_args$field_R_in, .d("field_R_in"))` | `500L` |
| `field_uniform` | `.g(mr_inference_args$field_uniform, .d("field_uniform"))` | `FALSE` |
| `field_M_cap` | `.g(mr_inference_args$field_M_cap, .d("field_M_cap"))` | `NULL` |
| `field_complement` | `.g(mr_inference_args$field_complement, .d("field_complement"))` | `TRUE` |
| `field_decompose` | `.g(mr_inference_args$field_decompose, .d("field_decompose"))` | `FALSE` |
| `field_scale_complement` | `.g(mr_inference_args$field_scale_complement, .d("field_scale_complement"))` | `c("selected", "none")` |
| `ij_residual` | `.g(mr_inference_args$ij_residual, .d("ij_residual"))` | `c("two_term", "winner", "winner_floor")` |
| `field_recovery` | `.g(mr_inference_args$field_recovery, .d("field_recovery"))` | `FALSE` |

Every value equals `fs_mr_inference()`'s own formal default, because it *is* that default.

---

### Gate F — results

M1 DGM, HR 1.75, n = 500, seeds `8316951 + sim_id`, sim_id 1–5, `z1_quantile` 0.25
(template default; super-population harm prevalence 0.1242, `k_inter` 1.32025).

#### Fa — byte-identical when nothing is asked for: **PASS**

`mr_inference_args` **omitted entirely**; both captures built explicitly, the pre-change one
against the pre-change *installed* package, timing elements stripped recursively.

| capture | `identical()` | `n_family` |
|---|---|---|
| `dina_1` | TRUE | 35 |
| `dina_2` | TRUE | 126 |
| `dina_3` | TRUE | 1745 |
| `dina_5` | TRUE | 37 |
| `grf_1` | TRUE | 772 |
| `grf_2` | TRUE | 770 |
| `grf_3` | TRUE | 778 |
| `grf_4` | TRUE | 783 |
| `grf_5` | TRUE | 785 |

Whole-list `identical(pre$Fa, post$Fa)`: **TRUE**. `dina_4` appears in neither capture: sim_id 4
is a non-detection on DINA, so no MR object exists — the same on both sides, and the key sets
are `identical()`.

**This is the criterion that makes the change add-only, and it holds on both engines.**

#### Fb — byte-identical against committed work: **PASS (substituted; reason stated)**

The literal check **cannot be evaluated**, for a reason that predates and is independent of
O-1. The seven committed DINA/M1 drivers wrote bundles under package **0.2.0** (and one under
0.2.1) on a **Mac Studio**; the tree is **0.3.5** on Linux. Reproducing
`sim_dina_maxeff_mr_m1_h10_knoise0_n500_batch_1_500.qmd`'s settings at sim_id 1–5 today gives:

| | sim_id 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| committed `detected` (v0.2.0) | 1 | 1 | 1 | 1 | 1 |
| reproduced `detected` (v0.3.5) | 0 | 0 | 1 | 0 | 0 |
| committed `sg_def` (1) | `{age >= 47} & {meno <= 0}` | | | | |
| reproduced `sg_def` (3) | `{er <= 4} & {nodes >= 2}` | | | | |

The identifier itself moved across that release history (the alignment repair, the `maxeff`
override semantics, admission resolution). No `identical()` against a v0.2.0 bundle is
achievable by any change, and a failure there would carry no information about O-1.

**Substitute, which isolates O-1 exactly:** run that driver's own settings (HR 1.00, n 500,
`sg_focus = "maxeff"`, `dina_select_statistic = "effect"`, `mr_draws = 5000`,
`include_complement = TRUE`, no `ci_method` — so the `"ij"` path with no field block, a
genuinely different route through the wrapper than Fa's) under the **pre-change** and
**post-change** installed packages and require `identical()`.

```
identical() on all 22 recorded columns: TRUE
differing columns: 0
```

#### Fc — now controllable: **PASS on both engines**

`field_decompose = TRUE, field_recovery = TRUE` through `mr_inference_args`, reported in the
Stage 0 report §3.3 format:

**DINA**

| | before | after |
|---|---|---|
| `field (harm)` | PRESENT | PRESENT |
| `field$complement` | PRESENT | PRESENT |
| `field-s (est2_s)` | PRESENT | PRESENT |
| `field$complement$decomp_fields` | **ABSENT** | **PRESENT** |
| `field$recovery` | **ABSENT** | **PRESENT** |
| `field$joint` | PRESENT | PRESENT |
| `field$joint_s` | PRESENT | PRESENT |
| `reselection$p_hat` | PRESENT | PRESENT |

**GRF** — identical pattern, the same two rows changing ABSENT → PRESENT.

Contents, sim_id 1:

| | DINA | GRF |
|---|---|---|
| `scale_sel` | 0.1379943 | 0.1305426 |
| `scale_win_mean` | 0.1373797 | 0.1347629 |
| `scale_win_cv` | 0.01190383 | 0.02460132 |
| `scale_ratio_c` | **1.004473** | **0.9686836** |
| `sens_H` | **0.7054431** | **0.2026579** |
| `ppv_H` | 0.6517603 | 0.2026967 |
| `sens_Hc` | 0.9445961 | 0.8601580 |
| `npv_Hc` | 0.9583106 | 0.8576637 |
| `q10 / q50 / q90` | 0.4754 / 0.6066 / 1.0000 | 0.0000 / 0.1447 / 0.6579 |
| `share_equal_1` | 0.3623624 | 0.0150000 |
| `n_used` | 999 | 1000 |

Standing invariants hold on both: `0 <= sens_H <= 1` TRUE, `0 <= ppv_H <= 1` TRUE,
`scale_ratio_c > 0` TRUE. (Over the 20-replicate Part V pilots these hold on every detected
replicate on both engines, together with `gamma` in `[0.025, 0.05]` for both the joint and
joint-s pairs, `fld_H_lo1s <= fld_H_est2`, `fld_Hc_est2 <= fld_Hc_up1s`, and finiteness.)

#### Fd — the inert three are now live: **PASS on both engines**

`field_complement = FALSE`, `field_scale_complement = "none"`, everything else as Fc:

| | DINA before | DINA after | GRF before | GRF after |
|---|---|---|---|---|
| `field$complement` | PRESENT | **ABSENT** | PRESENT | **ABSENT** |
| `field-s (est2_s)` | PRESENT | **ABSENT** | PRESENT | **ABSENT** |
| `field$joint` | PRESENT | **ABSENT** | PRESENT | **ABSENT** |
| `field$joint_s` | PRESENT | **ABSENT** | PRESENT | **ABSENT** |

```
pre-change:  Fd output == Fc output   DINA TRUE   GRF TRUE     (inert: asking changed nothing)
post-change: Fd output == Fc output   DINA FALSE  GRF FALSE    (control proven)
```

This is the provenance half of the fix: before, a campaign's meta would have recorded these as
set while they controlled nothing; now setting them changes the output, so the meta means what
it says.

#### Fe — the consistency branch is untouched: **PASS (substituted; reason stated)**

**By construction** the consistency branch cannot be affected: it does not call
`.fs_apply_mr()` at all, but `fs_mr_inference()` directly at `R/forestsearch_main.R:3395`.

The literal check against the committed `e1stud` bundle
(`fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n500_z1q60_nb20_e1stud`, v0.3.5, pop-os, R 4.6.1) is
**not evaluable**, and the blocker is upstream of MR entirely: the reproduction does not draw
the same trials. `n_true` — `sum(flag_harm == 1)` on the simulated trial, computed before any
search runs — differs (reproduced 161, 156, 140, 157, 175; committed 154, 172, 162, 171, 155),
while the super-population prevalence matches (0.3066 vs the recorded 0.307). The DGM source
files are unchanged since `e1stud` was built (`git log --since="2026-09-08 11:00"` over
`R/sim_aft_gbsg.R`, `R/setup_gbsg_dgm.R`, `R/simulate_from_dgm.R` and the DGM family returns
nothing), and `calibrate_k_inter()` is deterministic across sessions, so the most likely
account is that `forestsearch_version` does not pin a commit — several commits share `0.3.5`.
**This is recorded as an open side issue (§Side issues 2), not resolved here.**

**Substitute:** the same consistency configuration (effMaxSG, ε 0.20, z1q 0.60, `ci_method`
`"field"`, `field_decompose = TRUE`, `field_scale_complement = "selected"`, `ij_residual`
`"two_term"`, `mr_draws` 5000, `er_jcuts` 10, sim_id 1–5) run under the **pre-change** and
**post-change** installed packages:

```
columns compared: 38     identical() pre vs post: TRUE     differing columns: 0
```

covering `n_family`, `n_cons_qual`, the naive and MR H/Hc blocks, `fld_H_*`, `fld_Hc_*`,
`fld_Hc_est2_s` / `up1s_s`, `fld_Hc_scale_ratio`, the joint pair and `p_hat_H`.

---

## Part V — Verification for a future campaign

### V1. Does GRF's `effMaxSG` honor `effect_neighborhood` and share the inclusion-band logic?

**Yes on both counts — with two qualifications that matter for a campaign.**

`R/grf_subgroup_labels.R:355-378`, verbatim:

```r
.grf_frontier_select <- function(cand, dmin, rule = "effMaxSG", nbhd = 0.10,
                                 selection_rule = "neighborhood") {
  ...
  } else { # effMaxSG / effMinSG: restrict to the inclusion band first
    # Identical to the former hardcoded band under the default
    # selection_rule = "neighborhood": .compute_inclusion_band() applies
    # (1 - effect_neighborhood) * max(effect), which is the same comparison as
    # the previous effect >= emax * (1 - nbhd).
    in_band <- .compute_inclusion_band(
      hr_vec              = elig$effect,
      n_vec               = elig$size,
      selection_rule      = selection_rule,
      effect_neighborhood = nbhd)
    band <- elig[in_band == 1L, , drop = FALSE]
```

`forestsearch()` threads `effect_neighborhood` through as `nbhd = config$effect_neighborhood`
(`R/grf_main.R:293`). `.compute_inclusion_band()` (`R/subgroup_consistency_helpers.R:778`) is
the **one shared helper**, called by the consistency stages (`:601`, `:689`),
`forestsearch_helpers.R` (`:1252`, `:1258`), DINA (`dina_subgroup.R:519`, `:528`), MR
(`fs_mr_inference.R:147`) and GRF. So the band is the same band.

**Qualification 1 — frontier only.** The band exists only on `grf_selection = "frontier"`;
the policy-tree path never reaches `.grf_frontier_select()`, and `effect_neighborhood` is
documented "Used only when `grf_selection = 'frontier'`" (`R/grf_main.R:50-52`). Through
`forestsearch()` the default *is* frontier (`grf_selection = c("frontier", "tree")` →
`"frontier"`), which is also what the template sets — note this differs from the inner
`grf_subgroup_harm()`'s own roxygen, which names `"tree"` as its default.

**Qualification 2 — the band is applied differently from DINA's, and this is the one thing a
campaign must know.** DINA uses the band as a **sort key**
(`order(-in_band, -cand_n, -eff, idx)`, `R/dina_subgroup.R:519-535`), so an all-zero band is
harmless — the next key breaks the tie and a winner always exists. GRF uses it as a **filter**
with **no empty-band fallback**, deliberately; the source says so:

```r
    # No empty-band fallback here, deliberately.  The band CAN empty: when the
    # maximum effect is negative, (1 - nbhd) * emax exceeds emax, so even the
    # maximum fails its own test. ...  MR's .inband() does carry a "never empty"
    # fallback; reconciling the two is a decision pending on its own, and adding
    # a fallback here would pre-empt it and change behaviour beyond threading.
```

So on GRF an empty band yields a zero-row selection (no detection), while MR's re-selection
over the same band carries a "never empty" fallback (`R/fs_mr_inference.R:150-162`). The
identifier and MR can therefore disagree about an empty band on GRF in a way they cannot on
DINA. **That asymmetry is an open decision in the source, not something this task resolves.**
It did not bind in the pilot (GRF detected 20/20), but it is a property of the formula, not of
the data.

### V2. Admission per engine, extended

`.fs_resolve_admission(sg_focus, subgroup_method, hr.threshold = log(0.90),
hr.consistency = log(0.80), pconsistency.threshold = 0.90)`, measured:

| `subgroup_method` | `sg_focus` | effect floor (log) | consistency |
|---|---|---|---|
| consistency | `maxeff` | **NULL** | **NULL** |
| consistency | `maxeffCons` | −0.105361 | `c_cons` −0.223144, `p_star` 0.90 |
| consistency | **`effMaxSG`** | **−0.105361** | **`c_cons` −0.223144, `p_star` 0.90** |
| consistency | `effMinSG` | −0.105361 | `c_cons` −0.223144, `p_star` 0.90 |
| consistency | `maxSG` / `minSG` / `hr` | −0.105361 | `c_cons` −0.223144, `p_star` 0.90 |
| dina | `maxeff` / `maxeffCons` | −0.105361 | NULL |
| dina | **`effMaxSG`** | **−0.105361** | **NULL** |
| dina | `effMinSG` / `maxSG` / `minSG` / `hr` | −0.105361 | NULL |
| grf | `maxeff` / `maxeffCons` | −0.105361 | NULL |
| grf | **`effMaxSG`** | **−0.105361** | **NULL** |
| grf | `effMinSG` / `maxSG` / `minSG` / `hr` | −0.105361 | NULL |

**At the exact campaign setting — `effMaxSG`, ε 0.20 — DINA and GRF carry the effect floor
log(0.90) = −0.105361 and no consistency term; the consistency engine carries both.** ε is
**not** part of the admission set: it is the inclusion-band width applied *after* admission, so
this table is identical at ε 0.10 and ε 0.20. Consistency-engine `maxeff` remains the one cell
with neither floor.

### V3. Stem collision

Template stem (`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:409-412`):
`"%s_%s_fb_mr_field_m1_h%03d_knoise%d_n%d%s%s%s_%s"` over
`method_tag, focus_tag, 100*target_hr_harm, k_random_noise, n_sample, z1q_tag, nbhd_tag,
jcuts_tag, campaign_tag`. At `sg_focus = "effMaxSG"`, `fs_focus_tag()` returns `"effMaxSG"` on
**both** engines (unlike `maxeff`, which collapses to `"eff"` on DINA/GRF). At HR 1.75, n 500,
`knoise0`, z1q 0.25 (no tag), ε 0.20 (`_nb20`), J 10 (no tag):

```
dina_effMaxSG_fb_mr_field_m1_h175_knoise0_n500_nb20_<campaign>
grf_effMaxSG_fb_mr_field_m1_h175_knoise0_n500_nb20_<campaign>
```

**No collision is possible.** Every committed DINA/GRF bundle is pre-field and carries
`_mr_m1_` or `_fb_mr_m1_` — never `_fb_mr_field_`:

```
dina_eff_fb_mr_m1_h10_knoise0_n500_{combined_1_500, res_1_200, res_201_250, res_251_500}.rds
dina_maxcons_mr_m1_h10_knoise0_n500_{quickrun_res_1_500, res_1_20}.rds
dina_maxeffCons_mr_m1_h10_knoise0_n500_quickrun_res_1_500.rds
dina_maxeff_mr_m1_h10_knoise0_n500_quickrun_res_1_500.rds
grf_eff_fb_mr_m1_h10_knoise0_n500_{combined_1_500, quickrun_res_1_500, res_1_200, res_201_250, res_251_500}.rds
```

Four independent distinguishers separate them: the `field` token, the focus tag
(`effMaxSG` vs `eff`/`maxcons`/`maxeff`/`maxeffCons`), `h175` vs `h10`, and `_nb20`.
A glob for `^(dina|grf)_effMaxSG_fb_mr_field` and for `^(dina|grf)_.*_fb_mr_field_` both
return **nothing** against the committed results directory.

### V4. Recorder meaning, measured per engine

Canonical recorder = the committed `e1stud` bundle's **158** columns. The pilot computed the
**56** of them that are MR-side and engine-dependent (identification, classification, naive,
MR H/Hc including the winner-only and winner-floor variants, the whole field / field-s / joint
/ joint-s / p-hat blocks, `n_family`, `n_cons_qual`, `band_n`).

**Result, identical on both engines: 52 populated, 4 all-NA.**

All-NA, with the reason for each:

| column | why |
|---|---|
| `n_cons_qual` | **structurally NA without a consistency screen** — reads `grp.consistency$out_sg$result`, which neither engine produces |
| `band_n` | **structurally NA** — same source |
| `fld_H_kappa` | configuration: `field_uniform = FALSE` (the campaign does not set it) |
| `fld_H_note` | correct emptiness: the note is `NULL` when the field block is healthy |

`p_star` is **not** a recorder column; it is the admission set's consistency term and is `NULL`
on both engines (V2), by the branch's own design.

The remaining 102 canonical columns were not measured in this pilot and are engine-independent
or configuration-determined: `or_*` (oracle refit), `fb_*` (FB bootstrap, off under
`FS_S7_FB=none`), `betaHhat_*` / `nH_eval` / `nHc_eval` (the DGM eval frame), `err_msg`,
`covs`, the timing columns, and the `fld_H_q*` / `fld_Hc_q*` quantile and uniform blocks.

**The field / `_s` / joint / p-hat / recovery columns after the fix** — all populated on both
engines: `fld_H_est2`, `fld_H_lo1s`, `fld_H_se`, `fld_H_nout`, `fld_H_nin_mean`;
`fld_Hc_est2`, `fld_Hc_up1s`, `fld_Hc_se`, `fld_Hc_nout`, `fld_Hc_nfit`;
**`fld_Hc_scale_sel`, `fld_Hc_scale_win`, `fld_Hc_scale_cv`, `fld_Hc_scale_ratio`** (new);
`fld_Hc_est2_s`, `fld_Hc_up1s_s`, `fld_Hc_se_s`; `fld_joint_gamma`, `fld_joint_bonf_loH`,
`fld_joint_bonf_upHc`, `fld_joint_n`; `fld_joint_s_gamma`, `fld_joint_s_bonf_loH`,
`fld_joint_s_bonf_upHc`; `p_hat_H`, `p_hat_sum`, `p_hat_top1`; `n_family`.

The **nine recovery columns are new and are not among the 158** — they postdate `e1stud`
(added 2026-09-09), so a fresh render of the current template records **167**. Measured means
over detected replicates:

| | `sens_H` | `ppv_H` | `sens_Hc` | `npv_Hc` | `q10` | `q50` | `q90` | `share1` | `n_used` |
|---|---|---|---|---|---|---|---|---|---|
| DINA | 0.5961 | 0.6010 | 0.9212 | 0.9227 | 0.3418 | 0.5785 | 0.8936 | 0.1007 | 972.2 |
| GRF | 0.4466 | 0.4627 | 0.8773 | 0.8676 | 0.1320 | 0.4242 | 0.7943 | 0.0658 | 1000 |

Scale decomposition (means): DINA `scale_sel` 0.1413, `scale_win_mean` 0.1408, `cv` 0.02943,
`scale_ratio_c` 1.004; GRF 0.1425, 0.1408, 0.04595, 1.013.

### V5. Engine knobs

| knob | template | `forestsearch()` default | campaign-appropriate? |
|---|---|---|---|
| `dina_select_statistic` | `"effect"` | `c("effect", "dina")` → **`"effect"`** | **Yes.** Alignment requires ranking on the inferential effect; template pins what the package already defaults to. |
| `grf_select_statistic` | `"effect"` | `c("effect", "dr")` → **`"effect"`** | **Yes**, same reason. |
| `grf_selection` | `"frontier"` | `c("frontier", "tree")` → `"frontier"` | **Yes, and necessary** — the `effMaxSG` band exists only on the frontier path (V1). |
| `grf_depth` | `2L` | `2` | Yes — caps GRF at two-cut conjunctions, matching `maxk = 2L` on the consistency engine, so the engines are compared at the same conjunction depth. |
| `dmin.grf` | `0.0` | `0` | **A campaign must decide this one.** See below. |

**Both `*_select_statistic` resolve to `"effect"`** — confirmed by `match.arg()` on the
package formals, not by reading the template alone.

**`dmin.grf` is the one knob that is not obviously right.** It is GRF's own eligibility floor,
applied as `elig <- cand[cand$effect >= dmin, ]` *before* the band. Its units are **not** log
HR: `.grf_dr_candidates()` (`R/grf_subgroup_labels.R:271`) sets
`effect = mean(ctrl[S]) - mean(trt[S])` — a doubly-robust score difference on the outcome
scale. So `dmin.grf = 0.0` admits any non-negative DR contrast and is **not** comparable to the
admission effect floor log(0.90) = −0.105361 that MR's re-selection uses (V2). A GRF campaign
therefore runs with the identifier screening in DR units and MR re-selecting under a log-HR
floor. That is a units mismatch between the identifier's own screen and the admission set; it
is not created by this task and is not resolved here, but a campaign must state which floor it
intends.

Related and worth recording: GRF's candidate family is fixed by hardcoded grids —
`grid_probs = seq(0.1, 0.9, 0.1)` for depth 1 and `seq(0.2, 0.8, 0.2)` for depth 2 — so
`er_jcuts` (`FS_S7_ER_JCUTS`) is inert on GRF as it is on DINA. This is why GRF's family size is
nearly constant across replicates (V6).

### V6. Cost, 20 replicates per engine, `effMaxSG` ε 0.20, HR 1.75, n 500

Seeds `8316951 + sim_id`, sim_id 1–20; z1q 0.25 (prevalence 0.1242, `k_inter` 1.32025);
full campaign `mr_inference_args`; 20 workers, one replicate per worker.
**Data only — no projection is offered.**

| | DINA | GRF |
|---|---|---|
| **detection** | **16/20 = 80.0%** | **20/20 = 100.0%** |
| `n_family` mean | 336.5 | 776.9 |
| `n_family` median | 91 | 776 |
| `n_family` 10th pct | 25 | 770 |
| `n_family` 90th pct | 971 | 785 |
| `n_family` min | 15 | 766 |
| `n_family` max | 1745 | 787 |
| secs/replicate mean | 12.34 | 24.72 |
| secs/replicate median | 9.94 | 24.99 |
| secs/replicate min | 0.11 | 22.41 |
| secs/replicate max | 40.27 | 25.70 |
| total worker-seconds | 246.7 | 494.3 |
| wall (20 workers) | 0.7 min | 0.5 min |
| `n_sel` median (range) | 76 (60–127) | 88 (66–143) |

**The two engines have opposite cost profiles.** DINA's family spans two orders of magnitude
(15–1745) and its per-replicate cost tracks it (0.11 s for a non-detection to 40.27 s at the
largest family), so a DINA projection must be built from the family-size *distribution*. GRF's
family is nearly deterministic (766–787, a 2.7% spread) because the candidate grid is hardcoded
(V5), and its per-replicate cost is correspondingly flat (22.4–25.7 s) — GRF is roughly twice
DINA's mean cost per replicate but far more predictable. Neither profile resembles the FS walls.

DINA's `n_family` distribution at `effMaxSG` ε 0.20 is **identical** to the one measured at
`maxeff` in the Stage 0 report (mean 336.5, median 91, min 15, max 1745), which is the expected
behaviour and a useful consistency check: the family is DINA's qualifying candidate table,
fixed by the admission set, and the focus and band govern *selection within* it, not its size.

---

## Verification tallies

- **Full test suite:** `FAIL 0 | WARN 32 | SKIP 3 | PASS 5051` — **identical** to the tally
  recorded in `REPORT_field_recovery_2026-09-09.md`. No test was added, edited or skipped.
- **`R CMD check --as-cran`** (`rcmdcheck::rcmdcheck(args = "--as-cran")`, the certification
  surface, which builds the PDF manual): **0 errors | 1 warning | 2 notes** — **unchanged**
  from the baseline, and all three are the same pre-existing items, none in a file this task
  touched
  - Baseline for comparison: **0 errors | 1 warning | 2 notes**, recorded in
    `REPORT_field_recovery_2026-09-09.md` — WARNING `checking code files for non-ASCII
    characters` (`R/fs_bias_coverage.R`); NOTE `checking R code for possible problems`
    (`fs_plot_bias_coverage`'s nine NSE bindings); NOTE `checking HTML version of manual`
    (`no command 'tidy' found`, an environment note).
  - `quarto` and `pandoc` are not on this shell's PATH and must be prepended from
    `/usr/lib/rstudio/resources/app/bin/quarto/bin` (and `.../bin/tools/x86_64`) or the vignette
    rebuild fails during `R CMD build`, before any check result is produced.
- `R/fs_mr_inference_methods.R` is ASCII-clean (checked explicitly).

---

## Side issues, flagged not fixed

1. **`devtools::document()` aborts on `R/fs_bias_coverage.R:17`** — `` `r = se_mean / sd_emp` ``
   in an `@description` is parsed by roxygen as inline R. Falsified as pre-existing and
   unrelated (§F2). Same file as the pre-existing check WARNING.
2. **`forestsearch_version` does not pin a commit.** The `e1stud` bundle records `0.3.5`, as
   does HEAD, yet the same DGM call does not reproduce its trials (§Fe) while the DGM source
   files are untouched since it was built. Whatever the cause, a version string that several
   commits share cannot support an `identical()` comparison against committed work — which is
   what made both Fb and Fe unevaluable in their literal form. A commit hash in the meta would
   close this.
3. **GRF's empty-band asymmetry** (§V1): the identifier has no fallback, MR's re-selection does.
   The source names this as "a decision pending on its own." It did not bind in the pilot.
4. **`dmin.grf` units** (§V5): GRF's own eligibility floor is in DR-score units, the admission
   effect floor in log HR. A GRF campaign must state which floor it intends.
5. The **prevalence question from the Stage 0 report §6.1 is unchanged**: the template default
   `FS_S7_Z1Q = 0.25` gives 12.4%, while `summary_cert20.qmd` is a 31% (`z1q60`) study. Larry's
   decision; untouched here.

---

## What a DINA or GRF campaign still needs, and what it no longer needs

**No longer needed:** any `R/` change. Every field product now reaches both branches through
`mr_inference_args`, and the three formerly inert knobs now control rather than coincide, so a
run's meta means what it records.

**Still needed, all of them decisions rather than code:** the prevalence (`FS_S7_Z1Q`) and its
consequence for any `cert20`-derived comparator; the ε for `effMaxSG`; whether `dmin.grf` should
remain 0.0 given the units mismatch; and — for GRF specifically — whether the empty-band
asymmetry must be reconciled before the frontier path is used for a campaign at `effMaxSG`.
Whether DINA and GRF should *default* to the field constructions remains a separate decision and
was not taken here.

**No acceptance criteria were pre-registered, none are proposed, and this report makes no
recommendation. Report and wait.**
