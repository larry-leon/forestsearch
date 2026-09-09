# REPORT — Part D: the defaults become the recommended constructions

Date: 2026-09-08. Spec: `dev/tasks/TASK_cert20_2026-09-08.md`, Part D (classified *changes behaviour*; decided by Larry, F-2). Executor: Claude Code (Linux, unattended). **Result: Gate D PASS on all three sub-gates; Part D committed.**

## What changed

Three arguments of `fs_mr_inference()` now default to the finalized (F-3) constructions. `field_decompose` stays `FALSE` and **`ci_method` is untouched** — the open decision. The field block still runs only under `ci_method = "field"`; the gating line is unchanged, `R/fs_mr_inference.R`:

```r
  if (ci_method == "field") {
```

### D1 — `R/fs_mr_inference.R` signature

| argument | before | after |
|---|---|---|
| `return_reselection` | `FALSE` | `TRUE` |
| `field_complement` | `FALSE` | `TRUE` |
| `field_scale_complement` | `c("none", "selected")` | `c("selected", "none")` |
| `field_decompose` | `FALSE` | `FALSE` (unchanged) |
| `ci_method` | `c("ij", "wald", "field")` | unchanged (open decision) |

`match.arg(field_scale_complement)` therefore yields `"selected"`. Roxygen for the three now states that these are the recommended constructions and names the value that restores the prior behaviour; `devtools::document()` regenerated `man/fs_mr_inference.Rd` (NAMESPACE unchanged).

### D2 — `R/forestsearch_main.R` `.g_mr(...)` fallbacks

Before:

```r
        return_reselection = .g_mr(mr_inference_args$return_reselection, FALSE),
        field_complement = .g_mr(mr_inference_args$field_complement, FALSE),
        field_scale_complement = .g_mr(mr_inference_args$field_scale_complement, "none"),
        field_decompose = .g_mr(mr_inference_args$field_decompose, FALSE),
```

After:

```r
        return_reselection = .g_mr(mr_inference_args$return_reselection, TRUE),
        field_complement = .g_mr(mr_inference_args$field_complement, TRUE),
        field_scale_complement = .g_mr(mr_inference_args$field_scale_complement, "selected"),
        field_decompose = .g_mr(mr_inference_args$field_decompose, FALSE),
```

The surrounding comments were rewritten to say which value is now the gate's default and which restores the prior behaviour. Nothing else in the call changed: `ci_method` still falls back to `"ij"`, `include_complement` to `TRUE`, `ij_residual` to `"two_term"`.

### D3 — template and NEWS

`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`:

| knob | before | after |
|---|---|---|
| `FS_S7_FIELD_COMPLEMENT` | default `"FALSE"` | default `"TRUE"` |
| `FS_S7_FIELD_SCALEC` | default `"none"` | default `"selected"` |
| `return_reselection` | hard-coded `TRUE` | `mr_return_reselection`, `FS_S7_RETURN_RESEL` default `TRUE` |
| `FS_S7_FIELD_DECOMP` | default `"FALSE"` | unchanged |

Unset is now the recommendation for all three. `FS_S7_RETURN_RESEL` is **new**: the template previously hard-coded `return_reselection = TRUE`, so the pre-flip value was not reachable from the committed template by env alone and Gate D(iii) could not have been executed as written. The knob is default-inert (`unset` and `TRUE` both give `TRUE`) and is the only addition beyond a default flip in this part; flagged here rather than silently taken.

The knob echo gained `return_reselection=%s`; the status header now describes both complement knobs and says both default to the recommended constructions. `NEWS.md` gained one bullet under the development header, worded as the spec required.

### D4 — install and deparse

`devtools::install(dependencies = FALSE)` clean. Read back from the **installed** package:

```
  field_complement         = TRUE
  field_scale_complement   = c("selected", "none")
  return_reselection       = TRUE
  field_decompose          = FALSE
  ci_method                = c("ij", "wald", "field")
```

and, by walking the installed `forestsearch()` body to the single `fs_mr_inference(...)` call node (not by text grep, which the deparse line-wrap defeats):

```
  return_reselection       -> .g_mr(mr_inference_args$return_reselection, TRUE)
  field_complement         -> .g_mr(mr_inference_args$field_complement, TRUE)
  field_scale_complement   -> .g_mr(mr_inference_args$field_scale_complement, "selected")
  field_decompose          -> .g_mr(mr_inference_args$field_decompose, FALSE)
  ci_method                -> .g_mr(mr_inference_args$ci_method, "ij")
```

Installed matches source.

## Gate D — 5 replicates, effMaxSG ε 0.20, HR 1.50, n 500, z1q 0.60, J 10, seeds 8316951 + sim_id, sim_id 1–5

Three renders of the committed template, env only, `FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_WORKERS=5`, campaign tags `dflt_on` / `dflt_explicit` / `dflt_old`; walls 67 s / 68 s / 65 s. Timing columns excluded throughout: `fb_secs`, `fit_mr_secs`, `fld_H_secs`, `fld_H_uniform_secs`, `fld_Hc_secs`.

**(i) unset ≡ explicitly recommended — PASS.** `dflt_on` sets neither complement knob; `dflt_explicit` sets `FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_RETURN_RESEL=TRUE`. Both bundles carry the same 158 columns.

| check | value |
|---|---|
| non-timing columns `identical()` | **153 / 153** |
| `truth` `identical()` | TRUE |
| `meta` `identical()` (excluding `built_at`, `campaign_tag`) | TRUE |
| resolved `meta$field_complement` / `field_scale_complement` under *unset* | `TRUE` / `"selected"` |

**(ii) the recommended render ≡ committed `e1stud` rows 1–5 — PASS.** Comparator `results/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_e1stud_res_1_1000.rds` (same machine, `n_workers` 100 there vs 5 here; that is a `meta` field, not a results column).

| check | value |
|---|---|
| column sets equal (158 each, no side-only columns) | TRUE |
| complement / field-s / `joint_s` / scale / p̂ columns `identical()` | **58 / 58** |
| all common non-timing columns `identical()` | **153 / 153** |
| `truth` `identical()` | TRUE |

The 58 targeted columns are those matching `^fld_Hc_`, `^fld_joint`, `^p_hat_`, or `_scale`. The stronger 153/153 line is reported because it was free.

**(iii) explicitly OLD ≡ committed `p30sgnb20` rows 1–5 on pre-existing columns — PASS.** `dflt_old` sets `FS_S7_FIELD_COMPLEMENT=FALSE FS_S7_FIELD_SCALEC=none FS_S7_RETURN_RESEL=FALSE`. Comparator `results/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_p30sgnb20_res_1_1000.rds` (136 columns; it predates the `_s`, decompose and some p̂ recorder columns).

| check | value |
|---|---|
| common columns | 136 |
| of those, gated by the three flipped knobs (`^fld_Hc_`, `^fld_joint`, `^p_hat_`) | 37, excluded |
| pre-existing non-gated non-timing columns `identical()` | **95 / 95** |
| `truth` `identical()` | TRUE |
| the 37 gated columns are all-`NA` in the OLD render | **37 / 37** |

The 37 gated columns must be excluded rather than compared: `p30sgnb20` itself ran with `field_complement = TRUE` and `return_reselection = TRUE`, so it *has* values there, while an OLD render by construction does not. Their being uniformly `NA` in the OLD render — and not, say, partly populated — is the positive evidence that the knobs turned the blocks off cleanly. The harm-side and identification columns, which the complement block is documented not to touch, reproduce exactly.

**Verdict: PASS.** Part D committed. Gate bundles and their renders committed beside the campaign: `dflt_on`, `dflt_explicit`, `dflt_old` (`results/fs_effMaxSG_..._nb20_dflt_*_res_1_5.rds`, `dflt_*_batch_1_5.html`).

## Notes carried forward, not acted on

1. `.fs_apply_mr()` (`R/fs_mr_inference_methods.R`, the DINA and GRF branches) does **not** pass the three arguments at all, so it inherits the new `fs_mr_inference()` defaults automatically. That is consistent with the flip and needed no edit. Its practical effect on those branches is nil unless a caller also sets `ci_method = "field"` through `mr_inference_args`, since `.fs_apply_mr()` falls back to `"ij"` and the complement field block is gated on `"field"`; `return_reselection = TRUE` is add-only there.
2. `include_complement` keeps its `FALSE` signature default and its `TRUE` fallback in `forestsearch()`. Not in scope for this part; noted because its roxygen sentence "Default `FALSE`" sits immediately above the flipped block and is easy to misread as one of the three.
3. The `ci_method` default remains the open decision (which two-sided interval is primary in user-facing output). Untouched, as specified.
