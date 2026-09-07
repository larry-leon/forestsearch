# REPORT — Identifier variant `effMaxSG` on the prevalence-30% cells: Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_p30_effMaxSG_2026-09-06.md` (70bb63a3). Q1–Q4 at defaults (winner-only evaluated; winner-floor excluded from every table, figure and report line).
**Date:** 2026-09-06. No compute; no `R/` changes. Source at tip 27f78b12: `R/forestsearch_main.R`, `R/subgroup_consistency_helpers.R`, `R/subgroup_consistency_main.R`, `R/fs_mr_inference.R`, `R/fs_mr_inference_methods.R`, `R/forestsearch_helpers.R`; template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`.

---

## GATE 0: PASS — the gate's re-selection map under `effMaxSG` (and the field's `sel_one`) is the identifier's selection functional: same effect ordering, same neighbourhood rule and width, same size rule, consistency thresholds held fixed. One residual difference exists only on exact ties in size among in-band candidates (below); it is quantified in Stage 1's alignment-in-numbers check.

## 0a — The identifier's `effMaxSG` selection and the gate's re-selection, quoted

**Normalisation.** `effMaxSG` is the GLM-natural alias of the canonical `hrMaxSG` (`.normalize_sg_focus("effMaxSG")` → `"hrMaxSG"`, `R/forestsearch_main.R:1464–1471`); `fs_focus_tag("consistency", "effMaxSG")` → `"effMaxSG"`, so the stem tag is the alias as typed.

**The identifier (`sort_subgroups()`, `R/subgroup_consistency_helpers.R:540–630`).** Rows reaching it have already cleared the effect screen (`hr.subgroups$HR >= hr.threshold`, `subgroup_consistency_main.R:545`) and the consistency screen (`evaluate_subgroup_consistency()` returns NULL below `pconsistency.threshold`, per the `maxeffCons` comment at `:565–575`). For `hrMaxSG`:

```r
hr_vec <- as.numeric(result_new$hr)          # natural HR scale (effect_log_scale = FALSE)
N_vec  <- as.numeric(result_new$N); Pcons_vec <- as.numeric(result_new$Pcons); K_vec <- as.numeric(result_new$K)
in_band <- .compute_inclusion_band(hr_vec = hr_vec, n_vec = N_vec,
                                   selection_rule = selection_rule, effect_neighborhood = effect_neighborhood)
ord <- if (sg_focus == "hrMaxSG") order(-in_band, -N_vec, -Pcons_vec, -hr_vec, K_vec) else ...
```

with (`:778–800`) `in_nbhd = hr_vec >= (1 − effect_neighborhood) * max(hr_vec, na.rm = TRUE)` under `selection_rule = "neighborhood"`. So: **among consistency-qualifying candidates, those within 10% of the maximal HR form the band; the largest N in the band wins; ties in N break on Pcons, then hr, then K.** The documented rule (`forestsearch_main.R:515–517, 571`): "Among candidates with effect size within `effect_neighborhood` of the maximum, the largest N."

**The neighbourhood parameter.** `effect_neighborhood`, a formal of `forestsearch()` with default `0.10` (`:1247`; validated in [0, 1) at `:1510–1513` for `hrMaxSG`); `selection_rule` default `"neighborhood"`. The template pins both explicitly — `selection_rule <- "neighborhood"` (`:467`), `effect_neighborhood <- 0.10` (`:478`, "package default; … pinned rather than inherited so the value is visible; inert under maxeffCons") — and passes them to `forestsearch()` (`:866–867`). **Q2 default:** the identifier's own default 0.10, as found.

**Search-stage settings under `hrMaxSG`.** `stop_threshold` is reset to NULL for `hrMaxSG` (`:1628–1650`: "Neighborhood-based selection requires evaluating all candidates to determine the effect-size neighborhood") — the template already passes NULL (`:471`), so nothing changes; `max_subgroups_search` keeps its default `Inf` (`:1269`; the template does not set it; the `< 30` advisory at `:2203–2211` is inert); `use_twostage = TRUE` (the consistency-evaluation optimisation) and the consistency thresholds `hr.threshold 0.90`, `hr.consistency 0.80`, `pconsistency 0.90` are unchanged. The `maxeff` override block (`:1530–1545`: floors off, cap → Inf, minp → 0) applies **only** to `sg_focus = "maxeff"`, not to `hrMaxSG`; so `effMaxSG` runs the same screened family as `maxeffCons` and differs from it *only in the sort key* (`maxeffCons`: `setorder(-hr, K)`, `:559–563`).

**The gate's re-selection map (`fs_mr_inference()`).** `forestsearch()` derives the rule from the focus — `reselection = .g_mr(mr_inference_args$reselection, .fs_mr_reselection_from_focus(sg_focus, engine = "consistency"))` (`:3396–3399`) with `hrMaxSG → "effMaxSG"` (`R/fs_mr_inference_methods.R:85–110`) — and forwards `effect_neighborhood = effect_neighborhood` and `selection_rule = .g_mr(mr_inference_args$selection_rule, selection_rule)` (`:3401–3402`), i.e. **the identifier's own 0.10 and "neighborhood"**. The admission set is resolved once by `.fs_resolve_admission()` (`R/forestsearch_helpers.R:2333`) and carried to the gate; `.fs_admission_applies("effMaxSG", "consistency")` → `effect = TRUE, consistency = TRUE` — the same two floors as under `maxeffCons` (`t_g <- pmax(effect_floor, c_cons + z·σ_D)`, `fs_mr_inference.R:496–500`), so **the consistency thresholds are held fixed and the domain is unchanged between the two identifiers**. Per draw (`:537–541`):

```r
s <- .fs_mr_select(bs, .zcons(bs), sz, pass, reselection, effect_neighborhood, selection_rule, log_scale)
```

and in `.fs_mr_select()` (`:132–180`):

```r
.inband <- function() { eff <- exp(beta[passers]); sz <- sizes[passers]
  ib <- .compute_inclusion_band(hr_vec = eff, n_vec = sz, selection_rule = selection_rule, effect_neighborhood = nbhd) == 1L
  if (!any(ib)) ib <- rep(TRUE, length(passers)); passers[ib] }
pick <- switch(rule, ..., effMaxSG = { b <- .inband(); b[which.max(sizes[b])] }, ...)
```

— the **same** `.compute_inclusion_band()` helper (natural HR scale via `exp(beta)`, the same `(1 − 0.10)·max` floor, the same "neighborhood" rule), then the largest size in the band. The empty-band fallback cannot trigger (the argmax is always in its own band). **The field's `sel_one`** (`:743–748`) calls the same `.fs_mr_select(v, .zcons(v), sz, pass, reselection, effect_neighborhood, selection_rule, log_scale)` for every outer draw `v = w + ζ*_r` and every inner draw `v + ζ'_j`; `fast <- identical(reselection, "maxeff") && is.null(t_g)` is FALSE, so the per-draw path runs — exactly as it already does under `maxeffCons` (whose `t_g` is non-NULL).

**Same functional: yes, with one residual on exact ties.** Effect ordering (HR, natural scale), band (`≥ 0.9·max`), size rule (largest N), domain (screened family; floors fixed) all coincide. The identifier breaks *exact ties in N among in-band candidates* on `−Pcons, −hr, K`; the gate's `which.max(sizes[b])` takes the first such candidate in family order. Sizes are integers and do not change across draws, so the event "two in-band candidates with identical N" is possible; Pcons is not available per draw, so the gate cannot mirror that tie-break. This is a property of the existing gate under `effMaxSG` (not introduced here) and is measure-zero for the *effect* ordering; Stage 1b's check — the gate's map applied to the unperturbed effect vector reproduces the observed Ĥ — is the test in numbers, and the frequency of in-band size ties on the observed data will be reported. If it reproduces Ĥ on the smoke replicates, the correction runs as aligned; an R/ mirror of the Pcons tie-break would be a separate proposal for Larry.

## 0b — The template's focus lines and the proposed knob

Quoted: `sg_focus <- "maxeffCons"` (`:298`, with the comment that it and `subgroup_method` are "the identifier" and feed the stem); `method_tag <- if (identical(subgroup_method, "consistency")) "fs" else subgroup_method` (`:330`); `focus_tag <- forestsearch::fs_focus_tag(subgroup_method, sg_focus)` (`:351`) with the alias NOTE (`:352–355`); the stem `sprintf("%s_%s_fb_mr_field_m1_h%03d_knoise%d_n%d%s_%s", method_tag, focus_tag, …)` (`:371–373`), so the focus already names the stem (`fs_maxeffCons_…` → `fs_effMaxSG_…`); `sg_focus`/`focus_tag` are in the batch meta (`:1310–1311`) and the pooled meta (`:1422–1423`); `subgroup_method` is in the poolability gate, `sg_focus` is not (different foci already cannot pool because their stems differ).

Proposed (add-only; default reproduces `p30` exactly), replacing `:298`:

```r
# Identifier focus knob (TASK_p30_effMaxSG_2026-09-06): "maxeffCons" (the
# consistency-screened effect argmax; every committed campaign) or "effMaxSG"
# (the largest subgroup within effect_neighborhood = 0.10 of the maximal
# effect among the same screened candidates).  forestsearch() derives MR's
# re-selection rule from sg_focus, so the gate re-selects under the same
# functional (Stage 0 record); the stem follows focus_tag as it always has.
sg_focus <- .env_chr("FS_S7_FOCUS", "maxeffCons")
stopifnot(sg_focus %in% c("maxeffCons", "effMaxSG"))
```

The knobs echo line gains `focus=%s`; `sg_focus` joins the combine poolability gate (belt and braces; the stems already separate). Nothing else changes: `selection_rule`, `effect_neighborhood`, the thresholds and every MR knob are already explicit and forwarded.

## 0c — Cost anchors (`REPORT_p30_gate2_2026-09-06.md`)

`p30` at 100 workers: 21 / 23 / 24 / 36 min per 2,000-replicate cell (HR 1.00 / 1.50 / 1.75 at n = 500; HR 1.00 at n = 1000); fit+MR 47.5–55.4 s per replicate at n = 500, 92.5 s at n = 1000; complement fits 355–419 per replicate. Under `effMaxSG` the search evaluates the same screened family (no early stop already; cap Inf), and the field already takes the per-draw `sel_one` path, so the selection cost is unchanged; `.fs_mr_select`'s band computation per draw adds a vector max and a comparison (negligible). A larger Ĥ means a smaller complement and more distinct complements among near-tied candidates — the complement-fit count may rise. **Projection: ≈ 1.8–2.0 h for the four cells** (the p30 walls, 1 h 46 m, plus a small allowance), inside the 3 h ceiling.
