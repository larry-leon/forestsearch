# REPORT — `effMaxSG` band 0.20 at 31% prevalence: Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_p30sg_nb20_2026-09-07.md` (095d8c4d; supersedes `TASK_p30sg_jcuts20_2026-09-07.md`, not run). Decisions S-1–S-4 as stated (arm A = HR 1.50 n500 and HR 1.75 n500 only; ε = 0.20; no `maxSG` arm unless an aligned gate map is found, and then report and wait; 100 workers). Standing convention: winner-only and winner-floor excluded from every table, figure and report line.
**Date:** 2026-09-07. No compute beyond a two-replicate family probe (0c); no `R/` changes. Source at tip 095d8c4d (R/ unchanged since the installed build 8fbb3bdc, forestsearch 0.3.5): `R/forestsearch_main.R`, `R/fs_mr_inference.R`, `R/fs_mr_inference_methods.R`, `R/subgroup_consistency_helpers.R`, `R/forestsearch_helpers.R`, `R/get_FSdata_helpers.R`; template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at 22c5f854.

---

## GATE 0: PASS — all three knobs are document-level (the template already pins `effect_neighborhood` and `conf.cont_jcuts` and forwards both; `return_reselection` is an existing add-only pass-through of `forestsearch()` reached through `mr_inference_args`), and the p̂ pass-through needs no `R/` change (the gate already exposes `reselection$p_hat` under `return_reselection = TRUE`). `maxSG` (0d): the gate has an aligned re-selection map for the pure size rule; per S-3 and Larry's instruction it is reported here and **not run** — no arm C in this task.

## 0a — The band on both sides, the grid, and the J = 20 cuts

**ε is one value read by both sides.** The template pins `effect_neighborhood <- 0.10` (`:484`, "package default; band half-width for effMaxSG/effMinSG. Pinned rather than inherited so the value is visible") and `selection_rule <- "neighborhood"` (`:473`), and passes both to `forestsearch()` (`:872–873`: `selection_rule = selection_rule, # governs *MaxSG/*MinSG band (search, and MR's de-biasing family)`; `effect_neighborhood = effect_neighborhood, # band width; inert unless the focus is effMaxSG/effMinSG`). Inside `forestsearch()` the identifier's `sort_subgroups()` (`R/subgroup_consistency_helpers.R:540–630`) calls

```r
in_band <- .compute_inclusion_band(hr_vec = hr_vec, n_vec = N_vec,
                                   selection_rule = selection_rule, effect_neighborhood = effect_neighborhood)
ord <- if (sg_focus == "hrMaxSG") order(-in_band, -N_vec, -Pcons_vec, -hr_vec, K_vec) else ...
```

and the gate call (`R/forestsearch_main.R:3388–3402`) forwards the same formal — `effect_neighborhood = effect_neighborhood, selection_rule = .g_mr(mr_inference_args$selection_rule, selection_rule)` — which `fs_mr_inference()` hands to `.fs_mr_select(bs, .zcons(bs), sz, pass, reselection, effect_neighborhood, selection_rule, log_scale)` (`R/fs_mr_inference.R:599`, and the field's `sel_one` at `:808`), where (`:132–149`)

```r
.fs_mr_select <- function(beta, zcons, sizes, passers, rule, nbhd, selection_rule = "neighborhood", log_scale = TRUE) {
  .inband <- function() { eff <- if (log_scale) exp(beta[passers]) else beta[passers]; sz <- sizes[passers]
    ib <- .compute_inclusion_band(hr_vec = eff, n_vec = sz, selection_rule = selection_rule, effect_neighborhood = nbhd) == 1L
```

— the one helper (`R/subgroup_consistency_helpers.R:781–785`: `hr_floor <- (1 - effect_neighborhood) * hr_max; as.integer(!is.na(hr_vec) & hr_vec >= hr_floor)`), natural HR on both sides, multiplicative, width ε. **At ε = 0.20 the floor is 0.8·max hr — a constant log(0.8) = −0.223 log-HR below the maximizer (≈ three-quarters of an SE at n = 500), versus log(0.9) = −0.105 at ε = 0.10.** Nothing in `R/` keys on the value 0.10 except the formal's default (`forestsearch_main.R:1247`, validated in [0, 1) at `:1510–1513`).

**The grid.** `fs_conf.cont_jcuts <- list(er = 10)` (`:510`, "J-quantile grid on raw ER") is passed as `conf.cont_jcuts = fs_conf.cont_jcuts` (`:885`) into the consistency `method_args`; `cut_var_jq(x, J)` (`R/get_FSdata_helpers.R:41–48`) emits `x <= qj(x, k, J + 1)` for k = 1..J with `qj(x, k, J) = quantile(x, k/J)` (`:55`) — the cuts are the k/(J+1) **sample** quantiles of the analysis data (n = 500 or 1,000), so they vary by replicate around the population values. Population values (the super-population of 100,000 at seed 8316951, `z1_quantile = 0.60`; the raw GBSG n = 686 in brackets where different):

- J = 10: 0, 3, 9, 16, 29, 44, 68, 100, 174, 298 (GBSG: 0, 3, 9, 17, 30, 44, 69.8, 100, 173.5, 293.7 — the template's documented grid).
- J = 20: 0, 0, 1, 4, 7, 10, 14, 19, 25, 32, 40, **50**, **64**, 79, 94, 122, 159, 206, 288, 394 (GBSG: …, 40.8, 50.4, 64, 78.7, …).
- The planted boundary `er ≤ 59` (= q₀.₆₀, super-population share 0.6016) lies between the k = 12 cut (12/21 = 0.571 quantile: er ≤ 50, share 0.5715) and the k = 13 cut (0.619: er ≤ 64, share 0.6229). **The best representable sub-region `er ≤ 50` holds 0.950 of the truth's er-mass**; `er ≤ 59` is still not a candidate. On the J = 10 grid the neighbours are 44 (0.455) and 68 (0.636). Per-replicate sample cuts scatter around these: sim 2 (n = 500) has the k = 12/13 cuts at 58.1/68, sim 3 at 61/72.7 — on some replicates the sample grid lands within 1–2 units of 59.

## 0b — `return_reselection` and the p̂ recording

`forestsearch()` already forwards it as an add-only pass-through (`R/forestsearch_main.R:3410–3412`: `return_reselection = .g_mr(mr_inference_args$return_reselection, FALSE)`, "defaults reproduce prior output"), and the gate (`R/fs_mr_inference.R:957–961`) attaches

```r
if (isTRUE(return_reselection)) {
  p_hat <- tabulate(winner[!is.na(winner)], nbins = length(asm$names)) / draws
  names(p_hat) <- asm$names
  out$reselection <- list(winner = winner, p_hat = p_hat)
}
```

after `out` is fully built — the roxygen (`:290–297`) states "nothing in the arithmetic depends on this switch", and `sum(p_hat) == selection_rate`. The template builds `mr_inference_args` at `:562–566` (`ci_method, draws, include_complement, confirm_rule, field_uniform, field_complement, ij_residual`) and reads the gate's return as `g <- fs.est$mr_inference` (`:908`), recording `g$selected_label`, `g$n_selected` (`:920–921`). So the change is: add `return_reselection = TRUE` to that list, and record from `g$reselection$p_hat` the selected candidate's frequency (`p_hat[selected_label]`; labels coincide, as the Stage 1 alignment fits of p30sg used), the sum, and the three largest with their labels; plus `g$n_family` (the gate's kept family K) and, from the identifier's own table `fs.est$grp.consistency$out_sg$result` (columns `hr, N, Pcons, K` on the consistency-qualifying candidates), the count of consistency-qualifying candidates and the observed-effect band count `sum(hr >= (1 − ε)·max hr)` — the band the identifier's sort key saw. All new columns; no existing column's value changes. **Document-level; no `R/` change.**

**Naming.** The template's campaign guard admits alphanumerics only (`stopifnot(grepl("^[A-Za-z0-9]+$", campaign_tag))`, `:372`, a collision guard on the stem grammar), so the task's campaign names are written `p30sgnb20` (arm A) and `p30sgnb20j20` (arm B); the stem additionally gains `_nb20` and `_j20` tags when the knobs are non-default (empty at defaults, so every committed stem is unchanged): `fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_p30sgnb20`, `…_z1q60_nb20_j20_p30sgnb20j20`. Both knobs join the bundle meta and the combine poolability gate.

## 0c — Cost anchors

**p30sg (Gate 2 record, 100 workers, J = 10, ε = 0.10):** cell walls 26 / 30 / 31 / 39 min (HR 1.00 / 1.50 / 1.75 at n = 500; HR 1.00 at n = 1000); fit+MR 61–74 s/rep at n = 500 under load, 100 s at n = 1000; complement fits 457–578 per replicate.

**Family probe (standalone fits, HR 1.00 n500, sims 2 and 3, sequential, light load; the template's full configuration with `return_reselection = TRUE`):**

| sim | J | ε | fit+MR s | gate family K | consistency-qualifying | in band | selected (N) | p̂(Ĥ) / top-3 mass / Σp̂ | complement fits |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 10 | 0.10 | 25.0 | 1195 | 21 | 2 | !{er ≤ 14} & {er ≤ 49} (91) | 0.149 / 0.372 / 0.996 | 560 |
| 2 | 10 | 0.20 | 24.9 | 1195 | 21 | 8 | !{er ≤ 8} & {er ≤ 49} (134) | 0.160 / 0.319 / 0.996 | 705 |
| 2 | 20 | 0.10 | 35.2 | 2012 | 46 | 2 | !{er ≤ 9} & {er ≤ 25} (73) | 0.268 / 0.424 / 0.999 | 777 |
| 2 | 20 | 0.20 | 35.6 | 2012 | 46 | 6 | !{er ≤ 9} & {er ≤ 58} (142) | 0.095 / 0.334 / 0.999 | 959 |
| 3 | 10 | 0.10 | 26.3 | 1201 | 57 | 2 | {pgr ≤ 7} & !{size ≤ 27} (76) | 0.079 / 0.388 / 1.000 | 421 |
| 3 | 10 | 0.20 | 26.7 | 1201 | 57 | 7 | {er ≤ 3} (93) | 0.031 / 0.328 / 1.000 | 542 |
| 3 | 20 | 0.10 | 37.0 | 1913 | 90 | 3 | {er ≤ 4} & {size ≤ 35} (76) | 0.071 / 0.360 / 1.000 | 572 |
| 3 | 20 | 0.20 | 37.8 | 1913 | 90 | 14 | {er ≤ 4} (101) | 0.027 / 0.315 / 1.000 | 732 |

**Reading for the projection.** J = 20 enlarges the gate's family by ×1.6–1.7 (1195–1201 → 1913–2012) and the consistency-qualifying set by ×1.6–2.2, and the per-replicate cost by ×1.4 (25–27 → 35–37 s, light load); ε = 0.20 admits 6–8 candidates to the band instead of 2–3 and raises the complement-fit count by 25–30% at no measurable fit-time cost. Scaling the p30sg walls: **arm A ≈ 1.0–1.1 h** (30 + 31 min, plus the complement-fit allowance); **arm B ≈ 4.0–4.8 h** (n = 500 cells 26/30/31 min × 1.4 ≈ 36/42/43 min; n = 1000 cells 39 min × 1.4 ≈ 55 min each, the harm 1.5 n1000 cell taken as the null n1000 cell's cost plus the harm cells' complement allowance). Total ≈ 5–6 h, inside the 7 h ceiling (9 h hard timeout); Stage 1c measures under the smoke and Stage 2's driver re-projects arm B from arm A's realized walls before each arm B cell, deferring any cell whose projected finish would cross the ceiling.

## 0d — `maxSG` (decision S-3): the gate's map, quoted

The identifier's `maxSG` sorts the consistency-qualifying table by `(-N, -Pcons, K)` (`sort_subgroups()`, `R/subgroup_consistency_helpers.R:583–586`); `.fs_mr_reselection_from_focus()` maps `maxSG → "maxSG"` (`R/fs_mr_inference_methods.R:94`), `.fs_admission_applies("maxSG", "consistency")` returns `effect = TRUE, consistency = TRUE` (`R/forestsearch_helpers.R:2268–2287`, the same two floors as `effMaxSG`), and the gate's rule is `maxSG = passers[which.max(sizes[passers])]` (`R/fs_mr_inference.R:171`) — the largest N among the screened candidates, the same functional as the identifier's with the same residual on exact N ties (`-Pcons` tie-break vs first in family order) recorded for `effMaxSG` at p30sg Stage 0. `stop_threshold` is reset to NULL for `maxSG` as for `hrMaxSG` (`forestsearch_main.R:1630`), already the template's setting. **An aligned map exists.** Per S-3 and the instruction for this run, arm C is **not run**: it would need only the template's focus guard widened (`stopifnot(sg_focus %in% c("maxeffCons", "effMaxSG"))`, `:305`) and would be the four J = 10 cells under `FS_S7_FOCUS=maxSG`; reported here for Larry's decision.

## Proposed Stage 1 template change (document-level, add-only; defaults reproduce `p30sg` exactly)

1. `effect_neighborhood <- .env_num("FS_S7_NBHD", 0.10)` and `er_jcuts <- .env_int("FS_S7_ER_JCUTS", 10L)` read next to the focus knob (before the stem, which they tag when non-default); the later literal pin at `:484` becomes a pointer; `fs_conf.cont_jcuts <- list(er = er_jcuts)`.
2. `return_reselection = TRUE` in `mr_inference_args`; recorder columns `n_family, n_cons_qual, band_n, p_hat_H, p_hat_sum, p_hat_top1..3, p_hat_top_labels`; a render-side "Re-selection regime" callout after the subgroup-size callout (renders nothing on older bundles).
3. Knobs echo line gains `nbhd=`, `er_jcuts=`; both join the bundle meta and the combine poolability gate.

Knob-inert identity (Stage 1b): all non-p̂ columns shared with the `p30sg` bundle identical (≤ 1e-12) on the 5 replicates of HR 1.00 n500 at the p30sg seeds; the new columns are add-only and have no p30sg counterpart.
