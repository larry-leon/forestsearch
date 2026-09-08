# REPORT — Identifier comparison on the J = 10 harm cells (campaign `banddial`): Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_complement_variance_banddial_2026-09-07.md` (4b245516), Part B, **as instructed for this run**: the two J = 10 harm cells (HR 1.50 n500, HR 1.75 n500; 31% prevalence, `FS_S7_Z1Q=0.60`, seeds 8316951 + sim_id, sim_id 1–2,000) under **three** settings — `maxSG` and `minSG` at their own definitions, and `effMaxSG` at ε = 0.30 — six cells; `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `ij_residual = "two_term"`, `return_reselection = TRUE`, FB none, 100 workers; compute pre-approved at Gate 1 if the projection is ≤ 5 h wall (hard timeout 7 h), cells beyond the ceiling deferred and listed. Decisions T-2 (settings as instructed; no ε = 0.40) and T-5 (ceiling as instructed) recorded. Standing convention: winner-only and winner-floor excluded from every table, figure and report line. Parts A (fe5eeaf5) and C (b49be1a6) precede this record.
**Date:** 2026-09-07. No compute; no `R/` change. Tree at b49be1a6: `R/` last changed at 8fbb3bdc (2026-09-06, before every nb20 record); installed forestsearch 0.3.5, whose `sort_subgroups()`, `.fs_mr_select()` and `.fs_mr_reselection_from_focus()` bodies carry the lines quoted below (checked by `deparse()` of the installed functions).

---

## GATE 0: PASS for all three arms — every map quoted in the task document still holds on the current tree; no arm stopped.

## 0a — The sort keys and the gate maps, re-quoted from the current tree

**Identifier** (`sort_subgroups()`, `R/subgroup_consistency_helpers.R:578–584` and `:587–614`), over the consistency-qualifying table (every row cleared `pconsistency.threshold` = 0.90 and the effect screen):

```r
if (sg_focus == "maxSG") { data.table::setorder(result_new, -N, -Pcons, K); return(result_new) }
if (sg_focus == "minSG") { data.table::setorder(result_new,  N, -Pcons, K); return(result_new) }
if (sg_focus %in% c("hrMaxSG", "hrMinSG")) {          # effMaxSG / effMinSG normalize to these
  in_band <- .compute_inclusion_band(hr_vec = hr_vec, n_vec = N_vec,
                                     selection_rule = selection_rule, effect_neighborhood = effect_neighborhood)
  ord <- if (sg_focus == "hrMaxSG") order(-in_band, -N_vec, -Pcons_vec, -hr_vec, K_vec)
         else                        order(-in_band,  N_vec, -Pcons_vec, -hr_vec, K_vec)
```

with `.compute_inclusion_band()` (`:778–800`) under `selection_rule = "neighborhood"`: `hr_floor <- (1 - effect_neighborhood) * hr_max; as.integer(!is.na(hr_vec) & hr_vec >= hr_floor)` — natural HR, multiplicative, width ε. **At ε = 0.30 the floor is 0.7 × max hr**, log(0.7) = −0.357 log-HR below the maximizer (≈ 1.2 naive SEs at n = 500), against −0.223 at ε = 0.20 and −0.105 at ε = 0.10.

**Gate** (`.fs_mr_select()`, `R/fs_mr_inference.R:132–176`), per draw over the passers (the admission set: both floors, `.fs_admission_applies("maxSG" | "minSG" | "effMaxSG", "consistency")` = `c(effect = TRUE, consistency = TRUE)`, `R/forestsearch_helpers.R:2268–2285`):

```r
pick <- switch(rule,
  maxSG    = passers[which.max(sizes[passers])],
  minSG    = passers[which.min(sizes[passers])],
  effMaxSG = { b <- .inband(); b[which.max(sizes[b])] },   # .inband(): the same .compute_inclusion_band() on exp(beta[passers])
```

`.fs_mr_reselection_from_focus()` (`R/fs_mr_inference_methods.R:85–106`) maps `maxSG → "maxSG"`, `minSG → "minSG"`, `effMaxSG → "effMaxSG"` (identity), so `forestsearch()` derives the gate's rule from `sg_focus` with no template setting (`forestsearch_main.R:3399`). `fs_focus_tag("consistency", ·)` returns `maxSG` / `minSG` / `effMaxSG`, so the stems are `fs_maxSG_…`, `fs_minSG_…`, `fs_effMaxSG_…_nb30_…`.

**Alignment, functional by functional.** `maxSG`: largest N among the screened candidates on both sides; residual on exact-N ties only (identifier `-Pcons` then `K`; gate `which.max` = first in family order). `minSG`: smallest N, the same residual. `effMaxSG`: largest N within the band on both sides, the residual as at p30sg Stage 0 (no exact-N tie in the band has yet engaged it: nb20 Stage 1 3/3, Stage 1 of this task below). Note the floors: `n.min = NULL` in the template resolves to the adaptive floor max(60, ⌈0.10 n⌉) = **60 at n = 500** (`forestsearch()` docs, `:466–477`), with `d0.min = d1.min = 10` events, so `minSG` picks the smallest consistency-qualifying candidate of at least 60 patients.

**The empty-band fallback** (`.fs_mr_select()`, `:165`): `if (!any(ib)) ib <- rep(TRUE, length(passers))` — a draw whose perturbed effects empty the band falls back to all passers, so `effMaxSG` degrades to `maxSG` on that draw rather than losing it; on the identifier an all-zero band is harmless because `-in_band` is the leading key and `-N` breaks the tie — the same pick. Unchanged. (The identifier's other fallback, `all(is.na(hr_vec))` → `(-Pcons, -hr, K)`, cannot engage on these cells: every consistency-qualifying row carries a finite hr.)

**`stop_threshold`** (`forestsearch_main.R:1629–1663`): reset to NULL for `maxSG` / `minSG` as for `hrMaxSG`, without a warning because the template passes `stop_threshold = NULL` explicitly (`:508`) — already the setting under `effMaxSG`. `.validate_selection_rule()` (`subgroup_consistency_helpers.R:743–767`) accepts `selection_rule = "neighborhood"` for every focus and validates ε in [0, 1) — 0.30 admitted, 0.10 (the default, inert) under `maxSG` / `minSG`.

## 0b — Forwarding of ε = 0.30 and the tag guard

`effect_neighborhood <- .env_num("FS_S7_NBHD", 0.10)` is validated in the template (`:321–322`), passed to `forestsearch()` (`:919`), which forwards it to the gate as `nbhd` (nb20 Stage 0 0a; `forestsearch_main.R:3388–3402`) and to the identifier's sort — one value, both sides. The stem gains `_nb30` (`nbhd_tag`, non-default only). The campaign name `banddial` is alphanumeric: it passes the guard as it stood and as widened at Part C (`^[A-Za-z0-9_]+$`); no underscore is needed because the three settings are distinguished by the focus tag and the `_nb30` tag: `fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb30_banddial`, `fs_maxSG_…_z1q60_banddial`, `fs_minSG_…_z1q60_banddial` — three stems per HR, six in all, none colliding with a committed stem.

## 0c — What the template needs (Stage 1, document-level, add-only)

1. The focus guard `stopifnot(sg_focus %in% c("maxeffCons", "effMaxSG"))` (`:305` before the widening, `:313` after; nb20 side issue 3) widened to admit `maxSG` and `minSG`. Nothing else in the template keys on the value: `sg_focus` is read into the stem (via `fs_focus_tag()`), the knobs echo, the `forestsearch()` call, the bundle meta and the poolability gate (grep of every use, `:296–1515`).
2. The `band_n` recorder column (`:985–993`) counts the ε-band on the observed effects with the knob's value; under `maxSG` / `minSG` (ε at its inert default 0.10) it is informational — the sort key has no band term. The re-selection callout says so under those foci (one conditional phrase added). No column changes value under any focus already run.

## 0d — Cost anchors

nb20 arm A (the same two cells at ε = 0.20, J = 10, 100 workers): walls **30 / 31 min**, fit+MR 72–75 s/rep under load (36–38 s light), complement fits 580–615 per replicate. The band costs nothing (nb20 Stage 1c: ×1.0); what varies across the three settings is the number of distinct complements the field has to fit (`fld_Hc_nfit`) and the size of the selected complement. Prior: 6 × ~31 min ≈ **3.1 h**, plus a complement-fit allowance — inside the 5 h ceiling; Stage 1c measures under the smoke and the driver re-projects each cell from the realized mean wall, deferring (listed) any cell whose projected finish would cross the ceiling; every render caps its timeout at the 7 h hard stop.

## Proposed Stage 1 (per the instruction)

Knob-inert identity against nb20 arm A (campaign `bdinert`: `effMaxSG`, ε = 0.20, J = 10, HR 1.50 n500, sims 1–5 at the p30sg seeds; every shared non-timing column, the p̂ columns included since arm A carries them); smoke at the three settings on both cells (campaign `bdsmoke`, sims 1–5): band at ε = 0.30 ≥ band at ε = 0.20 on every replicate and |Ĥ| never smaller than arm A's; `maxSG` |Ĥ| ≥ ε = 0.30's; `minSG` |Ĥ| ≤ arm A's and ≥ 60; fields finite; γ in range; bound identities; p̂ valid; the nesting minSG ≤ ε 0.30 ≤ maxSG on the same replicates; alignment on one replicate per setting (standalone fits); projection at 100 workers. Gate 1 = compute go if all pass and the projection is ≤ 5 h.
