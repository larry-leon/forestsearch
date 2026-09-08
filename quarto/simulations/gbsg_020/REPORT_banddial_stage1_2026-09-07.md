# REPORT — Identifier comparison on the J = 10 harm cells (campaign `banddial`): Stage 1 (Knobs, identities, smoke, projection)

**Task:** `dev/tasks/TASK_complement_variance_banddial_2026-09-07.md` (4b245516), Part B as instructed (three settings: `maxSG`, `minSG`, `effMaxSG` ε 0.30; two J = 10 harm cells; six cells). Stage 0 record 3dd238fc (Gate 0 PASS, all three arms). Winner-only and winner-floor excluded from every table, figure and report line.
**Date:** 2026-09-07. Executor: Claude Code, unattended per the pre-authorization. No `R/` change.

---

## GATE 1: PASS — knob-inert identity **exact** against nb20 arm A (122/122 shared non-timing columns, the nine p̂ / band / family columns included, worst relative difference 0.0, `truth` identical); smoke at the three settings on both cells: ε = 0.30 admits at least as many candidates as ε = 0.20 on every replicate and |Ĥ| is never smaller than arm A's; `maxSG` never smaller than ε = 0.30's; `minSG` never larger than arm A's and always ≥ 60; the nesting minSG ≤ ε 0.30 ≤ maxSG holds on every replicate of both cells with the same family K; fields finite, γ in range, bound identities ≤ 1.7e-16, p̂ valid; alignment 4/4 (identifier top row = gate map, no exact-N ties, observed |Ĥ| = the top row's N). Projection ≈ 2.7 h at 100 workers, inside the 5 h ceiling (7 h hard timeout). **Compute go: all six cells; none deferred at the start** (the driver re-projects each cell from the realized mean wall and defers, listing it, any cell whose projected finish would cross the ceiling).

## 1a — Template change (document-level, add-only; defaults and every committed focus reproduce exactly)

- Focus guard widened: `stopifnot(sg_focus %in% c("maxeffCons", "effMaxSG", "maxSG", "minSG"))` (`:313`), with the comment recording the aligned maps (Stage 0 0a). Nothing else keys on the value.
- The re-selection callout's band-size line adds "; informational under this focus, whose sort key has no band term" when `sg_focus` is `maxSG` / `minSG`. No recorder column changes.
- Campaign tags: `bdinert` (identity), `bdsmoke` (smoke), `banddial` (Stage 2). Stems: `fs_effMaxSG_…_z1q60_nb30_banddial`, `fs_maxSG_…_z1q60_banddial`, `fs_minSG_…_z1q60_banddial`.

## 1b — Identities and smoke (all PASS; `smoke_banddial_check.R`, `align_banddial.R`, session scratchpad; 5 workers, sequential renders of 5 replicates at the p30sg seeds, 56–70 s each)

**Knob-inert identity** (`bdinert`: `effMaxSG`, ε = 0.20, J = 10, HR 1.50 n500, sims 1–5) against the committed nb20 arm A bundle rows: **122 / 122** shared non-timing columns identical (worst relative difference 0.0; 136 columns in the bundle, the excluded ones being the five timing columns and the FB placeholders), `truth` `identical()`, meta `sg_focus = effMaxSG`, `effect_neighborhood = 0.20`, `er_jcuts = 10`; |Ĥ| 161, 69, 113, 180, 87 and p̂(Ĥ) 0.013, 0.186, 0.064, 0.037, 0.159 on both sides.

**Smoke** (`bdsmoke`, sims 1–5; same data as nb20 arm A on every replicate: `n_true` identical; 5/5 detected under every setting, as arm A):

| Cell | setting | \|Ĥ\| per rep (arm A ε 0.20 beside) | relation to arm A | band_n (arm A) | cons-qual | K | fields finite | γ | bound ids | p̂(Ĥ) | top-3 | fit+MR s/rep (light) | complement fits | sens / spec (mean) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 | effMaxSG ε 0.30 | 210, 91, 188, 254, 118 (161, 69, 113, 180, 87) | ≥ every rep; band ≥ every rep | 35, 8, 26, 108, 7 (8, 1, 14, 28, 3) | 118, 227, 48, 337, 89 | 1233–1303 | yes | in range | 1.7e-16 | 0.009, 0.107, 0.004, 0.009, 0.144 | 0.21–0.41 | 36.8 | 718–883 | 0.70 / 0.82 |
| HR 1.50 | maxSG | 334, 445, 248, 457, 282 | ≥ ε 0.30 every rep | 4, 1, 2, 6, 2 (ε 0.10, informational) | same | same | yes | in range | 1.7e-16 | 0.100, 0.109, 0.182, 0.719, 0.025 | 0.27–0.81 | 27.8 | 571–886 | 0.80 / 0.34 |
| HR 1.50 | minSG | 62, 63, 64, 62, 63 | ≤ arm A every rep; ≥ 60 | same | same | same | yes | in range | 6.9e-17 | 0.011, 0.002, 0.020, 0.004, 0.092 | 0.50–0.93 | 26.1 | 79–175 | 0.22 / 0.92 |
| HR 1.75 | effMaxSG ε 0.30 | 210, 172, 203, 180, 118 (173, 127, 167, 153, 87) | ≥ every rep; band ≥ every rep | 49, 17, 32, 21, 5 (15, 3, 15, 3, 3) | 244, 388, 128, 508, 169 | 1233–1303 | yes | in range | 1.1e-16 | 0.023, 0.087, 0.003, 0.120, 0.223 | 0.22–0.46 | 37.9 | 670–822 | 0.92 / 0.92 |
| HR 1.75 | maxSG | 436, 454, 248, 457, 352 | ≥ ε 0.30 every rep | 2, 1, 3, 1, 2 (ε 0.10, informational) | same | same | yes | in range | 8.3e-17 | 0.281, 0.643, 0.107, 0.891, 0.056 | 0.30–0.93 | 28.5 | 358–668 | 0.87 / 0.27 |
| HR 1.75 | minSG | 62, 61, 62, 61, 63 | ≤ arm A every rep; ≥ 60 | same | same | same | yes | in range | 1.1e-16 | 0.004, 0.006, 0.121, 0.001, 0.051 | 0.57–0.97 | 26.9 | 62–93 | 0.21 / 0.92 |

**Nesting** on the same replicates: HR 1.50 — minSG (62–64) ≤ ε 0.30 (91–254) ≤ maxSG (248–457) on every replicate, same family K; HR 1.75 — minSG (61–63) ≤ ε 0.30 (118–210) ≤ maxSG (248–457) on every replicate, same family K.

Early read, not a check: `maxSG` returns half to nine-tenths of the sample (248–457 of 500 on both cells; the largest consistency-qualifying candidate is a broad `er ≤ c` or age rule at HR 1.50), with specificity 0.13–0.53 — the size rule without a band has no reason to stop at the planted region; `minSG` sits on the floor (62–64, the smallest qualifying two-factor rule of ≥ 60) with sensitivity 0.2; ε = 0.30 lands between ε = 0.20 and `maxSG` (91–254; band 5–108 candidates).

**Alignment in numbers** (standalone fits, HR 1.50 n500, `return_reselection = TRUE`; the candidate table from `fs$grp.consistency$out_sg$result`):

| sim | setting | consistency-qualifying | in band | identifier top row (N, Pcons) = gate map (N) | exact-N ties at the pick | observed \|Ĥ\| = top-row N | Σp̂ = selection rate | p̂(Ĥ) / argmax p̂ | selected rule |
|---|---|---|---|---|---|---|---|---|---|
| 1 | maxSG | 247 | (5 at ε 0.10) | row 1 (454, 0.940) = row 1 (454): yes | 1 | yes | 1.000 = 1.000 | 0.592 / 0.592 | `{er ≤ 289}` |
| 1 | minSG | 247 | — | row 1 (62, 0.990) = row 1 (62): yes | 1 | yes | 1.000 = 1.000 | 0.004 / 0.858 | `{age ≤ 46} & {size ≤ 25}` |
| 1 | effMaxSG ε 0.30 | 247 | 43 | row 1 (195, 1.000) = row 1 (195): yes | 1 | yes | 1.000 = 1.000 | 0.023 / 0.067 | `{age ≤ 53} & !{pgr ≤ 7}` |
| 4 | minSG | 48 | — | row 1 (64, 0.990) = row 1 (64): yes | 1 | yes | 0.999 = 1.000 | 0.117 / 0.486 | `{er ≤ 16} & {age ≤ 46}` |

On every fit the identifier's top row is the gate map's pick with no exact-N tie (the Stage 0 residual never engaged) and the observed Ĥ has that N; `settings$reselection` reads `maxSG` / `minSG` / `effMaxSG` respectively. Caveat as before: standalone fits reproduce the harness up to the ambient RNG, so their picks need not coincide with the smoke bundles' on the same sim ids; the alignment is a self-consistency of each fit. Under `maxSG` on sim 1 the pick is `er ≤ 289` — 454 of 500 patients (naive HR 1.01, de-biased 0.98): the size rule selects the broadest candidate that still clears the harm screen and the consistency floor.

## 1c — Projection (100 workers)

Light-load fit+MR per replicate against nb20 arm A's smoke (36.3 / 37.7 s at HR 1.50 / 1.75, realized walls 30 / 31 min): ε 0.30 **36.8 / 37.9 s (×1.0)** — as at nb20, the band costs nothing (complement fits 670–883, +15–25% on arm A's); `maxSG` **27.8 / 28.5 s (×0.75)**; `minSG` **26.1 / 26.9 s (×0.7)** (complement fits 62–175: the smallest Ĥ has the largest complement, and few distinct winners). Scaling arm A's walls: ε 0.30 ≈ 30 + 31 min; `maxSG` ≈ 23 + 24 min; `minSG` ≈ 22 + 23 min; **≈ 2.6 h of renders + combines and gates ≈ 2.7 h**, inside the 5 h ceiling (7 h hard timeout). Order: ε 0.30 (HR 1.50 → 1.75), `maxSG` (1.50 → 1.75), `minSG` (1.50 → 1.75); two seed-disjoint batches then combine per cell; fail-fast per cell; Gate 2 per cell after each combine (stop-on-failure per cell, the next proceeds); the driver (`stage2_banddial_driver.sh`, session scratchpad) starts with a 31 min per-cell projection and re-projects each subsequent cell from the realized mean wall, deferring (listed) any cell whose projected finish would cross the ceiling; every render caps its own timeout at the 7 h hard stop; every render is the committed template driven by `FS_S7_*` env only; the save guard is live.
