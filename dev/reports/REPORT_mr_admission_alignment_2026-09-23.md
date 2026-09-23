# REPORT — MR admission aligned with the screen's rounded rule; GBSG before / after

Task: `dev/tasks/TASK_mr_admission_alignment_2026-09-23.md` (commit `ba595f4b`).

## Result

- **The change bites, but only a little.** On the GBSG forest-search fit, the corrected HR moves from 1.2806 to
  1.2762 (−0.0044), and the field lower bound from 0.6365 to 0.6336 (−0.0030). The Bonferroni lower bound moves
  from 0.5099 to 0.5087, and the selection bias rises by 0.0037 on the log scale. The region, the unadjusted HR and
  the 1,744-candidate family are unchanged.
- **Draw level.** The admitted set differs on **2,963 / 5,000 draws (0.593)**, but the re-selected winner differs on
  only **26 / 5,000 (0.0052)**. On the observed fit, **one** candidate lies in the band [z_eff, z_exact):
  `{pgr <= 0} & {er <= 9}`, N 71, HR 1.683928, closed-form Pcons 0.898655. That matches the expectation.
- **GRF and DINA:** the MR output objects are identical before and after, every element (wall-clock timings excluded).

## Pins

- HEAD at the start: `ba595f4b` (the task-doc commit, on top of `ac860c6d`), on `feature/glm-extension`.
- forestsearch 0.3.5.9000. The installed build dates from 2026-09-23 02:32 UTC, but every fit here ran on the working
  tree via `devtools::load_all()`, so the installed build was not used.
- R 4.6.1 (2026-06-24), x86_64-pc-linux-gnu. `forestsearch()` ran on 102 multisession workers.
- **Seed: `seedit = 8316951` in all seven fits.** MR's `seed` defaults to `seedit`, so the main multiplier stream
  (`set.seed(8316951)`) and the field stream (`8316951 + 900000`) are the same draws before and after.
  - Verified: the FS fit's captured family and analysis frame are `identical()` before and after.
  - Re-running MR's re-selection loop on the captured inputs reproduces the stored before and after winner vectors
    exactly (Step 4). The draws are therefore shared.
- Untracked files already present at Step 0, none staged:
  - the two `quarto/simulations/actg175/binary_020/mr_or_harm/…_d5000/` directories;
  - `smoke_redes.html` and `smoke_relaunch.html`;
  - `quarto/simulations/gbsg_020/scripts_dinamr/logs/nullmr_findings.err`.
- Scratch payloads, not in the repo and not committed:
  - `~/Downloads/mr_admission_align_{before,after}_{fs,grf,dina}.rds`;
  - `~/Downloads/mr_admission_align_after_fs_d6.rds`;
  - `~/Downloads/mr_admission_align_step4.rds`.

## Amendment (Larry, 2026-09-23, mid-task; not in the committed task doc)

1. Option 1 confirmed. **Post-condition 8 is amended**: `R/forestsearch_main.R` is added to the permitted files, for
   the one-line pass-through of `pconsistency.digits` into `fs_mr_inference()`. The original list was written for a
   minimal change and contradicted the correct fix.
2. `fs_mr_inference()` gains `pconsistency.digits = NULL`. `NULL` falls back to the `subgroup.consistency()` default,
   so a direct caller is unaffected. The function is not exported, but it carries an Rd (`@keywords internal`), so the
   roxygen is updated and `man/fs_mr_inference.Rd` regenerated.
3. The gated comparison is the paired before / after at the default `digits = 2`, same seed; no additional
   conditions. An earlier draft of the amendment proposed an extra `digits = 6` fit, and that fit had already run. It
   is reported below as a supplementary check only; it is not a gate.

## Settings (read-only from `fs-glms-interpretable/quarto/gbsg/analysis_gbsg{,_grf,_dina}_mr.qmd`)

**Frame.**
- `survival::gbsg`: N 686, 246 hormonal / 440 none, 299 events.
- Derived columns: `time_months = rfstime / 30.4375`, `grade3 = (grade == "3")`, `id` = row number.
- `confounders = c("age", "meno", "size", "grade3", "nodes", "pgr", "er")`.

**Common to all three identifiers.**
- `is.RCT = TRUE`, `seedit = 8316951`, `est.scale = "hr"`, `use_lasso = use_grf = use_dina = FALSE`.
- Cuts: `max_n_confounders = 1000`, `conf_force = c("er <= 0", "pgr <= 0")`, `cut_type = "default"`,
  `cont.cutoff = 4`, `conf.cont_jcuts = list(er = 10, pgr = 10)`, `collapse_cuts = TRUE`.
- Selection: `sg_focus = "effMaxSG"`, `selection_rule = "neighborhood"`, `effect_neighborhood = 0.20`.
- Size floors: `n.min = 60`, `d0.min = d1.min = 10`, `maxk = 2`.
- Other: `m1.threshold = Inf`, `minp = 0.025`, `fs.splits = 1000`, `stop_threshold = NULL`, `hr.threshold = 1.00`.
- `mr_inference = TRUE` with `mr_inference_args`:
  - `ci_method = "field"`, `draws = 5000`, `include_complement = TRUE`, `confirm_rule = "point"`;
  - `field_uniform = FALSE`, `field_complement = TRUE`, `ij_residual = "two_term"`, `field_decompose = TRUE`;
  - `field_scale_complement = "selected"`, `field_recovery = FALSE`, `return_reselection = TRUE`.
  - Field Monte Carlo: package defaults, 1,000 outer × 500 inner.

**Per identifier.**
- **FS:** `subgroup_method = "consistency"`, `hr.consistency = 1.00`, `pconsistency.threshold = 0.90`,
  `consistency_method = "resample"`. `pconsistency.digits` is not set, so it is 2. Also
  `keep_declaration_field = TRUE`, `keep_field_matrix = FALSE`, `declaration_c0 = c(0.70, 0.75, 0.80, 0.85)`.
- **GRF:** `subgroup_method = "grf"`, `grf_selection = "frontier"`, `grf_select_statistic = "effect"`,
  `grf_depth = 2L`, `dmin.grf = 0.0`.
- **DINA:** `subgroup_method = "dina"`, `dina_select_statistic = "effect"`, `dina_args = list()`.

## The change

- `R/fs_mr_inference.R`, in the both-floors branch of the admission set:
  - before: `z <- qnorm((1 + p_star) / 2)`;
  - after: `z <- qnorm((1 + .fs_pcons_eff(p_star, digits)) / 2)`.
  - `t_g <- pmax(effect_floor, c_cons + z * sigma_D)` is otherwise unchanged. The effect-floor-only and unrestricted
    branches are untouched: `digits` is read only inside the both-floors branch.
- **Where `digits` comes from:** the argument is passed in directly (from the calling function's arguments), not
  read from `args_call_all`.
  - `fs_mr_inference()` had neither the argument nor a fit object, and `args_call_all` is not visible at that level.
  - `forestsearch()` has `pconsistency.digits` as an argument and now passes it (`R/forestsearch_main.R`, one line in
    the `fs_mr_inference()` call).
  - `NULL` falls back to `eval(formals(subgroup.consistency)$pconsistency.digits)`, the same fallback
    `fs_declaration_calibration()` uses.
- **Both consistency paths.** MR's floor is built by one expression, with no branch on `consistency_method`, so the
  change applies identically under `"resample"` and `"split"`. Under `"split"` the closed form was already the
  approximation; only the rounding changes.
- **Other callers.** `.fs_apply_mr()` (the GRF and DINA path) does not pass `digits`. Its admission set has no
  consistency floor, so `digits` is never read there.
- `R/fs_declaration_calibration.R` is unmodified. No selection rule, no screen and no part of the rounding design
  changed.

## Step 3 — GBSG before / after (digits = 2, seed 8316951)

| identifier | quantity | before | after | after − before |
|---|---|---|---|---|
| FS | region | `{er <= 0} & {pgr <= 26}`, N 75, 41 events | same | — |
| FS | unadjusted HR | 2.221839 | 2.221839 | 0 |
| FS | corrected estimate (HR) | 1.280592 | 1.276181 | −0.004411 |
| FS | field lower bound (H) | 0.636532 | 0.633562 | −0.002971 |
| FS | field-s upper bound (Hᶜ) | 0.803507 | 0.803570 | +0.000063 |
| FS | Bonferroni pair (H / Hᶜ) | 0.509882 / 0.839729 | 0.508684 / 0.839594 | −0.001198 / −0.000135 |
| FS | selection bias (log HR) | 0.539365 | 0.543074 | +0.003708 |
| FS | re-selection family size | 1,744 | 1,744 | 0 |
| FS | selection rate | 0.9928 | 0.9932 | +0.0004 |
| GRF | region | `{er <= 0}`, N 82, 45 events | same | — |
| GRF | unadjusted HR | 1.951393 | 1.951393 | 0 |
| GRF | corrected estimate | 1.110994 | 1.110994 | 0 |
| GRF | field lower bound (H) | 0.510471 | 0.510471 | 0 |
| GRF | field-s upper bound (Hᶜ) | 0.846291 | 0.846291 | 0 |
| GRF | Bonferroni pair | 0.403679 / 0.874401 | same | 0 / 0 |
| GRF | selection bias | 0.556622 | 0.556622 | 0 |
| GRF | re-selection family size | 858 | 858 | 0 |
| DINA | region | `{grade3 >= 1} & {pgr <= 20}`, N 113, 61 events | same | — |
| DINA | unadjusted HR | 1.314631 | 1.314631 | 0 |
| DINA | corrected estimate | 1.076093 | 1.076093 | 0 |
| DINA | field lower bound (H) | 0.637786 | 0.637786 | 0 |
| DINA | field-s upper bound (Hᶜ) | 0.842063 | 0.842063 | 0 |
| DINA | Bonferroni pair | 0.572467 / 0.886822 | same | 0 / 0 |
| DINA | selection bias | 0.192624 | 0.192624 | 0 |
| DINA | re-selection family size | 30 | 30 | 0 |

- The unadjusted HR is MR's `naive$est`. It matches `coxph()` on the region to printed precision.
- The direction is as the premise predicts. A lower admission floor admits more candidates, so the correction
  removes more selection bias (+0.0037) and every FS lower bound moves down. The complement's field-s upper bound is
  essentially unchanged.
- **Wall clock per fit:**
  - FS: 47.2 s before, 47.3 s after.
  - GRF: 27.2 s before and after.
  - DINA: 9.6 s before, 9.7 s after.

## Step 4 — draw-level diagnostic (FS, the 5,000 shared main-stream draws)

**Method.** The inputs to `fs_mr_inference()` were captured with a read-only `trace()` in the before fit. MR's
re-selection loop was then re-run on them under both floors. This is a reproduction check as well as a measurement:
the z_exact winners equal the before fit's `reselection$winner`, and the z_eff winners equal the after fit's, exactly.

**Thresholds.**
- z_exact = qnorm(0.95) = 1.644854.
- z_eff = qnorm((1 + 0.895) / 2) = 1.621082, with `pcons_eff` = 0.895.

**Results.**
- **Admitted set differs:** 2,963 / 5,000 draws (**0.5926**).
  - Median admitted count: 37 before, 38 after.
  - Mean extra admitted per draw: 1.36.
- **Re-selected winner differs:** 26 / 5,000 draws (**0.0052**).
  - No-winner draws: 36 before, 34 after. Two of the 26 are draws that had no winner before and gain one after.
  - The observed region wins on 26 draws both before and after.
  - Most frequent switch: `{er <= 0} & {pgr <= 114}` → `{er <= 9} & !{meno}` (3 draws). Every other switch occurs
    once.
- **Band on the observed fit:**
  - Candidates with T in [z_eff, z_exact): **1 of 1,744**.
  - It is `{pgr <= 0} & {er <= 9}` (internal `q2.1 & q4.1`): N 71, HR 1.683928, T 1.63837, closed-form Pcons 0.898655.
  - **This matches the expected candidate** exactly (HR and Pcons). Its rule was confirmed by membership against the
    raw covariates.
  - It also clears the effect floor. On the observed statistics the floor therefore admits 22 candidates at z_exact
    and 23 at z_eff.

**Reading.**
- The admitted set changes on most draws. This is because a 0.024 drop in z moves some candidate across the floor on
  most perturbations.
- The effMaxSG winner changes on only 0.5% of draws, because the new admits rarely win the size-within-neighborhood
  ranking.
- Selection bias is averaged over draws, so a 0.5% change in winners gives a 0.0037 shift on the log scale. Whether the
  Section 5 re-run is a refresh or a finding depends on how often the simulated fits carry candidates near the band.
  On this application the change is small.

## Supplementary (not a gate) — digits = 6 at the same seed

- `forestsearch(…, pconsistency.digits = 6)` ran on the modified package at seed 8316951, 47.1 s. `args_call_all`
  records 6.
- **The MR object is `identical()` to the baseline (before, exact-cutoff) MR object, every element, timings excluded.**
  - Region, family (1,744), corrected estimate, field bounds and Bonferroni pair are all the same to the last digit.
  - At 6 digits, z_exact − z_eff is about 1e-6, and no candidate falls in that band on any draw.
- The pass-through therefore carries a non-default value, and the whole behaviour change is the rounding.

## Gates

- **Gate 1 — PASS.** Seven FS quantities differ (table above), and the band is populated (Step 4).
- **Gate 2 — PASS.** For GRF and DINA, `identical()` holds on the full MR return object, before vs after, with
  wall-clock timings excluded. Every reported quantity is equal.
- **Gate 3 — PASS.** The effective-threshold expression `g - 0.5 * 10^(-digits)` appears in exactly one place in
  `R/`: `.fs_pcons_eff()`, `R/fs_declaration_calibration.R:192–193`. The other grep hits are comments. Line 213
  (`ceiling(round(p_k / step, 6)) * step`) is the settable search's grid start, which predates this task. MR calls
  the helper.
- **Gate 4 — PASS.** `devtools::test()`: 0 failures, 0 errors, 5,826 passes, 3 skips.
  - The skips predate this task.
    - Two are "primary fit did not identify".
    - The third is a multisession dev-load artifact: the workers load the installed 02:32 build, whose
      `forestsearch()` does not have `pconsistency.digits`.
  - **No existing test encoded MR's exact-cutoff admission, so none was updated.**
  - New file: `tests/testthat/test-mr-admission-rounded.R`, 4 tests, 9 expectations, all pass. It was written after
    the suite run started and was run separately with `load_all()` + `test_file()`.
    - MR's per-draw admitted set equals `{bs >= effect floor and zcons >= z_eff}` at digits 2. `.fs_mr_select` is
      stubbed to record each draw's admitted set.
    - The fixture's band is populated: some draw admits a candidate with zcons < z_exact.
    - `NULL` digits gives output identical to the `subgroup.consistency()` default.
    - At digits 6 the floor is z_eff(6), within 1e-5 of z_exact.
    - With no consistency floor (GRF / DINA shape), `digits` does not change the output.
  - Negative control: at digits 6, which is effectively the old rule, 13 of 239 draws violate the z_eff assertion. The
    test therefore discriminates between the old and new rules.
- **Gate 5 — PASS.** `rcmdcheck::rcmdcheck(args = "--as-cran")` returns 0 errors, 0 warnings and 1 NOTE: "Version
  contains large components (0.3.5.9000)", the deliberate dev-version marker. The PDF manual was built. `RSTUDIO_PANDOC`
  was set to the RStudio quarto tools directory.
  - The checked tarball includes the `R/`, `man/` and test changes. It was built before the NEWS.md entry was written.
  - The edited NEWS.md was then parsed with `tools:::.build_news_db_from_package_NEWS_md()`: 40 entries, 0 bad. The
    new entry is ASCII.
  - The version was not bumped.

## Files modified

`R/fs_mr_inference.R`, `R/forestsearch_main.R` (as amended), `tests/testthat/test-mr-admission-rounded.R` (new),
`man/fs_mr_inference.Rd`, `NEWS.md`, `dev/tasks/`, `dev/reports/`. Nothing was written to `fs-glms-interpretable`.

## Open items (not fixed here)

- The manuscript text still states the exact cutoff: Section 4.6's Eq. (8), Section 2.3.3's Step 3 and Section 4.4.
  Drafting that belongs to the manuscript chat.
- The consumers of the removed `pstar_implied` are still out of scope, as the task says.
- The Section 5 simulation re-run is a separate decision.
