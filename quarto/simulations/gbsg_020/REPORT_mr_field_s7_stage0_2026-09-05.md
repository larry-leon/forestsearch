# REPORT — MR (field) on Section 7 cells: Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_mr_field_section7_2026-09-05.md` (committed `bb516569`)
**Date:** 2026-09-05. Executor: Claude Code. Stage 0 only; stopped at Gate 0 per instruction.
**Decisions:** F1–F5 at defaults; F7 deferred to Gate 1.

---

## GATE 0: PASS — no stop condition triggered; four findings for Larry before Stage 1

- The wrapper needs **no arithmetic change** (0a).
- The seed schemes of the MR-only and FB documents are **identical**, line-for-line (0b).
- The naive/oracle cross-bundle identity **holds at machine precision** (max |diff| 3.6e-15, identical NA patterns, identical detection/membership); strict `identical()` is FALSE due to 0.2.0 vs 0.2.1 build vintages — see Finding 3 (0b).

Findings requiring Larry's acknowledgment (none are stop conditions; all affect Stage 1/2 planning):

1. **The committed FB range is sim_id 1–500, not 1–1000.** The task's Stage 2 text ("FB joined for the key cell on sim_id 1–1,000 (the committed range)") does not match the tree. The production FB bundle is `results/fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` (nb_boots = 300, 5 batches, sim_id 1–500). The only 1–1000 FB-stem file is `..._quickrun_res_1_1000.rds` (quickrun namespace, nb_boots ≠ production; the `fb_mr` `batch_1_1000.qmd` is committed in a `nb_boots <- 0L`, `quickrun <- TRUE` state — an MR-only quick run). Stage 2 should join FB on 1–500 only and record 501–2,000 as FB-absent.
2. **The h175 bundle exists under stem `h18`, not `h175`.** `rds_stem` uses `round(10 * target_hr_harm)` (h175 qmd, lines 164–166), so target_hr_harm = 1.75 → `h18`. The bundle is tracked: `results/fs_maxeffCons_mr_m1_h18_knoise0_n500_res_1_1000.rds`, and its save line appears verbatim in the committed render (`sim_fs_maxeffCons_mr_m1_h175_knoise0_n500_batch_1_1000.html`, matching mtime 2026-08-10 12:14). Note the collision hazard: target_hr_harm = 1.8 would also map to `h18` — the template's `campaign_tag` (task 1a) must guard this.
3. **Cross-bundle identity is exact-to-the-ulp, not bitwise.** The h10 MR-only bundle was built under forestsearch 0.2.0 (2026-08-10), the FB combined bundle under 0.2.1 (2026-08-18). All naive/oracle/MR-IJ numeric diffs are ≤ 3.6e-15 absolute (≤ 1.5e-15 relative), NA patterns identical, `detected`/`mr_ok`/`n_sel`/`n_harm`/`sens`/`spec`/`ppv`/`npv` identical — the FB join license holds. 19/500 rows spell the selected rule differently (`{meno == 0}` vs `!{meno}` — same member set, 0.2.0-era clause-spelling drift); Stage 1's identity checks must therefore compare at tolerance ≤ 1e-12 and compare memberships, not definition strings.
4. **The fit object has no `debias_gate` element.** The gate's return lands as `mr_inference` (plus `mr_harm_confirmed`); the task's parenthetical names an element that does not exist. Naming only — everything the task means by `debias_gate` is `fit$mr_inference`.

---

## 0a — Wrapper → gate: `ci_method = "field"` pass-through and the `field` element

**Forwarding.** `forestsearch()` (R/forestsearch_main.R) forwards `mr_inference_args` to `fs_mr_inference()` untouched, with `"ij"` as the only defaulting:

- `R/forestsearch_main.R:3368` — call opens: `fs_mr_inference(df = df.fs, candidates = fam, spec = gspec, ...`
- `R/forestsearch_main.R:3385` — `ci_method = .g_mr(mr_inference_args$ci_method, "ij"),` where `.g_mr <- function(a, b) if (is.null(a)) b else a` (line 3319). A supplied `ci_method = "field"` passes through **untouched**.
- `R/forestsearch_main.R:3384` — `include_complement = .g_mr(mr_inference_args$include_complement, TRUE),`
- `R/forestsearch_main.R:3386` — `seed = .g_mr(mr_inference_args$seed, seedit))`

**Return placement.** `R/forestsearch_main.R:3429–3430`:

```r
mr_inference = mr_out,
mr_harm_confirmed = if (!is.null(mr_out)) mr_out$harm_flag
```

**`field` element survival.** `fs_mr_inference()` (R/fs_mr_inference.R, signature lines 398–403: `ci_method = c("ij", "wald", "field")`, `field_R_out = 1000L`, `field_R_in = 500L`) attaches the field block at top level of its return, `R/fs_mr_inference.R:735`:

```r
if (!is.null(field)) out$field <- field
```

`mr_out` is stored unmodified, so the block survives as **`fit$mr_inference$field`** with elements (lines 687–705): `lambda_mean`, `lambda_sd`, `q05 q25 q50 q75 q95 q025 q975`, `n_out_used`, `n_in_used_mean`, `est2`, `lower_1s`, `lower_2s`, `upper_2s`, `se_field`, `lower_se`, `upper_se`, `R_out`, `R_in`, `seed_offset` (900000L), `timing_seconds`; degenerate case (`< 2` usable outer draws) yields a `note` + counts (lines 702–705).

**Complement (decision F3).** The field loop computes λ for the **winner only** — `lam[r] <- Zo[G, r] - mean(Zi[cbind(win, ii)])` (R/fs_mr_inference.R:672, 678) — and the `complement` list (lines 608–621) carries naive/debiased with IJ/Wald intervals only. **No `field` block exists for Ĥᶜ under `include_complement = TRUE`.** A complement field would be new arithmetic inside the field block, not a pass-through. → **F3 default applies: Ĥᶜ stays on IJ, with a note.**

**One add-only pass-through gap (task Protocol ¶ "add-only" clause).** The wrapper's call (lines 3368–3386) does **not** forward `field_R_out` / `field_R_in`; via `forestsearch()` the field always runs at the 1000/500 defaults. Not a stop (no arithmetic touched; defaults are the task's implied sizes). If Stage 1 wants these tunable from the template, the fix is an add-only pass-through in the Stage 1 record.

## 0b — Committed bundles, seed scheme, cross-bundle identity

**Paths (all under `quarto/simulations/gbsg_020/`, all git-tracked):**

| Cell | Bundle | rows / sim_id | meta |
|---|---|---|---|
| Key h10, MR-only | `results/fs_maxeffCons_mr_m1_h10_knoise0_n500_res_1_1000.rds` | 1000 / 1–1000 | n_sims 1000, nb_boots NULL, mr_draws 5000, seed_base 8316951, fs 0.2.0, built 2026-08-10 10:58, 115 workers, parallel_mode "sims" |
| Key h10, FB | `results/fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` | 500 / 1–500 | n_sims 500 (5 batches), **nb_boots 300**, mr_draws 5000, seed_base 8316951, fs 0.2.1 |
| h175 separation | `results/fs_maxeffCons_mr_m1_h18_knoise0_n500_res_1_1000.rds` (stem `h18`, Finding 2) | 1000 / 1–1000 | n_sims 1000, nb_boots NULL, mr_draws 5000, seed_base 8316951, fs 0.2.0, built 2026-08-10 12:14, 115 workers, parallel_mode "sims" |

The h175 bundle's schema is column-identical to the h10 MR-only bundle. Detections: h10 675/1000; h175 950/1000; FB 329/500 (FB `fb_H_est` finite on exactly the 329 detected rows). h175 truth block: marg_H = 1.769, cde_H = 2.036 (HR scale), consistent with the calibrated target_hr_harm = 1.75.

**Seed scheme — identical across the three documents** (`seed_base <- 8316951L` in all; qmd line numbers MR-only/h175 vs fb_mr):

| Draw | MR-only h10 & h175 (same lines) | fb_mr h10 |
|---|---|---|
| DGM calibration | `n_super = n_super, seed = seed_base)` :319 | :476 |
| Trial data | `simulate_from_dgm(..., seed = seed_base + sim_id)` :466 | :623 |
| Noise covariates | `set.seed(seed_base + sim_id + 1000000L)` :482 | :639 |
| Search + MR | `seedit = seed_base + sim_id` :502 | :659 |
| FB bootstrap | `nb_boots = nb_boots, seed = seed_base + sim_id` :595 | :760 |

**Cross-bundle naive/oracle identity on the 500 common sim_ids** (h10 MR-only vs FB combined): strict `identical()` FALSE on all 16 nv/or columns; measured differences all ≤ 3.6e-15 absolute / 1.5e-15 relative, zero NA-pattern differences. `detected`, `mr_ok`, `n_sel`, `n_harm`, `sens`, `spec`, `ppv`, `npv`, `status`, `mr_harm_flag`, `betaHhat_status`, `nH_eval`, `nHc_eval` identical. MR-IJ columns (`mr_H_est/lo/hi/se_ij`) likewise ≤ 3.6e-15. 19 rows differ in `sg_def`/`label` spelling only (e.g. sim 71: `{nodes <= 5} & {meno == 0}` vs `{nodes <= 5} & !{meno}`); membership metrics prove the same sets. **Verdict: the FB join is licensed** (Finding 3 sets the Stage 1 comparison convention: tolerance ≤ 1e-12, memberships not strings).

Check scripts (scratchpad, read-only; not committed): `stage0_bundle_checks.R`, `stage0_diff_characterize.R`.

## 0c — Recorder

Per-replicate columns, quoted from `sim_fs_maxeffCons_mr_m1_h10_knoise0_n500_batch_1_1000.qmd` `.na_record()` (lines 432–459) — identical in the h175 and fb_mr documents:

- Bookkeeping/detection/selection: `sim_id, detected, mr_ok, n_sel, status` (CONFIG-ERROR / NO-DETECTION / DETECTED trichotomy), `n_harm, n_true, label, sg_def, covs, err_msg`
- Truth attachment: `betaHhat_H, betaHhat_Hc` (+ `betaHhat_status, nH_eval, nHc_eval` appended by `fs_attach_betaHhat(results, eval_df, focus = "harm", ...)` at line 730)
- Cost: `fb_secs, fit_mr_secs, fb_err`
- Estimators × blocks (est/lo/hi/se): `or_H_*, or_Hc_*` (oracle, refit on the true harm subgroup, lines 634–638), `nv_H_*, nv_Hc_*` (naive, from `g$naive` / `g$complement$naive`), `fb_H_*, fb_Hc_*` (H2 = León Eq. 7, lines 625–631), `mr_H_est/lo/hi/se_ij, mr_Hc_est/lo/hi/se_ij` (from `g$debiased` / `g$complement$debiased`, lines 561–574)
- MR extras: `mr_harm_flag, ij_source`; classification: `sens, spec, ppv, npv`

MR block extraction (lines 545–575): `g <- fs.est$mr_inference`; detection is method-agnostic (`found <- !is.null(fs.est$sg.harm)`), `mr_ok` tracked separately so an MR failure never masks a detection. The MR call passes `mr_inference_args = list(ci_method = "ij", draws = mr_draws, include_complement = TRUE)` (line 513–514) — the template's knob delta is exactly `ci_method = "field"` plus the field-column recorder extension.

**Four coverage targets**, from `build_block()` (lines 1082–1116): per-replicate `Cov_oracle` (`or_<blk>_est`), `Cov_betaHhat` (`betaHhat_<blk>`, β(Ĥ) truth), and scalar-recycled `Cov_cde` (`truth$cde_*`, θ‡) and `Cov_marg` (`truth$marg_*`, θ†), via `.cover(target, lo, hi)` (lines 1057–1062, finite-only two-sided containment).

## 0d — Cost

`fit_mr_secs` = one `forestsearch()` fit **including** MR-IJ at draws = 5000 (comment, qmd line 539), measured per replicate while 115 workers ran concurrently under `parallel_mode = "sims"` (i.e., loaded-host per-replicate wall, the regime Stage 1c projects for):

| Bundle | n | mean | q0 / q25 / q50 / q75 / q90 / q99 / max (s) |
|---|---|---|---|
| h10 MR-only (1–1000) | 1000 | 28.3 | 5.98 / 17.4 / 32.5 / 36.6 / 38.9 / 43.8 / 45.5 |
| h175 (h18 stem, 1–1000) | 1000 | 35.0 | 7.69 / 33.7 / 36.2 / 38.8 / 40.9 / 45.3 / 49.3 |
| FB combined (1–500), fit_mr part | 500 | 7.4 | 2.31 / 4.37 / 6.47 / 7.42 / 8.70 / 36.1 / 40.2 |

(The FB document's smaller `fit_mr_secs` reflects a different concurrency regime — its wall went to the FB bootstrap, `fb_secs` mean 1973 s on the 329 detected rows; FB-side cost is irrelevant to Stage 2 since FB is joined, never recomputed.) The h175 render header: "Batch sim_id 1-1000: 950/1000 detections in 342.1 s" total wall at 115 workers. The field block will add its own `timing_seconds` on top of these — Stage 1c measures that at ~60 workers per task 1c.

---

## Deviations

- Task 0b names the h175 bundle path as `results/<stem>_res_1_1000.rds` with an `h175` stem; the actual tracked stem is `h18` (Finding 2). No file named `*h175*.rds` exists anywhere in the tree, history, or the parked dirs — resolved via the committed render's save line, not by re-running anything.
- Task 0a's "`debias_gate` element" does not exist under that name (Finding 4); verified against `R/` (no occurrence outside this task's wording) and quoted the actual element (`mr_inference`).
- Committed simulation documents and bundles were read-only throughout; no `R/` edits; nothing re-run.

**Stopped at Gate 0** per instruction. Stage 1 (template + identities + projection) awaits Larry's acknowledgment of Findings 1–3 and the F7 compute decision at Gate 1.
