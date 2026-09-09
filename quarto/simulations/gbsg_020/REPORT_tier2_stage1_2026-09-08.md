# REPORT — Tier 2 (Mac Studio), Stage 0 (installed-package quotes and machine record) and Stage 1 (5-replicate smoke, Gate 1 projection)

**Task:** `dev/tasks/TASK_tier2_mac_2026-09-08.md` (2a84f6c4). **Executor:** Claude Code on Mac-Studio-3, unattended (Larry offline until the morning). Runs in parallel with the Linux `cert20` campaign; outputs disjoint (campaign tag `tier2`, prevalence 12.4%, focus `maxeffCons`); **no `R/`, template or shared-document change of any kind.**
**Date:** 2026-09-08. Winner-only and winner-floor excluded from every table and line. **No recommendation change; report and wait.**

---

## GATE 1: COMPUTE GO — projection **≈ 2 h 45 m** at 12 workers (ceiling 9 h, hard timeout 11 h); **all four cells proceed, none deferred at the start.**

---

## Machine and HEAD

| item | value |
|---|---|
| host | `Mac-Studio-3.local` (Apple M4 Max, arm64, macOS 26.6.2) |
| physical / logical cores | 14 / 14 |
| RAM | 36 GiB (38,654,705,664 B) |
| R | 4.5.2 |
| repo HEAD at install | `2a84f6c4b013c557d0d4838b549f877452d788f4` — *Task spec as received: TASK_tier2_mac_2026-09-08.md …* (2026-09-08 20:46:06 −0700) |
| last commit touching `R/` | `c0f48a7c82a5326436f7d23d79cfe2b6f3df24b9` — *Merge remote-tracking branch 'origin/feature/glm-extension-mac' into feature/glm-extension* (2026-09-08 16:24:03 −0700) |
| `git pull` at start | `Already up to date.` (the Linux tree had not pushed tonight's commits at that moment) |
| install | `devtools::install(dependencies = FALSE, upgrade = FALSE)` → `* DONE (forestsearch)`, source install with `--install-tests`, library `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library` |
| installed version | **0.3.5** (`meta$forestsearch_version` on every bundle below) |
| BLAS environment for every render | `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` (single-threaded BLAS under a worker farm; Accelerate thread-safety note, `mac-render-memory-cap`) |

**Worker count chosen: 12** (template caps at `physical cores − 1` = 13; 12 leaves two physical cores and the parent process headroom). Measured peak across the probe renders: **19.3 GB total RSS over 13 R processes** (12 workers + parent), max single process 1.00 GB at n = 500 and 1.34 GB at n = 1500 — inside the standing rule (workers × per-worker RSS < 24 GB, no other R process running). 13 workers would project ≈ 21 GB; 12 is the headroom choice, and the dense stages here are memory-bandwidth-bound, so the extra worker buys little.

## Stage 0 — quotes from the **installed** package (0.3.5) and the committed template

**1. The `field_scale_complement` argument of `fs_mr_inference()`** — `formals()` of the installed function (the whole field block for context):

```
  field_R_out = 1000L
  field_R_in = 500L
  field_uniform = FALSE
  field_M_cap = NULL
  field_complement = FALSE
  field_decompose = FALSE
  field_scale_complement = c("none", "selected")
```

and the installed `man/fs_mr_inference.Rd` entry, verbatim:

> `\item{field_scale_complement}{\code{"none"} (default) or \code{"selected"}; consulted only when the complement field block runs.  Under \code{"selected"} the complement's field readings are studentized to the selected complement's own scale before differencing -- each outer reading \verb{zeta^c_{r, G_r}} is multiplied by \verb{s_sel / s_{G_r}} and each inner reading by \verb{s_sel / s_{G(v_r + zeta'_j)}}, with \code{s_g = sqrt(sum(Bc[, g]^2))} the candidate's complement influence-norm scale (the percentile-t root; PROPOSAL_complement_field_scale_2026-09-08_v2 s5, variant R1) -- and \code{field$complement} gains the \verb{_s} companio…`

**2. The `lam_cs` line of `.fs_mr_field_complement()`** — deparsed from the installed namespace (line numbers of the deparse, not of the source file), confirming **field-s is present at this HEAD**:

```r
  47:     lam_c <- rep(NA_real_, R_out)
  48:     lam_cs <- rep(NA_real_, R_out)
  ...
  64:         lam_c[r] <- Zo_c[G, r] - mean(Zi_c[cbind(wi[ok_in], ok_in)])
  65:         if (scale_on)
  66:             lam_cs[r] <- (s[sel]/s[G]) * Zo_c[G, r] - mean((s[sel]/s[wi[ok_in]]) *
  67:                 Zi_c[cbind(wi[ok_in], ok_in)])
  ...
  83:     if (scale_on) {
  84:         lfs <- lam_cs[ok_c]
  85:         qss <- stats::quantile(lfs, c(0.05, 0.25, 0.5, 0.75,
  86:             0.95, 0.025, 0.975), names = FALSE, type = 7)
  87:         est2s_w <- bdc - mean(lfs)
  88:         sd_cs <- stats::sd(lfs)
  89:     }
```

**3. The template's `FS_S7_FIELD_SCALEC` knob** — `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:596–603`, verbatim (unchanged; read-only in this task):

```r
# Studentized complement field knob (TASK_field_studentize_e1_2026-09-08):
# "none" (default) or "selected" -> fs_mr_inference(field_scale_complement).
# Under "selected" the gate attaches the field-s companions (est2_s,
# upper_1s_s, ... , se_field_s; joint_s) beside the unscaled complement field
# and the recorder captures the fld_Hc_*_s / fld_joint_s_* columns below.
# Add-beside: every pre-existing column is identical under either value.
mr_field_scalec <- .env_chr("FS_S7_FIELD_SCALEC", "none")
stopifnot(mr_field_scalec %in% c("none", "selected"))
```

forwarded at `:614` inside `mr_inference_args` as `field_scale_complement = mr_field_scalec`, beside `return_reselection = TRUE`.

## Stage 1 — the 5-replicate smoke (HR 1.75, n = 500, campaign tag `tier2smoke`)

Render: `FS_S7_FOCUS=maxeffCons FS_S7_HR=1.75 FS_S7_N=500 FS_S7_NSIMS=5 FS_S7_START=1 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=tier2smoke FS_S7_WORKERS=5` — **`FS_S7_Z1Q` deliberately unset**, so the M1 default 0.25 stands. Bundle `results/fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n500_tier2smoke_res_1_5.rds`, render `…_tier2smoke_batch_1_5.html`.

**Meta as set:** `field_complement TRUE`, `field_decompose TRUE`, `field_scale_complement "selected"`, `ij_residual "two_term"`, `fb_mode "none"`, `sg_focus "maxeffCons"`, `effect_neighborhood 0.1`, `er_jcuts 10` (J = 10 default), `harm_z1_quantile 0.25`, `seed_base 8316951`, `forestsearch_version "0.3.5"`, `hostname "Mac-Studio-3.local"`.

| check | result |
|---|---|
| rows / sim_id / duplicates | 5, sim_id 1–5, none |
| status | `DETECTED` on 5 / 5; `mr_ok` 5 / 5; CONFIG-ERROR rows 0 |
| **realized prevalence** | `n_true` mean 67.0 of n = 500 → **0.1340** (range 57–76); `meta$harm_prevalence_super` = **0.12418** — the M1 default 12.4 % configuration |
| finiteness on detected rows | **0 non-finite** in all 9 blocks: harm (6 cols), complement (7), `_s` (7), scale (4), joint (9), `joint_s` (9), IJ (8), naive (4), identification/β/p̂ (8) |
| interval invariants | harm `lo2s ≤ lo1s ≤ est2 ≤ hi2s` **TRUE** (5/5); complement `lo2s ≤ lo1s ≤ est2 ≤ up1s ≤ hi2s` **TRUE** (5/5); `_s` block **TRUE** (5/5); IJ `lo ≤ est ≤ hi` on both blocks **TRUE**; `se_field > 0`, `se_field_s > 0` **TRUE** |
| **γ in range** | `fld_joint_gamma` ∈ [0.02500, 0.02600] and `fld_joint_s_gamma` ∈ [0.02500, 0.02600] — both inside [0.025, 0.05] |
| achieved joint probability | min 0.9498 on all four probability columns (`joint`/`joint_s`, calibrated and Bonferroni) — **≥ 0.95 − 2/n_joint = 0.948** ✔; at the α/2 fallback on 0.800 (`joint`) and 0.600 (`joint_s`) of rows |
| joint wider than marginal | `bonf_loH ≤ fld_H_lo1s` and `bonf_upHc ≥ fld_Hc_up1s` **TRUE** on all rows |
| **ρᶜ recorded** | `fld_Hc_scale_ratio`: mean **1.0177**, SD 0.0152, range 1.0036–1.0422, share > 1 = 1.000 |
| R1 vs global rescale | `se_field_s / (ρᶜ · se_field)` mean 1.0001, range 0.9992–1.0010 |
| dominated-regime preview | `se_field_s / naive SE` mean **0.9994**; `se_field / naive SE` mean **0.9819** — ρᶜ ≈ 1 and field-s ≈ naive scale, as pre-registered for this regime |

**Pre-flight pairing against the comparator** (`fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n500_s7c_combined_1_2000.rds`, rows 1–5): **`n_true` `identical()` TRUE** (72 69 61 76 57 on both); `detected`, `n_sel` (66 69 69 99 66) and `label` identical; `betaHhat_H` max |difference| **0** on these rows. **`truth` is `identical()` FALSE but `all.equal()` TRUE** — `hr_causal` and `marg_Hc` agree to the last bit, `marg_H`, `cde_H`, `cde_Hc` differ by 1.3e−16, 6.5e−16, 3.8e−16 relative: the cross-machine BLAS-precision difference the task anticipated. **Per the task's cross-machine note, `truth` is therefore gated at the stated ~1e−8 tolerance, not by `identical()`**; the exact `identical()` gate is applied to `n_true` (which passes) and every within-bundle identity stays machine-local.

## Gate 1 — timing and projection

Per-round wall measured at the production worker count (12) by differencing two renders of the same cell, which nets out the ≈ 14 s quarto/DGM-build overhead:

| cell | n_sims = 12 (1 round) | n_sims = 60 (5 rounds) | per round (12 replicates) | implied overhead |
|---|---|---|---|---|
| HR 1.75, n = 500 | 27 s | 77 s | **12.5 s** | 14.5 s |
| HR 1.75, n = 1500 | 29 s | 88 s | **14.75 s** | 14.25 s |

Per-replicate MR cost from the 5-replicate smokes (5 workers, little contention): `fit_mr_secs` mean 9.45 s (HR 1.75 n500), 10.96 s (n1000), 11.15 s (n1500), 5.53 s (HR 1.00 n500 — 3/5 detected, undetected replicates are cheap); `fld_H_secs` 4.8–6.4 s, `fld_Hc_secs` 0.38–0.51 s.

Projection, 2,000 replicates per cell = 167 rounds of 12 at 12 workers, plus two batch overheads and one combine render per cell:

| # | cell | per round | projected wall |
|---|---|---|---|
| 1 | HR 1.75, n = 500 | 12.5 s | ≈ 36 m |
| 2 | HR 1.75, n = 1000 | 13.6 s (interpolated) | ≈ 39 m |
| 3 | HR 1.75, n = 1500 | 14.75 s | ≈ 42 m |
| 4 | HR 1.00, n = 500 | 12.5 s (conservative; the smoke is faster) | ≈ 36 m |
| | **campaign** | | **≈ 2 h 33 m + ≈ 12 m of combine renders ≈ 2 h 45 m** |

**≤ 9 h ceiling → Gate 1 GO on all four cells; none deferred, none dropped at the start.** Hard timeout 11 h; the defer order if a cell over-runs in flight is HR 1.75 n1500 first, then HR 1.00 n500.

## Artefacts committed with this record

`results/…_h175_knoise0_n500_tier2smoke_res_1_5.rds`, `…_h175_knoise0_n1000_tier2smoke_res_1_5.rds`, `…_h175_knoise0_n1500_tier2smoke_res_1_5.rds`, `…_h100_knoise0_n500_tier2smoke_res_1_5.rds` (the four smokes) and `…_tier2probe_res_1_12.rds` / `…_tier2probe_res_1_60.rds` at n = 500 and n = 1500 (the timing probes), each beside its render. The `tier2smoke` and `tier2probe` tags are outside the `tier2` combine glob by construction, so no probe row can ever enter a production pool.
