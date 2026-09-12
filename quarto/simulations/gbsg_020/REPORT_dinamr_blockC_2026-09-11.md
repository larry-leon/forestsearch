# REPORT — `dinamr` Block C, the deferred Block B cell, and the GRF cost probes

- **Date:** 2026-09-11. **Machine:** Mac Studio (darwin 25.6.0), R 4.5.2, `forestsearch` 0.3.5.
- **Task:** `dev/tasks/TASK_dinamr_blockC_grfprobe_2026-09-11.md` (committed as the first action of the session).
- **Predecessor:** `dev/tasks/TASK_dinamr_campaign_2026-09-10.md`, `REPORT_dinamr_2026-09-10.md`.
- **Commit only; not pushed. No `R/` change was made or needed.**

## Framing, restated

- **This does not certify DINA.** The fixed-family condition does not hold, so every coverage
  number in this record and in `summary_dinamr.qmd` is coverage of the
  **conditional-on-proposed-family estimand**, and every table says so.
- **Part B is cost measurement for GRF.** It is not a GRF campaign and not an adopted GRF
  configuration. No coverage claim, no acceptance criterion, no recommendation comes out of it.

---

# GATE 0 — Stage 0, verified from source

## Block C's cell definition in the predecessor

`dev/tasks/TASK_dinamr_campaign_2026-09-10.md:36`, the Part C block table:

```
| **C** | both prevalences | HR 1.00 at n = 500, 1000, 1500 — six cells |
```

with the defer order at line 38 ("Block C first (all six, 31% before 12.4%)"), and the grid's
provenance at line 10: "**Grid mirrors FS** (Larry, 2026-09-10): both prevalences, ε = 0.20, the
same HR and n structure as `cert20` / `tier2` / `p12ext`."

## What makes a planted-harm quantity undefined at HR 1.00 — and what does not

`FS_S7_HR` has exactly one meaning in this harness. Quoting the template
(`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`):

- **line 340** — `target_hr_harm  <- .env_num("FS_S7_HR", 1.0)   # calibrate k_inter to this Cox HR in the harm subgroup`
- **lines 344–346** — the harm rule, in the comment block that documents the prevalence knob:
  ```
  # z1_quantile is the engine's own formal: the harm rule is
  #   {er <= quantile(er, z1_quantile)} & {meno == 0}
  # (the same form at every value).
  ```
- **line 353** — `harm_z1_quantile <- .env_num("FS_S7_Z1Q", 0.25)`
- **line 488** — `dgm_model       <- "alt"`, a literal; there is no `.env_*` read for it anywhere
  in the document.
- **line 661** — `harm_col       <- "flag_harm"`
- **lines 703–708** — `k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm, ...)` then
  `dgm <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter, z1_quantile = harm_z1_quantile, ...)`

So the harm region `H` is planted by the **rule**, which is a function of `FS_S7_Z1Q` alone; `HR`
enters only as the calibration **target for the effect inside that region**.
`calibrate_k_inter()` (`R/sim_aft_gbsg.R:1004`) roots `dgm$hr_H_true` on `target_hr_harm`.

The harness does own a genuine global-null path — `model = "null"`, documented at
`R/sim_aft_gbsg.R:93` as "uniform treatment effect", and accepted by
`setup_gbsg_dgm()` at `R/setup_gbsg_dgm.R:89` — but `dgm_model` is pinned to `"alt"` at template
line 488 and is not overridable, so **that path is unreachable from this template.**

### The calibration, run

| `z1q` | target HR | `k_inter` | prevalence | overall causal HR | HR(H) | HR(Hᶜ) |
|---|---|---|---|---|---|---|
| 0.25 (12.4%) | 1.00 | **+0.567917** | 0.1242 | 0.6847 | 1.0005 | 0.6569 |
| 0.25 (12.4%) | 1.50 | +1.111548 | 0.1242 | 0.7041 | 1.5086 | 0.6569 |
| 0.60 (31%)   | 1.00 | **−8.579297** | 0.3065 | 0.7922 | 0.9999 | 0.7206 |
| 0.60 (31%)   | 1.50 | −19.497297 | 0.3065 | 0.8694 | 1.4990 | 0.7206 |

Confirmed on a live 2-replicate render at (12.4%, HR 1.00, n 500): the bundle's `truth` reads
`hr_causal 0.6847, marg_H 1.00045, marg_Hc 0.65689, cde_H 1.0000, cde_Hc 0.58478`.

### Verdict, and the STOP condition

**HR 1.00 is a planted region with HR 1.00 inside it, not a global null.** `k_inter` is not zero
at either prevalence; the planted contrast is HR(H) = 1.00 against HR(Hᶜ) = 0.657 / 0.721. The
treatment effect is heterogeneous at these cells. What is absent is **harm** (HR > 1), not the
region.

**The STOP condition is not met.** The kickoff stops the task only if the predecessor's Block C
definition is *ambiguous* between the two readings. It is not: the predecessor pins the cells
through `FS_S7_HR` (its own `campaign.sh` sets `FS_S7_HR=$H` and nothing else), that knob has a
single documented meaning at template line 340, and the harness offers no way to express a global
null on this path. The reading is forced by the source, not chosen.

## Which quantities are structurally undefined at these cells

**None of the classification set.** Because `flag_harm` is planted by the rule and not by the
effect size, it resolves at HR 1.00 exactly as at HR 1.50 (super-population prevalence 0.1242 and
0.3065, printed by the template's own "Harm rule:" line). Therefore:

| quantity | status at HR 1.00 | how it is reported |
|---|---|---|
| oracle θ̂(H) (`or_*`) | **defined** — a Cox fit on the planted region; its target is log 1.00 = 0 | normally, as at every other cell |
| `sens`, `spec`, `ppv`, `npv` | **defined** — cross-tabulated against the same planted H | normally; finite on every detected replicate of the smoke |
| `betaHhat_H`, `betaHhat_Hc` | **defined** — realized targets of the returned region | normally |
| the nine `fld_recov_*` columns | **defined and populated** (9 of 9 on the smoke) | normally |
| `p_hat_*`, `fld_Hc_scale_ratio` (ρᶜ) | **defined and populated** | normally |
| `n_cons_qual` | present, all-NA | **structural on DINA** (no consistency screen) — @sec-structna, never a failure |
| `band_n` | present, all-NA | **structural on DINA** — never a failure |
| `p_star` | not a recorder column | **structural on DINA** (admission-set term, NULL on this engine) — never a failure |

The three structural columns are structural on **every** DINA cell, at every HR; HR 1.00 adds
nothing to the list. **No quantity is dropped, NA-ed or suppressed on account of HR 1.00.**

## The "false-selection rate" label — recorded, and flagged

The kickoff asks that detection at a null be labelled a false-selection rate. That label is
**recorded as instructed and flagged**, because the Gate 0 finding does not support it as stated:
the planted region carries HR 1.00 and **DINA's effect floor is log(0.90)**, so the planted region
itself clears DINA's own admission floor at these cells. A replicate that returns H is making an
**admissible** selection under the criterion in force, not a false one. The document therefore
names the column `selection_rate` and carries the requested reading beside it.
**Which of the two names is right is Larry's call, not this report's** — it is recorded here and
in `summary_dinamr.qmd` @sec-null, and nothing was decided on it.

---

# GATE 1 — the projection and the go/no-go

## The ceiling and its source

**Ceiling 9 h wall for Part A; hard timeout 12 h.** Source: the kickoff's "Gate 1 — compute
go/no-go" bullet, "**Ceiling 9 h wall for Part A; hard timeout 12 h.**" The 9 h figure originates
in **Amendment 1** of the predecessor (`TASK_dinamr_campaign_2026-09-10.md:68–72`), which replaced
the body's 10 h with 9 h on the grounds that FS projections have run +16% (`cert20`) to −29%
(`p12ext`) against realized and DINA's ~100× family-size skew is less predictable.

## Which original Gate 1 probes covered a Block C cell

**Four of the ten did**, so no new 36-replicate probe was needed. `probe.sh` lines 30–33:

```
run ""    500  1.00 p124_h100_n500
run ""   1500  1.00 p124_h100_n1500
run 0.60  500  1.00 p31_h100_n500
run 0.60 1500  1.00 p31_h100_n1500
```

These include both corners the kickoff would otherwise have had me run — (12.4%, HR 1.00, n 500)
and (31%, HR 1.00, n 1500). All four bundles are on disk.

## The probe corners at HR 1.00, as measured

| prevalence | n | reps | selection rate | K q10 | K med | K q90 | K max | s q10 | s med | s q90 | s max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 12.4% | 500  | 36 | 0.7778 | 13.5 | 76.0  | 947.9  | 1195 | 0.0545 | 4.421 | 11.886 | 17.533 |
| 12.4% | 1500 | 36 | 0.3611 | 12.8 | 35.0  | 208.4  | 480  | 0.0790 | 0.113 | 6.417  | 9.754  |
| 31%   | 500  | 36 | 0.9444 | 31.9 | 359.5 | 1634.1 | 3053 | 3.2900 | 9.479 | 24.064 | 34.125 |
| 31%   | 1500 | 36 | 0.9444 | 50.7 | 258.5 | 1203.8 | 1699 | 4.2625 | 8.450 | 18.955 | 29.574 |

## Projecting from the distribution, not a mean

`projectC.R` (committed) bootstraps the 2,000-replicate total by resampling each corner's **36
measured per-replicate seconds** 4,000 times, so the projection carries an interval rather than
multiplying a mean. n = 1000 is not probed at HR 1.00 and is interpolated by **pooling the n 500
and n 1500 draws**, which interpolates the distribution rather than only its centre.

Calibration is on the campaign's **realized** Block A / Block B walls, reconstructed from the
committed bundles' mtimes (cell wall = batch 1 + batch 1001 + combine):

| block | HR | n | realized (h) | Gate 1 (h) | ratio |
|---|---|---|---|---|---|
| A | 1.50 | 500  | 0.518 | 0.412 | 1.256 |
| A | 1.50 | 1000 | 0.389 | 0.374 | 1.040 |
| A | 1.50 | 1500 | 0.330 | 0.275 | 1.196 |
| A | 1.75 | 500  | 0.572 | 0.412 | 1.388 |
| A | 1.75 | 1000 | 0.470 | 0.374 | 1.257 |
| A | 1.75 | 1500 | 0.429 | 0.275 | 1.556 |
| B | 1.50 | 500  | 1.395 | 1.168 | 1.195 |
| B | 1.50 | 1000 | 1.547 | 1.479 | 1.046 |
| B | 1.50 | 1500 | 1.627 | 1.617 | 1.006 |
| B | 1.75 | 500  | 1.581 | 1.168 | 1.354 |
| B | 1.75 | 1000 | 1.836 | 1.479 | 1.241 |

K over the **measured** corners (the HR 1.50 cells; the HR 1.75 cells were costed at the HR 1.50
corner and their ratio absorbs that assumption) is **1.090**, per-cell 1.006–1.256. Block A 1.275,
Block B 1.155, all eleven 1.184. The kickoff's prior ratios — A 1.27, B 0.92 — are
realized-over-**checkpoint** for B; realized-over-**Gate 1** is what is tabulated here, and A's
1.275 reproduces the quoted 1.27. **K = 1.275 (Block A's, the largest block-level value) is the
conservative multiplier used below.**

## Part A projection

| cell | basis | median (h) | 90% band | calibrated ×1.275 (h) |
|---|---|---|---|---|
| C 12.4% n 500  | measured | 0.2688 | 0.2611–0.2770 | 0.3427 |
| C 12.4% n 1000 | pooled draws | 0.1917 | 0.1843–0.1988 | 0.2443 |
| C 12.4% n 1500 | measured | 0.1146 | 0.1099–0.1195 | 0.1462 |
| C 31% n 500    | measured | 0.5696 | 0.5549–0.5840 | 0.7262 |
| C 31% n 1000   | pooled draws | 0.5328 | 0.5202–0.5454 | 0.6793 |
| C 31% n 1500   | measured | 0.4958 | 0.4849–0.5067 | 0.6322 |
| **Block C, six cells** | | **2.173** | 2.115–2.231 | **2.771** |

Reference to correct: the original Gate 1 put Block C at **2.17 h** — reproduced exactly on the
uncalibrated median, and revised up to **2.77 h** once Block A's realized calibration is applied.

**The deferred Block B cell (HR 1.75, n 1500, 31%)** is anchored on realized Block B walls two
ways, which agree to the second: (a) the n-profile within HR 1.75,
`wall(1.75,1500) = wall(1.75,1000) × wall(1.50,1500)/wall(1.50,1000)` = 6950 s; (b) the HR-profile
within n = 1500, `wall(1.50,1500) × wall(1.75,1000)/wall(1.50,1000)` = 6950 s. **1.931 h.**
References: the checkpoint's by-n re-projection 2.225 h (corrected down), the original Gate 1
1.617 h (corrected up).

| | h |
|---|---|
| deferred Block B cell | 1.931 |
| Block C, six cells (calibrated) | 2.771 |
| **PART A TOTAL** | **4.701** |
| ceiling | 9.000 |
| headroom | 4.299 (48%) |
| room left under the 12 h timeout for Part B's 1.5 h cap | 7.299 |

**GATE 1: GO — all seven cells run, none deferred.**

---

# Tooling

Everything is in `quarto/simulations/gbsg_020/scripts_dinamr/` and committed; nothing was left in
the session scratchpad. The committed drivers were reused verbatim (`campaign.sh`, `render.sh`,
`gate2.R`, `project.R`, `probe.sh` — `diff` clean against the copies actually executed).

| file | status | role |
|---|---|---|
| `blockC.cells` | **new** | the six Block C cells, in the kickoff's run order |
| `blockB_deferred.cells` | **new** | the one cell the Block A checkpoint deferred |
| `projectC.R` | **new** | Gate 1 for these seven: bootstrap over the per-replicate cost distribution, calibrated on realized walls |
| `grfprobe.sh` | **new** | the five Part B GRF corners, each under `/usr/bin/time -l` with a process-tree RSS sampler |
| `grfprobe.R` | **new** | Gate 3 plus the Part B cost readout |
| `grf_mechanism.R` | **new** | the empty-band diagnostic (below) |
| `gate2.R` | **one-line correction** | `fscomp()` sent Block B's n = 500 cell to `e1stud` unconditionally; `e1stud` was run only at HR 1.50 and 1.75, so at HR 1.00 that named a nonexistent file and the Amendment 3 assertion would have gone unevaluated on one of the six Block C cells. Now `cert20` there, which is on disk and criterion-matched. Regression: `gate2.R A` still reads **204 passes, 0 failures**. |

`stage1_checks.R` still carries the superseded bonf-vs-raw comparison, labelled as such; every
gate below used `gate2.R`'s corrected identity
(`log(fld_Hc_est2_s) + fld_Hc_lam_mean_s == log(fld_Hc_est2) + fld_Hc_lam_mean`).

One side issue, **not fixed**: the committed `render.sh` reads `$SP` for its log directory but
`campaign.sh` sets `SP` without exporting it, so a clean-environment run would `mkdir -p /logs`
and fail under `set -e`. The Part A run exported `SP` explicitly rather than editing either
script.

---

# PART A — results

*(filled in below as the cells land)*

---

# PART B — GRF cost probes

*(filled in below)*
