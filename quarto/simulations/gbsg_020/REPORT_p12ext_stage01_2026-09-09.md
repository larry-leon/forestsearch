# REPORT — campaign `p12ext`, Stage 0 and Stage 1 (smoke + Gate 1)

Date: 2026-09-10 (task document dated 2026-09-09). Executor: Claude Code, unattended.
Spec: `dev/tasks/TASK_p12ext_2026-09-09.md` (committed bbe33f66 before any work).
Host: `pop-os`, AMD Threadripper PRO 5995WX, 128 logical CPUs, 251 GB RAM. R 4.6.1.

## Stage 0 — record

### 0a. HEAD and installed package

- Repository HEAD at kickoff: **585d0949** (`feature/glm-extension`). `git diff --stat HEAD -- R/` empty: no working-tree `R/` change.
- After the task-document commit HEAD is **bbe33f66**; that commit touches `dev/tasks/` only, so the `R/` tree is unchanged from 585d0949.
- Installed `forestsearch` **0.3.5** at `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6`.
- Installed-vs-HEAD by `deparse()`:
  - `fs_mr_inference`: **MATCH** (17144 chars both sides)
  - `.fs_mr_field_complement`: **MATCH** (5611 chars both sides)
  - `.fs_mr_field_recovery`: **MATCH** (1969 chars both sides)
  - full sweep: **666/666 functions identical**.
- **No install was performed** — the installed package already matches HEAD. `load_all()` was not used at any point.

### 0b. Does `field_recovery = TRUE` reach the gate?

**Yes.** Verified from source, not from behaviour alone:

- `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:620` —
  `mr_field_recovery <- identical(.env_chr("FS_S7_FIELD_RECOV", "FALSE"), "TRUE")`
- template line 641 — `field_recovery = mr_field_recovery` inside `mr_inference_args`
- `R/forestsearch_main.R:3448` — `field_recovery = .g_mr(mr_inference_args$field_recovery, FALSE)`
- `R/fs_mr_inference.R:570` — formal `field_recovery = FALSE`; `:895` `rec_on <- isTRUE(field_recovery)`;
  `:1011` `field$recovery <- .fs_mr_field_recovery(kept, G_out, sel, Nall)`
- template `:1176-1187` — recorder reads `f$recovery` into the `fld_recov_*` columns.

Confirmed at run time: the smoke bundles carry `meta$field_recovery = TRUE` and all nine columns populated.

**Recorder columns filled (9):**

`fld_recov_sens_H`, `fld_recov_ppv_H`, `fld_recov_sens_Hc`, `fld_recov_npv_Hc`,
`fld_recov_q10`, `fld_recov_q50`, `fld_recov_q90`, `fld_recov_share1`, `fld_recov_n_used`.

`sens_H` / `ppv_H` / `sens_Hc` / `npv_Hc` are the membership-agreement cross-tab of each outer
draw's re-selection against the **observed** subgroup, averaged over used draws; `q10/q50/q90`
and `share1` describe the per-draw containment a_r/|Ĥ| behind that mean; `n_used` is the number
of outer draws contributing. The campaign runs **with** the knob.

### 0c. The committed `tier2` meta (quoted verbatim)

From `results/fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n500_tier2_combined_1_2000.rds`, `$meta`:

```
n_sample              : int 500
n_sims                : int 2000
nb_boots              : int 0
mr_draws              : int 5000
subgroup_method       : chr "consistency"
sim_id_start          : int 1
sim_id_end            : int 2000
seed_base             : int 8316951
n_batches             : int 2
sg_focus              : chr "maxeffCons"
focus_tag             : chr "maxeffCons"
ci_method             : chr "field"
field_complement      : logi TRUE
field_decompose       : logi TRUE
field_scale_complement: chr "selected"
ij_residual           : chr "two_term"
campaign_tag          : chr "tier2"
target_hr_harm        : num 1.75
harm_z1_quantile      : num 0.25
harm_prevalence_super : num 0.124
effect_neighborhood   : num 0.1
er_jcuts              : int 10
fb_mode_by_batch      : chr [1:2] "none" "none"
consistency_method    : chr "resample"
stop_threshold        : chr "NULL"
forestsearch_version  : chr "0.3.5"
```

`FS_S7_Z1Q` is unset (`harm_z1_quantile = 0.25`, the M1 default, `harm_prevalence_super = 0.124`);
`sg_focus = maxeffCons`; `er_jcuts = 10`; `seed_base = 8316951`; two batches of 1,000.
The five new cells are this configuration at new HR and n; the only meta additions are
`field_recovery = TRUE` (the knob did not exist when `tier2` ran) and `campaign_tag = "p12ext"`.

**Recorded for the timing section, not for configuration:** the four committed `tier2` cells were
produced on `Mac-Studio-3.local` at **12 workers** under R 4.5.2. Their per-replicate seconds
(15.7–19.1 s) are therefore **not** a timing reference for this host and are not used below.

## Stage 1 — 5-replicate smokes

Two renders, campaign tag `p12extsmoke`, every knob explicit and identical to the production
setting except `FS_S7_NSIMS=5`:

```
FS_S7_HR=1.50 FS_S7_N={500,1500} FS_S7_NSIMS=5 FS_S7_START=1 FS_S7_WORKERS=100
FS_S7_FOCUS=maxeffCons FS_S7_ER_JCUTS=10
FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE
FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_RETURN_RESEL=TRUE
FS_S7_FB=none FS_S7_CAMPAIGN=p12extsmoke
```

`FS_S7_Z1Q` was **not** set. Both renders exited 0.

Resolved meta on both smokes: `hr=1.50`, `z1q=0.25`, `prev_super=0.12418`, `focus=maxeffCons`,
`J=10`, `field_complement=TRUE`, `field_decompose=TRUE`, `field_scale_complement=selected`,
`field_recovery=TRUE`, `ij_residual=two_term`, `fb=none`, `workers=100`,
`forestsearch_version=0.3.5`, `seed_base=8316951`.

Truth targets (both n, HR 1.50): `marg_H = 1.509`, `marg_Hc = 0.6569`, `hr_causal = 0.7041`,
`cde_H = 1.671`, `cde_Hc = 0.585`.

### 1a. Completeness and finiteness

Both cells: 5/5 rows, `sim_id` 1–5, `detected = TRUE` 5/5, `mr_ok` 5/5, `status = DETECTED`,
`err_msg` all NA, `fld_H_note` / `fld_Hc_note` all NA, `ij_source = "ij"` on every row.

Non-finite counts on detected rows, by block (identical in both cells):

| block | cols | non-finite |
|---|---|---|
| naive | 8 | 0 |
| oracle | 8 | 0 |
| harm field | 18 | 0 |
| complement field | 9 | 0 |
| complement field-s | 9 | 0 |
| decompose (ρᶜ) | 4 | 0 |
| joint | 9 | 0 |
| joint-s | 9 | 0 |
| IJ two-term | 8 | 0 |
| β(Ĥ), β(Ĥᶜ) | 2 | 0 |
| p̂ | 3 | 0 |
| recovery | 9 | 0 |

`fld_H_lo2u` / `fld_H_hi2u` (the uniform band) are NA in both cells, as is
`fld_H_uniform_secs` — expected and correct: `FS_S7_UNIFORM` is unset, so `field_uniform = FALSE`
and the uniform block is not computed. It is excluded from the block above and from Gate 2.

### 1b. Interval invariants (detected rows, both cells)

`lo <= hi` and `lo <= est <= hi` hold 5/5 on: harm field two-sided, complement field two-sided,
complement field-s two-sided, IJ two-term harm, IJ two-term complement, naive harm, naive
complement. One-sided ordering holds: harm one-sided lower >= two-sided lower 5/5; complement
one-sided upper <= two-sided upper 5/5. Joint lower on Ĥ <= harm two-sided upper 5/5.

### 1c. Bound <-> quantile identities

The field's bounds are pivotal in the Λ* quantiles about `beta_deb`, and `est2 = to_eff(beta_deb −
lambda_mean)` (`R/fs_mr_inference.R:929-933`, `:1296-1313`), so the identity to check is

```
log(bound) = log(est2) + lambda_mean − q
```

Max absolute deviation on detected rows:

| identity | n = 500 | n = 1500 |
|---|---|---|
| harm lower two-sided vs q975 | 1.1e-16 | 5.6e-17 |
| harm upper two-sided vs q025 | 1.9e-16 | 1.1e-16 |
| harm lower one-sided vs q95 | 5.6e-17 | 1.1e-16 |
| complement lower two-sided vs q975 | 1.1e-16 | 1.1e-16 |
| complement upper two-sided vs q025 | 1.5e-16 | 1.1e-16 |
| complement upper one-sided vs q05 | 1.3e-16 | 5.6e-17 |
| complement lower one-sided vs q95 | 1.1e-16 | 1.7e-16 |
| harm SE-interval vs est2 ± z·se | 1.1e-16 | 1.7e-16 |
| complement field-s SE-interval vs est2_s ± z·se_s | 1.1e-16 | 2.2e-16 |
| IJ two-term harm lo/hi vs est ± z·se_ij | 1.1e-16 | 1.1e-16 |
| IJ two-term complement lo/hi vs est ± z·se_ij | 1.1e-16 | 1.1e-16 |

All at machine epsilon, well inside the 1e-12 tolerance.

### 1d. γ and the joint probability

| | n = 500 | n = 1500 |
|---|---|---|
| `fld_joint_gamma` | 0.025 (all 5) | 0.025 (all 5) |
| `fld_joint_s_gamma` | 0.025, 0.025, 0.025, 0.026, 0.025 | 0.025, 0.026, 0.025, 0.025, 0.025 |
| `fld_joint_prob` | 0.9508 0.9499 0.94835 0.95195 0.9509 | 0.9510 0.9510 0.9510 0.9499 0.9499 |
| `fld_joint_n` | 996 998 968 999 998 | 1000 1000 1000 998 998 |
| γ ∈ [0.025, 0.05] | TRUE | TRUE |
| joint prob >= 0.95 − 2/n_joint | TRUE | TRUE |
| `fld_joint_corr` | 0.076, 0.059, 0.088, 0.142, −0.113 | 0.049, 0.087, 0.046, −0.144, −0.043 |

### 1e. Realized prevalence against the M1 default (~0.124)

`n_true` is the planted harm region's size in the trial; `n_harm` = `n_sel` = |Ĥ|.

| | replicate values of `n_true`/n | mean |
|---|---|---|
| n = 500 | 0.1440 0.1380 0.1220 0.1520 0.1140 | **0.1340** |
| n = 1500 | 0.1373 0.1227 0.1200 0.1373 0.1127 | **0.1260** |

Super-population value in meta: **0.12418** on both. Five replicates give a binomial SE of
0.124·0.876/500 -> 0.0148 at n = 500 and 0.0085 at n = 1500 per replicate, so a 5-replicate mean
has SE ≈ 0.0066 (n = 500) and 0.0038 (n = 1500). The realized means sit +1.5 and +0.5 SE of
0.12418. Consistent with the M1 default; no shift.

|Ĥ| share for reference: n = 500 mean 0.1420; n = 1500 mean 0.1180.

### 1f. ρᶜ and the decompose columns

| | n = 500 | n = 1500 |
|---|---|---|
| `fld_Hc_scale_sel` (= naive robust SE) | 0.1347 0.1323 0.1383 0.1359 0.1422 | 0.0804 0.0775 0.0781 0.0810 0.0800 |
| `fld_Hc_scale_win` | 0.1338 0.1316 0.1360 0.1352 0.1417 | 0.0789 0.0762 0.0775 0.0795 0.0786 |
| `fld_Hc_scale_cv` | 0.0212 0.0242 0.0271 0.0327 0.0196 | 0.0291 0.0202 0.0199 0.0258 0.0192 |
| **ρᶜ** `fld_Hc_scale_ratio` | 1.007 1.006 1.017 1.005 1.003 | 1.019 1.017 1.009 1.019 1.018 |

All finite, all near 1.

### 1g. Recovery columns (first record at HR 1.50)

| column | n = 500 | n = 1500 |
|---|---|---|
| `fld_recov_sens_H` | 0.3099 0.4767 0.3596 0.2829 0.5864 | 0.5010 0.1798 0.4649 0.6956 0.6511 |
| `fld_recov_ppv_H` | 0.3560 0.4531 0.3558 0.2412 0.5804 | 0.4792 0.1587 0.4241 0.7808 0.7074 |
| `fld_recov_sens_Hc` | 0.8943 0.9099 0.8941 0.8669 0.9294 | 0.9343 0.8819 0.9231 0.9686 0.9599 |
| `fld_recov_npv_Hc` | 0.8658 0.9160 0.8970 0.8868 0.9361 | 0.9401 0.8938 0.9388 0.9526 0.9471 |
| `fld_recov_q10` | 0.0833 0 0 0 0 | 0.0313 0 0.1126 0.0442 0.3184 |
| `fld_recov_q50` | 0.1905 0.5217 0.2319 0.2537 0.5000 | 0.4688 0 0.3444 0.7427 0.6318 |
| `fld_recov_q90` | 1 0.9275 1 0.6866 1 | 1 0.7647 1 1 1 |
| `fld_recov_share1` | 0.1255 0.0231 0.1963 0.0591 0.3958 | 0.2160 0.0110 0.1330 0.3076 0.3487 |
| `fld_recov_n_used` | 996 998 968 999 998 | 1000 1000 1000 998 998 |

`0 <= sens_H, ppv_H <= 1` on every row in both cells.

### 1h. p̂ validity

`p_hat_H` n = 500: 0.1322 0.4736 0.1708 0.1188 0.3940; n = 1500: 0.1716 0.2628 0.1798 0.3092 0.3590.
`p_hat_sum` <= 1 on every row (0.961–1.000). `0 <= p_hat_H <= 1` on every row. Valid in both cells.

## Gate 1 — projection and decision

### The reference used

The smokes run 5 replicates on 100 declared workers, so they are a **light-load** measurement:
per-replicate `fit_mr_secs + fld_H_secs + fld_Hc_secs` means **33.8 s** (n = 500) and **54.7 s**
(n = 1500). Per the task these are not the projection basis.

The projection is taken from **committed `pop-os` bundles at 100 workers in the same
configuration family** — `sg_focus = maxeffCons`, M1 default prevalence (no `z1q` tag),
`er_jcuts = 10`, `effect_neighborhood = 0.10`, `field_complement = TRUE`, `fb = none`, 1,000
replicates per batch. These are loaded references by construction: the contention is already in
the numbers.

| loaded reference (2 batches each) | n | HR | s/replicate |
|---|---|---|---|
| `map1c` / `map1w` | 500 | 1.50 | 63.6, 63.6, 63.6, 64.0 |
| `s7c` / `s7w` | 500 | 1.00 | 58.6, 58.9, 58.5, 59.4 |
| `map1c` / `map1w` | 1000 | 1.00 | 107.8, 106.7, 108.0, 107.2 |
| `map1c` / `map1w` | 1500 | 1.50 | 211.5, 211.8, 210.3, 210.5 |

Implied contention against this campaign's own light-load smokes: **64.0/33.8 = 1.89×** at
n = 500 and **211.8/54.7 = 3.87×** at n = 1500 — bracketing the 2.0× / 3.1× figures quoted from
`REPORT_cert20`, and in the same direction (contention grows with n). Because three of the five
cells have a directly matching loaded reference, the factor is used only for the two that do not.

Two cells have no direct reference:

- **HR 1.50, n = 1000.** The n = 1000/n = 500 cost ratio at HR 1.00 is 107.8/58.6 = **1.84**;
  applied to HR 1.50 n = 500 (63.6) gives 117.0 s. The HR 1.50/HR 1.00 ratio at n = 500 is
  63.6/58.6 = 1.086, applied to HR 1.00 n = 1000 (108.0) gives 117.3 s. Both routes agree:
  **118 s** used.
- **HR 1.00, n = 1500.** The HR 1.00/HR 1.50 cost ratio at n = 1500 in the `cert20` cells (the
  only Linux pair at that n) is 237.4/265.1 = **0.896**; applied to HR 1.50 n = 1500 (211.8)
  gives 189.8 s. **195 s** used (rounded up).

### The projection

Wall per batch of 1,000 at 100 workers = 10 x (s/replicate); each cell is two batches plus a
combine render. Render/DGM overhead observed at ~1.5 min per render.

| # | cell | s/rep | batch wall | 2 batches | + combine/overhead | cell wall |
|---|---|---|---|---|---|---|
| 1 | HR 1.50 n500 | 64.0 | 10.7 min | 21.3 min | 4.5 min | **25.8 min** |
| 2 | HR 1.50 n1000 | 118 | 19.7 min | 39.3 min | 4.5 min | **43.8 min** |
| 3 | HR 1.50 n1500 | 211.8 | 35.3 min | 70.6 min | 4.5 min | **75.1 min** |
| 4 | HR 1.00 n1000 | 108.0 | 18.0 min | 36.0 min | 4.5 min | **40.5 min** |
| 5 | HR 1.00 n1500 | 195 | 32.5 min | 65.0 min | 4.5 min | **69.5 min** |
| | **total** | | | | | **254.7 min = 4.25 h** |

Adding the Stage 3 `summary_p12ext` render (~5–10 min on nine bundles) the total is **~4.4 h**.

This projection is conservative in one respect that can be checked against the record: for the
`map1c` n = 1500 HR 1.50 pair the two batch bundles were written 28.8 min apart, against the
10 x 211.5 s = 35.3 min the same arithmetic predicts — the realized wall ran **0.82x** the
projection. No correction is applied; the conservative figure stands.

The loaded references have `field_decompose = FALSE` and no field-s or recovery block, which this
campaign turns on. `HANDOFF_mr_field_linux_2026-09-08.md` §6 records these as costing nothing
measurable (fit+MR 36.1 vs 36.3 s; complement block 2.0 vs 1.9 s). The projection assumes that
holds; the driver re-projects from realized walls after every batch regardless.

### Decision

**Total projection 4.25 h (4.4 h with the report render) <= 8 h. Gate 1: GO.**
Hard timeout 10 h from Stage 2 start. Cell 5 (HR 1.00 n = 1500) is deferred first if the ceiling
is threatened; the driver re-projects each cell from realized walls of the same n before starting
it.
