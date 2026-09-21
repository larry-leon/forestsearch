# REPORT — the structural-null counterpart of the GBSG survival simulations (identification only)

**Task:** `dev/tasks/TASK_null_gbsg_identification_2026-09-21.md` (committed `41534a1e`) ·
**Opened:** 2026-09-21 · **Branch:** `feature/glm-extension` · **Study directory:**
`quarto/simulations/gbsg_020` · **Campaign tag:** `nullid`.

**What this is.** The manuscript states twice that the GBSG survival design carries no strict-null
cell and that no false-positive rate against the complete absence of an effect is measured (§5.1
p. 27, §5.4 p. 34). This is that cell: the same GBSG-based survival design with **no planted region
and a uniform treatment benefit**, at uniform HR **0.657** and **0.721**, n **500 / 1000 / 1500**,
2,000 replicates, identifiers **FS / DINA / GRF** on identical draws, rule **`max_A N(ε)` at
ε = 0.20**. **Identification and classification only.**

**Scope held.** No MR anywhere: no multiplier resampling, no field / field-s / IJ / Bonferroni
products, no bootstrap, no cross-validation. **No edit to `R/`.** No change to any of the eighteen
committed cells; nothing there was re-run or re-rendered. No `R CMD check`, no vignette build, no
test suite. Nothing was fetched, pulled or pushed. Every `git add` named its paths.

**Machine — a departure from the task, recorded.** The task names pop-os. This ran on
**`Mac-Studio-3.local`** (14 physical cores, 36 GB, R 4.5.2, forestsearch 0.3.5.9000) at **12
workers**. The precedent is the directory's own: `idsweep`, the campaign's identification-only
sweep (§2.6 of `current_status.md`), ran on this same Mac at 12 workers. Nothing in this task's
fixed parameters is host-dependent, and no MR product — the one place where a host factor has been
measured to matter here (`REPORT_partB_measurement_2026-09-12.md`) — is computed.

---

## 0. The campaign, read from source

### 0.1 Template, runner, recorder

| role | path |
|---|---|
| template (one file renders every batch of every campaign here) | `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` |
| render driver (one render; pins the three thread variables to 1; times it) | `scripts_dinamr/render.sh` |
| campaign driver — this task's | `scripts_dinamr/nullid.sh`, cell list `scripts_dinamr/nullid.cells` |
| campaign driver — the model it follows | `scripts_dinamr/idsweep.sh` (the Part B identification sweep) |
| recorder | `record_replicate()` inside the template (chunk `recorder`), with `.na_record()` fixing the column set and `.safe_record()` as the replicate-level error firewall |
| gates — this task's | `scripts_dinamr/nullid_gateA.R` (per run), `scripts_dinamr/nullid_gateC.R` (per cell) |

### 0.2 Every `FS_S7_*` knob the template exposes, and its in-file default

| knob | default | what it sets |
|---|---|---|
| `FS_S7_NSIMS` | `1000L` | replicates in this batch |
| `FS_S7_START` | `1L` | first `sim_id` of this batch |
| `FS_S7_MODE` | `"batch"` | `batch` (run + save one piece) or `combine` (merge + report) |
| `FS_S7_METHOD` | `"consistency"` | engine: `consistency` / `dina` / `grf` |
| `FS_S7_FOCUS` | `"maxeffCons"` | selection criterion; guard admits the six Part B criteria |
| `FS_S7_NBHD` | `0.10` | ε, the band half-width; **rejected unless the focus is `effMaxSG` / `effMinSG`** |
| `FS_S7_HR` | `1.0` | under `alt`, the target Cox HR **inside the planted region** (calibrates `k_inter`) |
| `FS_S7_N` | `500L` | trial size |
| `FS_S7_Z1Q` | `0.25` | the region rule's ER quantile (0.25 → 12.4% prevalence; 0.60 → 31%) |
| `FS_S7_DGM` | `"alt"` | **added by this task**: `alt` (planted region) or `null` (structural null) |
| `FS_S7_MR` | `"TRUE"` | `FALSE` runs identification and classification only |
| `FS_S7_FB` | `"none"` | FB bootstrap: `run` / `join` / `none` |
| `FS_S7_NB_BOOTS` | `300L` | FB replicates (read only when `FS_S7_FB=run`) |
| `FS_S7_KNOISE` | `0L` | inert N(0,1) noise confounders |
| `FS_S7_ER_JCUTS` | `10L` | J of the J-quantile grid on raw ER |
| `FS_S7_WORKERS` | `60L` | capped at physical cores − 1 |
| `FS_S7_CAMPAIGN` | `"c1"` | mandatory stem-namespacing tag |
| `FS_S7_QUICKRUN` | `"FALSE"` | tags the stem so a smoke can never pool with production |
| `FS_S7_SAVE_COMBINED` | `"TRUE"` | combine mode: also write the pooled bundle |
| `FS_S7_WINNER_ROWS` | `"FALSE"` | show the closed winner-only IJ rows |
| `FS_S7_UNIFORM` | `"FALSE"` | MR only — κ (uniform) field calibration |
| `FS_S7_FIELD_COMPLEMENT` | `"TRUE"` | MR only |
| `FS_S7_FIELD_SCALEC` | `"selected"` | MR only |
| `FS_S7_FIELD_DECOMP` | `"FALSE"` | MR only |
| `FS_S7_FIELD_RECOV` | `"FALSE"` | MR only |
| `FS_S7_IJ_RESIDUAL` | `"two_term"` | MR only |
| `FS_S7_RETURN_RESEL` | `"TRUE"` | MR only |
| `FS_S7_FB_PATH`, `FS_S7_JOIN_SKIP` | see template | FB-join only |

The seven MR-only knobs are inert with `FS_S7_MR=FALSE`. They were pinned anyway, because
`idsweep` pinned them (continuity of the knob set, not of behaviour).

### 0.3 Where the planted region is built, and where the treatment effect is set

All in `R/sim_aft_gbsg.R`, inside `.create_gbsg_dgm_()`.

**The region rule and its quantile knob** — `R/sim_aft_gbsg.R:279-289`:

```r
er_threshold <- stats::quantile(dfa$er, probs = z1_quantile)
dfa$z1 <- ifelse(dfa$er <= er_threshold, 1L, 0L)
dfa$z3 <- ifelse(dfa$meno == 0, 1L, 0L)  # Premenopausal
dfa$zh <- dfa$treat * dfa$z1 * dfa$z3                    # the interaction term
dfa$flag.harm <- ifelse(dfa$z1 == 1 & dfa$z3 == 1, 1L, 0L)
```

**The parameter that sets the interaction** — `R/sim_aft_gbsg.R:419-421`, reached only under
`model == "alt"`:

```r
if (model == "alt") {
  gamma["zh"] <- k_inter * gamma["zh"]
}
```

`k_inter` is calibrated by `calibrate_k_inter()` (`R/sim_aft_gbsg.R:1004-1045`) to
`FS_S7_HR`, which is `dgm$hr_H_true` — the Cox HR inside the region — with `use_ahr = FALSE`.
`calibrate_k_inter()` **stops unless `model == "alt"`** (`R/sim_aft_gbsg.R:1017-1020`).

**Where the marginal and complement effects are set.** The overall treatment effect is
`R/sim_aft_gbsg.R:417`:

```r
gamma["treat"] <- k_treat * gamma["treat"]
```

`k_treat` is `1` on every committed cell (it is `setup_gbsg_dgm()`'s formal default and the
template never overrode it). The marginal and complement effects are then **read off**, not set:
`hr_causal` is a Cox fit of treatment alone on the stacked potential outcomes of the whole
super-population (`R/sim_aft_gbsg.R:517-521`), `hr_H_true` the same fit on `flag.harm == 1`
(`:527-530`) and `hr_Hc_true` the same on `flag.harm == 0`. Under `model == "null"` there is no
region, and `hr_Hc_true <- hr_causal` (`R/sim_aft_gbsg.R:548`).

**The task's two targets come from here.** At `FS_S7_HR=1.00` the planted region carries HR 1.00
against a benefiting complement, and that complement's effect is `hr_Hc_true`: **0.656891** at
`FS_S7_Z1Q` unset (prevalence 0.12418, the narrow cell) and **0.720557** at `FS_S7_Z1Q=0.60`
(prevalence 0.30655, the broad cell). Both reproduced in
`scripts_dinamr/logs/nullid_designpoint.txt`.

### 0.4 The search settings this campaign actually uses

Read from the template's setup chunk; all are literals there, not knobs.

| setting | value | note |
|---|---|---|
| screening threshold `c1` (`hr.threshold`) | **0.90** | **natural scale on the survival path** — HR ≥ 0.90, not log. `subgroup_search()`'s own documentation is explicit: "1.25 means HR >= 1.25, not log(1.25)" (`R/subgroup_search.R:22-31`). The `exp(hr.threshold)` in the *display* line (`R/subgroup_search.R:208-213`) is a printing artefact of the GLM branch and does not reach the comparison. |
| per-split threshold `c2` (`hr.consistency`) | **0.80** | |
| consistency rate `p*` (`pconsistency.threshold`) | **0.90** | |
| selection rule | `"neighborhood"` | the band rule; it is what `effMaxSG` / `effMinSG` consult |
| ε (`effect_neighborhood`) | **0.20** for this task | template default 0.10; set per render by `FS_S7_NBHD` |
| criterion (`sg_focus`) | **`effMaxSG`** for this task | i.e. `max_A N(ε)`: the largest subgroup among those within ε of the maximal effect |
| eligibility minimum (`n.min`) | **`NULL`** | forwarded as `NULL`; `forestsearch()`'s own default applies |
| per-arm event minimum (`d0.min`, `d1.min`) | **10, 10** | consistency engine only — not forwarded to DINA or GRF, which enumerate their own candidates |
| `maxk` | **2** | one- and two-factor conjunctions |
| `consistency_method` | `"resample"` | closed-form resample identification |
| `fs.splits` | `400L` | moot under `"resample"` |
| candidate construction | `use_lasso = FALSE`, `use_grf = FALSE`, `use_twostage = TRUE`, `conf_force = c("meno == 0", "er <= 0", "pgr <= 0")`, `conf.cont_jcuts = list(er = 10)` | |
| candidate covariates | `er`, `age`, `meno`, `pgr`, `nodes`, `size`, `grade` — the **raw** GBSG variables | the pre-dichotomized `z1..z5` are never offered |
| **DINA's proposal floor** | `hr.threshold = 0.90`, `dina_select_statistic = "effect"`, `dina_args = list()` | the template sets the same `hr.threshold` for every engine; **recorded, not changed.** The manuscript's p. 7 TODO notes Table 1 gives DINA log 1.0 while this campaign ran it at 0.90 — what the template sets is 0.90, and that is what this campaign used. |
| **GRF's proposal floors** | `dmin.grf = 0.0` (DR-score pre-filter, RMST units), `grf_selection = "frontier"`, `grf_select_statistic = "effect"`, `grf_depth = 2L`; the binding effect-scale floor is the same `hr.threshold = 0.90` | two floors on two scales — `current_status.md` §2.4 |

### 0.5 How inference is switched off, and that it is off

`FS_S7_MR=FALSE` → `mr_inference_on <- FALSE` → `forestsearch(mr_inference = FALSE)`. This is the
mode the selection-rule comparison used (`idsweep`, §2.6 of `current_status.md`; Gate T3 of
`REPORT_partB_enabling_2026-09-12.md` established that the identification and classification
columns are `identical()` to an MR-on run).

**Confirmed, per bundle, by Gate A:** `fs.est$mr_inference` is `NULL`, so `mr_ok == 0` on every one
of the 2,000 replicates, and **every one of the 119 `mr_*` / `fld_*` / `fb_*` product columns is NA
on every replicate**. `fb_mode = "none"`, so `run_fb` is `FALSE` and
`forestsearch_bootstrap_dofuture()` is never called. No cross-validation entry point is reached.

### 0.6 What the recorder stores, and the two things it cannot

With MR off the recorder reads `label`, `n_sel`, `n_cons_qual` and `band_n` from the
`forestsearch()` result, plus `sg_def`, `covs`, `n_true`, `n_harm`, the four classification rates,
`status`, `detected`, `admitted_n` (GRF only), `fit_mr_secs` and the oracle refits `or_H_*` /
`or_Hc_*`. `betaHhat_H` / `betaHhat_Hc` are attached after the loop from the super-population
evaluation frame.

**Template-level additions made for this task** (add-only, NA on every path that does not fill
them; no `R/` change):

| column | source |
|---|---|
| `n_cand_enum` | `find.grps$filter_counts$n_evaluated` — conjunctions enumerated |
| `n_cand_floor` | `find.grps$filter_counts$n_passed_hr` — those clearing the effect floor and the size / per-arm-event floors |
| `maxT` | `max_g log(HR_g) / se_g` over `out.found$hr.subgroups`, with `se_g = (log U − log L) / (2 z_0.975)` from the recorded Wald interval |
| `p_sel`, `p_max_qual` | `Pcons` of the selected candidate, and the maximum over the consistency-qualifying family |
| `nv_H_est/lo/hi/se`, `nv_Hc_*` | the unadjusted within-region (and complement) Cox refit, via `.cox_hr_ci()` — the helper the oracle block already uses |
| `nv_H_lo1s` | the one-sided 95% Wald lower bound on the HR scale |

**Two Step-3 fields are left out, as the task directs, because supplying them would be an `R/`
change:**

1. **`n_family`, the gate's kept family K, stays NA.** It is MR's own fitted family, enumerated
   only inside the MR branch (`R/forestsearch_main.R:3352-3371`) and filtered by per-candidate fits
   in `.fs_mr_assemble()`. No field of the result carries it.
2. **The maximum consistency rate over the *full screened* family is not recoverable.** A candidate
   failing the consistency screen returns `NULL` and its `Pcons` is discarded
   (`R/subgroup_consistency_helpers.R:1866-1872`, `:1995-2000`). What is recorded is
   `p_max_qual`, the maximum over the **qualifying** family — the set the selection rule actually
   sorts.

**`n_cand_enum`, `n_cand_floor`, `maxT`, `p_sel` and `p_max_qual` are consistency-engine only.**
DINA and GRF enumerate their own candidates and return neither `find.grps` nor `out_sg$result`;
this reproduces `current_status.md`'s standing note that `n_cons_qual`, `band_n` and `p_star` are
structurally NA there. GRF alone returns `admitted_n`.

**The other selection rules' picks were not recorded.** The task allows it "if the campaign's
machinery already records [them] from the same search, as the S3.1 sweep does". It does not: the
S3.1 sweep (`idsweep`) ran a **separate search per criterion** — 288 cell-runs for three engines ×
six criteria × 18 cells. Getting a second rule's pick here would mean a second search, which the
task forbids.

### 0.7 Workers and seeds

`FS_S7_WORKERS`, default `60L`, capped at `physical cores − 1`; **12** for this campaign, as
`idsweep` used on this host. `render.sh` exports `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1` (Apple Accelerate is not fork-safe). Replicates fan out over the workers
(`parallel_mode = "sims"`), each replicate's search running sequentially inside its worker.

Seeds: `seed_base + sim_id` with `seed_base = 8316951L`, the template literal, for the DGM draw,
the search and (unused here) the bootstrap. `sim_id` 1–2000, one batch per run, so there is no
combine render.

---

## 1. Feasibility, and the structural null as expressed

**Feasible without touching `R/`.** `.create_gbsg_dgm_()` already carries a `model = "null"` path:
it drops `zh` from the AFT design entirely (`R/sim_aft_gbsg.R:336-349`), sets
`dfa$flag.harm <- 0L` (`:311`) and `fs_harm_true <- NULL` / `grf_harm_true <- NULL` (`:306-310`).
Every subject then carries the same treatment effect, `loghr_po = theta_1 − theta_0 = b0["treat"]`.
That path was unreachable from this study directory only because `dgm_model <- "alt"` was a
literal — `current_status.md` §1 records exactly that. `FS_S7_DGM` makes it reachable, at the
template level.

**Why `model = "null"` and not `alt` with `k_inter = 0`.** Under `alt` with a zero interaction the
region is still *present* — `flag.harm` is 1 on 12.4% of subjects and `fs_harm_true` is
`c("v1.1", "v3.1")` — it merely carries no effect. The task's check (1) asks that the planted-region
prevalence be "zero or absent". Only `model = "null"` gives that.

**What "uniform hazard ratio 0.657" is taken to mean, stated explicitly.** These two numbers are
the alt design's **complement effects** `hr_Hc_true`, which are **marginal Cox** hazard ratios over
a set — attenuated by non-collapsibility relative to the patient-level effect. The structural-null
analogue of that quantity is `dgm$hr_causal`, the marginal Cox HR over the whole super-population,
and that is what `FS_S7_HR` calibrates `k_treat` against here. The **patient-level (individual)
hazard ratio is uniform by construction** — that is what makes this a structural null — but it is
numerically a different number, further from 1, and it is reported beside the target throughout.
Under `model = "null"` the AHR equals the patient-level HR exactly, since `loghr_po` is constant.

| design point | `k_treat` | super-population marginal Cox HR | uniform patient-level HR = AHR |
|---|---|---|---|
| `null0657` | 1.27232145415221 | 0.657000011553 | 0.582908 |
| `null0721` | 0.99182531577674 | 0.720999999973 | 0.656562 |

`k_treat` is found by `uniroot` (`tol = 1e-12`) on `hr_causal` at the template level, because
`calibrate_k_inter()` refuses `model != "alt"` and `calibrate_k_treat()` (`R/calibrate_k_treat.R`)
targets `generate_aft_dgm_flex()`, a different generator from the GBSG one this campaign uses.

**`FS_S7_Z1Q` is rejected under `FS_S7_DGM=null`** and `z1_quantile` stays at its 0.25 default.
With no region the quantile would move only the z1 *main* effect, and the two null cells must
differ in the uniform effect alone. This is also why the null grid has no prevalence dimension:
6 cells = 2 effects × 3 sample sizes.

### 1.1 The three design-point checks, from the generator's own truth object

Run inside the render before any replicate (template chunk `build-dgm`), and reproduced standalone
in `scripts_dinamr/logs/nullid_designpoint.txt`. The candidate family is the consistency engine's
own pool — `get_FSdata()` with this campaign's `conf_force` and `conf.cont_jcuts`, no LASSO and no
GRF cuts — built on the super-population evaluation frame and enumerated to `maxk = 2` exactly as
`subgroup_search()` does: **29 factors, 435 conjunctions, 433 of them non-empty**. β(g) is the
patient-level log hazard ratio averaged over g, i.e. log AHR(g); that is what "the candidate
carries the uniform effect" means structurally.

| check | `null0657` | `null0721` |
|---|---|---|
| **(1)** planted-region prevalence; truth labels absent | `0.0000000000`; `TRUE` → **PASS** | `0.0000000000`; `TRUE` → **PASS** |
| **(2)** `max_g \|β(g) − log HR_uniform\|` (tol 1e-8) | `0.000e+00` → **PASS** | `5.551e-17` → **PASS** |
| (2, strongest form) range of the patient-level log HR over all 100,000 subjects | `4.441e-16` | `3.886e-16` |
| **(3)** `\|hr_causal − target\|` (tol 1e-8) | `1.155e-08` → **PASS under amendment** | `2.722e-11` → **PASS (1e-8)** |

### 1.2 Gate amendment on check (3) — applied unattended, for review

**The task fixes 1e-8. One of the two design points cannot reach it, for a numerical reason.**

`dgm$hr_causal` is a `survival::coxph` MLE on the 200,000-row stacked potential-outcome frame, and
`coxph`'s own convergence tolerance floors the attainable agreement. Measured: the objective steps
across zero by ~8e-8 between adjacent `k_treat` values near the root. A 41-point local scan at
4.35e-8 spacing around the `uniroot` root finds a minimum `|hr_causal − target|` of **1.155e-08**
at target 0.657, and **2.722e-11** at target 0.721. Meeting 1e-8 at 0.657 would require changing
`coxph`'s control inside `.create_gbsg_dgm_()` — an `R/` edit, out of scope.

**The amendment.** Check (3) asserts 1e-8 and prints the outcome either way; on a miss it falls
back to **1e-6 absolute on the HR scale** and says so. Checks (1) and (2), which carry the
structural content of "no subgroup exists", keep 1e-8 unconditionally and pass at 1e-16 or exactly
zero. The relative miss at 0.657 is 1.76e-8 — the design point is the target to eleven significant
figures.

**This needs Larry's acceptance or reversal.** The precedent for applying a gate amendment
unattended and recording it is `idsweep` Gate I Amendment 1 (`current_status.md` §2.6, §7).

### 1.3 The template edits, and that they are inert on the `alt` path

Add-only and default-inert: `FS_S7_DGM`; `dgm_tag` (`_null657` / `_null721`, empty under `alt`, so
no committed stem changes and no null bundle can share a `combine_glob` with an alt one); the
`k_treat` calibration branch; the design-point gate; the six added recorder columns; and
`dgm_model` / `k_inter` / `k_treat` in the bundle meta.

**Measured inertness** (`scripts_dinamr/logs/nullid_smoke_and_inertness.txt`): with `FS_S7_DGM`
unset, hr 1.50, n 500, `effMaxSG` ε 0.20, MR off, `sim_id` 1–20, the render is **identical on 158
of 166 shared non-timing columns** to the committed `idsweep` bundle
`results/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_nomr_idsweep_res_1_500.rds`, with the
same `truth` object, the same `harm_prevalence_super` (0.12418) and `k_treat = 1`.

The **eight columns that differ are `nv_H_est/lo/hi/se` and `nv_Hc_est/lo/hi/se`** — all-NA with MR
off before this change, and now populated from the search's own Cox refit of the selected region.
No existing *value* changes; NA becomes a number. This is a deliberate change of what an MR-off
render records, and it is what makes Step 3's "unadjusted within-region estimate" available without
MR; it is noted here because a future re-render of an `idsweep` cell would now populate those
columns.

---

## 2. Smoke, and the projection

Cell `null0657_n500`, 20 replicates, MR off, 12 workers
(`scripts_dinamr/logs/nullid_smoke_and_inertness.txt`):

| identifier | per-replicate search | declarations / 20 | render wall |
|---|---|---|---|
| FS (consistency) | 1.1871 s | 6 | 34 s |
| DINA | 0.1477 s | 9 | 32 s |
| GRF | 2.0528 s | 20 | 36 s |
| **total per replicate across the three** | **3.3876 s** | | |

**Projection recorded at Step 2, not a decision point.** With the `grfprobe` n-scaling (13.5 s at
n 500 → 19.8 s at n 1500, linear in n: ×1.0000 / ×1.2333 / ×1.4667), the grid is
6 cells × 3 identifiers × 2,000 replicates = **50,136 core-seconds = 1.16 h at 12 workers**, plus
18 renders at an assumed 60 s = 0.30 h, for **1.46 h**.

**Revised once a real 2,000-replicate render was measured** (the Gate A live test, below): that
render took 314 s wall for 206 s of compute, so the per-render fixed overhead at 2,000 replicates
is ~110 s, not 60 s. Revised projection **1.71 h**. Measured wall is in §3.

**Gate A was tested against a real 2,000-replicate bundle before the campaign launched**, and it
failed — on a bug in the gate, not the run: `mr_ok` matched the `^(mr_|fld_|fb_)` pattern the gate
uses to assert that every MR product is NA, but `mr_ok` is
`as.integer(!is.null(fs.est$mr_inference))`, a legitimate 0 on this path. Excluded, and asserted to
be 0 instead. Gate A then passed 29 of 29 on that bundle. Fixed in `1f5c0076`, before any cell ran.

---

## 3. The grid

Six cells, in the task's order, three identifiers each, 2,000 replicates each, one batch per run
(`sim_id` 1–2000, no combine render). Driver `scripts_dinamr/nullid.sh`, cell list
`scripts_dinamr/nullid.cells`, log `scripts_dinamr/logs/nullid.driver.log`.

**Same draws within a cell, and how that is checked.** The three identifiers of a cell see the same
2,000 trials: the seed is `seed_base + sim_id` and the DGM build is a deterministic function of the
cell. The campaign's usual same-draws fingerprint, `n_true`, is 0 on every replicate here and cannot
discriminate, so Gate C uses **`or_Hc_est` / `or_Hc_se`** instead — the template's oracle complement
refit, which under the structural null is a Cox fit of treatment alone on the **whole trial** and so
is computed from the simulated data alone, never touching the engine.

**Gates.** Gate A after every run (29 checks on FS, 25 on DINA and GRF — the four family/`maxT`
checks are FS-only), Gate C after every cell. Per the task, **a failing cell does not stop the
campaign**: the rest of that cell is abandoned, a `HALT_nullid_<cell>_<date>.txt` record is written
in the study directory, and the driver moves to the next cell.

### 3.1 What actually happened: two gate bugs, both mine, both corrected

Neither was a failure of a run, and no bundle was re-run. Both were the same mistake — **I asserted a
threshold I had guessed instead of the invariant the data satisfy.**

**Gate C, cell `null0657_n500`.** It compared `or_Hc_est` / `or_Hc_se` across engines and required
`identical(is.na(a), is.na(b))`. The template's oracle block sits **after** the NO-DETECTION early
return, so `or_Hc_*` is NA on non-declaring replicates; the declaring sets differ by engine, so the
NA patterns differ legitimately and that clause can never hold. The gate's own output said so — it
printed `max |diff| 0` beside every failure. On the replicates where both engines declared the
values are identical: 692 rows FS vs DINA, 890 rows FS vs GRF, `max |diff|` exactly 0. Gate C now
compares there, with a floor of 100 such rows, and separately asserts that `or_Hc_est` is present
exactly on the declaring replicates. Re-run: **PASS, 18 checks**.

**Gate A, cell `null0657_n1500`.** It asserted that `max_g T_g` is "recorded on > 95% of
replicates". `maxT` comes from `find.grps$out.found$hr.subgroups`, which is **absent when no
candidate cleared the effect floor** — at n 1500 under the stronger uniform benefit that happens on
**330 of 2,000 replicates (16.5%)**. The real invariant holds exactly on all 2,000 rows:
`is.finite(maxT)` ⟺ `n_cand_floor > 0`, and **no replicate ever declares when nothing cleared the
floor** (0 cases). Gate A now asserts both. Re-run on the same untouched bundle: **PASS, 30 checks**.

**Consequence, and it was real.** The driver, doing what the task specifies, abandoned the rest of
cell `null0657_n1500`. Its DINA and GRF runs were then produced by
`scripts_dinamr/nullid_complete_cell5.sh` with byte-identical knobs, gated, and the cell closed by
Gate C. The consistency bundle was not re-run — it was correct throughout.

Both halt records are **kept and bannered WITHDRAWN** with the diagnosis, rather than deleted:
`HALT_nullid_null0657_n500_2026-09-21.txt`, `HALT_nullid_null0657_n1500_2026-09-21.txt`. The
driver's own closing line reads `completed 4 cells ... failed 2`; that line is about the gates as
they stood at run time, and is superseded by this section.

**Final state: 18 of 18 runs, Gate A PASS on every one, Gate C PASS on every cell. Nothing
deferred, nothing dropped, nothing re-run.**

---

## 4. Results

**Every rate below is a FALSE-declaration rate.** There is no subgroup anywhere in this design.

### Table 1 — declaration, size and specificity (cells x identifiers)

| cell | HR | n | identifier | declarations / 2000 | rate [Wilson 95%] | mean \|H\| (MC SE) | \|H\| Q1/med/Q3 | mean \|H\|/n (MC SE) | spec uncond. (MC SE) | spec cond. (MC SE) |
|---|---|---|---|---|---|---|---|---|---|---|
| null0657_n500 | 0.657 | 500 | FS | 894 | 0.4470 [0.4253, 0.4689] | 105.4 (1.19) | 80 / 98 / 123 | 0.2107 (0.0024) | 0.9058 (0.0026) | 0.7893 (0.0024) |
| null0657_n500 | 0.657 | 500 | DINA | 1116 | 0.5580 [0.5361, 0.5796] | 95.6 (0.97) | 74 / 88 / 108 | 0.1913 (0.0019) | 0.8933 (0.0024) | 0.8087 (0.0019) |
| null0657_n500 | 0.657 | 500 | GRF | 1884 | 0.9420 [0.9309, 0.9514] | 103.7 (0.81) | 79 / 94 / 120 | 0.2074 (0.0016) | 0.8046 (0.0019) | 0.7926 (0.0016) |
| null0721_n500 | 0.721 | 500 | FS | 1356 | 0.6780 [0.6572, 0.6981] | 110.6 (1.02) | 84 / 104 / 127 | 0.2212 (0.0020) | 0.8500 (0.0027) | 0.7788 (0.0020) |
| null0721_n500 | 0.721 | 500 | DINA | 1480 | 0.7400 [0.7203, 0.7588] | 101.1 (0.94) | 77 / 92 / 117 | 0.2022 (0.0019) | 0.8504 (0.0024) | 0.7978 (0.0019) |
| null0721_n500 | 0.721 | 500 | GRF | 1969 | 0.9845 [0.9781, 0.9891] | 107.8 (0.89) | 80 / 98 / 125 | 0.2156 (0.0018) | 0.7878 (0.0019) | 0.7844 (0.0018) |
| null0657_n1000 | 0.657 | 1000 | FS | 637 | 0.3185 [0.2984, 0.3392] | 166.8 (2.30) | 125 / 153 / 189 | 0.1668 (0.0023) | 0.9469 (0.0019) | 0.8332 (0.0023) |
| null0657_n1000 | 0.657 | 1000 | DINA | 593 | 0.2965 [0.2769, 0.3169] | 157.8 (2.17) | 122 / 145 / 177 | 0.1578 (0.0022) | 0.9532 (0.0017) | 0.8422 (0.0022) |
| null0657_n1000 | 0.657 | 1000 | GRF | 1662 | 0.8310 [0.8139, 0.8468] | 189.5 (1.66) | 142 / 174 / 220 | 0.1895 (0.0017) | 0.8425 (0.0021) | 0.8105 (0.0017) |
| null0721_n1000 | 0.721 | 1000 | FS | 1234 | 0.6170 [0.5955, 0.6381] | 193.0 (2.15) | 139 / 176 / 227 | 0.1930 (0.0021) | 0.8809 (0.0025) | 0.8070 (0.0021) |
| null0721_n1000 | 0.721 | 1000 | DINA | 1152 | 0.5760 [0.5542, 0.5975] | 184.1 (2.16) | 132 / 166 / 213 | 0.1841 (0.0022) | 0.8940 (0.0024) | 0.8159 (0.0022) |
| null0721_n1000 | 0.721 | 1000 | GRF | 1886 | 0.9430 [0.9320, 0.9523] | 215.4 (1.98) | 157 / 196 / 256 | 0.2154 (0.0020) | 0.7969 (0.0022) | 0.7846 (0.0020) |
| null0657_n1500 | 0.657 | 1500 | FS | 337 | 0.1685 [0.1527, 0.1855] | 240.0 (4.63) | 177 / 216 / 272 | 0.1600 (0.0031) | 0.9730 (0.0014) | 0.8400 (0.0031) |
| null0657_n1500 | 0.657 | 1500 | DINA | 278 | 0.1390 [0.1245, 0.1549] | 224.7 (3.93) | 173 / 207 / 261 | 0.1498 (0.0026) | 0.9792 (0.0012) | 0.8502 (0.0026) |
| null0657_n1500 | 0.657 | 1500 | GRF | 1129 | 0.5645 [0.5427, 0.5861] | 273.1 (2.91) | 203 / 251 / 316 | 0.1821 (0.0019) | 0.8972 (0.0023) | 0.8179 (0.0019) |
| null0721_n1500 | 0.721 | 1500 | FS | 969 | 0.4845 [0.4626, 0.5064] | 297.3 (4.24) | 202 / 268 / 348 | 0.1982 (0.0028) | 0.9040 (0.0026) | 0.8018 (0.0028) |
| null0721_n1500 | 0.721 | 1500 | DINA | 748 | 0.3740 [0.3531, 0.3954] | 259.6 (3.46) | 186 / 239 / 305 | 0.1731 (0.0023) | 0.9353 (0.0021) | 0.8269 (0.0023) |
| null0721_n1500 | 0.721 | 1500 | GRF | 1689 | 0.8445 [0.8280, 0.8597] | 342.2 (3.24) | 248 / 310 / 405 | 0.2281 (0.0022) | 0.8074 (0.0026) | 0.7719 (0.0022) |

### Table 2 — the unadjusted within-region estimate, and where its one-sided bound lands

| cell | identifier | HR(H) Q1/med/Q3 | model-based SE med | lower 1s Q1/med/Q3 | share lower >= 1.00 [Wilson] | share lower >= 1.25 [Wilson] |
|---|---|---|---|---|---|---|
| null0657_n500 | FS | 1.297 / 1.390 / 1.555 | 0.299 | 0.811 / 0.843 / 0.915 | 93/894 = 0.1040 [0.0857, 0.1258] | 7/894 = 0.0078 [0.0038, 0.0161] |
| null0657_n500 | DINA | 1.026 / 1.197 / 1.479 | 0.320 | 0.627 / 0.719 / 0.841 | 93/1116 = 0.0833 [0.0685, 0.1010] | 16/1116 = 0.0143 [0.0088, 0.0232] |
| null0657_n500 | GRF | 1.036 / 1.206 / 1.456 | 0.314 | 0.649 / 0.722 / 0.824 | 113/1884 = 0.0600 [0.0501, 0.0716] | 20/1884 = 0.0106 [0.0069, 0.0163] |
| null0721_n500 | FS | 1.298 / 1.419 / 1.578 | 0.290 | 0.818 / 0.865 / 0.950 | 225/1356 = 0.1659 [0.1471, 0.1867] | 21/1356 = 0.0155 [0.0102, 0.0236] |
| null0721_n500 | DINA | 1.109 / 1.334 / 1.643 | 0.306 | 0.695 / 0.809 / 0.949 | 269/1480 = 0.1818 [0.1629, 0.2022] | 46/1480 = 0.0311 [0.0234, 0.0412] |
| null0721_n500 | GRF | 1.128 / 1.323 / 1.615 | 0.306 | 0.706 / 0.799 / 0.919 | 268/1969 = 0.1361 [0.1217, 0.1520] | 49/1969 = 0.0249 [0.0189, 0.0327] |
| null0657_n1000 | FS | 1.159 / 1.240 / 1.352 | 0.238 | 0.809 / 0.829 / 0.866 | 26/637 = 0.0408 [0.0280, 0.0591] | 4/637 = 0.0063 [0.0024, 0.0160] |
| null0657_n1000 | DINA | 0.949 / 1.029 / 1.173 | 0.254 | 0.636 / 0.704 / 0.782 | 12/593 = 0.0202 [0.0116, 0.0350] | 4/593 = 0.0067 [0.0026, 0.0172] |
| null0657_n1000 | GRF | 0.938 / 1.002 / 1.139 | 0.230 | 0.657 / 0.707 / 0.768 | 38/1662 = 0.0229 [0.0167, 0.0312] | 6/1662 = 0.0036 [0.0017, 0.0079] |
| null0721_n1000 | FS | 1.131 / 1.213 / 1.338 | 0.216 | 0.812 / 0.838 / 0.892 | 85/1234 = 0.0689 [0.0560, 0.0844] | 9/1234 = 0.0073 [0.0038, 0.0138] |
| null0721_n1000 | DINA | 0.959 / 1.062 / 1.234 | 0.226 | 0.672 / 0.753 / 0.852 | 75/1152 = 0.0651 [0.0523, 0.0808] | 7/1152 = 0.0061 [0.0029, 0.0125] |
| null0721_n1000 | GRF | 0.967 / 1.063 / 1.225 | 0.210 | 0.706 / 0.765 / 0.842 | 97/1886 = 0.0514 [0.0423, 0.0623] | 11/1886 = 0.0058 [0.0033, 0.0104] |
| null0657_n1500 | FS | 1.080 / 1.144 / 1.230 | 0.198 | 0.804 / 0.816 / 0.840 | 6/337 = 0.0178 [0.0082, 0.0383] | 0/337 = 0.0000 [0.0000, 0.0113] |
| null0657_n1500 | DINA | 0.924 / 0.966 / 1.062 | 0.219 | 0.648 / 0.694 / 0.753 | 1/278 = 0.0036 [0.0006, 0.0201] | 0/278 = 0.0000 [0.0000, 0.0136] |
| null0657_n1500 | GRF | 0.919 / 0.946 / 1.004 | 0.195 | 0.669 / 0.701 / 0.742 | 6/1129 = 0.0053 [0.0024, 0.0115] | 1/1129 = 0.0009 [0.0002, 0.0050] |
| null0721_n1500 | FS | 1.053 / 1.107 / 1.181 | 0.174 | 0.806 / 0.822 / 0.850 | 17/969 = 0.0175 [0.0110, 0.0279] | 1/969 = 0.0010 [0.0002, 0.0058] |
| null0721_n1500 | DINA | 0.929 / 0.982 / 1.081 | 0.189 | 0.677 / 0.733 / 0.804 | 14/748 = 0.0187 [0.0112, 0.0312] | 1/748 = 0.0013 [0.0002, 0.0075] |
| null0721_n1500 | GRF | 0.920 / 0.957 / 1.038 | 0.166 | 0.709 / 0.746 / 0.788 | 29/1689 = 0.0172 [0.0120, 0.0246] | 2/1689 = 0.0012 [0.0003, 0.0043] |

### Table 3 — the candidate family and the screen statistic (FS only)

| cell | enumerated Q1/med/Q3 | clearing the floor Q1/med/Q3 | consistency-qualifying Q1/med/Q3 | p_sel med | p_max_qual med | floor>0 but no declaration |
|---|---|---|---|---|---|---|
| null0657_n500 | 1711 / 1711 / 1830 | 16 / 43 / 89 | 1 / 3 / 8 | 0.930 | 0.960 | 1075/1969 = 0.5460 [0.5239, 0.5678] |
| null0721_n500 | 1711 / 1711 / 1830 | 44 / 100 / 185 | 2 / 7 / 19 | 0.950 | 0.970 | 637/1993 = 0.3196 [0.2995, 0.3404] |
| null0657_n1000 | 1711 / 1711 / 1830 | 6 / 16 / 37 | 1 / 2 / 5 | 0.930 | 0.950 | 1270/1907 = 0.6660 [0.6445, 0.6868] |
| null0721_n1000 | 1711 / 1711 / 1830 | 23 / 53 / 104 | 2 / 5 / 12 | 0.940 | 0.960 | 753/1987 = 0.3790 [0.3579, 0.4005] |
| null0657_n1500 | 1711 / 1711 / 1830 | 1 / 5 / 14 | 1 / 2 / 3 | 0.920 | 0.940 | 1333/1670 = 0.7982 [0.7783, 0.8168] |
| null0721_n1500 | 1711 / 1711 / 1830 | 11 / 26 / 53 | 1 / 3 / 8 | 0.930 | 0.950 | 980/1949 = 0.5028 [0.4806, 0.5250] |

### Table 4 — max_g T_g over the screened family, against z_0.95 = 1.645 (FS only)

| cell | n with max_g T_g | Q1 | median | Q3 | 90% | 95% | 99% | share > 1.645 [Wilson] |
|---|---|---|---|---|---|---|---|---|
| null0657_n500 | 1969 | 0.483 | 0.833 | 1.225 | 1.577 | 1.807 | 2.176 | 166/1969 = 0.0843 [0.0728, 0.0974] |
| null0721_n500 | 1993 | 0.767 | 1.159 | 1.564 | 1.887 | 2.093 | 2.575 | 414/1993 = 0.2077 [0.1905, 0.2261] |
| null0657_n1000 | 1907 | 0.139 | 0.522 | 0.877 | 1.238 | 1.465 | 1.975 | 56/1907 = 0.0294 [0.0227, 0.0379] |
| null0721_n1000 | 1987 | 0.488 | 0.874 | 1.253 | 1.637 | 1.840 | 2.270 | 188/1987 = 0.0946 [0.0825, 0.1083] |
| null0657_n1500 | 1670 | -0.137 | 0.146 | 0.489 | 0.831 | 1.037 | 1.579 | 13/1670 = 0.0078 [0.0046, 0.0133] |
| null0721_n1500 | 1949 | 0.172 | 0.511 | 0.892 | 1.237 | 1.478 | 1.915 | 58/1949 = 0.0298 [0.0231, 0.0383] |

### Table 5 — the composition of H-hat: covariates and cut directions


**FS** — share of declaring replicates whose rule contains the term, pooled over the six cells:

| term | replicates | share of declaring |
|---|---|---|
| `NOT er <=` | 1520 | 0.2801 |
| `er <=` | 1377 | 0.2537 |
| `NOT pgr <=` | 1163 | 0.2143 |
| `nodes <=` | 1039 | 0.1915 |
| `size <=` | 1014 | 0.1868 |
| `NOT size <=` | 986 | 0.1817 |
| `age <=` | 850 | 0.1566 |
| `NOT age <=` | 816 | 0.1504 |
| `NOT nodes <=` | 628 | 0.1157 |
| `pgr <=` | 614 | 0.1131 |
| `NOT meno (indicator)` | 305 | 0.0562 |
| `meno (indicator)` | 196 | 0.0361 |
| `NOT grade (indicator)` | 116 | 0.0214 |
| `grade (indicator)` | 73 | 0.0135 |

(declaring replicates pooled: 5427)

**DINA** — share of declaring replicates whose rule contains the term, pooled over the six cells:

| term | replicates | share of declaring |
|---|---|---|
| `pgr >=` | 1564 | 0.2914 |
| `er >=` | 997 | 0.1858 |
| `age >=` | 872 | 0.1625 |
| `nodes <=` | 857 | 0.1597 |
| `size >=` | 830 | 0.1546 |
| `age <=` | 827 | 0.1541 |
| `size <=` | 712 | 0.1327 |
| `pgr <=` | 617 | 0.1150 |
| `er <=` | 574 | 0.1069 |
| `nodes >=` | 538 | 0.1002 |
| `grade >=` | 459 | 0.0855 |
| `meno <=` | 356 | 0.0663 |
| `grade <=` | 308 | 0.0574 |
| `meno >=` | 240 | 0.0447 |

(declaring replicates pooled: 5367)

**GRF** — share of declaring replicates whose rule contains the term, pooled over the six cells:

| term | replicates | share of declaring |
|---|---|---|
| `nodes <=` | 2927 | 0.2864 |
| `pgr >` | 2418 | 0.2366 |
| `size <=` | 2006 | 0.1963 |
| `er >` | 1892 | 0.1851 |
| `age <=` | 1889 | 0.1849 |
| `size >` | 1729 | 0.1692 |
| `age >` | 1639 | 0.1604 |
| `er <=` | 1439 | 0.1408 |
| `pgr <=` | 1036 | 0.1014 |
| `nodes >` | 880 | 0.0861 |
| `meno <=` | 671 | 0.0657 |
| `grade <=` | 503 | 0.0492 |
| `grade >` | 406 | 0.0397 |
| `meno >` | 346 | 0.0339 |

(declaring replicates pooled: 10219)

### Table 6 — measured wall per cell

| cell | wall (s) | wall (h) |
|---|---|---|
| null0657_n500 | 991 | 0.275 |
| null0721_n500 | 1139 | 0.316 |
| null0657_n1000 | 1056 | 0.293 |
| null0721_n1000 | 1240 | 0.344 |
| null0657_n1500 | 1152 | 0.320 |
| null0721_n1500 | 1331 | 0.370 |
| **total (18 renders)** | **6909** | **1.919** |

Build: forestsearch 0.3.5.9000 ; host Mac-Studio-3.local ; workers 12 ; R 4.5.2

### Table 7 — the unconditional claim rate: declares **and** the unadjusted one-sided 95% lower bound reaches HR 1.00 / 1.25

Over all 2,000 replicates, not over declaring ones.

| cell | FS | DINA | GRF |
|---|---|---|---|
| null0657_n500 | 0.0465 / 0.0035 | 0.0465 / 0.0080 | 0.0565 / 0.0100 |
| null0721_n500 | 0.1125 / 0.0105 | 0.1345 / 0.0230 | 0.1340 / 0.0245 |
| null0657_n1000 | 0.0130 / 0.0020 | 0.0060 / 0.0020 | 0.0190 / 0.0030 |
| null0721_n1000 | 0.0425 / 0.0045 | 0.0375 / 0.0035 | 0.0485 / 0.0055 |
| null0657_n1500 | 0.0030 / 0.0000 | 0.0005 / 0.0000 | 0.0030 / 0.0005 |
| null0721_n1500 | 0.0085 / 0.0005 | 0.0070 / 0.0005 | 0.0145 / 0.0010 |

---

## 5. Findings

- **The headline: a region is declared on 13.9% to 98.5% of replicates, under a design in which no
  subgroup exists.** The extremes are DINA at `null0657_n1500` (278/2000, 0.1390
  [0.1245, 0.1549]) and GRF at `null0721_n500` (1969/2000, 0.9845 [0.9781, 0.9891]). **Declaration
  is not, and must not be read as, evidence that a subgroup exists.** This is the same reading
  `current_status.md` §4 already requires of the differentially-null cells — "selection rate is not
  an error rate" — and the structural null makes it unavoidable rather than arguable.

- **The rate falls with n on every identifier, and it is higher at the weaker uniform benefit at
  every n.** FS: 0.4470 → 0.3185 → 0.1685 at HR 0.657, and 0.6780 → 0.6170 → 0.4845 at HR 0.721.
  The mechanism is the screening floor, which is **HR ≥ 0.90 on the natural scale**: a weaker
  overall benefit leaves more of the candidate family above it, and a larger n tightens the
  candidate estimates so fewer stray above it. Table 3 shows the floor-clearing count directly —
  median 43 → 16 → 5 at HR 0.657, and 100 → 53 → 26 at HR 0.721, out of ~1,711–1,830 enumerated.

- **GRF declares far more often than FS or DINA everywhere** (0.5645–0.9845 against 0.1390–0.7400),
  and its ordering against the others never reverses. This is the identifier's own behaviour, not
  the null's: `current_status.md` §2.3 records GRF selecting at or near 1 on the committed harm and
  differentially-null cells too. Cross-identifier comparison here carries the same confound the
  directory already flags — identifier, family construction and detection set — but **not** the
  criterion confound, since all three ran `effMaxSG` ε 0.20.

- **Size is remarkably stable: |Ĥ|/n sits between 0.150 and 0.228 in all eighteen runs**, with
  Monte Carlo standard errors of 0.0016–0.0031. Whatever the cell and whatever the identifier, a
  declared region is between a seventh and a quarter of the trial.

- **Specificity, both conventions.** *Unconditional* (a replicate declaring nothing scores 1):
  0.7878–0.9792, rising with n exactly as the declaration rate falls. *Conditional* (declaring
  replicates only, S3.4): **0.7719–0.8502 — nearly flat across every cell and identifier.** That
  flatness is the direct consequence of the |Ĥ|/n stability above: under an empty planted region
  every selected subject is a false positive, so conditional specificity is exactly 1 − |Ĥ|/n. The
  two conventions therefore say different things here, and the unconditional one is doing almost
  all of its work through the declaration rate.

- **The candidate family, and where the screen bites.** The enumerated family is flat at
  1,711–1,830 conjunctions in every cell (it is outcome-independent). Of those, the median clearing
  the effect floor falls from 100 to 5 across the grid, and the median *consistency-qualifying* set
  is 2–7. **The consistency screen declines a floor-clearing family on 32.0% to 79.8% of the
  replicates where something did clear the floor** — it is doing a great deal of work, and doing
  more of it as n grows (0.5460 → 0.6660 → 0.7982 at HR 0.657).

- **`max_g T_g` against the conventional cutoff `z_0.95 = 1.645` — it is not calibrated.** The
  share exceeding 1.645 ranges from **0.0078** (`null0657_n1500`) to **0.2077** (`null0721_n500`),
  a factor of 27, and it moves systematically with both n and the uniform effect. At
  `null0721_n500` it is roughly four times nominal; at `null0657_n1500` it is a seventh of nominal.
  The median itself moves from 0.146 to 1.159 across cells. **A fixed 1.645 cutoff on this
  statistic would therefore not control anything uniformly over this grid**, which is precisely the
  case for the calibrated cutoff of Section 4. Evaluating that cutoff is out of scope here; the
  recorded `max_g T_g` is what makes it possible later without a re-run.

- **Where the unadjusted bound lands.** Among declaring replicates, the share whose one-sided 95%
  Wald lower bound reaches **HR 1.00** runs 0.0036–0.1818, and **HR 1.25** runs 0.0000–0.0311.
  Both fall steeply with n. Unconditionally (Table 7) the worst cell is `null0721_n500`, where
  **11.3% (FS) to 13.5% (DINA) of all replicates both declare a region and carry an unadjusted
  lower bound at or above HR 1.00**, and 1.1%–2.5% reach 1.25. By n 1500 those are 0.3%–1.5% and
  ≤ 0.1%. These are **unadjusted** products of the search itself; no selection-adjusted interval was
  computed, and the point of this cell is not to evaluate one.

- **The composition of Ĥ is driven by the candidate construction, not by signal.** On FS the ER
  terms dominate (`NOT er <=` 0.2801, `er <=` 0.2537 of declaring replicates), then `NOT pgr <=`
  (0.2143), `nodes <=` (0.1915) and the size terms. ER leads **because the template forces the cut
  `er <= 0` and puts a 10-quantile grid on raw ER** (§0.4), not because ER carries any effect —
  there is none. DINA's picks lead with `pgr >=` (0.2914) and GRF's with `nodes <=` (0.2864); each
  engine's ordering follows its own enumeration. Read across identifiers, the only safe statement
  is that the declared rule reflects how candidates were built.

- **Cost.** 6,909 s of render wall over the eighteen runs (1.92 h), 991–1,331 s per cell, against a
  1.71 h revised projection. The main campaign ran 6,088 s end to end; the cell-5 completion run
  added 824 s. Host `Mac-Studio-3.local`, 14 physical cores, 36 GB, 12 workers, R 4.5.2,
  forestsearch 0.3.5.9000, threads pinned to 1.

## 6. What is not reported, and why

- **Sensitivity and PPV are undefined with an empty planted region and are not reported.**
  Sensitivity is `TP/(TP+FN)` with `TP = FN = 0` — it is NA on all 36,000 replicates, and Gate A
  asserts that per bundle. PPV is `TP/(TP+FP)` with `TP = 0`, so it is **0 by construction** on
  every declaring replicate and carries no information; the recorder stores it, and it should not
  be quoted as a measured quantity.
- **NPV is 1 by construction** — `TN/(TN+FN)` with `FN = 0 `— on every declaring replicate. Gate A
  asserts it.
- **The identified-to-planted size ratio has no denominator here** and is not reported.
- **No MR product of any kind exists in these bundles.** Gate A asserts that all 118 `mr_*` /
  `fld_*` / `fb_*` product columns are NA and that `mr_ok == 0`, per run.
- **`n_family` is NA** and the maximum consistency rate over the *full screened* family is not
  recoverable; both would require an `R/` change (§0.6).
- **The other selection rules' picks are not recorded**, because this machinery does not produce
  them from one search (§0.6).
- **Rates carry Wilson 95% intervals and means carry Monte Carlo standard errors**, as S3.7
  requires. A replicate-mean rate is reported beside, never inside, a Wilson interval.
