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
