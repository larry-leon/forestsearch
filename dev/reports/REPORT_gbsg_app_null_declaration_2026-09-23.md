# REPORT — Realized declaration rate under a uniform marginal Cox HR of 0.75, GBSG application family

Task: `dev/tasks/TASK_gbsg_app_null_declaration_2026-09-23.md` (commit `d3cbcfae`). One cell: n = 686, 1,000 replicates, FS
only, no MR. Driver `quarto/simulations/gbsg_app_null/run_gbsg_app_null.R` (`e6477a22`); payload
`quarto/simulations/gbsg_app_null/results/gbsg_app_null_full.rds` / `.csv`, timing run `gbsg_app_null_timing.*`, logs in
`quarto/simulations/gbsg_app_null/logs/` (`5774bacf`).

## Result

- **Realized declaration rate: 329 / 1,000 = 0.329, Wilson 95% interval [0.301, 0.359].** Beside it,
  FŴ₀.₁₀(0.75) = 0.651 (read from the application payload, `declaration_calibration$table`, row c0 0.75). The gap is
  **−0.322**: the executed p\* = 0.90 screen declared in about half as many trials as the boundary quantity predicts.
- The gap is below 0.651, as expected. Its size comes with two qualifications, and they pull in opposite directions (see
  the per-candidate null and the allocation items below). Neither changes the order of magnitude.

## Pins

- HEAD at the run: `e6477a22` (driver commit), on `feature/glm-extension`. Task commit `d3cbcfae`.
- forestsearch 0.3.5.9000, installed build 2026-09-23 02:32:08 UTC. The last `R/` commit (`845fea56`, 2026-09-23 02:15 UTC)
  predates the build, so the installed package is HEAD's `R/`. The workers load the installed package, not `load_all()`.
- R 4.6.1 (2026-06-24), x86_64-pc-linux-gnu, 128 logical cores; run on 14 multisession workers.
- Seeds:
  - DGM `seed` 8316951 (the document's seed);
  - replicate b: `simulate_from_dgm` seed = `seedit` = 8316951 + b, b = 1…1000;
  - candidate-truth error draws: seed 20260923.
- Workers run under `.options.future = list(seed = TRUE)`, as the programme template does, so the replicate data are
  drawn under L'Ecuyer-CMRG. The DGM is built in the main process before any switch. Replicate 3 re-run sequentially under
  `RNGkind("L'Ecuyer-CMRG")` reproduces the parallel result exactly (same rule, N = 105, same censoring). Under
  Mersenne-Twister it draws a different trial.
- The timing run's 10 rows are identical to rows 1–10 of the full run.
- Pre-existing untracked files at Step 0: the two `actg175/binary_020/mr_or_harm/…_d5000/` directories,
  `smoke_redes.html`, `smoke_relaunch.html`, and `gbsg_020/scripts_dinamr/logs/nullmr_findings.err`. None were staged.

## The application's settings (Gate 1: PASS)

Read from `fs-glms-interpretable/quarto/gbsg/analysis_gbsg_mr.qmd`, read-only.

- Frame, L151–166: `survival::gbsg`, **N = 686** (asserted), 246 hormonal / 440 none, 299 events.
  - `time_months = rfstime / 30.4375`, `grade3 = ifelse(grade == "3", 1, 0)`, `id = row number`.
- Outcome `time_months`; event `status`; treatment `hormon`; id `id`. Seed `seedit = 8316951`.
- `confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er")`.
- Cuts: `cut_type = "default"`, `cont.cutoff = 4`, `conf.cont_jcuts = list(er = 10, pgr = 10)`,
  `conf_force = c("er <= 0", "pgr <= 0")`, `collapse_cuts = TRUE`, `max_n_confounders = 1000`.
- `maxk = 2`, `n.min = 60`, `d0.min = 10`, `d1.min = 10`, `fs.splits = 1000`.
- `hr.threshold` (c1) = 1.00, `hr.consistency` (c2) = 1.00, `pconsistency.threshold` (p\*) = 0.90.
- `sg_focus = "effMaxSG"` (normalized by `forestsearch()` to its alias `hrMaxSG`), `selection_rule = "neighborhood"`,
  `effect_neighborhood = 0.20`.
- `consistency_method = "resample"`, `is.RCT = TRUE`, `est.scale = "hr"`.
- `m1.threshold = Inf`, `minp = 0.025`, `stop_threshold = NULL`.
- `use_lasso = use_grf = use_dina = FALSE`, `subgroup_method = "consistency"`.
- The call does not set `use_twostage`, `outcome_type` or `effect_measure`, so the application ran at `forestsearch()`'s
  formal defaults: `use_twostage = TRUE`, `outcome_type = "survival"`, `effect_measure = NULL`.
  - The driver passes the first two explicitly, because `default_sim_params()` would otherwise substitute
    `use_twostage = FALSE`. `effect_measure` stays NULL, as in the application.
  - It also overrides that function's `use_lasso = TRUE`, `d0/d1.min = 12`, `fs.splits = 400` and `hr.threshold = 1.25`.
- **Verified at the receiving end.** Replicate 3's `args_call_all` shows every setting above as received, including
  `stop_threshold = NULL`.

## Calibration (Gate 2: PASS)

- `generate_aft_dgm_flex(model = "null", n_super = 5000, seed = 8316951)` on the application frame.
  - Continuous: age, size, nodes, pgr, er. Factors: meno, grade3.
  - `calibrate_k_treat(0.75, use_ahr = FALSE, tol_rel = 1)` gives **`k_treat` = 1.048469**.
- **Achieved overall marginal Cox HR = 0.750000** (within 1%: yes). `flag_harm` is identically 0 on all 5,000 rows of
  `df_super`.
- Other scales, for vocabulary only: the patient-level HR is uniform at exp(−0.3916) = **0.676**, and under `null` the
  AHR = CDE = 0.676. The calibration target is the marginal Cox HR, 0.75.

## The per-candidate null (finding)

- **Route.** `fs_oc_family_enumerate()` accepts `glm_dgm` objects only, so the driver reproduces its sections 2–4 on
  `df_super` with the package internals:
  - `get_FSdata()` cuts at population quantiles;
  - `dummy()` for both directions of each cut;
  - all combinations of up to two factors;
  - the empty / minp / rmin / size floors as population proportions at n = 686;
  - identical memberships collapsed.
  - The event floors `d0.min` / `d1.min` count sample events and have no population analogue, so they are not applied.
- **Truth.** Not `fs_betaHhat_table()`. Each candidate's marginal Cox HR is computed on the estimand the calibration
  itself uses (`calculate_hazard_ratios()`):
  - both potential outcomes per subject under a common extreme-value error, uncensored;
  - the error drawn 20 times per subject and stacked, to cut Monte Carlo noise.
  - The same construction gives an overall HR of 0.7476 against the calibrated 0.7500: noise of about 0.3%.
- **Family: 1,504 candidates.**
  - Counts: 66 cut columns, 2,211 enumerated; dropped 142 empty, 0 minp, 122 rmin, 329 size; 1,618 kept; 114 duplicates.
  - The application's own pre-reduction family on the observed sample is 1,744 (payload
    `family_size_prereduction`). The difference is sample- vs population-quantile cuts and sample floors.
- **β_true(g) on the HR scale: min 0.680, median 0.735, max 0.799.**
- **410 of 1,504 candidates (27%) have a true marginal HR above 0.75.** The highest are small high-receptor or older-age
  cells: `pgr > 340` at 0.799 (9% prevalence), `er > 293` at 0.796, `er > 174 & age > 61` at 0.788.
- **The maximum sits above 0.75.** Calibrating the overall marginal HR to 0.75 does not put every candidate inside FŴ's
  null set {β(g) ≤ log 0.75 for all g}. The spread comes from non-collapsibility of the Cox HR under a uniform AFT effect.
- **Direction of this effect.** Candidates nearer HR 1 are easier to declare. This design therefore declares at least as
  often as one with every candidate at or below 0.75, so on this account the realized rate overstates the boundary design,
  if anything.

## Event and censoring rates

- Simulated trials: mean event rate 0.426 (SD 0.019, range 0.364–0.483), censoring 0.574.
- GBSG: event rate 0.436 (299 / 686), censoring 0.564.
- They agree to about one point. The DGM's own fitted censoring model was used, with `analysis_time = Inf` (the
  `run_simulation_analysis()` default) and no administrative cut.

## Allocation (finding — material to the comparison)

- **The simulated trials are randomized exactly 1:1 (343 / 343).** `simulate_from_dgm()` defaults to `rand_ratio = 1`,
  and `run_simulation_analysis()` exposes no `rand_ratio`. GBSG is 246 / 440 (36% treated).
- At the same total events, a 36:64 split inflates each candidate's log-HR variance by about 1/(4 · 0.36 · 0.64) ≈ 1.09
  relative to 1:1, so SEs are about 4% larger. Wider SEs push more benefit-centred candidates across c1 = c2 = 1.00, so
  **the application-matched allocation would declare somewhat more often than 1:1.**
- This effect runs opposite to the per-candidate one. Its size is not measured here; see OPEN ITEMS.

## Among declaring replicates (329)

- **Declared subgroup size:** median 95; IQR 75–120; 10th / 90th percentiles 66 / 154; range 61–265. As a share of the
  trial, median 13.8% (10th / 90th percentiles 9.6% / 22.5%). The selections sit near the `n.min = 60` floor.
- **Rules:** 326 two-factor, 3 one-factor. There are **318 distinct rules** in 329 declarations, so no rule recurs more
  than 3 times. The most frequent is `{size <= 20} & !{meno}` (tumor ≤ 20 mm and premenopausal, 3 times); nine rules
  appear twice, and 308 appear once.
- **Most frequent factors in declared rules:**
  - `{grade3}` 28, `{size <= 20}` 24, `{age <= 46}` 19;
  - `!{meno}`, `!{nodes <= 3}` and `{age <= 53}` 16 each.
  - The receptor cuts are rare, although the per-candidate truth is highest in the high-receptor cells.
  - The declarations are dispersed across the family. That is the picture of multiplicity, not of a recurring real
    signal.
- Candidates surviving the c1 screen per replicate (`n_candidates_total`): median 63, IQR 28–126.

## Timing (Gate 3: PASS)

- **Timing run:** 10 replicates on 14 workers, 6.8 s wall-clock. Projected 1,000-replicate wall-clock **0.19 h**,
  against the 2 h gate and the 15–25 min estimate.
- **Full run:** 327 s (5.5 min) for the replicate loop. 08:29:37Z–08:37:06Z end to end, including the DGM calibration
  and the ~2 min family enumeration and truth.
- Per replicate: mean 4.3 s, range 2.9–10.6 s.
- **Errors: 0 of 1,000.** No replicate errored, and none raised the harness's `FS failed:` warning. The driver captures
  that warning explicitly, because the harness otherwise turns a failed fit into `any.H = 0` silently. No warnings of any
  kind were recorded. The denominator is 1,000.

## Deviations from the task text

- **§4 call path.** The driver calls `run_simulation_analysis(methods = list(FS = list()), fs_params = …)`, not the
  `run_fs = TRUE` path.
  - The `run_fs` path merges `fs_params` with `utils::modifyList()`, which drops `stop_threshold = NULL`.
    `forestsearch()` would then take its formal default (`pconsistency.threshold` = 0.90) and early-stop under `effMaxSG`,
    a setting the application did not use.
  - The `methods` path merges with `.modify_keep_null()`. The fit is otherwise identical.
- **§4 simulation.** The separate `simulate_from_dgm()` call in the task's snippet was not made.
  `run_simulation_analysis()` simulates internally with `seed = seed_base + sim_id`, `analysis_time = Inf` and
  `cens_adjust = 0`: the same trial the snippet would draw, apart from its `analysis_time` default of 48. The event and
  censoring rates reported above are those of the analysed trials (`p.cens`).
- **§4 seeds.** `seedit` is `8316951 + b` per replicate, the programme convention (template L1264), rather than the fixed
  document seed. Otherwise every replicate would share one split stream.
- **§6 catalogue.** `current_status.md` is generated from `status_curated.md` by `scripts_dinamr/current_status_regen.R`,
  and the directory's rule is to commit it alone as a child commit. The line therefore goes into `status_curated.md` (with
  the regen's hardcoded "Updated" line), followed by a separate regeneration commit: six commits instead of five.

## OPEN ITEMS

- **Allocation.** A second cell at `rand_ratio = 246/440` would put the simulated trials on the application's allocation.
  It needs either a `rand_ratio` pass-through in `run_simulation_analysis()` (an `R/` change) or a driver that calls
  `simulate_from_dgm()` + `forestsearch()` directly. It is out of scope here (one cell).
- **Per-candidate null.** 27% of the family sits above HR 0.75 (maximum 0.799). A cleaner comparison with FŴ would
  calibrate so that the **maximum** candidate HR, not the overall HR, is 0.75, which is a stricter interior null. With the
  allocation effect, this brackets the realized rate from both sides.
- **Harness side issue, not fixed.** `run_simulation_analysis()`'s legacy `run_fs` / `run_fs_grf` paths merge `fs_params`
  with `modifyList()`, silently dropping any explicit `NULL` (for example `stop_threshold = NULL`, `n.min = NULL`). The
  `methods` path already uses `.modify_keep_null()`. It is a one-line change per path if wanted.
- **Family comparison.** The population family (1,504) and the application's sample family (1,744) differ in
  construction. A candidate-for-candidate map between them was not attempted.
