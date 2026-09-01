# REPORT — the assembly skip (0.3.5), plus two read-only settlements: bootstrap outer parallelism and `stop_threshold`

**Task:** `dev/tasks/cc_task_assembly_skip_2026-09-02.md` · **Executed:** 2026-09-02 (session date 2026-09-01) by CC
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`

**Verdict up front: all gates pass, both settlements land.** The unadjusted per-candidate Cox fit runs on subset vectors, bit-identically against the committed 0.3.4 baseline (37/37 `identical()`, bootstrap payload included). The re-bucketed prize was **1.11 s** (against the 09-01-derived ≈1.03 s) and **1.01 s was realized** in the assembly bucket (F2 wall 4.82 → 3.58 s profiled; bootstrap replicate 5.49 → 4.18 s, B = 1000: 91.5 → 69.7 min sequential). Settlement A: **bootstrap outer parallelism is reproducible** — sequential and multisession/20 payloads are `identical()`; B = 1000 implied ≈15 min at 20 workers for the continuous configuration. Settlement B: **first-passer `stop_threshold` is sound for `maxeffCons`** — same winner as full evaluation, identical selected row, ordering established before the early-stop check.

---

## 1. Provenance and §2 verification

**Machine/repo:** pop-os, `~/Documents/GitHub/forestsearch`, branch `feature/glm-extension`, HEAD at task start `ef021f50`, porcelain clean. `packageVersion("forestsearch")` **0.3.4** ✓; R 4.6.1, 128 cores. *GATE:* the four medians-task commits `ff2f9ca2`, `4c923dfd`, `3921ffdd`, `ef021f50` in the log ✓. `vi.grf.min` default **`NULL`** ✓. Task doc committed alone: `dev/tasks/cc_task_assembly_skip_2026-09-02.md`, commit `4b826c98`.

### §2.1 `fit_cox_for_subgroup()` at HEAD (`ef021f50`, pre-edit numbering `R/subgroup_search.R:753–835`)

The assembly lines (L758–772):

```r
  # Create subgroup data
  data.x <- data.table::data.table(Y = yy, E = dd, Treat = tt, id.x = id.x)   # L759

  # Attach raw adjustment columns (row-aligned with the vectors) when present
  rhs <- "Treat"
  if (length(adj) > 0L && !is.null(df_clean)) {                               # L763
    adj_vars <- intersect(.fs_adjust_vars(adjust_covariates), names(df_clean))
    if (length(adj_vars) > 0L) {
      adj_cols <- as.data.frame(df_clean)[, adj_vars, drop = FALSE]
      data.x <- cbind(data.x, adj_cols)                                       # L767
    }
    rhs <- paste(c("Treat", adj), collapse = " + ")
  }

  df.x <- data.x[id.x == 1]                                                   # L772
```

Constructor: `data.table::data.table()`; columns `Y`, `E`, `Treat`, `id.x` bound from the passed vectors, plus (adjusted only) a `cbind` of raw covariate columns; subset `data.x[id.x == 1]`. The 0.3.4 state confirmed: unadjusted fast branch with direct CI arithmetic (L794–805), adjusted `summary.coxph` branch (L806–811), medians gone (NA return, L828–834). **Remaining consumers of `df.x`: exactly the two `coxph()` calls** — `survival::coxph(cox_fmla, data = df.x, robust = FALSE)` at L795 (unadjusted) and inside `summary(...)$conf.int` at L808 (adjusted). Nothing else reads `df.x` or `data.x`. *GATE ✓*

### §2.2 Assembly semantics

Pure column binding of the already-subset-aligned vectors plus a row subset (`data.x[id.x == 1]`, order-preserving): no NA handling, no sort, no filtering, no value-changing coercion. *GATE ✓* — the substitution's value-and-order identity premise holds.

### §2.3 Upstream NA guarantee

`subgroup.search()` cleans before the candidate loop (`R/subgroup_search.R:101–105` ← `prepare_search_data()` L322–332):

```r
prepare_search_data <- function(Y, Event, Treat, Z) {
  temp <- cbind(Y, Event, Treat, Z)
  temp <- na.exclude(temp)
  ...
```

`yy`, `dd`, `tt` are jointly row-wise NA-free by construction, so the environment-resolved `coxph()` (default `na.action`) cannot diverge from the frame-based fit. *GATE ✓*

### §2.4 The prize, re-sized from the committed record before any edit

Re-bucketing `medians_profile_F2_after_2026-09-02.Rprof` (0.3.4, committed at `ef021f50`) with the established `BUCKETS` + `FIT_SPLIT` (`assembly_rebucket_034_2026-09-02.log`): 481 samples × 0.010 s; fit bucket 3.69 s / 76.7 %; sub-split — solve 0.22 s (6.0 %), model-frame/matrix 0.56 s (15.2 %), vcov-summary 0, formula 0.06 s, **subset-assembly 1.11 s (30.1 % of fit)**, wrapper/other 1.74 s (47.2 %). **Predicted recovery: 1.11 s**, quoted beside the 09-01-derived ≈1.03 s — the two agree; no `summary.coxph`-style misattribution here, and the prize is far above the ≈0.3 s floor.

### §2.5 Transplants — authored nothing

`assembly_battery_2026-09-02.R` = copy of committed `medians_baseline_2026-09-02.R`; `diff` is exactly four lines: the env-var name (`FS_MEDIANS_OUT` → `FS_ASSEMBLY_OUT`, one code line + its header-comment mention) and the two output-filename lines (`medians_%s` → `assembly_%s` in `saveRDS` and the final `cat`). `assembly_profile_2026-09-02.R` = copy of committed `medians_profile_2026-09-02.R`; `diff` is exactly four lines: env-var name (code + header comment) and the two output-name lines (`.Rprof`, `.rds`). *GATE ✓* — committed as `934f31ed`. `assembly_compare_2026-09-02.R` = copy of committed `medians_compare_2026-09-02.R` with only the two `readRDS` input lines changed (baseline → `medians_postchange_2026-09-02.rds`, new → `assembly_postchange_2026-09-02.rds`).

**Baseline provenance (not re-run):** `medians_postchange_2026-09-02.rds`, committed at `ef021f50` — the 0.3.4 output of this identical fixture battery (script `medians_baseline_2026-09-02.R` at `4c923dfd`, run on installed 0.3.4 at HEAD `3921ffdd`), with within-stage self-consistency 3/3 (F2, E-ties, bootstrap) as recorded in `REPORT_medians_on_survivors_2026-09-02.md`.

---

## 2. The restructure as implemented (commit `f3975b99`)

Files: `R/subgroup_search.R`, `DESCRIPTION` (→ **0.3.5**), `NEWS.md`. `devtools::document()` clean (no NAMESPACE/man changes; roxygen never mentioned the frame). Post-edit anchors:

- **`idx <- id.x == 1`** computed once at the top (`R/subgroup_search.R:757`).
- **Unadjusted arm** (L776–796): three locals `Y <- yy[idx]; E <- dd[idx]; Treat <- tt[idx]` (L782–784), then the 0.3.3 fast block byte-for-byte except the fit line, now `survival::coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)` (L786) — literal formula, **no `data =`**, the three locals resolved from the fitting environment; coefficient named `Treat` as before; CI arithmetic, `try()` contract, and dimnames lines untouched.
- **Adjusted arm** (L798–821): the assembly lines relocated verbatim — `data.x` constructor (L802), the covariate `cbind` block (guard reduced from `length(adj) > 0L && !is.null(df_clean)` to the residual `!is.null(df_clean)`, since `length(adj) > 0L` is the branch condition), `df.x <- data.x[id.x == 1]` (L815), `cox_fmla` construction (L817), and the unchanged `summary(coxph(..., data = df.x))$conf.int` fit.
- Everything downstream — `trow` extraction, NA-medians return — unchanged; the medians site (`med_df`, L916 area) untouched, as gated.

---

## 3. Equality gates against the committed baseline

Battery run once at 0.3.5 (`FS_ASSEMBLY_OUT=postchange`) → `assembly_postchange_2026-09-02.rds`; compared via `assembly_compare_2026-09-02.R` under the battery's own (unchanged) volatile-field exclusion (`time_search`, `minutes_all`; `tmins_*` columns; data.table normalisation; 26 dropped paths of the recorded kinds).

**Equality matrix (`assembly_compare_2026-09-02.log`): 37/37 `identical() == TRUE`** — for each of F2, F5 (E-1: F5 exercises the relocated adjusted-arm assembly), F1, F3 (E-2: continuous untouched), Eties, Enamed, Efinite, Ezero, Edegen (E-3: E-degen re-proves the degenerate-fit warning/`try()` behaviour on the environment-resolved call — 2 coxph warnings at parity): `pruned`, `warnings`, `counters`, `sg.harm`; plus E-4, the 20-replicate gbsg bootstrap payload. **ALL EQUALITY GATES PASS.** Within-stage self-consistency at 0.3.5: F2, E-ties, bootstrap — 3/3 `TRUE`, as at both prior stages. `survfit` counts unchanged from 0.3.4 (F2 121, F5 111, E-zero 0, F1/F3 0).

---

## 4. Substitution proof and realized recovery

**Constructor counts** (`assembly_counts_2026-09-02.{R,log,rds}`; `trace()` on `data.table()` in the data.table namespace and `data.frame()` in base, innermost forestsearch caller attributed, one call each at 0.3.5):

| | `data.table()` from `fit_cox_for_subgroup` | medians-site `data.frame()` (`evaluate_combination_with_status`) |
|---|---|---|
| F2 unadjusted, 0.3.4 | 1,410 (structural: one construction per fitted candidate — cited, not re-installed) | 121 |
| F2 unadjusted, **0.3.5** | **0** (the single remaining `data.table()` is `format_search_results`'s result assembly, 1 call) | 121 — the separate constructor, counted and attributed |
| F5 adjusted, 0.3.5 | **1,410** — adjusted fits still construct | 111 |

**Profile** (transplanted profiler, F2, `Rprof(interval = 0.010)`; before = §2.4's re-bucketing of the committed 0.3.4 Rprof; after = `assembly_rebucket_035_2026-09-02.log`):

| | 0.3.4 | 0.3.5 |
|---|---|---|
| profiled wall | 4.82 s | **3.58 s** |
| fit bucket | 3.69 s / 76.7 % | 2.29 s / 64.3 % |
| — subset-assembly | **1.11 s** | **0.10 s** |
| — solve / model-frame / formula / wrapper | 0.22 / 0.56 / 0.06 / 1.74 s | 0.26 / 0.60 / 0.07 / 1.26 s |
| medians bucket | 0.22 s | 0.23 s (stable, as expected) |

Battery walls: F2 5.38 → 3.93 s, E-finite 4.91 → 3.54 s; F5 38.14 → 38.35 s (adjusted arm unchanged by design); F1 1.98 → 1.98 s, F3 1.65 → 1.87 s (continuous untouched; jitter). **Bootstrap:** 109.9 → 83.6 s / 20 reps, **5.49 → 4.18 s per replicate**; survival **B = 1000 projection 91.5 → 69.7 min** (sequential). **Realized against both predictions:** assembly-bucket recovery **1.01 s** vs §2.4's 1.11 s and the 09-01-derived ≈1.03 s — realized essentially in full (total profiled-wall Δ 1.24 s; the wrapper bucket also shed ~0.5 s, consistent with the removed `model.frame`-on-frame pathway shifting work). Residual gbsg cost: coxph wrapper/internals + model-frame (2.29 s fit bucket, 55 % wrapper) — the recorded deeper-Cox item, untouched.

---

## 5. Settlement A — bootstrap outer parallelism (read-only)

**Source mechanics.** Replicates are dispatched by `foreach::foreach(boot = seq_len(nb_boots), .options.future = list(seed = TRUE, globals = TRUE), .errorhandling = "pass") %dofuture% { ... }` (`R/bootstrap_analysis_dofuture.R:408–415`) — `seed = TRUE` is doFuture's parallel-safe per-element L'Ecuyer-CMRG stream pre-assignment. Bootstrap row indices do not use worker RNG at all: the B × N matrix is pre-generated **on the main process** — `boot_seed <- seed; boot_index_mat <- .make_boot_index_matrix(n, nb_boots, seed = boot_seed)` (`R/bootstrap_dofuture_main.R:382–387`, generator `set.seed(seed); sample.int(...)` at `R/bootstrap_calculations_helpers.R:26–34`) — and consumed per replicate as `in_boot <- boot_index_mat[boot, ]` (`bootstrap_analysis_dofuture.R:433–434`; the per-iteration `sample.int()` is a legacy fallback for direct callers only). Each replicate's `forestsearch()` re-run inherits the caller's `seedit` unchanged (`args_FS_boot <- args_FS_template`) and is forced to an inner sequential plan (`args_FS_boot$parallel_args$plan <- "sequential"`), so in-replicate RNG is re-seeded deterministically per replicate regardless of backend. The caller's outer plan enters at `setup_parallel_SGcons(parallel_args)` (`bootstrap_dofuture_main.R:366`) → `future::plan(multisession/multicore/callr/sequential, workers = ...)` (`R/subgroup_consistency_helpers.R:32–77`). One stale doc line noted for the record: `bootstrap_analysis_dofuture.R:200` claims per-iteration seeds "8316951 + boot * 100" — no such code exists.

**Measurement** (`assembly_settleA_2026-09-02.{R,log,rds}`; F1 continuous bootstrap configuration, `nb_boots = 40`, seed 8316952, at 0.3.5): sequential wall **91.5 s** (2.29 s/rep); `plan = "multisession", workers = 20` wall **35.7 s**; pruned payloads **`identical()` — TRUE**.

**Verdict (i): outer parallelism is reproducible.** Speedup 2.6× at nb = 40 (two waves of 20 with worker startup and serialization amortized poorly at this size); **implied B = 1000 at 20 workers ≈ 14.9 min** against 38.1 min sequential for this configuration — and it reframes every sequential projection in the record (e.g. the survival configuration's 69.7 min sequential at 0.3.5). No code change, no plan-default change; whether production runs adopt an outer parallel plan is Larry's call.

---

## 6. Settlement B — first-passer `stop_threshold` under `maxeffCons` (read-only)

**Measurement** (`assembly_settleB_2026-09-02.{R,log,rds}`; F3 anchor at 0.3.5, seed base as in the battery): as-is (default `stop_threshold = pconsistency.threshold = 0.90`, the `forestsearch()` formal at `R/forestsearch_main.R:1235`): **1 of 131** candidates evaluated, `early_stop_triggered = TRUE` at candidate 1, selection `{age <= 37} & !{cd40 <= 507}`. With `stop_threshold = NULL`: **131 of 131** evaluated, no early stop, **same selection**, and the selected row is **`identical()`** across the two runs — `Pcons 0.95, hr 87.91666667, N 66, E 66, g 1482, m 1, K 2`. The full run additionally carries 7 more consistency-passing rows (Pcons 0.90–0.96, hr 61.5–82.3, all below the winner's 87.9) that the early-stopped run never evaluates — extra survey, no bearing on selection.

**Ordering, from source.** The candidate pool is sorted **before** the evaluation loop: `found.hrs <- sort_subgroups_preview(found.hrs, sg_focus, ...)` at `R/subgroup_consistency_main.R:607`, where for `sg_focus == "maxeffCons"` the preview sort is `data.table::setorder(result_new, -HR, K)` (`R/subgroup_consistency_helpers.R:659–661`), i.e. effect-descending — the same key `sort_subgroups()` applies at selection (`(-hr, K)`). The early-stop check runs inside the Section-8 loop after each evaluation (`subgroup_consistency_main.R:829–832`: `pcons_m >= stop_threshold`). The in-source design note (`R/forestsearch_main.R:1586–1608`) states exactly this soundness condition — the selection key contains no Pcons term, and enumeration order equals the selection key, ties included, "so early stopping is sound only for … maxeffCons alone"; the guard at L1609+ blocks `stop_threshold` for every other consistency focus.

**Verdict: soundness confirmed** — matching winner, identical selected row, ordering established before the stop check, and the code's own soundness argument quoted. The bootstrap-level early-stop question stays with Larry; nothing here enables it.

---

## 7. Package health gates

- `devtools::test()` on 0.3.5: **0 failures | WARN 31 | 3 skipped | 4,864 passed** — parity with the recorded baselines; `tests/` untouched, `test-search-reproducibility.R` unmodified and passing. *GATE ✓* (`assembly_test_035_2026-09-02.log`)
- `devtools::check(document = FALSE, vignettes = FALSE, manual = FALSE)`: **0 errors | 0 warnings | 0 notes**, matching the recorded 0.3.4 result. *GATE ✓* (`assembly_check_035_2026-09-02.log`)
- Standing guards on 0.3.5 source: `fidelity_fs_oc_predict_2026-08-28.R` — **FIDELITY GATE: PASS (bit-identical)**; `prerefactor_reference_2026-08-29.R check` — **REFACTOR GUARD: PASS**.

---

## 8. Ten-line verdict

1. The unadjusted per-candidate Cox fit now runs on subset vectors — no frame — with the assembly relocated verbatim into the adjusted arm; `R/subgroup_search.R`, `DESCRIPTION` (0.3.5), `NEWS.md` only.
2. Every equality gate passes against the committed 0.3.4 baseline: 37/37 `identical()` across nine fixtures × four components plus the gbsg bootstrap payload; self-consistency 3/3.
3. The §2.4 re-bucketing put the prize at 1.11 s (09-01 figure ≈1.03 s — no misattribution this time); **1.01 s realized** in the assembly bucket, F2 profiled wall 4.82 → 3.58 s.
4. Constructor proof: zero `fit_cox_for_subgroup` frame constructions on the unadjusted path at 0.3.5 (1,410 → 0); F5 adjusted still constructs 1,410; the 121 medians-site `data.frame()` calls are separate and unchanged.
5. Bootstrap replicate 5.49 → 4.18 s; survival B = 1000 sequential projection 91.5 → 69.7 min.
6. Settlement A: replicate dispatch is `%dofuture%` with pre-assigned L'Ecuyer streams, main-process index matrix, and per-replicate deterministic `seedit` — and measured sequential vs multisession/20 payloads are `identical()`; B = 1000 ≈ 15 min at 20 workers (continuous config). Reproducible: verdict (i).
7. Settlement B: first-passer `stop_threshold` under `maxeffCons` selects identically to full evaluation (1 of 131 evaluated vs 131), selected row `identical()`, the (-HR, K) preview sort precedes the stop check — soundness confirmed; bootstrap-level use stays Larry's decision.
8. Tests 0 fail / WARN 31 / 4,864 pass; check 0/0/0; both standing guards bit-identical.
9. Residual gbsg cost is now coxph wrapper/internals + model-frame (recorded deeper-Cox item); stale doc line noted at `bootstrap_analysis_dofuture.R:200`.
10. Deviations: none of substance — the settlement-B driver needed two in-session script fixes before its committed form (a `modifyList()` NULL-deletion pitfall that silently kept the default `stop_threshold`, and the consistency-table location), both corrected and re-run; the committed script is the corrected one.

## Commits

1. `4b826c98` — task doc alone (`dev/tasks/cc_task_assembly_skip_2026-09-02.md`).
2. `934f31ed` — the two transplants: `assembly_battery_2026-09-02.R`, `assembly_profile_2026-09-02.R`.
3. `f3975b99` — the implementation: `R/subgroup_search.R`, `DESCRIPTION` → 0.3.5, `NEWS.md`.
4. (this commit) — artefacts + report: `assembly_postchange_2026-09-02.{rds,log}`, `assembly_compare_2026-09-02.{R,log}`, `assembly_rebucket_034_2026-09-02.log`, `assembly_rebucket_035_2026-09-02.log`, `assembly_profile_F2_after_2026-09-02.{Rprof,rds}`, `assembly_profile_after_2026-09-02.log`, `assembly_counts_2026-09-02.{R,log,rds}`, `assembly_settleA_2026-09-02.{R,log,rds}`, `assembly_settleB_2026-09-02.{R,log,rds}`, `assembly_check_035_2026-09-02.log`, `assembly_test_035_2026-09-02.log`, this report.

No push. Tree clean at close.
