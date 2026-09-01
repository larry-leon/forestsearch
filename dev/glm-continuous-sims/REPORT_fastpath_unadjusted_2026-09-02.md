# REPORT — unadjusted fast path for the per-candidate fit, under a bit-identity equality guard

**Task:** `dev/tasks/cc_task_fastpath_unadjusted_2026-09-02.md` (commit `c79f797d`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Category:** edits `R/` — adds routed code paths only; nothing moved, nothing deleted, no method change. Files changed: `R/glm_effect_estimators.R`, `R/subgroup_search.R`, `R/forestsearch_main.R`, `DESCRIPTION` (0.3.2 → 0.3.3), `NEWS.md`. Artefacts under `dev/glm-continuous-sims/`: `fastpath_baseline_2026-09-02.{R,rds,log}` (Stage A), `fastpath_postchange_2026-09-02.{rds,log}` (Stage C, same script), `fastpath_verify_2026-09-02.{R,rds,log}` (gates, routing, profile), `fastpath_profile_F1_2026-09-02.Rprof`, `fastpath_profile_F2_2026-09-02.Rprof`, this report. No push, no render, no payload or application change.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at task start b953705a · git status --porcelain: empty
Gate commits in log: fbd564de (dispatch fix, 0.3.1) · eb136a35 (0.3.2 orientation)
                     d7d5aae2 (re-profile task doc) · b953705a (re-profile artefacts)
packageVersion("forestsearch") 0.3.2 at start · R 4.6.1 · 128 cores
vi.grf.min default in force: NULL (formals(forestsearch)$vi.grf.min)
```

The task document arrived under its exact hyphenated name (`~/Downloads/cc_task_fastpath_unadjusted_2026-09-02.md`) → `dev/tasks/cc_task_fastpath_unadjusted_2026-09-02.md`, committed alone as **`c79f797d`**.

## 2. Design basis re-verified at HEAD (§2 — GATE passed)

`git diff 2a9ceeb6..HEAD -- R/ DESCRIPTION` is **empty** — the re-profile report's reading holds at HEAD verbatim. Line-anchored findings (pre-edit line numbers):

1. **Continuous call site and contract.** `R/subgroup_search.R:618` `glm_result <- fit_glm_for_subgroup(df_clean, id.x, estimator_fn)`. The body (L798–825): `df_sg <- df_clean[id.x == 1, , drop = FALSE]`; `if (nrow(df_sg) < 6L) return(NULL)`; `res <- tryCatch(estimator_fn(df_sg), error = function(e) NULL)`; `if (is.null(res) || is.na(res$estimate)) return(NULL)`; then `z_crit <- 1.96`, `lower <- res$estimate - z_crit * res$se` (upper analogous), returning `list(hr = res$estimate, lower = lower, upper = upper, med0 = NA_real_, med1 = NA_real_)`. **Failure contract: `NULL`** (three ways — tiny slice, estimator error, `NA` estimate).
2. **Continuous closure.** `.make_lm_estimator()` (`R/glm_effect_estimators.R:803`): per call, `fmla <- .build_adjusted_formula(outcome.name, treat.name, if (ps_adjust_method == "none") adjust_covariates else NULL)` (L826–827); the weighted branch `if (ps_adjust_method == "iptw" && "sw" %in% names(data_slice)) lm(..., weights = data_slice$sw) else fit <- stats::lm(fmla, data = data_slice)` (L830–834); `coef_val <- stats::coef(fit)[[treat.name]]`; `se_val <- sqrt(diag(stats::vcov(fit)))[[treat.name]]` (L835–836); error handler returns `list(estimate = NA_real_, se = NA_real_, converged = FALSE, ..., method_used = "ols_failed")`. All branch inputs are `force()`d at construction. `make_effect_estimator()`'s offset enforcement (L121–128) applies to the Poisson-rate measures only — the lm path takes no offset.
3. **Construction site.** `R/forestsearch_main.R:1721–1729` — `estimator_fn <- make_effect_estimator(... adverse_outcome = adverse_outcome, adjust_covariates = adjust_covariates)`, with `treat.name`, `outcome.name`, `adjust_covariates` in scope; `ps_adjust_method` is at this point the **unresolved formal** (possibly the length-3 default vector, first element `"none"`; resolution to `ps_adjust_resolved` happens at L2067). A second construction site at L2133 **rebuilds the closure when PS adjustment resolves active** (`ps_adjust_resolved != "none"`) — the rebuilt closure is untagged, so the fast tag cannot survive onto a weighted path.
4. **Survival fit.** Call site `R/subgroup_search.R:657–659` (passes `df_clean` and `adjust_covariates`). Body L722–748: **the assembly is pure column binding of already-clean vectors** — `data.x <- data.table::data.table(Y = yy, E = dd, Treat = tt, id.x = id.x)` (L728), optional adjustment columns cbind (adjusted arm only), `df.x <- data.x[id.x == 1]` (L741) — no NA handling, no sort, no coercion, no filtering (the vectors were cleaned upstream by `prepare_search_data()`, L322–332: `na.exclude(cbind(Y, Event, Treat, Z))`). Branch: `rhs <- "Treat"` unless `length(adj) > 0L` where `adj <- .fs_adjust_terms(adjust_covariates)` (L725, L731–739). Fit: `hr.cox <- try(summary(survival::coxph(cox_fmla, data = df.x, robust = FALSE))$conf.int, silent = TRUE)` (L744–748). **Downstream consumption:** `if (inherits(hr.cox, "try-error")) return(NULL)` (L750 — the error branch produces no row; the candidate is dropped at status 5); otherwise `trow <- if ("Treat" %in% rownames(hr.cox)) hr.cox["Treat", ] else hr.cox[1, ]` (L758) and the return list reads `trow[1]`, `trow[3]`, `trow[4]` (`exp(coef)`, `lower .95`, `upper .95`) as `hr`, `lower`, `upper` — named scalars whose names travel into the result row, so the fast branch must reproduce the full named 4-column row.
5. **Shared consumers of the closure** (where the fast closure also runs — intended, it is bit-identical): `fit_subgroup_effect()` (`R/Cox_estimation_helpers.R:210`, `tryCatch(estimator_fn(df_sg))`) and `fit_effect_models()` (L235–237), used by the bootstrap H/Hc fits (`R/bootstrap_analysis_dofuture.R:485, 514, 695, 709, 747, 761`; `R/bootstrap_dofuture_main.R:328` — note that site builds its **own** closure at L302, which stays slow); the literal-split consistency path (`R/subgroup_consistency_helpers.R:288–289`); display helpers (`R/summary_utility_functions.R:219`, with own closures at L829, L982); CV builds its own closures (`R/forestsearch_cross_validation.R:1451, 1507`). The consistency screen's closed form fits its **own** `lm` (`R/consistency_resample.R:292`) — not through the closure; it appears in §6's trace counts and is attributed there.
6. **Fit-eligibility floors.** Continuous: status 3 (event floors) is **skipped entirely** — `# Continuous: skip Status 3 entirely — d0.min/d1.min do not apply` (L609–610); the only pre-fit gates are `if (nx <= n.min) return(list(status = 4L, ...))` (L613, strict inequality) and `nrow(df_sg) < 6L` in the fit. **A single-arm slice can therefore reach a continuous fit** (E3 runnable). Survival: `meets_event_criteria(event_counts, d0.min, d1.min)` = `d0 >= d0.min && d1 >= d1.min` (L707–709) plus the same `n.min` gate — perfect separation (zero events in one arm) is **floor-unreachable** for `d0.min ≥ 1` (E5 finding).
7. **Treatment coding at the fit site.** `prepare_search_data()` (L322–332) builds `tt` as a column of `na.exclude(cbind(Y, Event, Treat, Z))` — a **numeric (double) matrix column** regardless of the input's integer/double class; `df_clean[[treat.name]]` retains the raw class (integer `treat_sim` on F1, numeric on F3). Both reach `lm.fit` as the same doubles: `model.matrix` emits a double matrix, and `cbind(1, t)` with a double `1` promotes an integer `t`. The fast path's matrix column is named `treat.name` exactly as `model.matrix` names a numeric treatment column (verified `identical()` against the slow fit in the pre-implementation probe, integer and double).

No substantive difference from the re-profile report's description — gate passed.

## 3. The design as implemented (§4)

**Routing predicate (continuous), `R/forestsearch_main.R:1731–1751`** (decided once, at closure construction):

```r
if (outcome_type == "continuous" &&
    identical(ps_adjust_method[[1L]], "none") &&
    length(.fs_adjust_terms(adjust_covariates)) == 0L) {
  estimator_fn <- .make_lm_estimator_fast(
    treat.name = treat.name, outcome.name = outcome.name,
    adverse_outcome = adverse_outcome)
  attr(estimator_fn, "fs_fast_unadjusted") <- TRUE
}
```

`ps_adjust_method[[1L]]` reads the unresolved formal safely (default vector's first element is `"none"`; a user-passed `"iptw"` makes the predicate FALSE — conservative, and the L2133 PS rebuild independently replaces the closure with an untagged slow one whenever PS adjustment actually resolves active, so the tag cannot leak onto a weighted path). The attribute is the only travelling signal; no signature changed.

**Continuous piece (a) — fast closure**, `R/glm_effect_estimators.R`: `.fs_lm_two_group()` (L896) is the two-group OLS core — `stats::lm.fit()` on the explicit `cbind(1, t)` matrix named `("(Intercept)", treat.name)`, then the slow chain's extraction replicated operation for operation from the R 4.6.1 sources, cited in the roxygen: `summary.lm`'s `rss <- sum(r^2)`, `resvar <- rss/rdf`, its "essentially perfect fit" warning with the exact condition and message, `p1 <- 1L:p`, `R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])`, the pivot handling `z$coefficients[Qr$pivot[p1]]` for `cov.unscaled`'s dimnames; `vcov.summary.lm`'s `sigma^2 * cov.unscaled` with `sigma <- sqrt(resvar)` (**not** `resvar` directly — `sqrt(x)^2` is not `x` bitwise); `stats:::.vcov.aliased()`'s `complete = TRUE` NA-expansion over aliased coefficients; finally `sqrt(diag(vc))[[treat.name]]`. A rank-deficient (single-arm) fit returns `estimate = NA, se = NA, converged = TRUE, method_used = "ols"` without error — exactly the slow chain's behaviour (E3 finding, §5). `.make_lm_estimator_fast()` (L940) wraps the core in the slow closure's exact contract: same signature, same sign-convention negation, same `n0`/`n1`, same `tryCatch` failure list (`"ols_failed"`, `converged = FALSE`).

**Continuous piece (b) — call-site vector routing**, `R/subgroup_search.R:622–627`: when the closure carries the tag, the candidate loop calls `fit_glm_for_subgroup_fast(yy, tt, id.x, estimator_fn)` (L884) — subsets the outcome/treatment **vectors** (no data-frame slice), reads `adverse_outcome` and `treat.name` from the tagged closure's environment (`force()`d there), and reproduces `fit_glm_for_subgroup()`'s contract exactly (< 6 rows → `NULL`; error or `NA` estimate → `NULL`; else the `hr/lower/upper/med0/med1` list with the same `z_crit <- 1.96` Wald arithmetic). The fit call site is unique — `evaluate_combination_with_status()` serves every search path (plain loop and batched both reach it via `search_combinations_parallel()`), so one branch routes all of them. The existing `fit_glm_for_subgroup()` call is verbatim the `else`.

**Survival — which tier and why**, `R/subgroup_search.R:754–785`: **the summary-skip only.** The unadjusted arm (`length(adj) == 0L`) fits the same `coxph(cox_fmla, data = df.x, robust = FALSE)` and computes the consumed `$conf.int` row directly, replicating survival 3.8.9 `summary.coxph()` operation for operation for the non-robust, `conf.int = 0.95`, `scale = 1` case (cited in the comment): `beta <- cox$coefficients`, `se <- sqrt(diag(cox$var))`, `z <- qnorm((1 + conf.int)/2, 0, 1)`, `tmp <- cbind(exp(beta), exp(-beta), exp(beta - z*se), exp(beta + z*se))` with dimnames `list(names(beta), c("exp(coef)", "exp(-coef)", "lower .95", "upper .95"))` (the `* scale` multiplications omitted — `x * 1` is exact in IEEE-754; the label `paste("lower .", round(100*0.95, 2), sep = "")` evaluates to `"lower .95"`, written literally). The `try(..., silent = TRUE)` wrapper and NULL-on-error contract are preserved identically. The **assembly is kept although §2.4 verified it is pure column binding**: `df.x` is consumed two lines later by the untouched medians computation (L788–791) — removing the assembly is the medians task's business (rank 2), not this one's, so the vector-fit tier would have saved nothing while medians remain. The adjusted arm keeps the original `summary(coxph(...))$conf.int` line verbatim.

**Housekeeping:** `DESCRIPTION` → 0.3.3; `NEWS.md` 0.3.3 entry; new internals documented `@keywords internal` + `@noRd`; `devtools::document()` clean (no Rd delta). `chol2inv` is base; `stats::`/`survival::` used with prefixes or existing imports — no NAMESPACE change.

**Pre-implementation probe** (scratchpad, before any edit): the prototype core was verified `identical()` against the installed slow closure on 8 lm cases (normal double/integer treat, `adverse_outcome = FALSE` negation, zero-variance y **including the warning**, single-arm double/integer, n = 6, 2-vs-78 unbalanced) and the conf.int replication against `summary.coxph` on 3 Cox cases (normal, heavy ties, monotone-likelihood degenerate) — **11/11 identical** before the design was committed to.

## 4. Stage A / Stage C mechanics and the equality gates (§3, §5)

**Volatile-field exclusion list, verbatim — the complete set, from evidence:**

- list fields: `find.grps$time_search`, `minutes_all` (wall-clock timings of the search and the whole call);
- bootstrap frame columns: `tmins_search`, `tmins_iteration` (per-replicate wall-clock);
- normalisation (not an exclusion): `data.table` → `data.frame` at three nodes (`grp.consistency$out_sg$result`, `$pareto_frontier`, `find.grps$out.found$hr.subgroups`) — their `.internal.selfref` external pointer differs across sessions; **all contents compared in full**.

Nothing else is excluded: no functions, environments or calls exist anywhere in the result object (verified by a recursive walk; the pruner would have flagged and recorded any). Every dropped path is recorded in the rds (32 unique, all four names above across the 15 saved objects). **Completeness proof:** F1, E1 and the 20-replicate bootstrap each ran twice within each stage's session — the two pruned sets `identical()` in both stages (`self-consistency F1/E1/bootstrap: TRUE` ×2), so the exclusion list is exactly sufficient, before and after the change.

**The equality matrix** (Stage A at 0.3.2, installed; Stage C at 0.3.3, installed; same script, same seeds, same order — `FS_FASTPATH_OUT` selects the output name only). Every fixture × every retained component, `identical()`:

| gate | fixtures | components compared | result |
|---|---|---:|---|
| **E-1** fast-path | F1 cont · F2 gbsg · F3 anchor | 25 + warnings + counters + sg.harm each | **all `identical()` TRUE** |
| **E-2** routing | F4 cont-adjusted · F5 surv-adjusted · F6 IPTW | same | **all `identical()` TRUE** |
| **E-3** edges | E1–E6 | same | **all `identical()` TRUE** |
| **E-4** bootstrap | F1, `nb_boots = 20`, seed 8316952, sequential | full pruned payload | **`identical()` TRUE** |

The 25 components are the full top-level result set (`grp.consistency`, `find.grps` incl. `out.found$hr.subgroups` at full precision, `df.est`, `df.predict`, `args_call_all`, `threshold_config`, `family_status`, `admission`, …). Selections and counts, both stages: F1 `!{wtkg <= 84} & !{cd40 <= 320}` (1 975 fits, 749 screened, dedup 852 → 749); F2 `{er <= 0} & {size <= 35}` (1 410 / 120 / 121 → 120); F3 `{age <= 37} & !{cd40 <= 507}` (4 702 / 1 / 146 → 131); F4 selects F1's subgroup via the adjusted fits (781 screened); F5 `{er <= 0} & {pgr <= 8}`; F6 selects nothing (see finding below).

**E-case findings:**

- **E1 (floor boundary, continuous):** the size gate is strict — `if (nx <= n.min) return(list(status = 4L, ...))` (`subgroup_search.R:613`) — so the task's literal `n.min = 75` would *exclude* the 75-subject slice, not fit it. `n.min = 74` was used so the slice sits exactly at the smallest admissible size, per the stated purpose; recorded as a deviation. Fit ran at n = 75, both arms; identical.
- **E2 (zero-variance slice):** the constant-y slice fit ran (3 fits); `.fs_lm_two_group()` replicates `summary.lm`'s "essentially perfect fit" **warning** (this R version warns, it does not stop) at the same point. Note the warning is muffled inside the search identically on both paths — the captured warning vector at the `forestsearch()` level is empty in both stages (the pre-implementation probe confirmed both paths emit it when the closure is called directly).
- **E3 (rank-deficient slice):** **reachable and exercised** — continuous has no per-arm floor (§2.6), so the `z1 = treat` construction put a single-arm slice through the fit: 2 of 4 fits returned `NULL` in both stages (the slow chain's `NA`-coefficient path: `vcov.lm`'s `complete = TRUE` aliased expansion gives `se = NA` without error, `converged = TRUE`, and the `NA` estimate maps to `NULL` at the fit contract). Identical.
- **E4 (heavy ties):** integer times 1:6, Efron handling through the same `coxph()` call on both branches; identical.
- **E5 (Cox failure contract):** true perfect separation is **floor-unreachable**: all events in one arm forces `d0 = 0` (or `d1 = 0`), and `meets_event_criteria()` requires `d0 >= d0.min && d1 >= d1.min` (`subgroup_search.R:707–709`) with floors ≥ 1 — demonstrated by the quoted gate, not weakened. The floor-passing degeneracy used instead: monotone likelihood (every treated event precedes every control event in the slice) — `coxph()` warned twice (identically captured in both stages), produced the identical extreme row through both branches, and the document's `try()` error branch remains verbatim for a genuine `coxph()` error.
- **E6 (boundary events):** treated events exactly at `d1.min = 5`; the `>=` gate admits it; boundary fit identical.

**A pre-existing defect surfaced by F6 (recorded, deliberately not fixed):** at 0.3.2 *and* 0.3.3, the IPTW continuous configuration (`ps_method = "logistic"`, `ps_adjust_method = "iptw"`) has **all 1 975 candidate fits fail** (`method_used = "ols_failed"`, search finds nothing): inside the closure, `stats::lm(fmla, data = data_slice, weights = data_slice$sw)` errors because `lm`'s `weights` argument is evaluated where `data_slice` is not visible (the formula's environment is `.build_adjusted_formula()`'s frame). Reproduced minimally outside the search. Fixing it here would change F6's results and violate the equality guard; it belongs to a follow-up task.

## 5. Routing proof (§6)

Trace counters on `stats::lm`, `survival::summary.coxph`, `fit_glm_for_subgroup[_fast]`, and the attribute as observed by `subgroup.search()`, with 0.3.3 installed:

| fixture | predicate / attr on search | fast fits | slow fits | `stats::lm` | `summary.coxph` |
|---|---|---:|---:|---:|---:|
| F1 cont | **TRUE** | **1 975** | 0 | 749 | 0 |
| F2 gbsg | FALSE (survival: branch is condition-local, not attribute) | 0 | 0 | 0 | **0** |
| F3 anchor | **TRUE** | **4 702** | 0 | 1 | 0 |
| F4 cont-adjusted | FALSE | 0 | **1 975** | 2 756 | 0 |
| F5 surv-adjusted | FALSE | 0 | 0 | 0 | **1 410** |
| F6 IPTW | FALSE | 0 | **1 975** | 1 975 | 0 |

- **Zero search-phase `lm`/`summary.coxph` on every fast route.** The residual `stats::lm` counts are fully attributed: F1's 749 = one per candidate reaching the consistency screen — the closed-form consistency's **own** fit, `MD = stats::lm(stats::reformulate(rhs, outcome.name), data = df)` at `R/consistency_resample.R:292`, not the estimator closure; F3's 1 = its single screen candidate (default `stop_threshold` early stop); F4's 2 756 = 1 975 search fits + 781 screen candidates; F6's 1 975 = the failing IPTW search fits. No unexplained counts.
- **The predicate routes both ways:** slow fits and `summary.coxph` are nonzero exactly on the routing fixtures (F4: 1 975 slow; F5: 1 410 `summary.coxph`), zero on the fast ones.
- F2's proof is the `summary.coxph` column: zero dispatches in the entire call (the survival fast branch keys on `length(adj) == 0L` inside the fit, not on the attribute — the closure does not exist on the survival path).

## 6. Realized recovery (§7)

Profiles under the re-profile script's exact procedure and bucket definitions (`Rprof` 0.010 s, `fit_glm_for_subgroup_fast` added to the fit bucket's name list):

**F1 continuous** (0.3.2 → 0.3.3): wall **4.74 → 1.91 s**; fit bucket **2.93 → 0.25 s** (61.9% → 13.2%) — **2.68 s recovered, matching the predicted ≈2.6–2.8 s**. The call is now the consistency evaluator's: 1.36 s, **71.6%** (the re-profile's rank 3 is the next lever). Everything else unchanged (evaluator 1.40 → 1.36 s, future 0.31 → 0.22 s).

**F2 gbsg** (0.3.2 → 0.3.3): wall **7.59 → 7.47 s**; fit bucket 3.97 → 3.77 s; medians 2.69 → 2.81 (noise). **Realized recovery ≈0.1–0.2 s — well short of the predicted ≈1.6 s**, and the shortfall is a finding about the *prediction*, not the change: re-parsing both raw profiles shows only **9 of 397** old fit-bucket samples were actually inside `summary.coxph`; the 09-01 report's "vcov/summary 40.3% (1.60 s)" sub-attribution was an artifact — `summary(coxph(...))` evaluates `coxph()` as `summary`'s *argument*, so the bare `summary` frame sat on the stack for the entire Cox fit (245 samples carried `summary` without `summary.coxph`) and the sub-split's priority order credited coxph's own wrapper time to the vcov/summary bin. The remaining fit-bucket samples sit where they always were: `coxph()`'s internals — assembly/subsetting (143 samples), `model.frame`/`model.matrix` (51), `coxph.fit` (33), wrapper (139) — i.e., §9's out-of-scope "deeper Cox surgery" item, plus the medians (rank 2, next task). The `summary.coxph` skip is correct, bit-identical, and worth its true ~0.1–0.2 s; recorded, no further edits here.

**Bootstrap (end-to-end):** Stage A 96.4 s / **4.82 s per replicate** (0.3.2); Stage C 48.8 s / **2.44 s per replicate** (0.3.3) — **2.06× faster**; **B = 1000: 83.7 → 40.7 min sequential**, matching the task's "toward ≈40 min".

## 7. Package health (§8 — GATES passed)

- `devtools::test()` on 0.3.3 source: **0 failures, 4 864 pass, 3 skips, WARN 31** — warning parity with the recorded 0.3.2 baseline (31), and **`test-search-reproducibility.R` passes unmodified** (the file the first dispatch attempt broke).
- `R CMD check` (`devtools::check(document = FALSE, vignettes = FALSE, manual = FALSE)`, `RSTUDIO_PANDOC` set): **`Status: OK` — 0 errors / 0 warnings / 0 notes** (8 m 52 s). Before-set from the record: 0.3.2's check at `eb136a35` was **0 errors / 0 warnings / 0 notes** (`REPORT_oc_signed_orientation_2026-08-31.md`), and `R/` was unchanged between that commit and this task's start.

## 8. Commits

1. `c79f797d` — task doc alone.
2. `b9d7b622` — Stage A script + baseline rds + log (before any `R/` edit).
3. `4dbf9f26` — implementation: `R/glm_effect_estimators.R`, `R/subgroup_search.R`, `R/forestsearch_main.R`, `DESCRIPTION` (0.3.3), `NEWS.md`.
4. (this commit) — Stage C rds + log, verify script + rds + log, two profile `.Rprof` files, this report.

Working tree clean at close. **No push.**

## 9. Verdict (ten lines)

1. The unadjusted fast path is in at 0.3.3: continuous `lm.fit` closure + vector call-site routing tagged at construction; survival direct conf.int replicating `summary.coxph` 3.8.9 operation for operation; adjusted/weighted/PS paths preserved verbatim as the `else`.
2. **Every equality gate passed with bit identity**: 12 fixtures × 25+ components + warnings + counters, and the 20-replicate bootstrap payload, all `identical()` across 0.3.2 → 0.3.3; the volatile exclusion is exactly four timing fields, completeness proven by within-stage double runs.
3. Routing proven: zero search-phase `lm`/`summary.coxph` on fast routes (F1/F2/F3), nonzero on slow routes (F4/F5/F6); every residual count attributed to source (`consistency_resample.R:292`; the screen's own fits).
4. Continuous recovery **on prediction**: fit bucket 2.93 → 0.25 s, single call 4.74 → 1.91 s.
5. Bootstrap replicate **4.82 → 2.44 s**; **B = 1000: 83.7 → 40.7 min** sequential.
6. Survival recovery ≈0.1–0.2 s, far below the predicted 1.6 s: the 09-01 sub-split's 40.3% "`summary.coxph`" was a stack-attribution artifact (bare `summary` frame during argument evaluation of `coxph()`); the survival fit cost is coxph internals + medians — rank 2 (medians) and the out-of-scope Cox surgery remain the survival levers.
7. Edge contracts verified bit-identical: floor-boundary fit, zero-variance warning path, rank-deficient `NA` propagation (reachable — continuous has no per-arm floor), heavy ties, monotone-likelihood degenerate; true perfect separation floor-unreachable by the quoted event gate.
8. Pre-existing defect recorded, untouched: the IPTW continuous closure fails every fit at 0.3.2 and 0.3.3 alike (`weights = data_slice$sw` resolution); fixing it would have violated the equality guard.
9. Tests 0 failures / WARN 31 parity with `test-search-reproducibility.R` unmodified; check `Status: OK` 0/0/0 against a recorded 0/0/0 baseline.
10. The continuous replicate is now the consistency evaluator's (71.6%) — the next lever is rank 3's count, then rank 2's medians on the survival path.
