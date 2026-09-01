# CC task — unadjusted fast path for the per-candidate fit, under a bit-identity equality guard

**File:** `dev/tasks/cc_task_fastpath_unadjusted_2026-09-02.md` · **Issued:** 2026-09-02 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1). If the file is absent under the exact name above, the hyphen-stripped stem `cc_task_fastpath_unadjusted_20260902.md` is the same file (as with the two prior tasks).
**Reads (in-repo):** `dev/glm-continuous-sims/REPORT_bootstrap_reprofile_2026-09-01.md` (the measurement this task implements against — its §4(a) is the design basis), `dev/glm-continuous-sims/bootstrap_reprofile_2026-09-01.R` (the fixture blocks, reused verbatim), `dev/glm-continuous-sims/REPORT_bootstrap_profile_2026-08-30.md` (background only).

**Why.** The re-profile at 0.3.2 established that the per-candidate fit bucket is 61.9% of a continuous replicate-configuration call (2.93 s of 4.74 s) and 52.4% of a gbsg call (3.97 s of 7.59 s), and that **the actual solve is 2–7% of that bucket** — the rest is data-frame subsetting/assembly (58.7% continuous, 63.9% anchor), `model.frame`/`model.matrix`, and on survival `summary.coxph` at 40.3% (1.60 s to read one confidence interval). It further established that an unadjusted fast path is routable by a construction-time predicate on both outcome types, leaving the adjusted/weighted paths untouched. This task implements that routing. The measured prize: up to ≈2.6–2.8 s per continuous replicate (B = 1000 from ≈84 toward ≈40 min) and ≈1.6–2.6 s per gbsg call.

**The governing rule, from the brief:** *a performance change that alters a result is a defect regardless of its speed.* The guard here is **bit identity** — `identical()` on the result components, not tolerance — plus routing proof and an edge-case smoke suite.

---

## ⚠ CATEGORY — this task edits `R/`

Called out per protocol, distinguishing the three kinds:

- **Adds routed code paths** (the substance): a fast closure variant for the continuous lm estimator; a fast branch inside the survival fit; a call-site branch in the search loop. When the routing predicate is FALSE, no existing line's behaviour changes — the slow paths are preserved verbatim as the `else`.
- **Moves existing code:** none. Nothing is relocated or deleted.
- **Changes the method:** no. Same estimators, same models, same selection. Outputs must be bit-identical where the predicate holds and trivially identical (untouched code) where it does not.

Files expected to change: `R/glm_effect_estimators.R`, `R/subgroup_search.R`, `R/forestsearch_main.R` (one attribute line at closure construction), `DESCRIPTION` (version → 0.3.3), `NEWS.md`. No driver, no application document, no payload, no test file is edited. If implementation genuinely requires touching any other `R/` file, **STOP** and report why.

**Compute:** Stage-A baselines ≈4–5 min; two install/test cycles; `R CMD check` (no manual, no vignettes) ≈5–10 min; Stage-C verification ≈4–5 min; post-change profile + 20-replicate bootstrap ≈4 min. **Estimate 60–90 minutes wall total.** No simulation campaign, no renders.

**Unattended.** Gates stop the task; never ask, never work around. On any equality-gate failure: report the component and the first differing values, leave the tree in a committed, documented state, and stop. Provenance gate is dirt-tolerant; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block: machine, repo path, branch, HEAD, `git status --porcelain`, `packageVersion("forestsearch")` (expect **0.3.2**), R version, cores. *GATE:* branch `feature/glm-extension`; commits `fbd564de` (dispatch fix), `eb136a35` (0.3.2 orientation), `d7d5aae2` (re-profile task doc), and the re-profile artefacts commit (`b953705a` expected) all in the log. Copy this document into `dev/tasks/` and commit it **alone**; report both filenames and the commit hash. Record the `vi.grf.min` default in force (expect `NULL`).

## 2. Re-verify the design basis at HEAD — GATE

The report describes source at `2a9ceeb6`; HEAD has moved (artefacts commit). `b953705a` touched only `dev/`, so `R/` should be unchanged — verify, don't assume. Re-read at current HEAD and quote line-anchored:

1. **Continuous call site:** `R/subgroup_search.R` ≈L618, `fit_glm_for_subgroup(df_clean, id.x, estimator_fn)`; the body of `fit_glm_for_subgroup()` (≈L798) with its `df_sg <- df_clean[id.x == 1, , drop = FALSE]` and its **exact return contract** — every field, and what it returns on fit failure (`NULL`? `NA` fields? quote the code).
2. **Continuous closure:** `.make_lm_estimator()` (`R/glm_effect_estimators.R` ≈L803–840): the formula build (≈L826), the `stats::lm()` call (≈L833), the coef and se extraction (≈L835–836), the branch conditions (`ps_adjust_method == "iptw"` weighted branch ≈L830; adjusted formula requires non-empty `.fs_adjust_terms(adjust_covariates)` under `"none"`), and the fact that all branch inputs are `force()`d at construction. Also `make_effect_estimator()`'s offset enforcement (≈L121–128, Poisson-rate path only).
3. **Closure construction site:** `R/forestsearch_main.R` ≈L1721–1729 — where `make_effect_estimator(...)` is called, with `ps_adjust_method` and `adjust_covariates` in scope.
4. **Survival fit:** `fit_cox_for_subgroup()` (`R/subgroup_search.R` ≈L722–748): the call site (≈L657–659), the data assembly (what exactly is built from `yy, dd, tt, id.x, df_clean` — quote it, and state whether it is **pure column binding of already-clean vectors** or performs any transformation: NA handling, sorting, coercion, filtering), the single branch `length(.fs_adjust_terms(adjust_covariates)) > 0L` (≈L732), the per-candidate `summary(survival::coxph(...))$conf.int` inside `try(..., silent = TRUE)` (≈L744–746), and the **exact downstream consumption**: which columns of `$conf.int` are read, in what order, and what happens when the `try()` errs (quote the error-branch code and the row it produces).
5. **Shared consumers of the continuous closure:** locate every caller of `estimator_fn` / the closure beyond the search loop (MR/consistency machinery, the bootstrap's H/Hc fits, CV). List them by file:line. This list defines where the fast closure will also run — which is intended, because it is bit-identical — and where the trace counters in §6 may legitimately show residual calls.
6. **Fit-eligibility floors:** where `n.min`, `d0.min`, `d1.min` are enforced relative to the fit — quote the gate. This determines whether a single-arm or zero-event slice can ever reach a fit (matters for §5's edge cases E3/E6).
7. **Treatment coding at the fit site:** establish from source the class and values of `df_clean[[treat.name]]` (and `tt` on the survival path) at the point of fitting — numeric 0/1, integer, factor? Quote the coercion (e.g., under `is.RCT`) that guarantees it. The fast path's model matrix must reproduce `model.matrix`'s column for this exact class.

*GATE:* if any of these differs **in substance** from the report's description (not mere line drift), STOP and report the difference — the design below assumes the report's reading.

## 3. Stage A — pre-change baselines, at HEAD, before any edit

One script, `dev/glm-continuous-sims/fastpath_baseline_2026-09-02.R`, run at the §1 HEAD with 0.3.2 installed. It produces `fastpath_baseline_2026-09-02.rds` holding, for every case below, the **full relevant result components** (see the volatile-field rule at the end of this section). Reuse the three fixture blocks **verbatim from `bootstrap_reprofile_2026-09-01.R`** — do not re-derive them.

**Core fixtures (predicate TRUE expected):**
- **F1 continuous** — MD40 fixture, seed 8316951+1, driver arguments, replicate configuration.
- **F2 survival** — gbsg effMaxSG fixture.
- **F3 anchor** — ACTG175 applied anchor, its document's arguments verbatim.

**Routing fixtures (predicate FALSE expected — these prove the untouched paths stay untouched):**
- **F4 continuous adjusted** — F1's call plus `adjust_covariates = c("age", "wtkg")` (verify both names are in F1's confounder set; if not, substitute two that are and record which).
- **F5 survival adjusted** — F2's call plus `adjust_covariates = "age"`.
- **F6 IPTW** — only if `ps_adjust_method` (or its equivalent) is reachable from the `forestsearch()` signature on the continuous path: F1's call with the IPTW setting and a minimal propensity configuration. Establish reachability from `formals(forestsearch)` and the §2.3 construction site; if not user-reachable, **record that finding and skip F6** — do not force it through internals.

**Edge fixtures (synthetic, fixed seeds, small — each runs in seconds):**
- **E1 floor-boundary subgroup (continuous).** `set.seed(20260902)`; n = 150; one binary confounder `z1 = rep(0:1, each = 75)`; treatment balanced within `z1`; `y = rnorm(150)` plus a small `z1` effect; a second noise confounder so the search has ≥2 candidates. Call `forestsearch()` with `n.min = 75` so the `z1 == 1` candidate sits exactly at the floor. Purpose: fits at minimum size, both arms present.
- **E2 zero-variance slice (continuous).** As E1 but `y` **constant** within `z1 == 1` (e.g., `y = ifelse(z1 == 1, 5, rnorm(150))`). Purpose: `sigma^2 = 0`, `se = 0` — the fast path must reproduce the slow path's exact zeros/values.
- **E3 rank-deficient slice (continuous).** Construct data where some enumerable candidate slice, *if it reaches a fit under the §2.6 floors*, contains a single treatment arm — `lm(y ~ treat)` then yields `NA` for the treat coefficient. If the floors provably gate this out before any fit (per §2.6), **demonstrate that with the gate's quoted condition instead of a run**, and record the edge as unreachable. Do not weaken floors to force it.
- **E4 heavy ties (survival).** Small synthetic survival set (n = 120, fixed seed) with many tied event times (integer times drawn from 1:6), two binary confounders, floors at defaults. Purpose: Efron handling flows through identically — the fast branch keeps `coxph()`, so this is a confirmation, not a risk.
- **E5 coxph failure contract (survival).** Small synthetic set engineered so at least one candidate slice passing the floors gives a degenerate Cox fit (e.g., perfect separation: all events in one arm; monotone likelihood). Purpose: the `try()` error/warning path — the fast branch must produce **the identical downstream row** the slow branch produces for that candidate. If the floors make even this unreachable, demonstrate with the quoted gate and record it.
- **E6 near-empty events at the floor (survival).** Events in one arm exactly at `d1.min` (or `d0.min`). Purpose: boundary fit succeeds and matches bitwise.

For each case F1–F6 and E1–E6 that runs: execute the full `forestsearch()` call, save the result object's comparison set, and also save a **20-replicate bootstrap** for F1 only (`forestsearch_bootstrap_dofuture()`, `nb_boots = 20`, seed 8316952, sequential — the §5 configuration of the re-profile task) as the end-to-end baseline.

**Volatile-field rule.** Before saving, identify from source the result components that legitimately vary run-to-run (timings, wall-clock fields, calls/matched calls, environments, RNG bookkeeping not consumed downstream). List them explicitly in the report, exclude exactly those, and keep **everything else** — selected subgroup definition and membership, `out.found` / `hr.subgroups` (all columns, full precision), consistency outputs, per-candidate frames, bootstrap payload frames. The exclusion list is part of the deliverable; an over-broad exclusion is a gate failure in spirit.

Commit the script and the baseline `.rds` (with a short log) as the second commit.

## 4. Implement — the design

Three components. Slow paths remain verbatim; every new path is selected by a predicate that is FALSE by default construction only when the configuration requires the old path.

### 4.1 The routing predicate

**Continuous (decided once, at closure construction, `forestsearch_main.R` ≈L1721):**

```
fast_ok <- outcome_type == "continuous" &&
           ps_adjust_method == "none" &&
           length(.fs_adjust_terms(adjust_covariates)) == 0L
```

using the in-scope names verified in §2.3 (adjust spelling to source; the three conditions are the spec). When TRUE, tag the closure: `attr(estimator_fn, "fs_fast_unadjusted") <- TRUE`. The attribute is the only signal that travels — no signature changes anywhere.

**Survival (decided at the fit, inside `fit_cox_for_subgroup()`):** the existing branch condition at ≈L732 already computes it — `length(.fs_adjust_terms(adjust_covariates)) == 0L` selects the treatment-only model. The fast branch lives inside that arm.

**Scope guard:** binary/count outcomes are **out of scope** (§9): the closure fast variant is built only under `outcome_type == "continuous"`; no glm-family fast path of any kind.

### 4.2 Continuous — two pieces

**(a) Fast closure** (in `R/glm_effect_estimators.R`, beside `.make_lm_estimator()`; internal, `@noRd`). Same signature as the slow closure (takes the data slice), so every §2.5 consumer works unchanged. Body: extract `y <- data_slice[[outcome.name]]` and the treatment column per §2.7's verified class; build `X` as exactly the two-column matrix `model.matrix(~ treat)` would produce for that class (for numeric 0/1: `cbind("(Intercept)" = 1, treat)` with the treat column named as `model.matrix` names it — verify the name against one slow-path fit); then:

```
z      <- stats::lm.fit(X, y)
rss    <- sum(z$residuals^2)
rdf    <- length(y) - z$rank
sigma2 <- rss / rdf
p1     <- seq_len(z$rank)
R      <- chol2inv(z$qr$qr[p1, p1, drop = FALSE])
se_all <- sqrt(diag(R) * sigma2)
```

with coefficient and se taken at the treat position **respecting `z$qr$pivot`**, exactly as `summary.lm()` does — copy `summary.lm`'s pivot handling for the extraction, cite the lines copied in a comment. Return exactly the slow closure's output structure (fields, names, types — §2.2's contract). Rank-deficiency (`z$rank < ncol(X)`): return precisely what the slow path returns for that case, established from the E3 finding — if E3 is unreachable, implement the same `NA` propagation `lm` + the extraction lines would give, and say so in the report.

`lm.fit` on the identical model matrix is the same QR the slow path runs, and the se ops above are `summary.lm`'s own — this is why bit identity is achievable, and the §6 gate is where that claim is checked, not assumed.

**(b) Call-site vector routing** (in `R/subgroup_search.R`). This is where the measured 58.7% lives — the per-candidate `df_clean[id.x == 1, , drop = FALSE]` over the full frame. At the search's candidate loop: when `isTRUE(attr(estimator_fn, "fs_fast_unadjusted"))`, hoist **once** before the loop `yvec <- df_clean[[outcome.name]]` and `tvec <- df_clean[[treat.name]]` (names per §2.3/§2.7 scope), and per candidate call a new internal `fit_glm_for_subgroup_fast(yvec, tvec, id.x)` that subsets the two vectors and runs the §4.2(a) operations directly (no data-frame slice at all), returning the identical contract of `fit_glm_for_subgroup()` including its failure behaviour. The existing `fit_glm_for_subgroup()` and its call remain verbatim as the `else`. If the loop exists in more than one place (plain loop and batched path — §2 will show), route **both**, identically.

### 4.3 Survival — skip `summary.coxph`, conditionally skip assembly

Inside the unadjusted arm of `fit_cox_for_subgroup()`:

- **Always:** replace `summary(coxph(...))$conf.int` with fitting `fit <- survival::coxph(<same formula/data as now>, robust = FALSE)` and computing the consumed quantities directly: `co <- fit$coefficients[[1]]`, `se <- sqrt(fit$var[1, 1])`, `z975 <- stats::qnorm(0.975)` (or the exact z `summary.coxph` uses for its default `conf.int = 0.95` — read `survival:::summary.coxph` and replicate its conf.int arithmetic **operation for operation** for the one-coefficient case, citing the lines), producing exactly the `$conf.int` entries the downstream consumes per §2.4, in the same structure. Preserve the `try(..., silent = TRUE)` contract exactly: same wrapper, same error handling, same row on failure.
- **Only if §2.4 verified the assembly is pure column binding of already-clean vectors** (no NA handling, no sort, no coercion, no filtering): additionally skip the per-candidate table assembly and fit on subset vectors directly (`survival::coxph(survival::Surv(yy[idx], dd[idx]) ~ tt[idx], robust = FALSE)` with `idx <- id.x == 1`), preserving row order. If the assembly does **anything** beyond binding, keep it and take only the summary skip — state which tier was implemented and why, quoting the assembly line that decided it.

**Medians are untouched.** The `survfit` medians computation (≈L761–763) stays exactly where it is — that is the next task, not this one.

### 4.4 Housekeeping in the same commit

`DESCRIPTION` Version → **0.3.3**. `NEWS.md`: one entry under 0.3.3 — performance: unadjusted fast paths for the per-candidate fit (continuous lm, survival Cox), routed by construction-time predicate, verified bit-identical; adjusted/weighted paths unchanged. Roxygen for the new internals: `@keywords internal`, `@noRd`, markdown-on conventions (literal `% < > &`). Tidyverse style; no `:::` into other packages in package code (reading `survival:::summary.coxph` to *replicate* its arithmetic is fine; *calling* it is not); confirm `stats` and `survival` imports cover `lm.fit`, `chol2inv`, `qnorm`, `coxph` as used. `devtools::document()` runs clean.

Commit the implementation (R files + DESCRIPTION + NEWS + generated Rd if any) as the third commit. Then reinstall.

## 5. Stage C — the equality gates, post-change

Rerun the **same** Stage-A script against 0.3.3 (same seeds, same order), saving `fastpath_postchange_2026-09-02.rds`. Then compare, component by component, with `identical()` under the §3 volatile-field exclusion:

- *GATE E-1 (fast-path fixtures F1, F2, F3):* every retained component `identical()` TRUE. Selected subgroup, `out.found` / `hr.subgroups` full precision, consistency outputs, per-candidate frames.
- *GATE E-2 (routing fixtures F4, F5, F6-if-run):* every retained component `identical()` TRUE — these ran untouched code; any difference means the change leaked.
- *GATE E-3 (edge fixtures E1–E6 as runnable):* `identical()` TRUE per case; unreachable edges documented per §3 with the quoted gate.
- *GATE E-4 (bootstrap):* the F1 20-replicate bootstrap payload `identical()` TRUE against Stage A, volatile fields excluded.

On any failure: report the case, the component, and the **first differing values verbatim** (both sides, full precision), then STOP. No tolerance fallback, no "close enough", no partial acceptance.

## 6. Routing proof and the smoke verdict

Equality proves harmlessness; this section proves the fast path actually **ran**.

- With 0.3.3 installed, wrap F1 and F2 single calls with `trace()` counters on `stats::lm` (F1) and on `summary.coxph` via `trace(summary, signature/where appropriate — implement however cleanly counts `summary.coxph` dispatches)` (F2). Expected: **zero** search-phase calls on the fast routes. Any nonzero count must be attributed by source line to a §2.5 consumer or a display site and explained; unexplained counts are a STOP.
- Repeat on F4 and F5: counts must be **nonzero** (slow path active), confirming the predicate routes both ways.
- Print, for each of F1–F6, the predicate's evaluated value and the attribute's presence — a one-line routing table for the report.

## 7. Realized recovery — measure, don't extrapolate

- Re-profile **one** single call each for F1 and F2 under `Rprof(interval = 0.010, line.profiling = TRUE)`, the re-profile script's procedure and bucket definitions verbatim. Report the bucket table beside the 09-01 report's: fit-bucket seconds and share, before → after, and total wall before → after.
- Rerun the F1 bootstrap timing (already produced in §5's Stage-C run): per-replicate mean, against 5.02 s pre-change, and the B = 1000 projection against 83.7 min.
- One line comparing realized recovery to the report's predicted ≈2.6–2.8 s (continuous) and ≈1.6–2.6 s (gbsg, tier-dependent). If realized falls well short, say where the remaining samples sit — that is a finding for the record, not a trigger for further edits in this task.

## 8. Package health gates

- `devtools::test()`: full suite passes, with `test-search-reproducibility.R` called out explicitly (the file the first dispatch attempt broke — it must pass **unmodified**). *GATE.*
- `R CMD check` with `--no-manual --no-vignettes --no-build-vignettes` (or `devtools::check(document = FALSE, vignettes = FALSE, manual = FALSE)`): no new ERRORs, no new WARNINGs, and NOTEs limited to those present before the change (record the before set by checking once at Stage A if not already known from the record — a single line suffices). *GATE.*

## 9. Out of scope — explicit

- **The medians move (rank 2).** Next task, after this report lands.
- **Binary/count fast paths.** Closed-form or matrix-level shortcuts for glm families cannot be certified bit-identical the way lm/QR can, and the OR work is a separate workstream with its own open decisions.
- **Deeper Cox surgery** (`coxph.fit`/`agreg.fit` direct calls, pre-sorting): risk outweighs the residual ≈0.7 s; revisit only if §7 shows an unexpected remainder.
- **`stop_threshold` semantics, batch-size defaults, `parallel_args` defaults:** recorded items, untouched here.
- **No renders, no push, no payload or application changes.**

## 10. Report

`dev/glm-continuous-sims/REPORT_fastpath_unadjusted_2026-09-02.md`:

1. Provenance; the §2 re-verification with quotes at HEAD, including the return contracts, the assembly finding (pure binding or not, and therefore which §4.3 tier), the treatment-coding finding, and the §2.5 consumer list.
2. The design as implemented: files and functions touched, the predicate as written, the attribute mechanism, which survival tier, the `summary.lm`/`summary.coxph` lines replicated (cited).
3. Stage A/C mechanics: the volatile-field exclusion list, verbatim; the equality matrix — every fixture × every retained component, `identical()` result; the E-case findings including any unreachable edges with their quoted gates.
4. The routing table and trace counts, with attributed residual calls if any.
5. §7's before/after bucket tables, the bootstrap per-replicate and B = 1000 line, and realized-vs-predicted.
6. §8's test-suite and check results (the NOTE set before/after).
7. Ten-line verdict.

Commits, in order: (1) this task doc alone; (2) Stage-A script + baseline artefacts; (3) the implementation + version/NEWS/docs; (4) Stage-C/verification artefacts + post-change profile outputs + the report. Explicit paths at every stage; working tree clean at close. **No push. No render. Nothing outside the files named in the CATEGORY block plus `dev/` artefacts.**
