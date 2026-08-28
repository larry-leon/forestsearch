# REPORT — OC-wrapper build: `fs_oc_predict()` and `fs_oc_family_enumerate()`

**Task:** `dev/tasks/cc_task_oc_wrapper_build_2026-08-28.md` (commit `47769eb3`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Evidence base:** `dev/glm-continuous-sims/REPORT_oc_wrapper_inventory_2026-08-28.md` (A3 S1–S6, A5, C2), read before writing code.
**Category (as the task states it):** ADDS files to `R/`; edits no existing `R/` file; changes no behaviour, no method. Exports exactly two new functions.

---

## 0. Provenance (§1 — GATE passed)

```
hostname                 pop-os
pwd                      /home/larryleon/Documents/GitHub/forestsearch
branch                   feature/glm-extension
HEAD at start            c016df02   (the inventory report commit; descends from itself)
git status --porcelain   (empty — clean tree at start)
git log --oneline -4     c016df02 report — OC-wrapper step 0 ...
                         86342569 task doc — OC-wrapper step 0 ...
                         6afb3df3 MRCT literature search
                         27dd09e2 docs(reports): record the 0.2.2-vs-0.2.0 ...
packageVersion           forestsearch 0.2.2
```

`ls ~/Downloads/cc_task_oc_wrapper_build*` matched **exactly one** file:
`cc_task_oc_wrapper_build_20260828.md` (hyphens stripped from the date stem, as the task anticipated). Copied to `dev/tasks/cc_task_oc_wrapper_build_2026-08-28.md` and committed as the first action: **`47769eb3`**.

**Parity baseline** — `devtools::test()` before any edit:

```
[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4514 ]
```

---

## 1. §2 — the `sigma_D` question, from source (GATE: closed → `"split"` only)

`R/consistency_resample.R` read in full; call sites in `R/subgroup_consistency_helpers.R` (L1453, L1474, L1825, L1829) and the dispatch in `R/subgroup_consistency_main.R` (L388–398).

**(1) Rate formula, comparison scale, `sigma_D` as implemented.**
`consistency_resample()` L437–447:

```r
  delta <- pieces$beta_hat - c_cmp
  ...
  if (method %in% c("closed", "both")) {
    out$rate_closed <- max(0, 2 * stats::pnorm(delta / pieces$sigma_D) - 1)
  }
```

with `c_cmp` (L430–436): the caller's `comparison_threshold` if supplied, else `log(thr_nat)` for ratio measures, else `thr_nat` for identity measures (MD, RD, IRD). Every package call site passes `method = "closed"`.
`sigma_D`, GLM path (`.consistency_glm_pieces`, L307–318):

```r
  beta_hat <- unname(cf[[treat.name]])
  dfb <- tryCatch(.dfbeta_glm(fit), error = function(e) NULL)
  ...
  dfb <- as.numeric(dfb)
  sigma_D <- sqrt(sum(dfb^2))
```

and identically on the Cox path (L79–86, `residuals(fit, type = "dfbeta")`). `.dfbeta_glm()` (L168–176) is the exact leave-one-out influence `A x_i e_i / (1 - h_i)` with `A = (X'WX)^{-1}` — so `sigma_D` is the **HC-type sandwich standard error** of the treatment coefficient of the subgroup fit (`lm(Y ~ Treat [+ adj])` on the subgroup for MD). The file header says the same (L15–16: "the split-perturbation SD is sigma_D = sqrt(sum_i dfbeta[i, treat]^2), the robust (sandwich) SE of the treatment coefficient").

**(2) Is `sigma_D = k · SE(β̂_g)` for a constant or simply-derivable `k`?** **No.** `sigma_D` is a function of the realized dfbeta vector: the per-subject residuals `e_i`, leverages `h_i` and the realized design in the subgroup. The document's `se_g` (`worked-scenario` L1774–1783; A3 S1(b)) is `seQ1000 * sqrt(2) * sqrt(piQ / Pg)`, a **scale-table** quantity built from `fs_dgm_scale(dgm)`'s `V_eff[Q]` and `P_g[Q]`. The population limit of `sigma_D` on a candidate `g` would be `sqrt(Var(Y | T=1, g) / (n P_g p) + Var(Y | T=0, g) / (n P_g (1-p)))`, whose treated-arm variance in a *mixed* candidate carries a purity term `bint² · PQg(1 − PQg)` on top of the residual variance and the covariate-driven variance — a quantity `fs_dgm_scale()` does not tabulate per candidate (it reports `V_eff` for `Q`, `Qc`, `S`), and whose ratio to `se_g` is neither constant nor a simple factor of `(n, Pg)`. Writing `sigma_D ≈ se_g` (asymptotic equivalence under equal arm variances) or deriving a per-candidate arm-variance formula would each be **an approximation chosen by me**, i.e. a method decision. Per the task's gate, none was made.

**(3) `fs.splits` on the `"resample"` path.** Not used. `evaluate_subgroup_consistency()` L1450–1466: `rr <- consistency_resample(df.x, ..., method = "closed", ...)`; `p.consistency <- round(rr$rate_closed, pconsistency.digits)`. `n.splits` enters only the `.consistency_via_splits()` fallback taken when `rate_closed` is `NA` (non-convergent fit).

**(4) RNG on the `"resample"` path.** None. `method = "closed"` never reaches the `set.seed`/`rbinom`/`rnorm` block (L448–458, `method %in% c("mc", "both")` only). Two-stage screening is bypassed entirely: `evaluate_consistency_twostage()` L1810–1840, comment verbatim — "the two-stage split machinery exists only to limit the number of refits via early stopping; the resampling approximation returns the rate from a single fit, so Stage 1/2 are bypassed entirely" — with `cox_resample || glm_resample` short-circuiting before any split is drawn; `n.splits.screen`/`twostage_args` are consumed only on the fall-through.

**Gate decision:** (2) yields no exact population expression → **`"split"` gate implemented; `consistency_method = "resample"` left unimplemented with an informative `stop()`** (`R/fs_oc_predict.R`, the first thing the function does after `match.arg`). The error text states the formula, the definition of `sigma_D`, the file, and why a stand-in is not a port.

**What `"resample"` would need (for Larry):** a per-candidate population value for `sigma_D(g; n)`. Candidates: (a) an extension of `fs_dgm_scale()` that returns the arm-specific conditional-outcome variances on every candidate `g` (`V_arm0(g)`, `V_arm1(g)` are already computed per region inside `fs_dgm_scale()` L155–160 — the API just needs the family's memberships passed as `regions`, which `fs_oc_family` carries in `memb`), giving `sigma_D²(g) ≈ (V_arm1(g)/p + V_arm0(g)/(1-p)) / (n P_g)`; or (b) the cruder identification `sigma_D ≡ se_g`. Then the gate would be `rate_g = 2Φ((Bhat_g − c2)/sigma_D(g)) − 1 ≥ pconsistency.threshold`, i.e. `Bhat_g ≥ c2 + sigma_D(g) · Φ^{-1}((1 + 0.90)/2)`, combined with `Bhat_g ≥ c1`. Which of (a)/(b), and whether the finite-sample `1/(1−h)` inflation is modelled, is the method decision.

---

## 2. §3 — `R/fs_oc_family.R`: `fs_oc_family_enumerate(dgm, forestsearch_args, n, max_M = 2000L, verbose = FALSE)`

Signature as suggested (no renames needed; the arguments are read from `forestsearch_args` under their `forestsearch()` names).

**Established from source, not guessed:**
- Harm indicator: `dgm$df_super$flag_harm` (`R/generate_glm_dgm.R` L313, L435 `in_Q <- df_super$flag_harm == 1`; `fs_dgm_scale()` default `harm_col = "flag_harm"`).
- `n.min` rule: `forestsearch_main.R` L1295–1325 — the supplied value (formal default 60); `n.min = NULL` opts into `max(60L, ceiling(n.min.frac * N))`. Implemented verbatim with `N = n`. The size floor is applied as the task specifies, `Pg >= n.min / n` (the document's `Pg >= 0.12`). Note for the record: the search's sample rule is strict, `nx <= n.min → drop` (`subgroup_search.R` L592), i.e. `nx > n.min`; at the population level the two coincide except on the boundary.
- `minp` (formal default 0.025), applied per constituent factor as `meets_prevalence_threshold()` does (`all(colMeans(x) >= minp)`).
- `rmin` is **not** a `forestsearch()` formal: the search receives `subgroup.search()`'s default `rmin = 5` (`forestsearch_main.R` L2922–2925 comment: "rmin is not a forestsearch formal ... The default (subgroup.search rmin = 5) is preserved for every other focus"), relaxed to 0 together with `minp → 0` under `sg_focus = "maxeff"` only. Both mirrored. The rule itself (`extract_idx_flagredundancy()`, L500–514: each added factor must shrink the membership by more than `rmin` subjects) is a count in subjects of the analysis sample; on the population frame it is applied as a proportion, `rmin / n`, which is its meaning at trial size `n`. This is the one place a sample rule needed a population reading; it is stated in the roxygen.
- Status 0 of the search (`has_positive_variance()`: every entry of `x'x` positive) drops combinations with an empty pairwise intersection — complement pairs. Kept, counted as `empty`.
- Cuts: `get_FSdata()` on `df_super` with `use_lasso = FALSE, use_grf = FALSE, dina_cuts = NULL` and the eleven cut-related entries; `get_FSdata()` validates an outcome and an event column, so two constant numeric columns `.fs_oc_y`, `.fs_oc_event` are appended for the contract (they feed nothing with LASSO off). Both directions via `dummy()` (`acm.disjctif`, `qK.0`/`qK.1`). Combinations via `generate_combination_indices()` / `get_covs_in()` / `get_subgroup_membership()`, the composition of `forestsearch_main.R` L3264–3277.
- Scale: `fs_dgm_scale(dgm)`; `piQ = P_g[Q]`, `tauQc = |m_tau[Qc]|`, `bint = |m_tau[Q]| − tauQc`, `seQ1000 = sqrt(V_eff[Q] / (1000 · piQ))`, `se_g = seQ1000 · sqrt(1000/n) · sqrt(piQ/Pg)` (at `n = 500`, `sqrt(1000/500)` is exactly `sqrt(2)`). The document's common-sign assertion on `m_tau` is enforced.

**Return:** the nine fields (`lab, Pg, PQg, beta_g, se_g, sens_g, spec_g, ovl, M`) plus `PQ` (needed by `Enpv`; carried explicitly rather than back-derived so the port stays bit-exact), `memb` (population membership matrix, for tests and for any later `regions` extension), `scale`, `n`, `args_used` (with the resolved `n.min`, `minp`, `rmin`, `maxk`), `cuts`, `counts`. Class `fs_oc_family`; `print()` reports M, prevalence range, floors, drops.

**Size guard:** after the floors and duplicate collapse, before any M×M allocation: `stop()` with M, `max_M`, the MB per M×M matrix and the O(M³) flop count.

**Realized M — the MD40 fixture** (`build_md40_dgm.R`'s construction, rebuilt in-session: gate `effect_Q`/`prevalence_Q` diff 0; `fs_dgm_scale` regions identical to the tracked n=500 payload's `$scale$regions`), under the driver's arguments (`confounders.name` = the 12 analysis variables, `conf.cont_jcuts = list(age = 10, preanti = 10)`, `n.min = 60`, `maxk = 2`, `sg_focus = "maxeffCons"`), `n = 500`:

| stage | count |
|---|---|
| cut expressions from `get_FSdata()` | 36 (10 age, 7 preanti after `collapse_cuts` merged the repeated `preanti <= 0`, 4 wtkg/karnof-2/cd40-4/cd80-4, 6 binary) |
| cut columns (both directions) | 72 |
| combinations enumerated (≤ 2) | 2628 |
| dropped — empty intersection (status 0) | 119 |
| dropped — `minp` | 0 |
| dropped — `rmin` (5/500) | 99 |
| dropped — size `Pg < 60/500` | 631 |
| kept | 1779 |
| collapsed — identical membership | 178 |
| **M** | **1601** |

Enumeration time 8.6 s. Prevalence range 0.120–0.920; purity range 0.000–1.000; `P(Q) = 0.3446`. The population J-grid differs from the driver's *sample* grid as expected: on `df_super` the 8/11 quantile of `preanti` is 672 (the driver documents 677.9 from a sample), and the `age` grid contains 33. This M is 100× the document's 16 — the `fs_sym_root()` on a 1601×1601 matrix is what `max_M` exists for; the default `max_M = 2000` admits it.

On the small synthetic test DGM (3 covariates, J = 4): 18 columns → 171 → M = 104.

---

## 3. §4 — `R/fs_oc_predict.R`: `fs_oc_predict(dgm = NULL, forestsearch_args = list(), n, c1 = NULL, c2 = NULL, family = NULL, consistency_method = c("resample", "split"), draws = 2e5, seed = NULL, ...)`

- `c1`/`c2` default to `forestsearch_args$effect.threshold` / `$consistency.threshold`; explicit values override; missing both is an error naming what to pass. `n` is a required explicit argument.
- `family = NULL` → `fs_oc_family_enumerate(dgm, forestsearch_args, n, ...)`; a supplied `fs_oc_family` is used as-is after a structural check (`.fs_oc_check_family`). No other constructor.
- `consistency_method`: `"resample"` first, as the package's formal; it `stop()`s (§1). `"split"`: `(W1 + W2 >= 2*c1) & (W1 >= c2) & (W2 >= c2)`, verbatim.
- Body = `worked-predictions` L1826–1856 / `predict_scn()` L1988–2002 with the family construction removed, operations in the same order: `Rho <- ovl / sqrt(outer(Pg, Pg)); Sg <- Rho * outer(se_g, se_g)` → `set.seed(seed)` → `fs_sym_root(Sg, scale = 2)` → `W1`, `W2` (each `matrix(rnorm(R*M), R, M) %*% t(root_half) + rep(beta_g, each = R)`) → `Bhat` → `pass` → `det` → `Bmask`/`max.col(ties = "first")` → `P1`, `p_sel`, `det_rate`, `sel_c` → `EnH = n * sum(sel_c * Pg)`, `Esens`, `Espec`, `Eppv`, `Enpv`, `EbetaH`, `noise`/`Enaive_bias`, `mass_below`.
- Return: the nine quantities, MC SEs (`sqrt(p(1−p)/draws)`, the document's `.mc_se_prop`) for `det_rate`, `P1`, `p_sel`; `sel_c`; `M`, `lab`; `settings` (`n, c1, c2, consistency_method, draws`); `seed`; the family. Class `fs_oc_predict`; `print()`.

**Not implemented, and what each would need:**
- *Grid:* a wrapper looping `fs_oc_predict()` over `(c1, c2, n)` with a shared family per `n` and common random numbers (the `predict_scn` pattern: re-seed per cell). Cost is one M×M root per `n`, reusable across `(c1, c2)`.
- *Inversion:* `uniroot` on `det_rate(c1)` (or on any returned functional) at fixed draws — monotone in `c1` under common random numbers (test 5 checks the ladder), so a bracketed root is well-posed; the c2 ceiling of the document (`pnorm((m_Q − c2)/(seQ1000·√2))²`) bounds attainable targets.
- *Null-DGM path:* `flag_harm` is all-zero under `model = "null"` and `fs_dgm_scale()` has no `Q` row; the document reuses the alternative's covariance root with all means collapsed to `tauQc` (A3 S10). The enumerator currently refuses an empty `Q` (`stop`). A null path needs a decision on where the covariance scale comes from without `Q`.
- *Binary/OR path:* `fs_dgm_scale()` rejects ratio measures (needs a delta-method layer); the gate would also move to the log scale (`c = log(threshold)`). Nothing in the wrapper assumes MD except through the scale table.

---

## 4. §5 — fidelity harness (GATE passed, bit-identical)

`dev/glm-continuous-sims/fidelity_fs_oc_predict_2026-08-28.R`: extracts the `anchor`, `worked-scenario` and `worked-predictions` chunks **by label from the .qmd** and evaluates them (in the document's directory, so the payload read is the document's own), wraps the chunk's nine fields as an `fs_oc_family`, calls `fs_oc_predict(family = ., n = 500, c1 = 30, c2 = 10, consistency_method = "split", draws = 2e5, seed = 20260825)`, and asserts `identical()` per quantity. Output verbatim:

```
Document chunks evaluated: M = 16 ; Rdraw = 2e+05
quantity     identical                  document             fs_oc_predict
det_rate     TRUE            0.96355999999999997       0.96355999999999997
EnH          TRUE             132.84322367055503        132.84322367055503
Esens        TRUE            0.32855277357396867       0.32855277357396867
Espec        TRUE            0.76736777301263281       0.76736777301263281
Eppv         TRUE             0.3592225203795214        0.3592225203795214
Enpv         TRUE            0.69495433543964846       0.69495433543964846
EbetaH       TRUE             31.192664686668142        31.192664686668142
Enaive_bias  TRUE             30.935708656675445        30.935708656675445
mass_below   TRUE            0.46759931919133213       0.46759931919133213
P1           TRUE      [16] 0.58729500000000001 ... [16] 0.58729500000000001 ...
p_sel        TRUE      [16] 0.00059000000000000003 ... [16] 0.00059000000000000003 ...
sel_c        TRUE      [16] 0.00061231267383453033 ... [16] 0.00061231267383453033 ...
det_rate_se  TRUE         0.00041899956085895859    0.00041899956085895859

identical() on the nine-quantity numeric vector: TRUE
FIDELITY GATE: PASS (bit-identical)
```

RNG consumption order is the chunk's (one `set.seed`, then `W1` in full, then `W2`); no reordering was needed.

---

## 5. §6 — tests: `tests/testthat/test-fs-oc-predict.R`

Fixture: a 2000-row synthetic `glm_dgm` (age, preanti continuous; V binary factor; `flag_harm = age > 34 & preanti <= 745`), plus the hand-built M = 3 family. Draw counts 4000–5000.

1. Interface contract — all nine quantities finite, `det_rate ∈ [0,1]`, `sum(sel_c) == 1`, SE form, print.
2. Family construction — `Pg == colMeans(memb)`, `ovl` symmetric with `diag == Pg`, `PQg`/`sens`/`spec` in [0,1] and equal to direct recomputation from `memb`, every `Pg >= 60/500`, no duplicate memberships, `PQ == mean(flag_harm)`; a second test pins `beta_g`/`se_g` to `fs_dgm_scale()` algebra.
3. `conf.cont_jcuts` — J = 10 gives 10 age cuts, J = 4 gives 4; cut values equal `round(quantile(df_super$age, k/(J+1)), 1)`.
4. Determinism — identical families; identical predictions at the same seed; different at a different seed.
5. Override semantics — explicit `c1`/`c2` override the args (recorded in `settings`, change the rate); `det_rate` monotone non-increasing on the ladder `c1 ∈ {20,30,40,50,60}` under common draws; `EnH` scales with `n`; missing thresholds error.
6. Size guard — `max_M = 5` below the realized M errors from both entry points.
7. Gate parity — **not applicable** (§1 gate closed); replaced by a test that `"resample"` (default and explicit) stops with the `sigma_D` explanation.
Plus: `n.min = NULL` resolves to `max(60, ceiling(0.2·500)) = 100` and tightens the family; a malformed family is rejected.

New file alone: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 63 ]`.

**v5 §9 — a test shown to fail against its defect.** With the line `if (!is.null(seed)) set.seed(seed)` replaced by a comment (defect: seed ignored) the suite gives `[ FAIL 3 | ... | PASS 60 ]`: test 4's two `expect_identical` on the two same-seed predictions (L149, L150) and test 5's override comparison that relies on common draws (L180). Restored; `grep` confirms the line is back and the suite passes. Checked.

---

## 6. §7 — close-out

`devtools::document()` — wrote `NAMESPACE` (+`export(fs_oc_family_enumerate)`, `+export(fs_oc_predict)`, `+S3method(print, fs_oc_family)`, `+S3method(print, fs_oc_predict)` — nothing else), `man/fs_oc_family_enumerate.Rd`, `man/fs_oc_predict.Rd`.

`devtools::test()` after the edits:

```
[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4577 ]
```

Parity against the §1 baseline `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4514 ]`: WARN 31 = 31 (parity), SKIP 3 = 3, FAIL 0 = 0, PASS 4514 + 63 = 4577 — the 63 new expectations and nothing else.

`R CMD check` (`devtools::check(document = FALSE, args = "--no-manual")`; a first attempt failed at vignette build because `pandoc` is not on `Rscript`'s path in this shell — an environment matter, fixed by `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64`):

```
── R CMD check results ───────────────────────────────── forestsearch 0.2.3 ────
Duration: 10m 36.6s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

Version `0.2.2 → 0.2.3` (`DESCRIPTION`). `NEWS.md` entry names both exports, states the family is a population enumeration under the `forestsearch()` cut specification, and that `consistency_method` mirrors the package's (with `"resample"` refused and why).

Commit (explicit paths only: `R/fs_oc_family.R`, `R/fs_oc_predict.R`, `tests/testthat/test-fs-oc-predict.R`, `man/fs_oc_family_enumerate.Rd`, `man/fs_oc_predict.Rd`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`, the fidelity script, this report): **`__COMMIT__`**. `git status --porcelain` after the commit: __DIRT__. **No push.**

`devtools::install()`: `packageVersion("forestsearch")` → __INSTALLED__.

No existing `R/` file was edited (`git diff --stat HEAD~2 -- R/` shows only the two additions).

---

## 7. Verdict (ten lines)

1. `fs_oc_family_enumerate()` and `fs_oc_predict()` are added; no existing `R/` file, function or method changed.
2. The §2 gate closed: `sigma_D` is the subgroup fit's sandwich SE (`sqrt(sum(dfbeta²))`), a sample statistic with no scale-table expression; no stand-in was chosen.
3. Hence `consistency_method = "split"` is implemented and `"resample"` stops with the reason; §6 test 7 is replaced by that refusal.
4. The §5 fidelity gate passed bit-identically on every quantity, vectors included, against the document's own chunks at seed 20260825, 2×10⁵ draws.
5. The family is a true population enumeration: cuts at `df_super` quantiles via `get_FSdata()`, both directions, ≤ `maxk` combinations, the three structural floors, identical-membership collapse.
6. On the MD40 fixture under the driver's arguments the realized family is M = 1601 (from 2628 combinations), two orders above the document's 16 — the size guard is not decorative.
7. `rmin` (not a `forestsearch()` formal; `subgroup.search` default 5) is applied as `rmin/n` in population proportion — the one sample rule that needed a population reading, stated in the roxygen.
8. Tests: 63 new expectations, all pass; the seed test demonstrably fails when the seed handling is removed.
9. Test parity and check status are recorded above verbatim; version 0.2.3; NEWS written; one commit, explicit paths, no push.
10. Next decisions are Larry's: the `sigma_D` population rule for `"resample"`, then grid / inversion / null / binary paths on this seam.
