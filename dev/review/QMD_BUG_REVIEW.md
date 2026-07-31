# QMD bug review — findings

Targets (maintainer copies; the `dev/review/` copies are byte-identical, verified
by `diff`):

* `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` — 2507 lines,
  69 chunks. Unmodified vs `HEAD` (last touched in `9f53a25`).
* `quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd`
  — 1329 lines, 11 chunks. Modified in the working tree (`n_sims` 100→50,
  `nb_boots` 300→20, `mr_draws` 5000→500, plus the `sg_focus` move, the
  `focus_tag` fold, and the provenance block).

Nothing was edited. Everything below is report-only.

---

## 1. Verdict

**Neither file will error on render.** No unresolved symbol, no `sprintf` arity
mismatch, no dead field path, no argument that does not match a real formal.
Both setup chunks were executed in a clean session; every extraction path in
both documents was resolved by introspection on a real `forestsearch` fit and a
real bootstrap object built for this review.

Scope limit, stated plainly: the four heavy fits and three `NB = 1000`
bootstraps in the multimethod file were **not** run end-to-end (hours of
compute). What was executed is listed per check in §3. The probe fit used the
same engine, `sg_focus`, and `mr_inference` settings as the main analysis, so
the object *shapes* are the real ones; the *values* are not.

**Nothing found changes a computed number.** Every finding is either a wrong
number *displayed in prose*, a wrong label, or neither.

### Would error — none.

### Would silently produce a wrong displayed value

| # | File | Line(s) | What |
|---|---|---|---|
| **F1** | sim | 98 / 172 | Batch writes `..._res_1_50.rds`; the filename, title and every "1_100" reference say otherwise. Breaks a later `combine`. |
| **F2** | multimethod | 612, 617, 1177, 1351 | Four hard-coded parameter values in prose contradict the pinned code values. |
| **F3** | multimethod | 753, 1286 | `stop_threshold = 0.95` is silently reset to `NULL` by the package and raises a render-visible warning; the adjacent comment is wrong. |
| **F4** | multimethod | 1058, 1572, 1583, 1646, 1801, 1811, 1878 | Seven surviving "tier"/"tiers" — the Tier→FB/MR sweep missed them. |
| **F5** | multimethod | 2220–2237 | The Computational Timing table omits the (D) GRF fit and its bootstrap, but "Total" includes them. |
| **F6** | multimethod | 1124–1156 | "GRF-Cut Stability Across Bootstraps" renders a single `(no GRF cut) — 100%` row because the main fit sets `use_grf = FALSE`. |
| **F7** | multimethod | 1940–2022 | ~60 lines of CV narrative (incl. an inline "50 × 10 = 500") render in full while every CV chunk is skipped. |

### Neither

| # | File | Line(s) | What |
|---|---|---|---|
| N1 | multimethod | 1524, 1755 | `mr_in_replicates` not pinned on the DINA/GRF bootstraps (§2 — confirmed no-op). |
| N2 | sim | 515–526 | Twelve unguarded `rec$… <- g$…` assignments; provably safe against the current package, but asymmetric with the guarded FB block. |
| N3 | sim | 566–567 | The comment justifying the FB guard mis-states the failure mode ("silently dropped" — it actually errors). |
| N4 | sim | 399–402 | `mr_ok`, `n_sel`, `label`, `covs` are recorded and persisted but never read. |
| N5 | multimethod | 119–121, 1520–1522 | Two comments describe behaviour the code does not have. |

---

## 2. `mr_in_replicates` — quantified

**Call counts across both files.**

| Function | Calls | Carry `mr_in_replicates = FALSE` | Do not |
|---|---|---|---|
| `forestsearch_bootstrap_dofuture()` | 4 | multimethod L914 (main), sim L548 | multimethod **L1524** (DINA), **L1755** (GRF) |
| `forestsearch_tenfold()` | 1 | multimethod L1916 | — |
| `forestsearch_Kfold()` | 1 | multimethod L2038 | — |
| **Total** | **6** | **4** | **2** |

The brief said "three call sites" in the multimethod; that is right for that
file (L914, L1916, L2038). The sim file has a fourth, at L548.

**The no-op characterisation is correct.** Verified, not read off comments:

1. `formals()` on all three functions returns `mr_in_replicates = FALSE`
   (`bootstrap_dofuture_main.R:193`, `forestsearch_cross_validation.R:295`,
   `:784`). Omitting the argument and passing `FALSE` are identical.
2. Each entry point strips the flag from the replicate/fold call regardless of
   what the fit carries:
   * `bootstrap_analysis_dofuture.R:597` — `if (!isTRUE(mr_in_replicates)) args_FS_boot$mr_inference <- FALSE`
   * `forestsearch_cross_validation.R:402` and `:866` — `cv_args$mr_inference <- FALSE`

   Note it assigns `FALSE`, not `NULL`, so this is not itself a `[[<-`-deletes
   hazard.
3. On the real bootstrap object built for this review (fit carried
   `mr_inference = TRUE`; bootstrap run with `mr_in_replicates = FALSE`),
   `fs_bc$mr_replicates` is `NULL`.
4. Neither document reads `$mr_replicates` anywhere. The only occurrence in
   either file is sim L556, inside a comment.

So the two unpinned DINA/GRF bootstraps are a **documentation-completeness gap
only**: the three pinned sites carry a five-line comment explaining the pin, and
a reader comparing the DINA/GRF blocks against the main block sees an
unexplained asymmetry. No estimate changes.

---

## 3. Per-check results

### 3.1 Symbol resolution in definition order

**Sim — PASS (executed).** The setup chunk was extracted and evaluated in a
clean `--vanilla` session with `library()` and `source()` stubbed. Stubbed:
`forestsearch`, `survival`, `data.table`, `gt`, `ggplot2`, `foreach`,
`doFuture`, `future`, and `source("betaHhat_truth.R")`. It ran to completion,
no `object … not found`.

`focus_tag` (L161) now resolves **before** `rds_stem` consumes it (L167) — the
reordering is correct.

Resolved knob values every downstream chunk consumes:

```
parallel_mode  "boots"      n_sims        50L        nb_boots      20L
mr_draws       500L         diag_mode     FALSE      k_random_noise 0
run_mode       "batch"      sim_id_start  1L         sim_id_end    50L
subgroup_method "consistency"  sg_focus   "maxcons"  focus_tag     "maxcons"
method_tag     "fs"         target_hr_harm 1.0       n_sample      500L
rds_stem       "fs_maxcons_fb_mr_m1_h10_knoise0_n500"
rds_path       "fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_50.rds"
combine_glob   "fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_*.rds"
run_fb         TRUE         run_smoke     TRUE       save_rds      TRUE
consistency_method "resample"   selection_rule "neighborhood"
stop_threshold NULL         effect_neighborhood 0.10
hr_threshold   0.90         hr_consistency 0.80      pconsistency  0.90
fs_splits      400L         maxk          2L         n_min         NULL
d0_min/d1_min  10L/10L      use_lasso/grf/dina FALSE  use_twostage TRUE
seed_base      8316951L     n_workers     115        z975          1.95996
dgm_model      "alt"        analysis_time 84         n_super       100000L
```

**Both files — PASS (executed).** A forward-reference detector walked every
chunk in document order (plus every inline `` `r ` `` expression and every
chunk-option expression), accumulating assigned names and diffing reads against
`codetools::findGlobals` and the loaded package namespaces.

* Sim: 6 flags, all false positives — `build_eval_frame`,
  `betaHhat_theta_dagger_check`, `rule_covs`, `attach_betaHhat` (all four
  confirmed defined in `betaHhat_truth.R`, which exists in the render
  directory); `s` (the `%dofuture%` loop variable); gt/ggplot NSE column names.
* Multimethod: 18 flags, all false positives — `$`-field names harvested by
  `all.vars` from the `eval =` guards (`sg.harm`, `mr_inference`), `within()`
  columns, and gt NSE column names.

**No genuine forward reference in either file.**

### 3.2 Rename completeness

**`t1`/`t2` as substring — PASS.** Full-file token scan
(`[A-Za-z_.][A-Za-z0-9_.]*t[12][A-Za-z0-9_.]*`):

* multimethod: one hit, `ggplot2` (L303).
* sim: `ggplot2` (L76) and `_t1_t2_m1_` at **L141**. That one is inside the
  deliberate back-compatibility note — *"the stem changed from
  `<engine>_t1_t2_m1_...` to `<engine>_<focus>_fb_mr_m1_...`"* — so it is a
  citation of the old name, not a survival. Correct as written. (Confirmed:
  the 30-odd `fs_t1_t2_*.rds` / `dina_t1_t2_*.rds` files in the sim directory
  do not match the new `combine_glob`, exactly as that note claims.)

No `t2m`-class survivals.

**`fb_*` / `mr_*` assigned but never read — 4 found (sim).** Accounting for the
dynamic `r[[paste0(k, "_", suffix, "_est")]]` construction and `.get_se()`, the
`.na_record()` columns split as:

* read: `sim_id`, `detected`, `n_harm`, `n_true`, `betaHhat_H/Hc`, `fb_secs`,
  `fit_mr_secs`, `fb_err`, `mr_harm_flag`, `ij_source`, `sens/spec/ppv/npv`,
  and all 24 `{or,nv,fb,mr}_{H,Hc}_{est,lo,hi}` + 8 SE columns.
* `sg_def` — read, but by `attach_betaHhat()` in `betaHhat_truth.R`
  (`results$sg_def`), not in the `.qmd`.
* **never read: `mr_ok` (L399), `n_sel` (L399), `label` (L401), `covs`
  (L401)** — assigned at L508, L512, L511, L505.

Consequence (N4): benign — they travel into the saved bundle. Two are worth a
look anyway. `covs`'s comment says *"variable names -> id figures"*, but this
document has no identification-figure chunk (11 chunks, none). And `mr_ok`
exists precisely so *"an MR failure never masks a true detection"* (L498–L499),
yet no table reports it: an MR failure is currently invisible in the output.
`ij_source` is reported; `mr_ok` is not.

**`fb_*` / `mr_*` read but never assigned — none.** No silent `NA` column from
this cause.

**`gate` / `Tier` in any identifier, comment, heading, `cat()` string or table
label — 7 survivals, all in the multimethod (F4).** Word-boundary scan for
`\btiers?\b`:

| Line | Text |
|---|---|
| 1058 | "descriptive columns and naive effect are **tier**-independent" |
| 1572 | "For DINA the two **tiers** are **not** expected to agree…" |
| 1583 | "so it agrees across **tiers** regardless" |
| 1646 | "descriptive columns and naive effect are **tier**-independent" |
| 1801 | "As with DINA, the two **tiers** are **not** expected to agree…" |
| 1811 | "so it agrees across **tiers** regardless" |
| 1878 | "descriptive columns and naive effect are **tier**-independent" |

All seven are narrative prose, all present in the committed version, so the
sweep missed them rather than the working tree reintroducing them. Under the
new vocabulary these read as "FB/MR-independent" and "FB and MR are not expected
to agree". They change displayed wording, not numbers.

Legitimate substring holders confirmed and correctly left alone: `frontier`
(×20 multimethod, ×3 sim), `delegate`/`delegates`/`delegated` (L1186, L1439,
L1675), `gatekeeper` (L1209), `ungated` (L631), `Negate` (sim L1266), and the
ordinary English "not to gate anything" (sim L679). No `T1`/`T2` uppercase.

### 3.3 `sprintf` / `cat` arity — PASS (executed)

An AST walker parsed every chunk and inline expression in both files, extracted
every `sprintf()`/`gettextf()` call, and compared conversion specifications
(handling `%%` as literal and `%*d` as consuming an extra argument) against the
number of arguments supplied.

| File | Calls | Literal format strings | Arity mismatches |
|---|---|---|---|
| multimethod | 47 | 47 | **0** |
| sim | 57 | 57 | **0** |

Every format string is a literal, so all 104 were statically checkable — none
had to be taken on trust. The checker was validated against a synthetic file
containing a too-few-args case, a too-many-args case, a `%%` case and a
non-literal `fmt`; it caught all of them.

`cat()` itself takes no format string, so the only arity risk is the
`cat(sprintf(...))` composition, which the above covers.

### 3.4 Field paths against real objects — PASS (executed)

A real fit and a real bootstrap were built on GBSG for this check
(`subgroup_method = "consistency"`, `sg_focus = "effMaxSG"`,
`mr_inference = TRUE`, `nb_boots = 8`; selected
`{er <= 0} & {size <= 35}`). Every path below was resolved on those objects, not
read from source.

**Multimethod — all resolve.**

* `fs$args_call_all` — 85 entries against 85 formals: **no formal missing, no
  extra**. `fs_param_provenance()` therefore never emits "(not captured)".
  The 19 `NULL`-valued formals are *recorded* as `NULL` rather than dropped
  (this is `af9cc0e` working), so `.fs_fmt_val()` renders them "NULL".
* All 36 `.fs_selection_knobs` are real formals of `forestsearch()`.
* `fs$mr_inference$…`: `naive$est/lower/upper`, `debiased$est/lower/upper/
  se_ij/se_wald/ij_source`, `complement$naive$est`, `complement$debiased$est`,
  `complement$n`, `complement$selection_bias`, `n_family`, `selection_bias`,
  `timing_seconds`, `ci_method`, `harm_flag`, `selected_label`, `n_selected`,
  `settings$reselection` (= `"effMaxSG"`), `settings$confirm_rule`,
  `settings$t_confirm` — all present and non-`NULL`.
* `fs$mr_harm_confirmed`, `fs$sg.harm`, `fs$df.est$treat.recommend`,
  `fs$find.grps$out.found$hr.subgroups` (L1016),
  `fs$grp.consistency$out_sg$result$hr` (L2189),
  `fs$grp.consistency$sg.harm.id` — all present.
* `fs_bc$H_estimates` / `$Hc_estimates` (a `data.table`): `H0`, `H0_lower`,
  `H0_upper`, `H2`, `H2_lower`, `H2_upper`, `sdH2` all present.
* Payload assembly executed verbatim against the probe objects: `P$table` built
  with all 13 columns (`method, region, n, pct, naive_{est,lo,hi},
  fb_{est,lo,hi}, mr_{est,lo,hi}`), **no all-`NA` column**;
  `meta$timings$fs` fully populated (`label`, `n_family`, `selection_bias`,
  `boot_sec`, `mr_sec`, `speedup`). `.boot_sec()` on a missing key returns
  `NA_real_` without erroring (`list[["missing"]]` is `NULL` in R, not an
  error — verified).
* CV/LOO reads verified against source rather than by running CV:
  `summarize_cv_results()` documents and returns `identification_summary`,
  `grf_cut_summary`, `no_subgroup_decomposition`, `pconsistency_distribution`,
  `fold_numeric_summary`, `original_agreement` — all six the qmd reads exist;
  `cv_summary_tables()` returns `combined_table`
  (`R/cv_summary_tables.R:100`) and `metrics_table`.
* `summarize_bootstrap_results()` returns `table`, `diagnostics_table_gt`,
  `plots$H_distribution`/`$Hc_distribution`,
  `subgroup_summary$original_agreement`/`$factor_presence`/`$grf_cut_freq` —
  every slot the qmd reads.
* `sg_tables()` returns `tab_estimates`, `sg10_out` — both read.

Two paths resolve to `NULL` on the probe fit, neither harmful:
`fs$grf_cuts` and `fs$dina_cuts` (the main fit sets `use_grf = use_dina =
FALSE`). `fs$grf_cuts` is passed to `summarize_cv_results(original_grf_cuts =)`
at L1953, inside a `run_cv = FALSE` chunk. See F6 for the visible consequence.

**Sim — all resolve.** `boot$H_estimates$H2`, `$H2_lower`, `$H2_upper`,
`$sdH2`, and the `Hc_estimates` siblings; `fs.est$mr_inference$debiased$est/
lower/upper/se_ij/ij_source`, `$naive$*`, `$complement$*`, `$selected_label`,
`$n_selected`, `$harm_flag`.

**Additional static check (not requested).** Every named argument at every call
site in both documents was matched against the real `formals()` of the function
called (allowing R's partial matching). Of the 27 package entry points used,
only `dina()` has `...`; the rest were fully checked. **Zero mismatched
argument names.** No renamed or removed argument is being passed anywhere.

### 3.5 The new, never-executed code — PASS (executed)

**(1) The provenance block.** Extracted verbatim and executed against the
resolved setup values. It builds a 20-element `meta` with no error. All seven
named symbols resolve at that point (`sg_focus`, `focus_tag`,
`consistency_method`, `stop_threshold`, `selection_rule`,
`effect_neighborhood`, `n_workers`).
`paste(deparse(stop_threshold), collapse = " ")` yields the character string
`"NULL"` — confirmed, not an error.
`utils::packageVersion("forestsearch")` is reachable: `"0.2.0"`.

**Combine-path propagation — confirmed against the real pre-provenance batch.**
`fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds` on disk has a 9-element
`meta` with none of the provenance fields. Reading it through the combine
block's `%||% NA_character_` idiom:

```
sg_focus  focus_tag  consistency_method  stop_threshold  forestsearch_version
   NA         NA             NA                NA                NA          (all character)
```

Clean `NA_character_`, no error. The custom `%||%` (sim L298) also treats
`all(is.na(x))` as missing, which does not misfire here: `stop_threshold` is
stored as the *string* `"NULL"`, and `all(is.na("NULL"))` is `FALSE`.

**(2) `focus_tag`'s new branches — every legal value exercised.** Each
`sg_focus` was pushed through the `switch()` and, independently, through the
package's own `.normalize_sg_focus()`:

| `sg_focus` | `focus_tag` | package-normalized rule |
|---|---|---|
| `hr`, `eff`, `maxcons` | `maxcons` | `hr` |
| `effMaxSG`, `hrMaxSG` | `effMaxSG` | `hrMaxSG` |
| `effMinSG`, `hrMinSG` | `effMinSG` | `hrMinSG` |
| `maxSG` | `maxSG` | `maxSG` |
| `minSG` | `minSG` | `minSG` |
| `maxeff` | `maxeff` | `maxeff` |
| `maxeffCons` | `maxeffCons` | `maxeffCons` |

The map is a bijection onto the package's normalized rules: every alias pair
folds to one tag, and **no two distinct rules share a tag**. `maxSG`/`minSG`
stay distinct from `effMaxSG`/`effMinSG` — confirmed, they fall through the
`switch()` default to themselves. No stem collision is possible.

**(3) Agreement key vector.** L714 reads exactly

```r
c("n_sample", "nb_boots", "mr_draws", "subgroup_method", "seed_base")
```

No provenance field was added. Correct — verified empirically in F1 below, where
a pre-provenance bundle passes the key check on its own.

### 3.6 Chunk options and eval guards — PASS

**Sim.** One guard: `#| eval: !expr run_fb` on `paired-panel` (L1239).
`run_fb` is defined at L198, 1041 lines earlier, and is `TRUE` under the current
config. No downstream chunk reads anything that chunk creates.

**Multimethod.** 25 guarded chunks. Every guard variable is defined earlier:
`fs` (L682) for the four `debias-MR-*` chunks; `fs_dina` (L1449) for five;
`fs_grf` (L1686) for five; `run_cv`/`run_loo` (L139–140) for eleven;
`eval = FALSE` for one.

`run_cv = FALSE`, `run_loo = FALSE` — objects those chunks would create are
`fs_kfold`, `metrics_tables`, `cv_diag`, `fs_OOB`, `cv_out`, `tables`. Every
read of them is inside a chunk carrying the same guard, with exactly one
exception: the forest plot at L2120–2121, which uses
`if (exists("fs_OOB")) fs_OOB else NULL` and the same for `fs_kfold`.
`timings$kfold` and `timings$oob` are pre-seeded to `NA_real_` at L147–148, so
the timing table is well-formed either way. **Nothing downstream breaks.**

The `FALSE` branches do leave a cosmetic problem — see F7.

**`#| message: false`** on the three bootstrap chunks (L906, L1511, L1746):
confirmed it suppresses only `message()` output. All three chunks' own reported
output goes through `cat()` (L935, L1539, L1770), which is unaffected. The
"Low events in NEW subgroup" condition the comments cite is re-surfaced properly
by `summarize_bootstrap_events()` in the un-suppressed `event-summary` chunk
(L1097). No chunk relies on a message for its output.

### 3.7 NULL semantics — one asymmetry, no live defect

An AST scan found every `$<-` / `[[<-` whose right-hand side is an extraction or
an opaque call.

**Multimethod — clean.** All 12 hits are `timings$x <- (proc.time() - t0)["elapsed"]`,
which cannot be `NULL`.

**Sim — 12 unguarded assignments (N2), at L515–517 and L519–526:**

```r
rec$nv_H_est  <- g$naive$est          rec$nv_Hc_est <- g$complement$naive$est
rec$nv_H_lo   <- g$naive$lower        rec$nv_Hc_lo  <- g$complement$naive$lower
rec$nv_H_hi   <- g$naive$upper        rec$nv_Hc_hi  <- g$complement$naive$upper
rec$mr_H_est  <- g$debiased$est       rec$mr_Hc_est <- g$complement$debiased$est
rec$mr_H_lo   <- g$debiased$lower     rec$mr_Hc_lo  <- g$complement$debiased$lower
rec$mr_H_hi   <- g$debiased$upper     rec$mr_Hc_hi  <- g$complement$debiased$upper
```

None is `%||%`-guarded, unlike the FB block 50 lines later (L568–L578), whose
comment explains exactly this hazard.

**This is not currently a defect.** `fs_mr_inference.R:432–451` constructs
`naive` and `debiased` unconditionally in the returned list, so neither can be
`NULL` when `g` is non-`NULL`. For the complement,
`fs_mr_inference.R:415–427` builds `naive` and `debiased` together or returns
`list(note = "complement subgroup could not be fit")` with no `debiased` at
all — and the guard at L518 tests `!is.null(g$complement$debiased)`, which
correctly excludes that case. `NaN`/`NA` values store fine; only `NULL` deletes.

The other flagged assignments (`rec$sens <- cl["sens"]`, `rec$or_H_est <-
oH["est"]`, …) index *named atomic vectors*, which return `NA` for a missing
name, never `NULL`. Safe.

**N3 — the comment at L566–567 mis-states the failure mode.** It says a short
record is *"silently dropped by the run-loop's rbind/.errorhandling = 'remove'
-> a hidden, non-random loss."* Verified: `rbind()` on data frames with
differing columns **errors** (`numbers of columns of arguments do not match`),
and `.errorhandling = "remove"` covers errors raised inside the loop body, not
inside `.combine`. Under the current `parallel_mode = "boots"` the loop is
`do.call(rbind, lapply(...))`, which errors outright. The guard is right; its
stated rationale is not, which could mislead someone deciding whether the MR
block needs the same treatment.

---

## 4. Findings the brief did not anticipate

### F1 — the batch writes `_res_1_50.rds`, not `_res_1_100.rds` (sim)

`n_sims <- 50L` (L98, changed from `100L` in the working tree) with
`sim_id_start <- 1L` gives `sim_id_end = 50`, so L172 resolves

```
rds_path = "fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_50.rds"
```

while the file is named `..._batch_1_100.qmd`, its rendered HTML is
`..._batch_1_100.html`, and the batch already on disk is
`..._res_1_100.rds` (written 2026-07-31 04:17 by the committed `n_sims = 100`
version). Both match `combine_glob = "..._res_*.rds"`.

Executed against the real on-disk bundle plus the batch this render would
write:

```
agreement check   -> STOP: combine mode: batches disagree on meta$nb_boots (300, 20).
duplicated sim_id -> 50   (would STOP the pool on the next check)
```

So re-rendering as-is leaves two mutually-incompatible batches in the same glob,
and the next `run_mode = "combine"` halts. It halts loudly rather than pooling
silently — the guards work — but the render gives no warning that it has just
created that state. The three per-batch knobs the header calls "the main knobs"
(`n_sims`, `nb_boots`, `mr_draws`) are all count-invariant meta, and two of them
are in the agreement key.

This is not a code defect; it is the working-tree downscale interacting with a
directory of prior output. Worth resolving before the next render.

### F2 — four hard-coded parameter values in prose contradict the code

Verified against the resolved setup values:

| Line | Prose claims | Code actually uses |
|---|---|---|
| 612 | `hr.threshold` = **1.25** | `fs_hr_threshold <- 1.0` (L158) |
| 617 | `d0.min`, `d1.min` = **12** | `d0.min = 10`, `d1.min = 10` (L738–739) |
| 1177 | DINA sections mirror "harm threshold `hr.threshold = 1.25`" | `dina_screen_hr_threshold <- 0.9` (L220), `fs_dina_hr_threshold <- 0.90` (L231) |
| 1351 | standalone DINA harm floor "`m_diff = log(1.25)`" | `m_diff_dina <- log(fs_dina_hr_threshold)` = `log(0.90)` (L1198) |

The rest of that table is correct (`hr.consistency` 1.0, `pconsistency.threshold`
0.90, `maxk` 2, `n.min` 60), and the claim at L1213 that (A) and (C) select the
same cut holds (`0.9` and `0.90` are the same value).

Two nearby prose numbers I checked and found **correct**: the ITT Cox HR at L547
and L2182 ("0.69 (0.54, 0.89)" — recomputed as 0.695 (0.544, 0.888)) and the
N/censoring at L475–477 (N = 686, 43.6% events → ~56% censoring).

The per-fit provenance tables (which read `args_call_all`) are unaffected; only
the hand-written prose drifted. But a reader hits L612 nine hundred lines before
the provenance table that contradicts it.

### F3 — `stop_threshold = 0.95` is inert, warns at render, and the comment is wrong

L753 and L1286 both pass `stop_threshold = 0.95` with the comment
`# default = pconsistency.threshold (0.90 here), NOT 0.95`.

Both fits use an `effMaxSG` focus. Observed on the probe fit, which reproduces
this exactly:

```
Warning: stop_threshold = 0.95 reset to NULL for sg_focus = 'hrMaxSG'.
  Neighborhood-based selection requires evaluating all candidates …
  To suppress this warning, pass stop_threshold = NULL explicitly.
```

and `fs$args_call_all$stop_threshold` comes back as `NULL`, not `0.95`.

Three consequences: the YAML sets `warning: true`, so this warning renders
twice in the output; the visible code says `0.95` while the provenance table
says `NULL`; and the comment's claim about the default is wrong on both counts
(the resolved value is `NULL`, and it is not `pconsistency.threshold`).

The file already documents an analogous "captured before the override" problem
for `sg_focus = "maxeff"` at L390–399, and that caveat correctly does **not**
fire here (it is keyed on `maxeff`). This is a *different* override that the
package now re-syncs properly — so the provenance table is right and the
hand-written comment is the only thing wrong.

A second render-visible warning from the same fits, this one substantive rather
than cosmetic:

```
Warning: max_subgroups_search = 50 truncated the candidate pool from 123 to 50.
  Candidates ranked below 50 under the sg_focus = 'hrMaxSG' preview ordering
  were not evaluated.
```

`fs_max_subgroups_search <- 50` (L160) is a deliberate setting, so this is
informational — but it does mean the main analysis searches a truncated pool,
and nothing in the surrounding prose says so.

### F5 — the Computational Timing table silently omits the GRF track

L2220–2237 lists ten components: GRF (preliminary), ForestSearch, Bootstrap,
FS (DINA screening), DINA (standalone), FS (subgroup_method=dina),
FS DINA bootstrap, kfold, oob, Total.

Missing: `timings$fs_grf_select` (set at L1722) and `timings$fs_grf_bootstrap`
(set at L1768 or L1776). Section (D) is a full fit plus an `NB = 1000`
bootstrap — plausibly the second-largest cost in the document — and its time is
inside "Total" but attributed to no row. The label/value vectors are both length
10, so nothing errors; the table just under-accounts.

### F6 — the GRF-cut stability section is empty by construction

The main fit sets `use_grf = FALSE` and `use_lasso = FALSE` (L698–699), so
candidate cuts come only from `conf_force` + `conf.cont_jcuts` + default cuts.
Confirmed on the probe fit built with the same settings:

```
   Rank      grf_cut     N Percent
1:    1 (no GRF cut)     8     100
```

So the "GRF-Cut Stability Across Bootstraps" table (L1138–1155) renders as a
single 100% "(no GRF cut)" row, under 13 lines of narrative (L1126–1136)
explaining how to read a dominant GRF row versus a frequent no-cut row — whose
own diagnostic advice ("relax `dmin.grf` / `frac.tau`") does not apply, because
GRF was never asked to contribute. Same for the "Methodology Overview" bullet at
L282 ("Using LASSO and/or Generalized Random Forests"). The numbers are right;
the framing implies a mechanism that is switched off.

### F7 — the CV narrative renders in full while every CV chunk is skipped

With `run_cv = FALSE` the `cv-skipped-note` chunk prints one line, and then
L1940–2022 renders ~60 lines of prose ("**Which subgroups emerged across
resamples?**", "**When a fold failed to identify a subgroup, why?**", …) each
introducing a table that does not appear. L1958–1960 includes an inline
expression that always evaluates:

```
`r paste0(Ksims, " × 10 = ", Ksims * 10)`   ->   "50 × 10 = 500"
```

so the rendered document asserts 500 sim × fold pairs for a run that did none.
Same shape for the LOO section at L2025–2032, though it is much shorter.

### N5 — two comments describe behaviour the code does not have

* **L119–121.** `# Detect available cores (limited to 2 cores for CRAN checks)`
  sits above `n_cores <- max(1L, floor(0.80 * detectCores(logical = FALSE)))`,
  which does no such limiting. (On this machine that resolves to 102.) The
  package's `.default_parallel_workers()` does the CRAN capping; this line does
  not call it.
* **L1520–1522.** `# Diagnostic run: sequential plan + details = TRUE …
  Switch plan back to "multisession" for production speed` sits directly above
  a call that already passes `plan = "multisession"`. Half-applied edit.

Neither affects output.

---

## Appendix — what was executed

| Check | Method |
|---|---|
| Setup-chunk resolution | Chunk extracted, evaluated in `R --vanilla` with `library()`/`source()` stubbed |
| Forward references | AST walk in document order; `codetools::findGlobals` vs accumulated assignments + loaded namespaces; chunk options and inline `` `r ` `` included |
| `sprintf` arity | AST walk; conversion-spec counter handling `%%` and `%*d`; validated against a synthetic file with known mismatches |
| Field paths | Real GBSG `forestsearch()` fit (`effMaxSG`, `mr_inference = TRUE`) + real `forestsearch_bootstrap_dofuture()` object (`nb_boots = 8`); every path evaluated |
| Payload assembly | Export-chunk helpers run verbatim against those objects |
| Combine path | Real on-disk `..._res_1_100.rds` + the batch this render would write |
| `focus_tag` | All 11 legal `sg_focus` values through the `switch()` and `.normalize_sg_focus()` |
| `mr_in_replicates` | `formals()` on all three entry points; strip sites read; `$mr_replicates` checked on the real bootstrap object |
| Argument names | Every named argument at every call site matched against real `formals()` |
| NULL assignment | AST scan for `$<-`/`[[<-` with extraction or opaque-call RHS |
| R semantics | `list[["missing"]]` → `NULL`; `paste(deparse(NULL), collapse=" ")` → `"NULL"`; `rbind` on ragged data frames → error |

R 4.6.1, forestsearch 0.2.0.
