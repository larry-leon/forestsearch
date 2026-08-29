# REPORT — `fs_family_report()`: what is and is not deterministic about the candidate family

**Task:** `dev/tasks/cc_task_fs_family_report_2026-08-30.md` (commit `424d27ba`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Ran after:** `cc_task_oc_wrapper_confs_2026-08-30.md` (`7e4a5dfd`, `58024dfe`, `49f6e1ce`) and `cc_task_vi_grf_default_2026-08-30.md` (`5332f00f`, `06728783`, `8d6cc64c` — that task stopped at its §2 gate, outcome (B), and made no edit), not concurrently with either.
**Category:** **additive only.** One new file `R/fs_family_report.R` (one export, one `print()` method), one new test file, `man/fs_family_report.Rd`, `NAMESPACE` (+1 export, +1 `S3method`), `DESCRIPTION` (0.2.6 → 0.2.7), `NEWS.md`. **No existing `R/` file was edited.** No default, driver, document or OC-wrapper file changed. No simulations, no renders, no push.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 8d6cc64c · git status --porcelain: empty
git log -6: 8d6cc64c 06728783 5332f00f 49f6e1ce 58024dfe 7e4a5dfd   (both prior tasks' commits present)
packageVersion forestsearch 0.2.6
```

`ls ~/Downloads/cc_task_fs_family_report_2026*`: exactly one match, `cc_task_fs_family_report_20260830.md` (hyphens stripped) → `dev/tasks/cc_task_fs_family_report_2026-08-30.md`, committed **`424d27ba`**.

**Parity baseline** `devtools::test()` on the unchanged tree: `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4737 ]` — the run made at the start of the `vi.grf.min` task on the byte-identical tree (that task edited nothing; HEAD moved only by its two report commits), reused here rather than re-run.

---

## 2. §2 — the stage inventory, verified from source

Every row of chat's table checked against the governing line. **Corrections and additions are marked ▶.**

| stage | governed by | verified status | governing line(s) |
|---|---|---|---|
| cut construction | `conf.cont_jcuts`, `cut_type`, `cont.cutoff`, `conf.cont_medians`, `conf.cont_medians_force`, `conf_force`, `defaultcut_names`, `exclude_cuts`, `collapse_cuts`, `collapse_cuts_args` | **data-dependent, not disableable** — continuous cuts at sample quantiles | `get_fsdata.R` L273–277 (`is.continuous(x, cutoff = cont.cutoff)`), L302–310 (`get_conf_force_jq()`), L396–420 (`get_conf_force()` default grid), L476–481 (medians), L486 `confs <- c(conf.categorical, conf.cont_Medcuts)`; `fs_oc_family.R` L6–8 documents the same fact from the other side ("cuts land at POPULATION quantiles" is what makes the OC family deterministic) |
| LASSO screen | `use_lasso` | disableable | `get_fsdata.R` L362 `if (use_lasso && cut_type == "default")`, L386 |
| DINA / GRF subgroup paths | `use_dina`, `use_grf`, `subgroup_method` | disableable | `forestsearch_main.R` L2462–2492 (GRF cut generation, `if (use_grf …)`), L2590–2620 (DINA); `.fs_family_status()` (`forestsearch_helpers.R` L2037–2049) already names the three switches |
| GRF VI ordering | `vi.grf.min` | ordering only at ≤ 0; `NULL` skips Section 5 | `forestsearch_main.R` L2752 `if (!is.null(vi.grf.min))`; L2814–2818 `vi_max <- max(…); if (vi_max > 0) { vi_ratio <- vi / vi_max; which(vi_ratio > vi.grf.min) }` — `grf::variable_importance()` is non-negative, so nothing is dropped at ≤ 0; the sort at L2812–2813 (`order(vi.cs, decreasing = TRUE)`) is what remains |
| confounder cap | `max_n_confounders` | applied only inside the VI block → inert when `vi.grf.min = NULL` ▶ **and also inert whenever the cap ≥ the number of cut columns** (the report says so when `data` is supplied) | L2820, inside `if (vi_max > 0)` inside `if (!is.null(vi.grf.min))` |
| dummy expansion | — | deterministic | L2838–2846 `dummy(df.confounders)` |
| combination enumeration | `maxk` | deterministic given the columns | `subgroup_search.R` (`generate_combination_indices`, `calculate_max_combinations`) |
| per-factor prevalence floor | `minp` | data-dependent ▶ **disabled at `minp = 0`, and `sg_focus = "maxeff"` sets it to 0** | `subgroup_search.R` L572, L688–690 `all(colMeans(x) >= minp)`; `forestsearch_main.R` L1484–1489 |
| redundancy | `rmin` | data-dependent, column-order sensitive ▶ **not disableable: `rmin` is not a `forestsearch()` formal (the search receives `subgroup.search()`'s default 5; `maxeff` sets 0 via `search_overrides`), and at 0 it still flags an added factor that removes no subject** | `subgroup_search.R` L500–514 `extract_idx_flagredundancy()` (walks `seq_len(ncol(x))`; `if (nx.last - nx.this <= rmin) flag.redundant <- TRUE`); `forestsearch_main.R` L2922–2925 |
| subgroup size floor | `n.min`, `n.min.frac` | data-dependent, not disableable | `subgroup_search.R` L612–613 `nx <- sum(id.x); if (nx <= n.min) return(status 4)`; `forestsearch_main.R` §1A2 (`n.min = NULL` → `max(60, ceiling(n.min.frac * N))`) |
| arm-count floors | `d0.min`, `d1.min` | ▶ **skipped for continuous *and count*** (chat: "continuous"); binary → events per arm; survival → events per arm; other GLM types → per-arm sample size | `subgroup_search.R` L586–609: `is_continuous <- outcome_type %in% c("continuous", "count")`; L593 "Continuous: d0.min/d1.min not meaningful — skip"; L598 (binary events), L605 (fallback n), L609 "skip Status 3 entirely" |
| effect screen | `effect.threshold` / `hr.threshold` | data-dependent ▶ **disabled under `sg_focus = "maxeff"`** (`disable_effect_floor`) | `forestsearch_main.R` L2915–2919 `disable_effect_floor = !.admit_applies[["effect"]]`; `.fs_admission_applies()` L2276–2277 |
| near-duplicate removal | `sg_focus` | data-dependent, **no off switch**; `maxeff` → exact-membership dedup; every other focus incl. `maxeffCons` → statistics-keyed | `subgroup_consistency_main.R` L556–590 (`if (identical(sg_focus, "maxeff")) .maxeff_membership_dedup(...) else remove_near_duplicate_subgroups(...)`); helpers L1102–1130: key = columns 2:10 = (K, n, E, d1, m1, m0, HR, L, U) rounded to `tolerance = 0.001`, `!duplicated(dup_key)`; L2036–2060: subject-membership key, fewest-cut representative |
| ▶ **candidate cap** (missing from chat's table) | `max_subgroups_search` | data-dependent when finite (truncates the *preview-ordered* pool before consistency; the order is sample-fitted); inert at `Inf`; `maxeff` forces `Inf` | `forestsearch_main.R` L3030 `stop_Kgroups = max_subgroups_search`, L1482–1488; `subgroup_consistency_main.R` L605–625 (`found.hrs[seq_len(min(n_pool, stop_Kgroups)), ]`, with the warning "truncated the candidate pool") |
| ▶ **m1 filter** (missing from chat's table) | `m1.threshold` | data-dependent when finite; inert at the default `Inf` | `subgroup_consistency_main.R` L516–531 (`hr.subgroups$m1 <= m1.threshold`) |
| consistency screen | `consistency_method`, `pconsistency.threshold`, `consistency.threshold` / `hr.consistency`, `fs.splits`, `use_twostage`, `twostage_args` | data-dependent — the gate ▶ **disabled under `maxeff` (`pconsistency.threshold → 0`)** | `forestsearch_main.R` L3025–3045; L1476–1489 |
| early stopping | `stop_threshold`, `sg_focus` | meaningful for `maxeffCons` only; reset to `NULL` for `hr` / `hrMaxSG` / `hrMinSG` / `maxSG` / `minSG` (aliases `eff`, `maxcons` → `hr`; `effMaxSG` → `hrMaxSG`; `effMinSG` → `hrMinSG`) and by the `maxeff` block; truncates what is evaluated, the prefix winner is the global winner | `forestsearch_main.R` L1554–1574 (the argument), L1575–1609 (the reset, with `.sync_args_call_all`), L1480–1486; `.FS_SG_FOCUS_ALIASES` (`forestsearch_helpers.R` L905–908) |
| winner selection | `sg_focus`, `selection_rule`, `effect_neighborhood` | deterministic given the qualifying set | `forestsearch_main.R` L1453–1459 (validation), L2154–2163 |
| time cap | `max.minutes` | **inert** — forwarded and never compared | `forestsearch_main.R` L2911 → `subgroup_search.R` L84, L157, L241 (pass-through only; L25 roxygen: "Currently inert"); `run_simulation_analysis.R` L63 says the same |

**Status-classification changes to chat's table (prominent):**

1. ▶ **Arm-count floors are skipped for `count` outcomes as well as `continuous`** (`subgroup_search.R` L587–588). Chat's row said continuous only. The function's row keys on `outcome_type %in% c("continuous", "count")`.
2. ▶ **Two stages were missing**: the `max_subgroups_search` candidate cap (data-dependent when finite — and its default `Inf` is *not* the engine's own default of 200 for `stop_Kgroups`; `forestsearch()` passes `Inf`, so it is inert unless the caller sets it) and the `m1.threshold` filter (inert at `Inf`). Both sit between the effect screen and the consistency stage. The report has 19 rows, not 17.
3. ▶ **`minp` is disableable** (set 0, or `maxeff`), so it is not intrinsic; **`rmin` and `n.min` are** (no formal for `rmin`; `rmin = 0` still prunes; `n.min` has no off value). Chat's table did not separate the three floors by disableability; the function does — the intrinsic set is exactly {cut construction, `rmin`, `n.min`, near-duplicate removal}.
4. ▶ **`maxeff` relaxes more than chat listed**: effect floor off, consistency screen pass-all, `minp` 0, `rmin` 0, `max_subgroups_search` `Inf`, `use_twostage` `FALSE`, `stop_threshold` `NULL` (L1476–1489, L2919, L2925). Each affected row carries a "(maxeff: …)" marker.

Everything else in chat's table verified as stated.

**One engine finding, reported not corrected (out of scope):** `forestsearch_main.R` L3031–3036 says in its own words that **`args_call_all` is not re-synced after the `maxeff` override block** (the locals are passed to the engine instead), whereas the `hr`-family `stop_threshold` reset at L1607–1609 *does* re-sync. So a fitted `maxeff` object's `args_call_all$stop_threshold` / `$pconsistency.threshold` / `$minp` carry the user's stale values. `fs_family_report()` therefore re-resolves from the raw arguments rather than trusting a fitted object's recorded values, and the drift-guard test uses an `hr`-family focus, where the engine's record is faithful.

---

## 3. The function

`R/fs_family_report.R`, exported: **`fs_family_report(x, data = NULL, outcome_type = NULL)`**; `print.fs_family_report()` registered as an S3 method.

- `x`: a named list of `forestsearch()` arguments (unspecified ones take the formal defaults, evaluated in the argument environment so that `stop_threshold = pconsistency.threshold` resolves), or a fitted `forestsearch` object (reads `args_call_all`; `outcome_type` from the object). Anything else, or a list with a name that is not a `forestsearch()` formal, errors with a clear message. `outcome_type` is required when `x` does not carry it.
- `data`: optional; when supplied, `get_FSdata()` is called exactly as `fs_oc_family_enumerate()` calls it (LASSO/GRF off, two constant columns for the outcome/event contract) to count the cut columns, and `calculate_max_combinations(L, maxk)` gives the combinations. Nothing is fitted; the search is not run.
- Resolution mirrored: `.normalize_sg_focus()`; the `maxeff` overrides; the `hr`-family `stop_threshold` reset; `rmin` from `formals(subgroup.search)`; `n.min = NULL` → `max(60, ceiling(n.min.frac * N))` (with `data`) or a note (without); `max_n_confounders` inert when `vi.grf.min = NULL` (or, with `data`, when the cap ≥ cut columns); `d0.min`/`d1.min` inert for continuous and count.
- Return: `data.frame` of class `c("fs_family_report", "data.frame")`, 19 rows, columns `stage`, `arguments`, `values`, `status` ∈ {`deterministic`, `disabled`, `inert`, `data-dependent`, `data-dependent (not disableable)`}, `note`. Attributes `verdict`, `status_counts`, `data_supplied`, `outcome_type`, and with `data` `n_cut_columns`, `n_combinations`.
- `print()`: verdict line; the counts line when `data` was supplied; the table (`stage`, `status`, `values` truncated to 60 characters); a footer listing every `data-dependent (not disableable)` row with its note.
- Roxygen states that it reports and changes nothing, and that no combination of arguments makes the family deterministic while cuts are placed at sample quantiles. `@examples` are runnable (no search is run).

---

## 4. Tests — `tests/testthat/test-fs-family-report.R`

| # | test | expectations |
|---|---|---|
| 1 | shape and contract: statuses from the allowed set, every `arguments` name a real formal, no duplicate stage, attributes, prints (verdict line and footer), `data = NULL` and `data` supplied, three error paths | 18 |
| 2 | **drift guard**: `forestsearch()` on the 300-row synthetic continuous fixture with `sg_focus = "eff"` (→ `hr`) and `stop_threshold` at its default, `vi.grf.min = NULL`; asserts `fit$args_call_all$stop_threshold` is `NULL`, the report's resolved value equals it, the early-stopping row is `disabled`, the cap row `inert` and the VI row `disabled` on the fitted object | 7 |
| 3 | **coverage guard**: `setdiff(names(formals(forestsearch)), classified ∪ .FR_OUT_OF_SCOPE)` is empty, the two sets are disjoint, and every out-of-scope name is a real formal. `.FR_OUT_OF_SCOPE` (47 names: data/column names, estimand, supplied screening results and engine tuning, execution/printing/plotting, MR) lives in the test file | 3 |
| 4 | status sensitivity: `vi.grf.min` `NULL`/numeric flips the VI and cap rows; `use_lasso`, `use_grf` flip theirs; `maxeff` disables floor/screen/effect/cap and switches dedup to exact-membership; `outcome_type` `count`/`binary` moves the per-arm row; `stop_threshold` under `maxeffCons` vs `eff`; then, under **ten** argument combinations (drivers' set; `vi.grf.min = NULL`; all front ends off; `maxeff`; `eff`; `hrMaxSG`; survival; binary; floors zeroed/`Inf`; `conf_force` values with `conf.cont_jcuts = NULL`), the cut-construction, near-duplicate, `n.min` and `rmin` rows are `data-dependent (not disableable)` and the verdict says "intrinsic to the method" | 76 |
| 5 | fitted-object path: the same fit's table `identical()` to the table from its `args_call_all` | 3 |

Run in the repository (`devtools::test(filter = "fs-family-report")`): **107 expectations, 0 failures, 0 warnings, 0 skips.**

**Injected defect (v5 §9).** A copy of the function with the near-duplicate row's status hard-coded to `"disabled"` (one token changed: `DDN` → `DIS` on that `add()` call), run against the same test file: **20 failures**, all in test 4 — `row_status(r, "near-duplicate removal")` ≠ `"data-dependent (not disableable)"` and `status_counts[["data-dependent (not disableable)"]] >= 4` false, once per combination. Tests 1, 2, 3, 5 unaffected, as expected (they do not assert that row's status). The copy was not committed.

---

## 5. Worked `print()` on the drivers' own argument set

The ACTG175 MD40 drivers' arguments (`sim_fs_maxeffCons_mr_md40_knoise0_n{500,700}_batch_1_1000.qmd`: 13 confounders incl. `str2`, `conf.cont_jcuts = list(age = 10, preanti = 10)`, `n.min = 60`, `d0.min = d1.min = 12`, `maxk = 2`, `vi.grf.min = -0.2`, `sg_focus = "maxeffCons"`, `selection_rule = "neighborhood"`, `effect_neighborhood = 0.10`, `stop_threshold = NULL`, `consistency_method = "resample"`, `effect.threshold = 30`, `consistency.threshold = 10`, `pconsistency.threshold = 0.90`, `fs.splits = 400`, `use_twostage = TRUE`, LASSO/DINA/GRF off), with `data =` the ACTG175 arms 1/3 frame (n = 1,083):

```
Candidate family is data-dependent: 8 of 19 stages vary with the sample; 4 are disableable, 4 are intrinsic to the method.
On the supplied data: 76 cut columns, 2,926 combinations of up to maxk.

 stage                            status
 cut construction                 data-dependent (not disableable)
 LASSO screen                     disabled
 DINA / GRF subgroup paths        disabled
 GRF variable-importance ordering data-dependent
 confounder cap                   inert
 dummy expansion                  deterministic
 combination enumeration          deterministic
 per-factor prevalence floor      data-dependent
 redundancy (rmin)                data-dependent (not disableable)
 subgroup size floor              data-dependent (not disableable)
 per-arm floors                   inert
 effect screen                    data-dependent
 near-duplicate removal           data-dependent (not disableable)
 candidate cap                    inert
 m1 filter                        inert
 consistency screen               data-dependent
 early stopping                   disabled
 winner selection                 deterministic
 time cap                         inert
 (values column: conf.cont_jcuts = list(age = 10, preanti = 10); cut_type = "default"; cont.cutoff = 4; ... | use_lasso = FALSE | use_dina = FALSE; use_grf = FALSE; subgroup_method = "consistency" | vi.grf.min = -0.2 | max_n_confounders = 1000 | ... | maxk = 2L (76 cut columns -> 2,926 combinations) | minp = 0.025 | rmin = 5 (subgroup.search() default) | n.min = 60L; n.min.frac = 0.1 | d0.min = 12L; d1.min = 12L | effect.threshold = 30 | sg_focus = "maxeffCons" -> statistics-keyed dedup | max_subgroups_search = Inf | m1.threshold = Inf | consistency_method = "resample"; pconsistency.threshold = 0.9; consistency.threshold = 10; fs.splits = 400L; use_twostage = TRUE; twostage_args = list() | stop_threshold = NULL | sg_focus = "maxeffCons" (canonical "maxeffCons"); selection_rule = "neighborhood"; effect_neighborhood = 0.1 | max.minutes = 3)

Intrinsic to the method -- no argument switches these off:
  * cut construction: continuous cuts are placed at SAMPLE quantiles (get_FSdata(): get_conf_force() / cut_var_jq(), medians for cut_type = "median"); only a full set of cut VALUES supplied through conf_force with the defaults excluded would fix them, and the floors and near-duplicate removal below still vary
  * redundancy (rmin): each added factor must shrink the SAMPLE membership by more than rmin subjects, walking the columns in their current order (extract_idx_flagredundancy()); rmin = 0 (maxeff) still drops exact-membership prefixes, so it cannot be switched off, and rmin is not a forestsearch() formal
  * subgroup size floor: sum(id.x) > n.min on the SAMPLE (subgroup.search() status 4); a floor on a sample count with no off switch
  * near-duplicate removal: candidates whose (K, n, E, d1, m1, m0, HR, L, U) agree to 0.001 collapse to the first (remove_near_duplicate_subgroups(), subgroup_consistency_main.R L579); statistics are sample-fitted; NO argument disables it
```

The same arguments with `vi.grf.min = NULL` — the "fixed family" belief the function exists to correct — give: **"Candidate family is data-dependent: 7 of 19 stages vary with the sample; 3 are disableable, 4 are intrinsic to the method."** The VI row becomes `disabled`; the four intrinsic rows do not move. (The 76 cut columns on the ACTG175 frame are 2 × 38 cut expressions, against 2 × 37 = 74 for the corrected OC-wrapper family enumerated on `df_super` (`oc_wrapper_grid_corrected_2026-08-30.rds`): the frame adds `wtkg <= 74`, and **15 of the 30 continuous thresholds differ** between the two — `age <= 40` vs `39`; `preanti <= 56 / 214 / 432 / 678 / 873 / 1071` vs `39 / 188 / 402 / 672 / 875 / 1055`; `wtkg <= 74, 82` vs `83`; `cd40 <= 348 / 338 / 421` vs `349 / 339 / 420`; `cd80 <= 988 / 649 / 1208` vs `989 / 645 / 1210` — the same arguments, two frames, two cut grids. That is the sample-quantile point in one line, and it is why the cut-construction row can never be `deterministic`.)

---

## 6. Close-out

`devtools::document()`: wrote `NAMESPACE` (+`export(fs_family_report)`, +`S3method(print,fs_family_report)`) and `man/fs_family_report.Rd`; nothing else regenerated.

`devtools::test()` after the addition: `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]` — parity against the baseline `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4737 ]`: WARN 31 = 31, SKIP 3 = 3, PASS +107 = exactly the new file's 107 expectations; no existing test changed behaviour.

`R CMD check` (`devtools::check(document = FALSE, args = "--no-manual")`, `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64` on `PATH`, pandoc 3.8.3): 

```
── R CMD check results ───────────────────────────────── forestsearch 0.2.7 ────
Duration: 10m 26.6s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
Status: OK
```


`DESCRIPTION` 0.2.6 → **0.2.7**; `NEWS.md` entry states that it reports and changes nothing. Staged by explicit path (`R/fs_family_report.R`, `tests/testthat/test-fs-family-report.R`, `man/fs_family_report.Rd`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`, this report); `git status --porcelain` clean after. **No push.** `devtools::install()`: `Execution halted`.

Commits: `424d27ba` task doc; `d5ef8c7c` function, tests, man, NAMESPACE, DESCRIPTION, NEWS and this report; the hash line in the follow-up commit.

---

## 7. Verdict (ten lines)

1. Chat's 17-row stage table verified from source; **four changes**: arm-count floors are skipped for count *and* continuous; the `max_subgroups_search` cap and the `m1.threshold` filter were missing (19 rows now); `minp` is disableable while `rmin` and `n.min` are not; `maxeff` relaxes six things, not two.
2. The intrinsic set — what no argument switches off — is exactly {cut construction at sample quantiles, `rmin`, `n.min`, statistics-keyed near-duplicate removal}; every other data-dependent stage has a switch.
3. `fs_family_report(x, data, outcome_type)` returns a 19-row classed data frame with a verdict attribute; `print()` leads with the verdict and ends with the intrinsic rows — the footer is the deliverable.
4. It mirrors the engine's resolution (`sg_focus` aliases, `maxeff` overrides, `stop_threshold` reset, `max_n_confounders` gating, `d0.min`/`d1.min` by outcome type) and a drift-guard test holds it there.
5. Engine finding, not fixed: `args_call_all` is not re-synced after the `maxeff` override block (the engine's own comment, L3031–3036); the report re-resolves from raw arguments for that reason.
6. Five tests, 107 expectations pass; the coverage guard forces a decision on any new `forestsearch()` formal; an injected `disabled` on the dedup row fails 20 expectations.
7. On the drivers' arguments: 8 of 19 stages vary with the sample (4 disableable, 4 intrinsic); with `vi.grf.min = NULL`, 7 of 19 (3 + 4); 76 cut columns / 2,926 combinations on the ACTG175 frame.
8. Additive only: one new `R/` file, one export, one `print()` method, one test file, `man/`, `NAMESPACE`, `DESCRIPTION` 0.2.7, `NEWS.md`; **no existing `R/` file edited**; no default, driver or document touched.
9. `devtools::test()` `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]` (+107, warning parity 31 = 31); `R CMD check` 0 / 0 / 0 at 0.2.7; installed 0.2.7.
10. Out of scope and untouched: any `screening = "none"` setter, family passing, the `args_call_all` re-sync.
