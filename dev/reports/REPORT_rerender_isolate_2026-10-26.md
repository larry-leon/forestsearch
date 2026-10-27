# REPORT — isolate the sources of the re-render differences (build ladder + bisect); fresh gbsg effMaxSG run

**Task:** `dev/tasks/cc_task_isolate_rerender_2026-10-26.md` (committed alone, `0c00959b`) · **Date:** 2026-10-26 · **Machine:** pop-os, Linux, 128 physical cores (documents resolve `n_cores` = 102)
**Record:** `dev/reports/rerender_isolate_2026-10-26/` · **Fresh runs:** `quarto/applications/actg175/rerun_2026-10-26/{C,B}_<rung>/`, `quarto/applications/gbsg/rerun_2026-10-26/G_<rung>/`
**Legacy files:** none written. The MD5 listing of every file under `quarto/applications/{gbsg,actg175}` outside the `rerun_2026-10-26/` directories (131 files) matched `legacy_tree_md5_before.txt` after all 12 ladder renders and at the end (`gate_log.txt`). The tracked LOO cache `quarto/GuoHe/cv_out_analysis_gbsg_survival_effMaxSG_hrMaxSG_neighborhood.rds` also stayed unchanged (MD5 `4844ae28…`, checked after every render).

## Verdict

1. **Did anything move at the `lean` rung? No.** None of the 130 key quantities (C 28, G 29, B 73) changes between `depr` and `lean`. The `depr`→`lean` output-line diffs are 0 (C), 0 (B) and 10 wall-clock lines (G). C, B and G payloads are leaf-identical between `depr` and `lean`, apart from G's run-directory cache path.
2. Every legacy record reproduces on this machine under its own era build, on every key quantity. That holds for B and G too, even though both legacy renders were made on a Mac (0.2.0, 11 workers): the platform bucket is empty. The legacy→era line diffs are worker counts, batch lines, timings and one `data.table` attach message.
3. **C:** all 16 changed quantities first depart at `cran040`. The bisect names **`7713942e` "fix(mr): admit on the screen's rounded consistency threshold" (2026-09-23)**. Its effect is β̃ 21.79 → 22.05, the field one-sided lower bound −38.59 → −38.31, and the Bonferroni pair −51.45 → −51.16. The search, Ĥ, the naive estimates and p̂(Ĥ) are unchanged.
4. **G (fresh effMaxSG run):** the same subgroup `{er <= 0} & {size <= 35}` and the same MR row 1.44 (0.62, 3.36), LOO 1.56 (0.81, 3.02) and naive 2.54. Six quantities move at `cran040`: the MR point estimate 1.442 → 1.444, FB 1.94 (0.89, 4.21) → 1.93 (0.89, 4.20), and the candidate count 124 → 121 (to-evaluate 123 → 120).
5. G's headline takes **three** values along history: 1.442 (era) → 1.449 (every intermediate commit tested, 2026-08-29 … 09-18) → 1.444 (0.4.0). The G bisect was specified as two-valued, so every step exited 125. It hit its 90-min cap without a verdict.
6. **B:** 11 quantities first depart at `cran040`. These are MR-gate values in the third decimal (OR 1.794 → 1.795) and the FS-dina bootstrap/T1 row (0.85 (0.23, 3.11) → 0.87 (0.25, 2.97)). Package history only: no `vi.grf.min` change applies at this rung.
7. **B:** 8 quantities (plus 3 clause-order-only reprints) first depart at `depr`, the deprecation bucket. The screened family goes 145 → 160 (DINA-screen 168 → 184), near-duplicates 9 → 23, and the consistency pool 136 → 137. The selected subgroup is the same (clauses reprinted in a different order). Five MR-gate values move at both `cran040` and `depr` (e.g. upper 5.082 → 5.097 → 5.095).
8. The previous task's differences are therefore all explained: C entirely by `7713942e`; B by package history plus the deprecation; G by package history. None of them by the lean fits.
9. G's fresh run is committed at `quarto/applications/gbsg/rerun_2026-10-26/G_lean/` (HTML and payload). There is no legacy G payload in the repo (`_payloads/analysis_gbsg_survival_effMaxSG/` does not exist), so the leaf comparison is fresh-vs-era (supplementary).
10. The lean build is reinstalled from the tree and passes `test-lean-fits.R` (344 / 0 failures, invoked with the installed package attached; see deviations). All five task worktrees are removed.

## 1. Provenance

| item | value |
|---|---|
| branch / start HEAD | `master` / `7518c955`; task commit `0c00959b` |
| GATE `merge-base --is-ancestor ef058d81 HEAD` | PASS |
| installed at start | forestsearch 0.4.0.9000, Built `R 4.6.1; ; 2026-10-26 17:45:59 UTC; unix` (as expected) |
| toolchain | quarto 1.10.18 · R 4.6.1 · survival 3.8.12 · grf 2.6.1 · speff2trial 1.0.5 |
| final install | 0.4.0.9000, Built `2026-10-27 01:56:21 UTC`, `.fs_logit_two_group` present (`final_install.txt`) |

§0: the previous task's in-place outputs had already been discarded (the working tree was clean at start). The only leftover was the **empty** untracked directory `quarto/applications/actg175/_payloads/analysis_actg175_binary_multimethod_fixed_family/`. `git log --all` on that path prints nothing, so it was removed with `rmdir`. Nothing was left to copy, so `dev/reports/rerender_2026-10-26/inplace_outputs/` was not created; the previous record and `~/Downloads` hold its HTMLs. The §0.3 gate passed.

## 2. Inventory summary (`inventory.txt`)

| doc | legacy era commit | version at era | legacy platform | settings | reads / writes |
|---|---|---|---|---|---|
| C `analysis_actg175_continuous_intervals.qmd` | `f3f188f9` 2026-09-16 | 0.3.5 | this machine | fs.splits 500, plan sequential, MR field draws 5000, no bootstrap | reads nothing; writes `<qmd dir>/_payloads/analysis_actg175_continuous_oc_intervals/…_effmaxsg_payload.rds` |
| B `analysis_actg175_binary_multimethod_fixed_family.qmd` | `43b051b6` 2026-08-16 | 0.2.0 | macOS Tahoe 26.6, aarch64, R 4.5.2, 11 workers | NB 1000, mr_draws 5000, fs.splits 500, multisession n_cores, run_cv/run_loo FALSE, `vi.grf.min = -0.2` ×2 | reads nothing; writes `<qmd dir>/_payloads/<stem>/<stem>_payload.rds` (no legacy payload exists) |
| G `analysis_gbsg_survival_effMaxSG.qmd` | `43b051b6` 2026-08-16 | 0.2.0 | macOS Tahoe 26.6, aarch64, R 4.5.2, 11 workers | NB 1000, mr_draws 5000, fs.splits 1000, multisession n_cores, LOO + 50×10-fold CV on, GH off (focus ≠ maxeff) | LOO cache read/write via `GH_DIR` (default `../../../quarto/GuoHe`, parent-relative, escapes a copy); payload `<qmd dir>/_payloads/<stem>/…` (no legacy payload exists) |

No document source changed since its era commit. No absolute or `here::here()` write paths; the `~/GitHub/fs-post-selection-manuscript_jasa` strings are comment examples.

### Render wall-clocks (s; `timing.txt`)

| doc | era | cran040 | depr | lean |
|---|---|---|---|---|
| C (cap 20 m) | 89 | 89 | 90 | 91 |
| G (cap 3 h) | 647 | 276 | 284 | 198 |
| B (cap 30 m) | 706 | 518 | 489 | 394 |

G's first render (647 s) was under the 45-min guard, so the full G ladder ran. A G render under 10 min (cran040, 276 s) qualified the G bisect.

## 3. Ladder tables with attribution

Full tables: `ladder_C.md/.csv`, `ladder_G.md/.csv`, `ladder_B.md/.csv` (columns: legacy, era, cran040, depr, lean, path, first departs, cause, changes at). In the `path` column each symbol is one rung in the order era, cran040, depr, lean, compared with legacy: `=` same, `x` changed, `~` clause order only.

| doc | quantities | none | first departs at era | at cran040 | at depr | at lean |
|---|---|---|---|---|---|---|
| C | 28 | 12 | 0 | 16 | 0 | **0** |
| G | 29 | 23 | 0 | 6 | 0 | **0** |
| B | 73 | 51 | 0 | 11 | 8 (+3 clause order only) | **0** |

**C (all `=xxx`, first departs `cran040`, bisect `7713942e`):**

| quantity | legacy = era | cran040 = depr = lean |
|---|---|---|
| H MR (IJ two-term) β̃ / IJ SE / CI | 21.79 / 50.60 / (−77.38, 120.96) | 22.05 / 50.78 / (−77.48, 121.58) |
| H MR (field) est₂ / λ-SD / one-sided lower | 16.74 / 32.31 / −38.59 | 17.16 / 32.47 / −38.31 |
| Hc MR (IJ two-term) β̃ᶜ | −32.74 (15.06) | −32.72 (15.05) |
| Hc field unstud. / field-s upper | −17.37 / −17.66 | −17.38 / −17.64 |
| Bonferroni pair (unstud.) / (field-s) | −51.45, −15.38 / −51.45, −15.82 | −51.16, −15.37 / −51.16, −15.82 |
| joint corr / aligned draws | −0.038 / 983 | −0.042 / 985 |
| selection rate; field Λ*-mean; complement fits | 0.972; +5.049; 1149 | 0.974; +4.894; 1147 |
| top re-selected share | 0.077 | 0.078 |

Unchanged: N, ITT, Ĥ, n(Ĥ) = 79, T_obs, Pcons, the band line, both naive rows, p̂(Ĥ) = 0.064, the 2nd/3rd re-selected, n_family 4935.

**G (all `=xxx`, first departs `cran040`):**

| quantity | legacy = era | cran040 = depr = lean |
|---|---|---|
| passed effect threshold / candidates found | 124 / 124 | 121 / 121 |
| to evaluate (on frontier 14, in band 2) | 123 | 120 |
| consistency: evaluated (passed 9, frontier 3, band 2, m = 1) | 123 | 120 |
| MR harm confirmation HR (MR point estimate) | 1.442 | 1.444 |
| summary: Full bootstrap | 1.94 (0.89, 4.21) | 1.93 (0.89, 4.20) |

Unchanged: the combinations evaluated (2211), Pcons 0.99, the rank-1 row, the selected subgroup, max HR 2.54, the naive H/Hc rows, the Pareto frontier (3 members, selected row), FB success 976 (97.6%), MR family 1744, the scale cross-check, 10-fold CV sens/PPV/exact/partial, the LOO line, agreement and KfA rows, and the summary Naive / MR 1.44 (0.62, 3.36) / LOO 1.56 (0.81, 3.02).

**B, first departs `cran040` (package history; legacy = era for all):** the FS MR harm OR 1.794 → 1.795; the Pareto selected row's Split/FSBC CIs (0.50, 25.61) / 3.30 (0.46, 23.60) → (1.23, 10.48) / 3.42 (1.17, 10.02); the SG_CIs H MR-gate upper 5.08 → 5.10 → 5.09 (depr); the gate de-biased OR T2 1.794 (0.633, 5.082) → 1.795 (0.632, 5.097) → (0.632, 5.095); the complement upper 1.026 → 1.027; selection bias 0.6915 → 0.6919 → 0.6922; IJ SE 0.531 → 0.533 → 0.532; complement selection bias −0.0499 → −0.0501; FS-dina bootstrap success 978 → 977; the FS-dina SG_CIs H bootstrap 0.85 (0.23, 3.11) → 0.87 (0.25, 2.97); the FS-dina gate T1 0.851 (0.233, 3.106) → 0.867 (0.253, 2.969).

**B, first departs `depr` (deprecation; the candidate family):** FS passed effect screen 145 → 160 (and the table count and `family_status` "screened survivors"); near-duplicates removed 9 → 23; consistency pool 136 → 137; on frontier 10 → 11; DINA-screen passed 168 → 184 and pool 165 → 166. The selected subgroup is reprinted as `!{cd40 <= 380} & !{wtkg <= 86}` (same clauses, three quantities). Consistency passed (5), the selected H/Hc naive rows and FS bootstrap rows are unchanged.

### Output-line diff counts between adjacent rungs (`rungdiff_counts.csv`)

| doc | legacy→era | era→cran040 | cran040→depr | depr→lean | legacy→lean |
|---|---|---|---|---|---|
| C (60 lines) | 0 / 0 | −18 / +18 | −1 / +1 (wall-clock) | **0 / 0** | −19 / +19 |
| G (447 → 445) | −23 / +18 (workers, batches, timing) | −44 / +47 | −9 / +9 | −10 / +10 (**all wall-clock**) | −54 / +52 |
| B (770 → 1619) | −25 / +25 (workers, batches, attach message) | −23 / +860 (new print blocks) | −59 / +71 | **0 / 0** | −101 / +950 |

### Payloads (supplementary, read-only; `payload_CB_rungs.txt`, `payload_G_compare.txt`)

- **C:** legacy payload → era_C: **0 of 187 leaves differ**. era → cran040: 61. cran040 → depr: 0. depr → lean: 0.
- **B:** era → cran040: 9 of 108 leaves. cran040 → depr: 7. depr → lean: 0.
- **G:** cran040 → lean and depr → lean: only `extras.loo.cache` (the run-directory path). era_G → lean: 138 leaves differ. Most are the provenance table shifting by one row: `pconsistency.digits` is a new formal, from `7713942e`. The substantive ones are `extras.mr.debiased.est` 1.442118 → 1.443947 (lower 0.6193 → 0.6206, upper 3.3584 → 3.3596, se 0.43130 → 0.43084, ij_draws 4964 → 4966, selection_bias 0.55565 → 0.55457), FB H2 1.9353 (0.8895, 4.2105) → 1.9340 (0.8897, 4.2039), and four 10-fold-CV identification-summary rows (fold counts). 15 leaves exist only in lean.

## 4. Bisect results

### C — conclusive (16 min, 9 steps; `bisect_C_steps.txt`, `bisect_C_result.txt`, `bisect_C_gitlog.txt`)

The pattern warranted it (era 21.79 = legacy; cran040 22.05). Run in worktree `fs_worktrees/bisect_C`: `git bisect start a0025013 f3f188f9`, 404 commits.

**First bad commit: `7713942e07f6268457320160eae64adf8017b8b0`**, 2026-09-23 16:26:15 −0700, *fix(mr): admit on the screen's rounded consistency threshold*.
> MR's consistency floor used qnorm((1 + p_star)/2), the exact cutoff, while the search admits on round(Pcons, pconsistency.digits) >= p_star. Build the floor from .fs_pcons_eff(p_star, digits) instead. fs_mr_inference() gains pconsistency.digits (NULL -> subgroup.consistency() default); forestsearch() passes its own. Forest search only: GRF/DINA carry no consistency floor.

`git show --stat`: `R/forestsearch_main.R | 1 +`, `R/fs_mr_inference.R | 23 +++++-`, `tests/testthat/test-mr-admission-rounded.R | 108 +` (new). NEWS.md entry, added by the follow-up commit `b3463489`: "MR admission at the screen's rounded threshold … **Corrected estimates, field bounds and the Bonferroni pair can move for forest search.** … GRF and DINA carry no consistency screen and are unaffected" (full text in `bisect_C_result.txt`). This matches the ladder: only MR quantities move in C. In B, the FS gate moves at cran040 while the GRF-path rows do not.

Steps (good = 21.79, bad = 22.05): c4baf796 (09-21) good · 9224dfe9 (09-24) bad · 85ecbbfc (09-23) good · ba595f4b good · 1b11967c bad · d637ff07 bad · b3463489 bad · 6fabaed5 bad · 7713942e bad.

### G — inconclusive (90-min cap; `bisect_G_steps.txt`, `bisect_G_gitlog.txt`)

The pattern warranted it (era 1.442 = legacy; cran040 1.444), and a single G render took ≤ 10 min. Run in worktree `fs_worktrees/bisect_G`: `git bisect start a0025013 43b051b6`, 1062 commits. All 14 completed steps returned **1.449**, which is neither value, so each exited 125 as specified. Step 15 (`c1752b99`) was killed by the cap mid-render. The tested commits span 2026-08-29 (`d9ccc24c`) to 2026-09-18 (`4e1b6dd4`).

What the steps do establish:
- G's MR estimate moves **twice**. The first move, 1.442 → 1.449, lies in `(43b051b6, d9ccc24c]`. Eleven commits touch `R/`/`DESCRIPTION` there (`16e6bd96` 08-17 … `202ef15b` 08-29 "vi.grf.min default -0.2 -> NULL (0.3.0)"; listed in `bisect_G_brackets.txt`).
- The second move, 1.449 → 1.444, lies in `(4e1b6dd4, a0025013]`, a range that contains `7713942e`. That commit is the plausible candidate, given C's result and its NEWS entry naming the GBSG application, but it was **not tested for G**.

A G bisect would need a three-valued script (e.g. bad = "≠ 1.442", then a second bisect from 1.449 → 1.444) to finish; out of scope here.

## 5. G — fresh vs legacy

The fresh run (`G_lean`, 198 s on 102 workers) against the legacy HTML (Mac, 0.2.0, 11 workers):

| | legacy | fresh (lean) |
|---|---|---|
| selected subgroup | `{er <= 0} & {size <= 35}` (n 61, Pcons 0.99) | same |
| Naive HR | 2.54 (1.25, 5.17) | same |
| Full bootstrap | 1.94 (0.89, 4.21), 976/1000 success | 1.93 (0.89, 4.20), 976/1000 |
| Multiplier resampling | 1.44 (0.62, 3.36); harm-confirmation HR 1.442 | 1.44 (0.62, 3.36); 1.444 |
| Leave-one-out | 1.56 (0.81, 3.02), n 64; sens_H 0.92, ppv_H 0.88, exact 0.90 | same |
| 10-fold CV | sens_H 86.1, PPV_H 67.8, exact 0.0, partial 432 (86.4%) | same |
| candidates passing the effect screen | 124 | 121 |

Payload: no legacy payload exists (see §3). The fresh payload is `quarto/applications/gbsg/rerun_2026-10-26/G_lean/_payloads/analysis_gbsg_survival_effMaxSG/analysis_gbsg_survival_effMaxSG_payload.rds`. Its fresh LOO result is beside it as `cv_out_…rds`.

## 6. Deviations (verbatim)

1. §0.1/§0.2: the in-place outputs were already discarded (clean tree at start), so no copies were made and `inplace_outputs/` was not created. The leftover untracked payload directory was empty; `git log --all -- <path>` printed nothing; it was removed with `rmdir`.
2. §4.1, G: the LOO cache path `gh_dir <- Sys.getenv("GH_DIR", unset = "../../../quarto/GuoHe")` is parent-relative and escapes a copy directory. From a run directory it resolves to the nonexistent `quarto/applications/quarto/GuoHe`, so the render would have failed at `saveRDS`. It is not absolute or `here::here()`, and the document exposes it as an environment knob. G was rendered with `GH_DIR=<run dir>/_payloads` exported (no source edit). As in the legacy render, the cache read missed and LOO was computed; the cache file was written inside the run directory and is committed with its payloads. The tracked `quarto/GuoHe` cache was MD5-checked after every render and is unchanged. The decision was recorded in `inventory.txt` before any render.
3. §3, installs of git rungs used `quick = TRUE` as specified, which skips the build step. `Packaged` is therefore NA in their `build.txt`; the worktree SHA is recorded instead. CRAN 0.4.0 was downloaded once per ladder (3 times, the same permitted package): Version 0.4.0, Repository CRAN, Packaged `2026-10-10 04:12:05 UTC`.
4. §3, `era_B`: B's era commit is G's (`43b051b6`). B used the `era_G` worktree and was reinstalled from it before its render.
5. §5, G: the bisect script follows the task's 0 / 1 / 125 rule. G's three-valued history made every step 125, and the 90-min cap stopped it (see §4).
6. §7.1: the literal `Sys.setenv(NOT_CRAN = "true"); testthat::test_file("tests/testthat/test-lean-fits.R")` in a bare `Rscript` ran 339 expectations with 1 error, `could not find function "forestsearch"`. Nothing attaches the package in that call, so it cannot pass for any build. The installed build was then tested with `testthat::test_file(…, package = "forestsearch", load_package = "installed")`: **344 expectations, 0 failed, 0 errors, 0 skipped, 0 warnings**, the same count as `REPORT_lean_fits_2026-10-26.md`. Both runs are in `final_install.txt`.
7. §7.2: `git worktree list` still shows two worktrees that predate this task and were not created by it (`.claude/worktrees/mr-terminology`, `.claude/worktrees/replication-check`). They were left alone. All five task worktrees (`era_C`, `depr`, `bisect_C`, `era_G`, `bisect_G`) were removed and pruned, along with `~/Documents/GitHub/fs_worktrees/`.
8. §5: bisect step renders under `~/fs_bisect_{C,G}/step_<n>/` are outside the `quarto/` project, so the `_quarto.yml` format settings did not apply. This changes formatting only: C's steps reproduced the ladder's 21.79 / 22.05 exactly. The step directories are left on disk, outside the repo.
9. §6: C and B key-quantity specs are evaluated from the previous task's `compare.R` parse tree (not copied); G's 29 specs are new (`ladder.R`). Added beyond the task: a `changes at` column (every rung where a value changes from the previous rung) and supplementary payload comparisons across rungs (`payload_CB_rungs.*`). The first `ladder.R` run errored on a names bug (`file.exists()` drops names) before producing any output; it was fixed and rerun.
10. Renders inside `quarto/` update the gitignored project cache `quarto/.quarto/`.
11. Spec conflict: §7.4 places the report at `dev/reports/REPORT_rerender_isolate_2026-10-26.md`, but the closing block's first command allows only `dev/reports/rerender_isolate_2026-10-26/`. The report was committed at the path §7.4 names. After commit (b), the first closing command therefore prints exactly one line, this report's path, and nothing else. Same for the commit (b) pre-check `git diff --cached --name-only` against that filter.

## 7. Files

Record (`dev/reports/rerender_isolate_2026-10-26/`): `inventory.txt`, `provenance.txt`, `legacy_tree_md5_before.txt`, `gate_log.txt`, `timing.txt`, `install_log.txt`, `final_install.txt`, `install_rung.sh`, `run_rung.sh`, `ladder.R`, `ladder_{C,G,B}.{csv,md}`, `rungdiff_counts.csv`, `outputs_<D>_<rung>.txt`, `bisect_{C,G}.sh`, `bisect_{C,G}_{steps,gitlog,start}.txt`, `bisect_{C,G}_run.log`, `bisect_C_result.txt`, `bisect_G_brackets.txt`, `payload_G_compare.{R,txt}`, `payload_CB_rungs.{R,txt}`, `closing_block_pre_b.txt`.

## Closing block (run after commit (a), before the record was staged; re-run and printed after commit (b))

```
$ git diff --name-only 7518c955 HEAD | grep -v -E '^(dev/(tasks/cc_task_isolate_rerender_2026-10-26\.md|reports/rerender_isolate_2026-10-26/)|quarto/applications/(actg175|gbsg)/rerun_2026-10-26/)'
$ git status --porcelain | grep -v '^?? '
$ find quarto/applications/gbsg quarto/applications/actg175 -type f -not -path '*/rerun_2026-10-26/*' -exec md5sum {} + | sort -k2 | diff - dev/reports/rerender_isolate_2026-10-26/legacy_tree_md5_before.txt && echo LEGACY-UNCHANGED
LEGACY-UNCHANGED
$ git worktree list
/home/larryleon/Documents/GitHub/forestsearch                                      60a29cd5 [master]
/home/larryleon/Documents/GitHub/forestsearch/.claude/worktrees/mr-terminology     82580bf4 [feature/mr-terminology]
/home/larryleon/Documents/GitHub/forestsearch/.claude/worktrees/replication-check  5bad3df7 [wt/replication-check]
$ Rscript -e 'cat(as.character(packageVersion("forestsearch")), exists(".fs_logit_two_group", asNamespace("forestsearch")), "\n")'
0.4.0.9000 TRUE
```
