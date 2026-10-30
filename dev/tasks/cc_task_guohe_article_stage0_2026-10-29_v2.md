# CC task — Guo & He (S4) article, Stage 0: reproduction gate and timing ladder — v2

**File:** `dev/tasks/cc_task_guohe_article_stage0_2026-10-29_v2.md` · **Issued:** 2026-10-30 by chat, commissioned by Larry (new pkgdown article summarizing Section S4 of León & Anderson 2026b, arXiv:2610.28379).
**Repo:** `forestsearch` · **Branch:** `feature/articles-guohe` (exists; created by v1) · **Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone). Filenames may arrive with hyphens stripped; locate by pattern.

**v2 supersedes v1** (`dev/tasks/cc_task_guohe_article_stage0_2026-10-29.md`, committed `e751c215` on `feature/articles-guohe`). v1 stopped at its §1.5: `devtools::install(upgrade = "never")` is rejected by devtools 2.5.2 (`upgrade` must be `TRUE`, `FALSE` or `NA`); nothing after it ran. v1 stays in `dev/tasks/` as the record of that stop. **This document is the complete, authoritative task.** Only §1 (resume state, the corrected install and its gate) and the commit list in §6 differ from v1; §2–§5 and §7 are unchanged.

**What this task answers.** Two things the article cannot be written without:

- **R** — that the exported `fs_mr_inference()` and `guohe_algorithm3()` of the installed 0.4.1 build, called on Guo & He's Section 5.2 nested family, reproduce the committed per-replicate records that Section S4 of the manuscript was built from (`quarto/GuoHe/`, September 2026, package 0.3.5 with `forestsearch:::fs_mr_inference`). If they do, the article's code path *is* the manuscript's.
- **T** — the per-replicate cost of each component at the manuscript's settings and at two reduced rungs, and the cost of the truth curve β(c) at reduced large-sample sizes, so Larry can choose the article's Monte Carlo size from measured numbers.

**What this task decides: nothing.** It produces a report with numbers and a package-facts record. The article's cells, replicate counts and settings are Larry's decision afterwards.

---

## ⚠ CATEGORY

**No package code is modified.** No edit to any file under `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `_pkgdown.yml`, or anything under `quarto/`. The committed scripts and bundles in `quarto/GuoHe/` are **read-only inputs**: sourced and read, never edited, never re-run as campaigns. No `R CMD check`, no test suite, no render.

**One install** from the branch root in §1, so that the installed package is the 0.4.1 tree at HEAD (library write only; no repo file changes).

**Writes:** this document under `dev/tasks/`; the scripts, `results/` and `logs/` under `dev/reports/guohe_article_stage0_2026-10-29/` (the directory v1 created, with its two logs kept); the report `dev/reports/REPORT_guohe_article_stage0_2026-10-29.md`.

**Compute (its own line):** single core throughout (`future::plan("sequential")`; `guohe_algorithm3(parallel = FALSE)`). §3: 4 replicates × (one `guohe_algorithm3()` at B = 2000 + one `fs_mr_inference()` at 5,000 draws, 1000/500 field + 151 `coxph()` fits) — bound from the record: 39 s + 2.4 s + ≈ 0.5 s per replicate (`quarto/GuoHe/REPORT_guohe_fixedr_t7_2026-09-09.md`, `REPORT_mr_field_stage1_2026-09-05.md`), so ≈ 3 min. §4: 3 replicates × three rungs of each component ≈ 3 min; two truth curves at n_big = 2e5 and 5e4 (31 `coxph()` fits each on ≤ 150k rows) ≈ 1 min. **Caps:** §3 15 min, §4 25 min, in the driver. Nothing at 2,000 replicates; nothing parallel.

**Unattended.** Gates stop the task; never ask, never work around. The provenance gate is dirt-tolerant: pre-existing untracked files are never staged and never a stop unless they sit under this task's own paths — and this task's own folder, left by v1, is expected.

**Scripts are CC's.** Chat could not dry-run (no R in its container). Write the scripts to the specification in §3–§4, keep them in the record directory, and list every assumption you had to make in the report's "Script notes". A gate written into a script must test exactly the condition stated here — never a stricter one.

---

## 1. Resume state, install, first commit — GATE

1. `git rev-parse --abbrev-ref HEAD` must be `feature/articles-guohe`; `git log --oneline -2` must show `e751c215` (v1's task-doc commit) directly on top of `c2967362`. **MASTER_SHA = `c2967362e7ee9a6361cb1b07fdf83c0098e712db`** (from v1's §1.1). *GATE:* branch or history differ → stop and report.
2. `git status --porcelain` recorded verbatim. Expected: only `?? dev/reports/guohe_article_stage0_2026-10-29/` (v1's `logs/install.log` and `logs/s0_env.log`). Anything else is pre-existing dirt: tolerated, never staged.
3. Copy this document to `dev/tasks/cc_task_guohe_article_stage0_2026-10-29_v2.md`, `git add` it by path, commit it alone. Report both filenames (as received, as committed).
4. Record the task start: `Rscript -e 'cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"), "\n")'` → write it to `dev/reports/guohe_article_stage0_2026-10-29/results/keyvalues.txt` as `task_start_utc = …` (create `results/` if absent).
5. **Install, corrected:** from the repo root,
   `Rscript -e 'devtools::install(".", upgrade = FALSE, build_vignettes = FALSE)'` → `logs/install_v2.log`.
   If devtools rejects its *own* arguments again (an error raised by `devtools::install()` before any build starts, as in v1), run `R CMD INSTALL .` from the repo root instead → `logs/install_v2_rcmd.log`, and say so in the report's Script notes. An error raised *by the package build or install itself* is a stop.
6. **Installed-build gate**, in a fresh `Rscript`: `pd <- utils::packageDescription("forestsearch")`; record `pd$Version`, `pd$Built`, `pd$Packaged` (`"<absent>"` if `NULL`), the library path (`dirname(dirname(dirname(attr(pd, "file"))))`), `packageVersion("survival")`, `R.version.string`, `Sys.info()[c("sysname","machine")]`, `parallelly::availableCores()` — all into `results/keyvalues.txt`.
   *GATE:* `pd$Version == "0.4.1"`, **and** the date-time inside `pd$Built` (the third `;`-separated field, e.g. `2026-10-30 16:12:03 UTC`) is later than `task_start_utc`. The second condition is the machine check that the installed build is this task's build of HEAD, whichever install route ran. Either fails → stop and report (quote `pd$Built` and `task_start_utc`).
7. *GATE:* `"fs_mr_inference" %in% getNamespaceExports("forestsearch")` and `"guohe_algorithm3" %in% getNamespaceExports("forestsearch")` both `TRUE`.

## 2. Source verification — quote, do not paraphrase

Quote from the current tree, with line numbers, into the report:

- `R/fs_mr_inference.R`: the full signature of `fs_mr_inference` (the `function(` block, ≈ lines 584–609) and the `admission` handling lines that treat `NULL` floors (≈ 661–679).
- `R/guohe_algorithm3.R`: the full signature of `guohe_algorithm3` (≈ 232–248).
- `quarto/GuoHe/mr_vs_guohe_sim.R`: `mv_gh52_base()` (≈ 79), `mv_gh52_regen()` (≈ 93–98), `mv_cand_idx_52()` (≈ 107–110), `mv_spec52` (≈ 114–117), the executing block of `mv_mr()` (≈ 132–151), and the three seed constants `MV_DRAWS`, `MV_MULTIPLIER`, `MV_SEED_MR` (≈ 60–62).
- `quarto/GuoHe/guohe_sec52_sim.R`: `gh52_candidates()` (≈ 116–123), `gh52_subgroup_fits()` (≈ 132–146), `gh52_naive()` (≈ 149–162), and in `gh52_one_rep()` the `boot_seed` line and the `guohe_algorithm3()` call with its arguments (≈ 199–211).
- `quarto/GuoHe/guohe_sec52_truth.R`: `gh52_sim_data()` (≈ 68–82) and the signature of `gh52_truth_curve()` and of `gh52_truth_at()`.
- `quarto/GuoHe/mr_field_complement_vs_guohe_run.R`: the direct `fs_mr_inference()` argument block (≈ 104–131) and the lines that name the emitted complement/joint columns (≈ 141–234) — these are the column names §3 compares against; take them from the file, not from this document.
- The N1 comparison standard as the complement driver states it (search for `MFC_DISCRETE` and the `all.equal` tolerance, ≈ 310–340).

*GATE:* `mv_mr()` must call `forestsearch:::fs_mr_inference` with `reselection = "maxeff"` and `admission = list(effect_floor = NULL, consistency = NULL)`; `gh52_one_rep()` must call `guohe_algorithm3()` with `orient = +1`, `min_events = 5L`, `level = 0.05`, `diagnostics = FALSE` and `seed = boot_seed` where `boot_seed = seed + 500000L`. If either differs, stop and report — the reproduction below would compare unlike things.

## 3. R — reproduction of the committed S4 records on the installed 0.4.1 (GATE)

Script `s1_reproduce.R` in `dev/reports/guohe_article_stage0_2026-10-29/`. At the top: `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`; `library(forestsearch); library(survival)`; `future::plan("sequential")`; `source("quarto/GuoHe/mr_vs_guohe_sim.R")` (path relative to the repo root; it sources the Sec 5.2 DGM and family helpers). Nothing under `quarto/` is written.

Cells and replicates: `t7_beta2_00` and `t7_beta2_03`, m ∈ {1, 2} — four replicates. For each (id, m), with `base <- mv_gh52_base(id)`:

1. **Data and family.** `cand <- mv_gh52_regen(id, m)`; `k <- length(cand$cuts)`; `fits <- gh52_subgroup_fits(cand$df, cand)`; `tru <- readRDS(sprintf("quarto/GuoHe/guohe_sec52_truth_beta2_%02d.rds", b2*10))`; `nv <- gh52_naive(fits, cand, tru)`.
   Compare with `rep <- readRDS(sprintf("quarto/GuoHe/guohe_repro_%s.rds", id))$results[m, ]`: `nv$c_hat` vs `rep$c_hat_naive` (**identical**), `nv$point` vs `rep$naive_point`, `nv$lower` vs `rep$naive_lower`, `nv$gamma_s` vs `rep$gamma_s_naive`, `k` vs `rep$n_cand` (**identical**). Also assert `readRDS(...)$seed_base == base`.
2. **Guo & He at r = 1/12, B = 2000** — the manuscript's reported column. Call the **exported** `guohe_algorithm3(data = cand$df, outcome = "survival", treatment = "treat", candidates = cand$names, time = "time", event = "event", orient = +1, B = 2000L, r = 1/12, level = 0.05, seed = base + m + 500000L, min_events = 5L, diagnostics = FALSE)`. Compare: `fit$selected` → cutpoint `cand$cuts[match(fit$selected, cand$names)]` vs `rep$c_hat_gh` (**identical**); `log(fit$bound_one_sided)` vs `rep$gamma_s - rep$r2_dist`; `log(fit$debiased)` vs `rep$r2_bias + rep$gamma_s` (r2 is the second entry of the grid `c(1/3, 1/12, 1/21, 1/30)` — confirm from `GH52_R_GRID` in the sourced file). Record `fit$lean_cox`, `fit$fast`, `fit$B`, `fit$n_estimable`.
3. **MR, field and complement, at the manuscript's settings.** The selected label is `sel_lab <- cand$names[which.max(replace(fits$est, !is.finite(fits$est), -Inf))]`. Call the **exported** `fs_mr_inference(df = cand$df, candidates = mv_cand_idx_52(cand), spec = mv_spec52, selected_members = mv_cand_idx_52(cand)[[sel_lab]], admission = list(effect_floor = NULL, consistency = NULL), reselection = "maxeff", draws = 5000L, multiplier = "poisson", ci_method = "field", seed = base + m + 700000L, return_reselection = TRUE, include_complement = TRUE, field_R_out = 1000L, field_R_in = 500L, field_uniform = FALSE, field_complement = TRUE, field_scale_complement = "selected")`.
   Compare with `mf <- readRDS(sprintf("quarto/GuoHe/mr_field_vs_guohe_%s.rds", id))$results[m, ]`: `mr$selected_index` vs `mf$sel` (**identical**); `log(mr$debiased$est)` vs `mf$mr_est`; `mr$debiased$se_ij` vs `mf$mr_se_ij`; `log(mr$debiased$lower_1s)` vs `mf$mr_lower_1s`; `mr$selection_bias` vs `mf$mr_bias_sel`; `mr$fixed_bias` vs `mf$mr_bias_fix`; `log(mr$field$est2)` vs `mf$fld_est2`; `log(mr$field$lower_1s)` vs `mf$fld_lower_1s`; `mr$field$lambda_mean` vs `mf$fld_lambda_mean`; `unname(mr$reselection$p_hat[mr$selected_index])` vs `mf$p_hat_H`.
   Compare with `mc <- readRDS(sprintf("quarto/GuoHe/mr_field_complement_vs_guohe_%s.rds", id))$results[m, ]`, using the column names the driver emits (§2): the complement two-term estimate (`log(mr$complement$debiased$est)`), the complement field-s upper bound at 0.95 (`log(mr$field$complement$upper_1s_s)`), and the Bonferroni pair from `mr$field$joint_s` (`bonf_lower_H`, `bonf_upper_Hc`, logged). Record `mr$field$joint_s$bonf_joint_prob`, `$gamma`, `$corr`.
4. **Standard (the record's own N1):** selection keys and integer columns `identical()`; floats `isTRUE(all.equal(x, y, tolerance = 1e-8))`. For every float also record `identical(x, y)` and `abs(x - y)`.

*GATE (stop):* any selection key not `identical()`, or any float failing `all.equal` at **1e-6**. *Not a stop, record as a finding:* a float that fails at 1e-8 but passes at 1e-6 — report the largest such difference and which quantity; the article's path is then "numerically equal, not bit-reproducible" and the report says so.

Write every comparison to `results/reproduction.csv` (id, m, quantity, stored, recomputed, abs_diff, identical, all_equal_1e8, all_equal_1e6) and the key = value record to `results/keyvalues.txt`.

**Package facts for the article** (same script, after the four replicates, on the last `mr` and `fit`): `names(mr)`, `names(mr$debiased)`, `names(mr$field)`, `names(mr$field$complement)`, `names(mr$field$joint_s)`, `names(mr$complement)`, `names(mr$reselection)`, `names(fit)`; and the defaults `formals(fs_mr_inference)[c("draws","multiplier","ci_method","include_complement","field_R_out","field_R_in","field_complement","field_scale_complement","reselection")]` and `formals(guohe_algorithm3)[c("orient","B","r","level","min_events","fast","parallel")]`, deparsed. Save to `results/package_facts.md`.

## 4. T — timing ladder (single core; cell `t7_beta2_03`, seeds base + 101, 102, 103)

Script `s2_timing.R`. Three fresh replicates (data seeds `base + 101L` … `base + 103L`, with `base <- mv_gh52_base("t7_beta2_03")`, i.e. `mv_gh52_regen("t7_beta2_03", 101L)` and so on; MR seed = data seed + 700000L; GH seed = data seed + 500000L). Per replicate, `proc.time()[["elapsed"]]` around each component:

| component | rung 1 (manuscript) | rung 2 | rung 3 |
|---|---|---|---|
| data + family + 151 naive `coxph()` fits (`gh52_subgroup_fits`) | once | — | — |
| `guohe_algorithm3()` at r = 1/12 | B = 2000 | B = 500 | B = 200 |
| `fs_mr_inference()` as in §3.3 (`include_complement = TRUE`) | draws 5000, 1000/500 | draws 2000, 500/250 | draws 1000, 200/100 |

Record per rung the seconds (mean and max over the three replicates) and, for the MR rungs, the bounds themselves (`log(mr$field$lower_1s)`, `log(mr$debiased$lower_1s)`, `log(mr$field$complement$upper_1s_s)`) so the report can show how far the reduced rungs move a bound on one data set (informational; no gate).

**Truth curve** (β₂ = 0.3): time `gh52_truth_curve(0.3, n_big = 2e5, c_step = 1, seed = 20260721L + 30L, verbose = FALSE)` and the same at `n_big = 5e4`. Its internal gates (anchor β(30) = β₂, monotonicity) are the function's own; if one fails at 5e4, record it and continue (the article may need the larger size); if it fails at 2e5, stop. Compare each with the committed cache `quarto/GuoHe/guohe_sec52_truth_beta2_03.rds` through `gh52_truth_at()` at c = 30, 31, …, 60: max |Δβ(c)| and the value at c = 30 and c = 60. Record the cache's `n_big` and `c_step` (read from the object; report its `str()` top level).

**Projection table** (pure arithmetic from the measured seconds; `results/projection.md`): wall-clock for an article of cells {β₂ = 0, 0.3} × n_rep ∈ {50, 100, 200} × rung ∈ {1, 2, 3}, single core, plus the truth-curve cost once per β₂ > 0 cell; and the same at 2 workers assuming replicate-level parallelism scales linearly (labelled as an assumption). Add the manuscript's own shape for reference: 6 cells × 2000 at rung 1 (from the same per-replicate figure; labelled "not an article size").

## 5. Report — `dev/reports/REPORT_guohe_article_stage0_2026-10-29.md`

1. **Provenance** (first lines): MASTER_SHA, branch, v1's stop (`e751c215`; the devtools argument) and this document's commit, the install route taken, installed version / `Built` / `Packaged`, `task_start_utc`, R and survival versions, machine, cores; the §1 and §2 gate outcomes; the quotes of §2.
2. **R, the reproduction**: `results/reproduction.csv` transplanted as a table (four replicates × every quantity), then one line per gate: selection keys identical (count), floats at 1e-8 (count), floats at 1e-6 (count), and the largest absolute difference with its quantity. State plainly whether the 0.4.1 exported functions reproduce the S4 records bit-for-bit, to 1e-8, or to 1e-6.
3. **T, the ladder**: the component table with measured seconds; the bound movement across rungs on the three data sets; the truth-curve timings, gate outcomes and the comparison with the cache; `results/projection.md` transplanted.
4. **Package facts**: `results/package_facts.md` transplanted.
5. **Script notes**: every assumption made in writing the scripts; every deviation from this specification with the reason (file, what, why); the install route. "None" is an acceptable entry only if true.
6. **OPEN ITEMS** (never a stop): anything the article will have to decide that this task surfaced (an element name that differs from the one this document used, a default that moved, a cache field absent).
7. **Ten-line verdict**, numbers only. **No recommendation of an article size** — that is Larry's decision from the projection table.

Stop-point figures are quoted verbatim from `results/keyvalues.txt`; a paraphrased figure is ABSENT.

## 6. Commit and post-conditions — GATE

Commit by explicit path, in this order, no push, no fetch, no merge:
1. (done by v1) `dev/tasks/cc_task_guohe_article_stage0_2026-10-29.md` — `e751c215`;
2. (§1.3) `dev/tasks/cc_task_guohe_article_stage0_2026-10-29_v2.md` alone;
3. `dev/reports/guohe_article_stage0_2026-10-29/` (scripts, `results/`, `logs/` — including v1's `install.log` and `s0_env.log`) and `dev/reports/REPORT_guohe_article_stage0_2026-10-29.md`.

Before each commit `git diff --cached --name-only` must list only paths under `dev/`. After the final commit, these must hold and be printed in the report's closing block:

```
git diff --name-only c2967362 HEAD | grep -v '^dev/'              # must print nothing
git status --porcelain | grep -v '^?? ' | grep -v '^.. dev/'      # must print nothing
git status --porcelain | grep '^?? ' | grep -v ' dev/'             # only pre-existing dirt listed in §1.2
```

If the first line prints anything, do not fix it — stop and report what it prints. Copy the report and `results/` to `~/Downloads/guohe_article_stage0_2026-10-29/`.

## 7. Out of scope

- No article text, no `vignettes/articles/` file, no `_pkgdown.yml` entry — the article is the next task, after Larry sizes it.
- No re-run of any September campaign; no change to any `quarto/GuoHe/` file or bundle.
- No `R/` edit (the exported functions are used as a user would use them). No check, no tests, no render, no push, no merge.
- No parallel plan: the measurements are single-core by design so the projection is a rate, not a machine.
