# REPORT — A-002 magnitude check, v6, 2026-11-10

**Verdict: MATERIAL.** The H2 point estimate for H differs between the arms by −0.449 on the HR scale
(B − A), which is ≥ 0.05. For Hc the difference is +0.015.

## 1. Header

| Item | Value |
|---|---|
| `master` SHA | `217656e2b7142ad62932d12da91698752f174a96` |
| `fs-a002-check` before this report | `cdffe0da74160a4b73fc36a902f11aed87d25e08` (v6 task commit); `R/` identical to `master` (0 diff lines) |
| Arm A tree | the repo checkout (code as is) |
| Arm B tree | v1 patched copy `/tmp/a002_2026-11-10/fs_patched`; differs from the repo in `R/forestsearch_main.R` only (row-wise `treat.recommend`, v1 report §3) |
| Version (DESCRIPTION) | 0.4.1.9002 |
| R / future / doFuture / pkgload | 4.6.2 / 1.80.0 / 1.4.0 / 1.5.3 |
| B (`nb_boots`) | 500 |
| Bootstrap seed | 8316951 |
| Parallel route | `parallel_args = list(plan = "multicore", workers = 100L, show_message = FALSE)`; the log reports "Bootstrap will use: 100 workers with 'multicore' plan" in both arms (128 cores) |
| Wall clock, arm A | 42.6 s in all (fit 2.5 s, bootstrap 37.3 s), well under the 30 min cap |
| Wall clock, arm B | 41.5 s in all (fit 2.5 s, bootstrap 36.5 s), well under the 30 min cap |
| Run | `Rscript --vanilla v6/app_v6.R <tree> <A|B>` from a terminal, arm A then arm B, each under `timeout 1800` (`v6/run_arms.sh`) |

Changes from the document (`quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd`), as
carried in `v6/app_v6.R`, a copy of `v5/app_v5.R`. Marks: `# A002:` (v1), `# A002v3:`,
`# A002v5:`, `# A002v6:`.

- v1: `library(forestsearch)` → `pkgload::load_all(<tree>)`; the `forestsearch()` call's
  `parallel_args` plan `"multisession"` → `"sequential"`; the bootstrap chunk stops after its
  timing line.
- v3: the tree path and arm label are command-line arguments; `mr_inference = FALSE` in
  `forestsearch()` (`mr_inference_args` left as is); v1's fit-only mode removed; bootstrap
  `seed = 8316951L`.
- v5: the v3 cluster block is removed. A fork gate is added (§2). The bootstrap takes the
  multicore `parallel_args` above.
- v6, the only edits to the v5 script: the bootstrap's `nb_boots` 20L → 500L, and the output
  path `v5/<arm>.rds` → `v6/<arm>.rds`. The document's `NB` variable is still 20L from v5; the
  bootstrap call does not read it, because it passes `nb_boots` as a literal. Left as is, since
  the task names only `nb_boots` and the output paths. The fork-gate log line still reads "A002v5".

## 2. Gate checklist

| Gate | Result | Evidence |
|---|---|---|
| Step 1: `git status --porcelain` empty on `fs-a002-check` | PASS | 0 lines |
| Step 1: task copied to `dev/tasks/` and committed | PASS | `cdffe0da dev: add A-002 magnitude task v6`; `cmp` with the `~/Downloads` copy identical |
| Step 3: fork gate, arm A | PASS | `parallelly::supportsMulticore()` TRUE; a `multicore` future (pid 523579 ≠ main 523453) returns FALSE for `get_dfpred(df.predict, sg.harm` in `deparse(forestsearch::forestsearch)` |
| Step 3: fork gate, arm B | PASS | `supportsMulticore()` TRUE; the future (pid 525643 ≠ main 525615) returns TRUE |
| Step 4: arm A ≤ 30 min | PASS | 42.6 s, exit 0 |
| Step 4: arm B ≤ 30 min | PASS | 41.5 s, exit 0 |
| Step 5: both fits select `{er <= 0} & {size <= 35}` | PASS | both `unname(sg.harm)` = `{er <= 0}`, `{size <= 35}` |
| Step 5: per replicate, same found / not-found, `M.1`, `M.2` | PASS | 500 of 500 replicates equal, with `boot_id` aligned |
| Step 7: `git diff --name-status master..fs-a002-check` only added files under `dev/tasks/`, `dev/reports/`; `git status --porcelain` empty | run after this report's commit | result in the reply |

## 3. HR-scale results (B = 500, seed 8316951)

Naive H0 is the Cox HR on the original data (the bootstrap's `H_estimates$H0` / `Hc_estimates$H0`).
The H2 columns are the bias-corrected estimate and its 95% CI (`H2`, `H2_lower`, `H2_upper`).
Each arm's H2 uses the 489 replicates of 500 that found a subgroup.

| Set | Quantity | Arm A | Arm B | B − A |
|---|---|---:|---:|---:|
| H  | naive H0 | 2.537 | 2.537 | 0.000 |
| H  | H2 estimate | 1.912 | 1.463 | −0.449 |
| H  | H2 95% CI | (0.918, 3.980) | (0.654, 3.271) | (−0.264, −0.709) |
| Hc | naive H0 | 0.608 | 0.608 | 0.000 |
| Hc | H2 estimate | 0.633 | 0.648 | +0.015 |
| Hc | H2 95% CI | (0.435, 0.920) | (0.393, 1.068) | (−0.043, +0.148) |

Unrounded B − A on the H2 estimate: H −0.44910080, Hc +0.01459718. Naive H0 is bit-identical
across the arms.

Where the difference arises: in all 489 replicates that found a subgroup, the per-replicate
bias-adjustment terms (`H_biasadj_1`, `H_biasadj_2`, `Hc_biasadj_1`, `Hc_biasadj_2`) differ between
the arms (|B − A| > 1e-12). Every replicate selected the same subgroup in both arms.

## 4. Replicates with a subgroup found

489 of 500 in each arm, the same 489. The 11 that found none, in both arms, are replicates 6, 19,
43, 285, 301, 318, 327, 386, 414, 456 and 474. Among the 489 (identical in both arms):

- 226 distinct subgroups were selected; 27 were single-factor.
- `M.1 = {er <= 0}` in 239.
- The original subgroup `{er <= 0} & {size <= 35}` was re-selected in 65 (13.0% of 500).
- The most frequent selections:

| Subgroup | Replicates |
|---|---:|
| {er <= 0} & {size <= 35} | 65 |
| {er <= 0} & {nodes <= 7} | 30 |
| {er <= 0} | 27 |
| {er <= 0} & !{size <= 20} | 17 |
| {er <= 0} & {pgr <= 45} | 9 |
| {er <= 0} & {pgr <= 26} | 8 |
| {er <= 0} & {pgr <= 48} | 7 |
| {pgr <= 0} & {size <= 35} | 7 |

## 5. Verdict

**MATERIAL**: |B − A| on the H2 estimate is 0.449 for H (≥ 0.05) and 0.015 for Hc.

## 6. Notes

- Neither log carries a warning or an error. v5's `nb_boots < 100` notice does not apply at B = 500.
- The Step 5 selection check compares `unname(sg.harm)`: `sg.harm` carries the names `q1.1` and
  `q27.1` (the v5 report, §6).
- The first 20 replicates reproduce v5's 20 exactly: the same found / not-found and the same
  `M.1`, `M.2` (`boot_id` 1–20). The bootstrap index matrix is drawn from the seed alone.
- Nothing under `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `quarto/`, `DESCRIPTION` or
  `NAMESPACE` changed; no install. The repo receives this report and the v6 task file only.
  Scratch: `/tmp/a002_2026-11-10/v6/` (`app_v6.R`, `run_arms.sh`, `compare_v6.R`,
  `compare_v6.rds`, `A.rds`, `B.rds`, `runA/armA.log`, `runB/armB.log`).
