# REPORT — A-002 magnitude check, v5, 2026-11-10

**Verdict: MATERIAL.** The H2 point estimate for H differs between the arms by −0.409 on the HR scale
(B − A), which is ≥ 0.05. For Hc the difference is +0.034.

## 1. Header

| Item | Value |
|---|---|
| `master` SHA | `217656e2b7142ad62932d12da91698752f174a96` |
| `fs-a002-check` before this report | `99c2afa9b5f4b8c24a1134b6bf3fe1e7505f0283` (v5 task commit); `R/` identical to `master` (0 diff lines) |
| Arm A tree | the repo checkout at `99c2afa9` (code as is) |
| Arm B tree | v1 patched copy `/tmp/a002_2026-11-10/fs_patched`; differs from the repo in `R/forestsearch_main.R` only (row-wise `treat.recommend`, v1 report §3) |
| Version (DESCRIPTION) | 0.4.1.9002 |
| R / future / doFuture / pkgload | 4.6.2 / 1.80.0 / 1.4.0 / 1.5.3 |
| B (`nb_boots`) | 20 |
| Bootstrap seed | 8316951 |
| Parallel route | `parallel_args = list(plan = "multicore", workers = 100L, show_message = FALSE)`; the log reports "Bootstrap will use: 100 workers with 'multicore' plan" in both arms (128 cores) |
| Wall clock, arm A | 13.7 s in all (fit 2.4 s, bootstrap 8.6 s), well under the 30 min cap |
| Wall clock, arm B | 13.7 s in all (fit 2.4 s, bootstrap 8.6 s), well under the 30 min cap |
| Run | `Rscript --vanilla v5/app_v5.R <tree> <A|B>` from a terminal, arm A then arm B, each under `timeout 1800` (`v5/run_arms.sh`) |

Changes from the document (`quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd`), as
carried in `v5/app_v5.R`, a copy of `v3/app_v3.R`. Marks: `# A002:` (v1), `# A002v3:`, `# A002v5:`.

- v1: `library(forestsearch)` → `pkgload::load_all(<tree>)`; the `forestsearch()` call's
  `parallel_args` plan `"multisession"` → `"sequential"`; the bootstrap chunk stops after its
  timing line.
- v3: the tree path and arm label are command-line arguments; `mr_inference = FALSE` in
  `forestsearch()` (`mr_inference_args` left as is); v1's fit-only mode removed; bootstrap
  `seed = 8316951L`.
- v5: the v3 cluster block is removed (`makeCluster`, `clusterEvalQ`, `plan(cluster)`, the worker
  gate, and the Step 2 stop). A fork gate is added (§2). The bootstrap call takes
  `nb_boots = 20L` and the multicore `parallel_args` above (`NB` 1000L in the document → 20L).
  Output goes to `v5/<arm>.rds`.

## 2. Gate checklist

| Gate | Result | Evidence |
|---|---|---|
| Step 1: `git status --porcelain` empty on `fs-a002-check` | PASS | 0 lines |
| Step 1: task copied to `dev/tasks/` and committed | PASS | `99c2afa9 dev: add A-002 magnitude task v5`; `cmp` with the `~/Downloads` copy identical |
| Step 3: fork gate, arm A | PASS | `parallelly::supportsMulticore()` TRUE; a `multicore` future (pid 518272 ≠ main 518237) returns FALSE for `get_dfpred(df.predict, sg.harm` in `deparse(forestsearch::forestsearch)` |
| Step 3: fork gate, arm B | PASS | `supportsMulticore()` TRUE; the future (pid 518596 ≠ main 518388) returns TRUE |
| Step 4: arm A ≤ 30 min | PASS | 13.7 s, exit 0 |
| Step 4: arm B ≤ 30 min | PASS | 13.7 s, exit 0 |
| Step 5: both fits select `{er <= 0} & {size <= 35}` | PASS | both `sg.harm` = `c(q1.1 = "{er <= 0}", q27.1 = "{size <= 35}")`, `identical()` across arms |
| Step 5: per replicate, same found / not-found, `M.1`, `M.2` | PASS | 20 of 20 replicates equal, with `boot_id` aligned (§4) |
| Step 7: `git diff --name-status master..fs-a002-check` only added files under `dev/tasks/`, `dev/reports/`; `git status --porcelain` empty | run after this report's commit | result in the reply |

## 3. HR-scale results (B = 20, seed 8316951)

Naive H0 is the Cox HR on the original data (the bootstrap's `H_estimates$H0` / `Hc_estimates$H0`).
The H2 columns are the bias-corrected estimate and its 95% CI (`H2`, `H2_lower`, `H2_upper`).
Each arm's H2 uses the 18 replicates of 20 that found a subgroup.

| Set | Quantity | Arm A | Arm B | B − A |
|---|---|---:|---:|---:|
| H  | naive H0 | 2.537 | 2.537 | 0.000 |
| H  | H2 estimate | 1.977 | 1.568 | −0.409 |
| H  | H2 95% CI | (0.006, 630.6) | (0.002, 1099.8) | (−0.004, +469.1) |
| Hc | naive H0 | 0.608 | 0.608 | 0.000 |
| Hc | H2 estimate | 0.624 | 0.658 | +0.034 |
| Hc | H2 95% CI | (0.307, 1.268) | (0.037, 11.577) | (−0.270, +10.309) |

Unrounded B − A on the H2 estimate: H −0.40892845, Hc +0.03350199. Naive H0 is bit-identical
across the arms.

Where the difference arises: in all 18 replicates that found a subgroup, the per-replicate
bias-adjustment terms (`H_biasadj_1`, `H_biasadj_2`, `Hc_biasadj_1`, `Hc_biasadj_2`) differ between
the arms (|B − A| > 1e-12). Every replicate selected the same subgroup in both arms.

## 4. Replicates with a subgroup found

18 of 20 in each arm. Replicates 6 and 19 found none in either arm. Selections (`M.1`, `M.2`) were
identical in both arms:

| boot_id | found | M.1 | M.2 |
|---:|:-:|---|---|
| 1 | 1 | !{er <= 115} | {pgr <= 79} |
| 2 | 1 | {er <= 0} | {pgr <= 8} |
| 3 | 1 | {er <= 0} | {size <= 35} |
| 4 | 1 | {er <= 0} | {pgr <= 26} |
| 5 | 1 | {er <= 0} | {pgr <= 46} |
| 6 | 0 | — | — |
| 7 | 1 | {er <= 0} | !{size <= 20} |
| 8 | 1 | {er <= 0} | {pgr <= 170} |
| 9 | 1 | {er <= 0} | — |
| 10 | 1 | {er <= 11} | !{meno} |
| 11 | 1 | {er <= 0} | {size <= 35} |
| 12 | 1 | {er <= 0} | — |
| 13 | 1 | !{er <= 194} | !{nodes <= 1} |
| 14 | 1 | {er <= 0} | {size <= 35} |
| 15 | 1 | !{er <= 107} | {pgr <= 119} |
| 16 | 1 | !{er <= 91} | {pgr <= 84} |
| 17 | 1 | {er <= 0} | {pgr <= 41} |
| 18 | 1 | {er <= 0} | {size <= 35} |
| 19 | 0 | — | — |
| 20 | 1 | {er <= 0} | {nodes <= 7} |

The original subgroup `{er <= 0} & {size <= 35}` was re-selected in replicates 3, 11, 14 and 18
(4 of 20).

## 5. Verdict

**MATERIAL**: |B − A| on the H2 estimate is 0.409 for H (≥ 0.05) and 0.034 for Hc.

## 6. Notes

- The task says it supersedes v1–v4. No v4 task reached this branch or `~/Downloads`; v5 is
  complete on its own and was carried out as written.
- The Step 5 selection check was first scripted as `identical(sort(sg.harm), target)`. That
  returned FALSE because `sg.harm` carries the names `q1.1`, `q27.1`. The values matched in both
  arms, and the two arms' `sg.harm` were `identical()` to each other. The check was corrected to
  compare `unname(sg.harm)`, and the gate passed. No run was repeated.
- Both arms emitted the same single warning, the package's own notice that `nb_boots < 100` may
  give unreliable confidence intervals. With B = 20 the H CIs are very wide in both arms.
- Nothing under `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `quarto/`, `DESCRIPTION` or
  `NAMESPACE` changed; no install. The repo receives this report and the v5 task file only.
  Scratch: `/tmp/a002_2026-11-10/v5/` (`app_v5.R`, `run_arms.sh`, `compare_v5.R`,
  `compare_v5.rds`, `A.rds`, `B.rds`, `runA/armA.log`, `runB/armB.log`).
