# TASK — A-002 magnitude check, v6, 2026-11-10

Supersedes v1–v5 (v6 = v5 with B = 500). Same as v3 except the parallel route: v3 stopped because
`forestsearch_bootstrap_dofuture()` always replaces the caller's plan, and admits only
multisession / multicore / callr / sequential (`subgroup_consistency_helpers.R:43`). This task
uses `plan = "multicore"`: forked workers are copies of the main R process, so they run
whatever code `pkgload::load_all()` loaded there. No cluster, no install. Complete on its own.
No change to `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `quarto/`, `DESCRIPTION`,
`NAMESPACE`.

## Purpose

Same GBSG forest-search application, same seed, B = 500, run twice: arm A = repo code as is,
arm B = the v1 patched copy `/tmp/a002_2026-11-10/fs_patched` (row-wise `treat.recommend`).
**MATERIAL** if the H2 point estimate of H or of Hc differs by ≥ 0.05 on the HR scale;
otherwise **IMMATERIAL**.

## Rules

- Gates are stop-on-failure, never stop-to-ask. Do not push.
- Run from a terminal with `Rscript --vanilla` (not RStudio, where multicore is disabled).
- Scratch: `/tmp/a002_2026-11-10/v6/`. Nothing written into the repo except the report.
- Arms run one after the other, each ≤ 30 min wall clock; on the cap, stop and report.

## Steps

1. On `fs-a002-check`, `git status --porcelain` empty (**gate**). Copy this file to
   `dev/tasks/` and commit it on `fs-a002-check`.
2. Copy `/tmp/a002_2026-11-10/v5/app_v5.R` to `v6/app_v6.R` (it already carries the multicore route and fork gate; change only `nb_boots` to 500L and the output paths). For reference, v5 was made from v3 by removing the v3 cluster block
   (`makeCluster`, `clusterEvalQ`, `plan(cluster)`, its worker gate and the stop before the
   bootstrap call). Call the bootstrap with
   `nb_boots = 500L, seed = 8316951L,
   parallel_args = list(plan = "multicore", workers = 100L, show_message = FALSE)`.
   Keep everything else from v3 (`mr_inference = FALSE`, tree path as a command-line argument
   loaded with `pkgload::load_all()`, save to `v6/<arm>.rds`).
3. **Fork gate**, each arm, after `load_all()` and before the bootstrap:
   `parallelly::supportsMulticore()` is TRUE, and under
   `future::plan(future::multicore, workers = 2)` a future returns
   `grepl("get_dfpred(df.predict, sg.harm", paste(deparse(forestsearch::forestsearch), collapse = ""), fixed = TRUE)`
   equal to FALSE for arm A and TRUE for arm B. Otherwise stop.
4. Run arm A, then arm B.
5. **Gate:** both fits select `{er <= 0} & {size <= 35}`; per replicate, the same found /
   not-found indicator and the same `M.1`, `M.2` in both arms.
6. Write `dev/reports/REPORT_a002_magnitude_v6_2026-11-10.md`: header (SHA, B, seed, wall
   clock per arm, changes from the document), gate checklist, one HR-scale table (H and Hc:
   naive H0, H2 estimate, H2 95% CI; arm A, arm B, B − A), replicates with a subgroup found,
   one-line verdict. No recommendation.
7. Commit the report on `fs-a002-check`. **Gate:** `git diff --name-status
   master..fs-a002-check` lists only added files under `dev/tasks/` and `dev/reports/`;
   `git status --porcelain` empty. Copy the report to `~/Downloads`; print the table and
   verdict.
