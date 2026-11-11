# TASK — A-002 magnitude check, v3, 2026-11-10

Supersedes v1 and v2. Reuses v1's completed Steps 0–3 (branch `fs-a002-check`, patched copy
`/tmp/a002_2026-11-10/fs_patched`, script `/tmp/a002_2026-11-10/app.R`). Complete on its own.
No change to `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `quarto/`, `DESCRIPTION`,
`NAMESPACE`. No package install.

## Purpose

Same GBSG forest-search application, same seed, B = 500, run twice: arm A = repo code as is,
arm B = the patched copy (row-wise `treat.recommend`). **MATERIAL** if the H2 point estimate
of H or of Hc differs by ≥ 0.05 on the HR scale; otherwise **IMMATERIAL**.

## Rules

- Gates are stop-on-failure, never stop-to-ask. Do not push.
- Scratch: `/tmp/a002_2026-11-10/v3/`. Nothing written into the repo except the report.
- Arms run one after the other, each ≤ 30 min wall clock; on the cap, stop and report.

## Steps

1. Check out `fs-a002-check`; `git status --porcelain` empty (**gate**). Copy this file to
   `dev/tasks/` and commit it on `fs-a002-check`.
2. Copy `/tmp/a002_2026-11-10/app.R` to `v3/app_v3.R`, with these marked changes:
   - `mr_inference = FALSE` in the `forestsearch()` call (MR runs after selection and does not
     enter the bootstrap H/Hc estimates).
   - The arm's tree path is a command-line argument, loaded with `pkgload::load_all()`.
   - Before the bootstrap: `cl <- parallel::makeCluster(100)`;
     `parallel::clusterEvalQ(cl, pkgload::load_all(<tree>, quiet = TRUE))`;
     `future::plan(future::cluster, workers = cl)`; stop the cluster on exit. Plain
     multisession workers would run the installed package, not the arm's code.
   - Call `forestsearch_bootstrap_dofuture()` with `nb_boots = 500L`, `seed = 8316951L`, and
     `parallel_args` that leave this cluster plan in place (read the bootstrap's plan
     handling to choose them; if every value resets the plan, stop and report the lines).
   - Save the bootstrap object to `v3/<arm>.rds`.
3. **Worker gate**, each arm, before the bootstrap: every worker's
   `deparse(forestsearch::forestsearch)` contains `get_dfpred(df.predict, sg.harm` for arm B
   and does not for arm A.
4. Run arm A, then arm B.
5. **Gate:** both fits select `{er <= 0} & {size <= 35}`; per replicate, the same found /
   not-found indicator and the same `M.1`, `M.2` in both arms.
6. Write `dev/reports/REPORT_a002_magnitude_v3_2026-11-10.md`: header (SHA, B, seed, wall
   clock per arm, the changes from the document), gate checklist, one HR-scale table (H and Hc:
   naive H0, H2 estimate, H2 95% CI; arm A, arm B, B − A), replicates with a subgroup found,
   and the one-line verdict. No recommendation.
7. Commit the report on `fs-a002-check`. **Gate:** `git diff --name-status
   master..fs-a002-check` lists only added files under `dev/tasks/` and `dev/reports/`;
   `git status --porcelain` empty. Copy the report to `~/Downloads`; print the table and
   verdict.
