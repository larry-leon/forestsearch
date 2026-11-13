# TASK — re-run the manuscript_rev1 simulations and analyses under forestsearch 0.5.0 (check only), 2026-11-12

Session name: `fs-rev1-rerun-check-0.5.0`

Host `pop-os`. Runs unattended overnight. Read this whole file before starting.

## Purpose

The `manuscript_rev1` results in fs-post-selection were produced on 2026-11-10 and 2026-11-11
with forestsearch 0.4.1.9003. The package is now 0.5.0 on `cran-prep`, after many code fixes.
Re-run the same simulations and analyses, unchanged, with 0.5.0 and report what moved.

This is a check of the package. **It changes no manuscript result:** nothing is written into
fs-post-selection, and no rev1 payload is replaced.

## Rules

- **Lean run.** The main session does everything. No verifier, polish or report-drafting
  agents. One short report, no detail file.
- **No `R CMD check`, no `devtools::check()`, no vignette build, no test suite.** None is part
  of this task.
- **fs-post-selection is read-only.** No branch, no commit, no file written, no file changed.
- **fs-glms-interpretable:** not touched.
- **forestsearch:** stay on `cran-prep`. The only commits are this task file and the report
  (path-limited commits; do not stage anything else). No file under `R/`, `tests/`, `man/`,
  `vignettes/` or `quarto/` changes.
- **Packages.** No install into the user or site library. The 0.5.0 build goes into
  `/tmp/chk050/lib` only.
- **No push, pull, merge or tag.**
- **Gates** stop on failure; they never stop to ask. There are only three (marked **Gate**).
  Everything else is recorded, not gated.
- **A failing item does not stop the night.** If one item errors, record the error text and
  the call that raised it, then go on to the next item.
- Do not re-verify what the rev1 reports already established (for example that per-simulation
  seeds depend on the simulation id only). Use their record.

## Sources to read first (all in fs-post-selection, read-only)

- `manuscript_rev1/REVISION_LOG.md` — the authoritative list of `re-done` payloads.
- `dev/reports/REPORT_a002_table_s4_2026-11-10.md` — Table S4 run.
- `dev/reports/REPORT_rev1_setup_inventory_2026-11-10.md` — generators and settings.
- `dev/reports/REPORT_rev1_reruns_2026-11-10.md` — S1, S2, A1, A2 runs.
- The Mac timing report from branch `rev1-mac-timing` (find it by name in `dev/reports/`).

## Items

Expected set, to be confirmed against `REVISION_LOG.md`:

| Item | Payload in `manuscript_rev1/_payloads/` | What it is |
|---|---|---|
| A1 | `analysis_gbsg_survival_multimethod_payload.rds` | GBSG survival application |
| A2 | `analysis_actg175_binary_multimethod_fixed_family_payload.rds` | ACTG175 binary application |
| S4 | `fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` | Table S4 simulation, n = 500, sims 1–500 |
| S1 | `fs_maxeffCons_fb_mr_m1_h10_knoise0_n1000_combined_1_500.rds` | n = 1000, h10, sims 1–500 |
| S2 | `fs_maxeffCons_fb_mr_m1_h15_knoise0_n1000_combined_1_400.rds` | n = 1000, h15, sims 1–400 |

- If `REVISION_LOG.md` marks any other payload `re-done`, add it after S2 and treat it the
  same way. If it marks one of the five as not re-done, say so and still run it.
- **Run order:** A1, A2, S4, S1, S2, then any added item.

## Step 0 — start

1. Copy this file into forestsearch `dev/tasks/` and commit it (path-limited).
2. Record, for each of the three repositories: branch, HEAD, `git status --porcelain`.
3. **Gate 1:** forestsearch is on `cran-prep` and `DESCRIPTION` at HEAD says
   `Version: 0.5.0`.

## Step 1 — package route

1. `git archive` forestsearch `cran-prep` HEAD into `/tmp/chk050/src` and `R CMD INSTALL` it
   into `/tmp/chk050/lib`.
2. Every run script's first R line: `.libPaths(c("/tmp/chk050/lib", .libPaths()))`.
3. **Gate 2:** in the main process and in one multisession worker (2 workers),
   `packageVersion("forestsearch")` is 0.5.0 and `find.package("forestsearch")` is under
   `/tmp/chk050/lib`.

## Step 2 — generators

For each item, use **exactly the generator that produced the payload now in
`manuscript_rev1/_payloads/`**: the path, the git blob and the edits that the rev1 report for
that item records (including any fallback it records as taken). Where two reports name
different blobs for the same generator, use the one from the report that produced the current
payload, and record both.

- Copy each generator (and, for the simulations, its combiner) to `/tmp/chk050/<item>/` and
  run it there.
- Allowed edits, each marked `# CHK050:` and listed in the report:
  - the `.libPaths()` line;
  - worker count at most 100 for the simulations (for A1 and A2, the generator's own formula
    capped at 100);
  - output paths, so that everything is written under `/tmp/chk050/<item>/`;
  - the same batch settings (`n_sims`, `sim_id_start`) the rev1 run used.
- Nothing else changes: seeds, DGM, sample sizes, `nb_boots`, `mr_draws`, thresholds and all
  analysis settings are those of the rev1 run.

## Step 3 — run

1. **A1, then A2.** Record each wall clock.
2. **Smoke.** Run sims 1–2 of S4, of S1 and of S2. An item whose smoke errors is recorded and
   skipped. Record the three smoke wall clocks.
3. **S4, S1, S2 in full,** each combined with its own combiner into a file with the original
   payload name. Record each wall clock, the worker count, and the wall clock the rev1 report
   recorded for the same item.
4. **Safety cap: 12 hours** from the start of Step 3, as a safety stop only. On the cap, stop
   the running item, combine the simulations that finished, list their ids, and go to Step 4
   with what completed.
5. Copy every new payload to `~/fs050_check/payloads/` under its original name (this directory
   is outside all repositories and is kept). Record each file's size and md5.

## Step 4 — compare (records, not gates)

Reference for each item: the payload in `manuscript_rev1/_payloads/` at fs-post-selection
HEAD. For A1 and A2, if `REVISION_LOG.md` says the current payload came from the Mac Studio,
compare **also** with the Linux payload committed by the rev1-reruns task (read it with
`git show <commit>:<path>` into scratch), so that a platform difference is not mistaken for a
code difference.

**Applications (A1, A2).** For every estimate field (naive, FB and MR estimates and limits,
for each method, for H and Hc) and every selected-subgroup definition:

- identical or not at full precision;
- the maximum absolute difference per field;
- every displayed table cell (as the rev1 fragment formats it), new beside reference, with a
  `same` / `changed` flag.

**Simulations (S4, S1, S2).** Pair by simulation id.

- Per simulation, over every column the two payloads share except timing columns:
  - non-numeric columns (detection flag, selected subgroup definition and the like): the
    number of simulations that differ;
  - numeric columns: the number of simulations that differ and the maximum absolute
    difference.
- Table level: for each rev1 fragment that reads the payload (the revision log lists them),
  run a scratch copy with only the payload path changed, and tabulate every displayed cell,
  new beside reference, with a `same` / `changed` flag.

**Classify each item** as exactly one of:

- `IDENTICAL` — every compared value identical;
- `SELECTION SAME` — detection and the selected subgroup agree in every simulation (or, for an
  application, for every method), and only estimates differ;
- `SELECTION DIFFERS` — detection or the selected subgroup differs somewhere. List the
  simulation ids (first 20 and the count) or the methods;
- `FAILED` — the item did not run; give the error.

For every item that is not `IDENTICAL`, state which fields moved (naive, oracle, FB, MR,
classification summaries). State facts only. No recommendation, and no view on which version
is right.

## Step 5 — report and close

1. Write `dev/reports/REPORT_rev1_rerun_check_0.5.0_2026-11-12.md` in forestsearch:
   - **top table:** one row per item — classification, simulations run, wall clock now, wall
     clock recorded by rev1, number of displayed cells changed, largest displayed change;
   - provenance: the three repositories' branch and HEAD, the package SHA, each generator's
     path and blob, every `# CHK050:` edit;
   - one comparison table per item (displayed cells, new beside reference);
   - the per-field difference counts;
   - any error text;
   - plain markdown, tables for results, bullets for the rest.
2. Commit the report on `cran-prep` (path-limited) and copy it to `~/Downloads`.
3. **Gate 3 (post-conditions):**
   - fs-post-selection and fs-glms-interpretable: branch, HEAD and `git status --porcelain`
     are the same as in Step 0;
   - forestsearch: since the Step 0 commit, the only files changed are this task file and the
     report;
   - nothing was pushed, merged or tagged.
4. Delete `/tmp/chk050/lib` and `/tmp/chk050/src`. Keep `~/fs050_check/payloads/`.
5. Print the top table and the Gate 3 result as the closing message.
