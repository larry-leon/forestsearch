# CC task — isolate the sources of the re-render differences (build ladder + bisect), and a fresh gbsg effMaxSG run — no legacy file is overwritten

**File:** `dev/tasks/cc_task_isolate_rerender_2026-10-26.md` · **Issued:** 2026-10-26 by chat, commissioned by Larry: isolate the reason(s) for the differences found by `dev/reports/REPORT_rerender_compare_2026-10-26.md`; re-run `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` as a **separate fresh run** compared against the legacy (current) document; **never overwrite any legacy file — in particular no `_payloads` file** (the gbsg payload is being incorporated in the manuscript).
**Repo:** `forestsearch`. Larry may merge `feature/vi-deprecation-lean-fits` into `master` before or after this task; §1 accepts either branch. CC never fetches, pulls, pushes or merges. The one permitted download is the CRAN 0.4.0 package (§3).
**Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone).

---

## CATEGORY — measurement only. No `R/` change, no document-source edit, no in-place render. Installs: several, all transient (§3), ending with the lean build reinstalled and verified.

- **Why.** The previous task compared each document with a legacy record that was rendered under an older package (C: 0.3.5, 2026-09-16, this machine; B: 0.2.0, 2026-08-16, Mac, 14 workers), so its differences mix package history, platform and this workstream's two changes. This task renders each document under a **ladder of builds on this machine** — the legacy record's own era, CRAN 0.4.0, the deprecation-only commit, the lean-fit HEAD — so every changed quantity is attributed to the rung where it first moves; for the continuous document (95 s a render) a `git bisect` then names the commit.
- **No legacy file is touched; every fresh run gets its own dated directory beside the legacy document.** C and B render into `quarto/applications/actg175/rerun_2026-10-26/<doc>_<rung>/`, G into `quarto/applications/gbsg/rerun_2026-10-26/G_<rung>/` — each holding that run's HTML, its `_payloads/`, the render log and the build it ran under. The document is rendered as a copy placed in that directory, so every relative write lands there. Before any render, the MD5 of every file (tracked or not) under `quarto/applications/gbsg/` and `quarto/applications/actg175/`, excluding the new `rerun_2026-10-26/` directories, is recorded; after every render and at the end that listing must be identical (§4, §7 GATES). The `lean` rung's HTML and every rung's payloads, logs and build stamps are committed; the diagnostic rungs' HTMLs stay on disk, gitignored. The previous task's in-place renders (three modified tracked files, one new untracked payload directory) are **reverted** in §0 — those were not accepted — after their copies are kept in the record.
- **Where the outputs live.** `quarto/applications/actg175/rerun_2026-10-26/{C,B}_<rung>/` and `quarto/applications/gbsg/rerun_2026-10-26/G_<rung>/`, one directory per run, beside the legacy document and never inside its `_payloads/`.
- **Compute:** C ≈ 2 min a render × 4 rungs + a bisect of ≈ 8 renders; B ≈ 7 min × 4; the gbsg document unknown until §2 prints its settings (cap 3 h a render; up to 4 renders, plus a bisect only if one render takes ≤ 10 min). Sequential, nothing else running during renders. Expect 1.5–3 h.
- **Unattended.** Gates stop the task; never ask, never work around.

---

## 0. Restore the working tree — first, before anything else

1. If the previous task's in-place outputs are still in the working tree, copy them into `dev/reports/rerender_2026-10-26/inplace_outputs/` as `new_C.html`, `new_B.html`, `new_C_payload.rds`, `new_B_payload.rds` (gitignored there; kept on disk). If Larry has already discarded them in GitHub Desktop, say so and move on — the record directory and `~/Downloads` already hold the HTMLs.
2. `git checkout -- quarto/applications/actg175/analysis_actg175_continuous_intervals.html quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.html "quarto/applications/actg175/_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_effmaxsg_payload.rds"` and remove the untracked `quarto/applications/actg175/_payloads/analysis_actg175_binary_multimethod_fixed_family/` directory (a by-product of the previous B render; no legacy file existed there — confirm with `git log --all -- <that path>` printing nothing before removing).
3. **GATE:** `git status --porcelain | grep -v '^?? '` prints nothing; the only untracked entries are pre-existing dirt (`quarto/extreme_subgroups/fixed_random/`) and anything the previous task's record directory holds.
4. **Snapshot:** `find quarto/applications/gbsg quarto/applications/actg175 -type f -not -path '*/rerun_2026-10-26/*' -exec md5sum {} + | sort -k2 > <record dir>/legacy_tree_md5_before.txt`. This file is the "never overwritten" proof; §4 and §7 compare against it (the `-not -path` clause is part of every later listing too, so the new rerun directories never enter the comparison).

## 1. Provenance — GATE

`git rev-parse --abbrev-ref HEAD` (master or the feature branch) and HEAD; **GATE** `git merge-base --is-ancestor ef058d81 HEAD` (the previous task's record is in history, so the lean fits and the deprecation are too). Installed `packageDescription("forestsearch")[c("Version","Built")]` — expected 0.4.0.9000 built `2026-10-26 17:45:59 UTC`; record whatever it is. `quarto --version`, R, survival, grf, speff2trial. Copy this document to `dev/tasks/`, commit alone. Record directory: `dev/reports/rerender_isolate_2026-10-26/`.

## 2. Inventory — print before any build or render

Documents: **C** `quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd`; **B** `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd`; **G** `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd`.

For each, to `inventory.txt`:
- every file it **reads** (`readRDS`, `read.csv`, `source`, `load`, `here::here`, `file.path` to a sibling, `_payloads/…`) and every file it **writes** (`saveRDS`, `write*`, `ggsave`, `qs::qsave`, …), with the literal path expressions; whether any path is absolute, `here::here()`-based or repo-root-based (those would escape a copy directory — see §4);
- `nb_boots`, `mr_inference_args$draws`/`R_out`/`R_in`, `fs.splits`, `parallel_args` (plan, workers), every `eval` flag, every `vi.grf.min`/`max_n_confounders` line;
- the legacy record: the sibling `.html` and every payload under its `_payloads/` (list with sizes and MD5), `git log -1 --format='%H %ad %s' --date=short -- <html>` (the era commit), the `forestsearch` version string in the legacy HTML (provenance line and/or sessionInfo) and its platform line, and `DESCRIPTION`'s `Version:` at the era commit (`git show <era>:DESCRIPTION | grep Version`).

## 3. The build ladder

Rungs, each installed into the normal library with `devtools::install(<path>, dependencies = FALSE, upgrade = FALSE, quick = TRUE)` (or `install.packages()` for CRAN), **one at a time, sequentially**, and verified by `packageDescription("forestsearch")[c("Version","Packaged","Built")]` plus `exists(".fs_logit_two_group", asNamespace("forestsearch"))` (TRUE only for the lean rung) printed into the log before each render:

| rung | source | how |
|---|---|---|
| `era_<doc>` | the era commit of that document's legacy HTML (§2) | `git worktree add ~/Documents/GitHub/fs_worktrees/era_<doc> <sha>` then install from it. If `DESCRIPTION` at that commit needs a dependency that is not installed, or the install fails, the rung is **skipped** with the error quoted — never install dependencies. |
| `cran040` | CRAN 0.4.0 | `install.packages("forestsearch", repos = "https://cloud.r-project.org")`; verify Version 0.4.0, Repository CRAN, Packaged 2026-10-10 |
| `depr` | commit `60a6496b` (deprecation closing block; no lean fits) | worktree + install |
| `lean` | the main tree at HEAD | `devtools::install(".")` from the repo — this is also the **final** state (§7) |

Worktrees live outside the repo and are removed at the end (`git worktree remove --force`, `git worktree prune`). `git worktree` is local; it is not a fetch.

Render order (installs interleaved): C on `era_C` → `cran040` → `depr` → `lean`; then the C bisect (§5) if warranted; then G on `era_G` → `cran040` → `depr` → `lean`; then B on `era_B` → `cran040` → `depr` → `lean`; then the G bisect if warranted (§5); then the final lean reinstall (§7). Where two rungs are the same commit (e.g. `era_G` = the 0.4.0 submission), render once and say so. **G's compute guard:** after G's first render, if it took more than 45 min, drop the `depr` rung for G (render `cran040` and `lean` only); if more than 90 min, render only `lean` and report G's ladder as partial — and in either case no G bisect.

## 4. Render procedure — a dedicated dated directory per run, beside the legacy document

Rerun roots: `quarto/applications/actg175/rerun_2026-10-26/` (C and B) and `quarto/applications/gbsg/rerun_2026-10-26/` (G). For document D and rung R the run directory is `<rerun root>/<D>_<R>/` (e.g. `quarto/applications/gbsg/rerun_2026-10-26/G_lean/`). Steps:

1. `mkdir -p <run dir>`; copy the `.qmd` there, plus every file §2 found it **reads** (preserving relative layout, as read-only inputs; a payload it writes and later reads back in the same render needs nothing copied). If §2 found an absolute, `here::here()` or repo-root-based **write** path, do **not** render that document: report it under deviations and continue with the others — rendering it would risk the legacy file.
2. Write `<run dir>/build.txt`: the rung name, the installed `packageDescription("forestsearch")[c("Version","Packaged","Built")]`, the worktree SHA where applicable, and the date.
3. `cd <run dir> && timeout <cap> quarto render <file>.qmd > render.log 2>&1`; wall-clock recorded in `build.txt`. Every relative write (`_payloads/…`) lands inside the run directory.
4. Remove the copied `.qmd` from the run directory afterwards (the source at HEAD is the record of what was rendered; `build.txt` names the SHA). Keep the HTML, `_payloads/`, `render.log`, `build.txt`.
5. Caps: C 20 min, B 30 min, G 3 h. A render that errors or hits its cap is reported with the log tail and the ladder continues.

Bisect steps (§5) are diagnostics, not records: they render under `~/fs_bisect_<D>/step_<n>/` outside the repo, and only their extracted values are committed.

After **every** render: `find quarto/applications/gbsg quarto/applications/actg175 -type f -not -path '*/rerun_2026-10-26/*' -exec md5sum {} + | sort -k2 | diff - <record dir>/legacy_tree_md5_before.txt` must be empty (**GATE** — stop at once if it is not, and report which file changed; do not restore it silently).

What is committed from the rerun roots (§7): every run directory's `_payloads/`, `render.log` and `build.txt`; the `<D>_lean/*.html`; and a `.gitignore` at each rerun root containing `*.html` and `!*_lean/*.html` (verify with `git check-ignore -v` that the lean HTML is tracked and the diagnostic HTMLs are ignored). The diagnostic rungs' HTMLs stay on disk.

## 5. Bisect — C certainly; G only if cheap

Run only if the ladder shows the pattern "era rung reproduces the legacy value, `cran040` does not" for the document's headline quantity — for C: Ĥ's `MR (IJ two-term)` point estimate β̃ (legacy 21.79; previous task's new value 22.05); for G: the analogous MR/field point estimate, chosen from §2's legacy HTML. If the era rung does **not** reproduce the legacy value, the legacy record is not reproducible on this machine under its own era build (platform/environment, or a payload-dependence) — report that and skip the bisect.

Procedure (in a dedicated worktree, so the main tree never moves): `git worktree add ~/Documents/GitHub/fs_worktrees/bisect_<D> <era sha>`; inside it `git bisect start a0025013 <era sha>` (bad = `a0025013`, master at the start of this workstream, whose `R/` files are byte-identical to CRAN 0.4.0 per the Stage-0b MD5 record; good = the era commit); `git bisect run <record dir>/bisect_<D>.sh`, where the script: installs the worktree (`dependencies = FALSE`; exit 125 on failure = skip), renders a fresh copy of D under a per-step directory with a 20-min cap (exit 125 on failure), extracts the headline quantity from the HTML (reusing `dev/reports/rerender_2026-10-26/extract_lib.R`), exits 0 if it equals the legacy value at printed precision and 1 if it equals the new value (125 if neither — print both). Record every step's SHA and value in `bisect_<D>_steps.txt`. The result: the first bad commit's SHA, date, message, and `git show --stat` file list, plus its `NEWS.md` lines if any. Cap the bisect at 60 min for C, 90 min for G, and run G's only when a single G render took ≤ 10 min.

## 6. Comparison and attribution

Reuse `dev/reports/rerender_2026-10-26/extract_lib.R` and the key-quantity definitions of that task's `compare.R` (add G's, read from its own result chunks). For each document, one ladder table `ladder_<D>.csv` / markdown: quantity × {legacy record, era, cran040, depr, lean}, with a final column naming the **first rung at which the value departs from the legacy value**, mapped to a cause:

- departs already at `era` → "not reproducible under the era build on this machine" (platform/environment or stale record; for B, where the legacy is a Mac render, this is the platform bucket);
- first departs at `cran040` → "package history between the era commit and 0.4.0" (and the bisect's commit, for C/G);
- first departs at `depr` → "deprecation (candidate family)";
- first departs at `lean` → "lean fits" — **this bucket is expected to be empty**; anything in it is the headline finding.

Also the line-diff counts between adjacent rungs. For G, the legacy payload under `quarto/applications/gbsg/_payloads/` is compared leaf by leaf with the fresh `lean` payload (`payload_G_compare.txt`) — read-only.

## 7. Finish — lean reinstall, cleanup, record, post-conditions — GATE

1. `devtools::install(".", dependencies = FALSE, upgrade = FALSE)` from the main tree; **GATE:** `packageVersion() == "0.4.0.9000"`, `exists(".fs_logit_two_group", asNamespace("forestsearch"))` is TRUE, and `Sys.setenv(NOT_CRAN = "true"); testthat::test_file("tests/testthat/test-lean-fits.R")` passes (the lean build is back, proven). Record `Built`.
2. Remove all worktrees; `git worktree list` shows only the main tree.
3. **GATE:** the `find … -not -path '*/rerun_2026-10-26/*' … md5sum` listing equals `legacy_tree_md5_before.txt` exactly.
4. Report `dev/reports/REPORT_rerender_isolate_2026-10-26.md`: provenance; inventory summary (settings, era commits and versions, read/write paths, render wall-clocks per rung); the ladder tables with attribution; the bisect result(s) with the commit shown; G's fresh-vs-legacy comparison (quantities and payload); deviations verbatim; a ten-line verdict whose first line answers "did anything move at the `lean` rung?"; closing block with:
   ```
   git diff --name-only <start HEAD> HEAD | grep -v -E '^(dev/(tasks/cc_task_isolate_rerender_2026-10-26\.md|reports/rerender_isolate_2026-10-26/)|quarto/applications/(actg175|gbsg)/rerun_2026-10-26/)'   # must print nothing
   git status --porcelain | grep -v '^?? '                                                      # must print nothing
   find quarto/applications/gbsg quarto/applications/actg175 -type f -not -path '*/rerun_2026-10-26/*' -exec md5sum {} + | sort -k2 | diff - dev/reports/rerender_isolate_2026-10-26/legacy_tree_md5_before.txt && echo LEGACY-UNCHANGED
   git worktree list
   Rscript -e 'cat(as.character(packageVersion("forestsearch")), exists(".fs_logit_two_group", asNamespace("forestsearch")), "\n")'
   ```
5. Commit by explicit path, in two commits: (a) the rerun roots — each run directory's `_payloads/`, `render.log`, `build.txt`, the `<D>_lean/*.html`, and the two `.gitignore` files ("rerun 2026-10-26: fresh runs under the build ladder, legacy untouched"); (b) the record — the report, `inventory.txt`, `legacy_tree_md5_before.txt`, the ladder CSVs, `bisect_*.txt`/`.sh`, `payload_G_compare.txt`, logs. Before each, `git diff --cached --name-only` must be a subset of the allowed paths (§7 closing block, first command).
6. Copy the report, the ladder CSVs, the bisect step files and every `<D>_lean/*.html` to `~/Downloads`; print `ls -la`.

## 8. Out of scope

- No source edits (B keeps its two `vi.grf.min = -0.2` lines: under `cran040` they reproduce the screened family, under `depr`/`lean` the deprecated one — that contrast is the point).
- No in-place render, no record acceptance, no payload update — Larry decides after reading the ladder.
- No bisect for B (seven-minute multi-method renders; its remaining differences are the Mac-vs-Linux and 0.2.0-history buckets, which the ladder brackets).
- No push, merge, check or `R/` change.
