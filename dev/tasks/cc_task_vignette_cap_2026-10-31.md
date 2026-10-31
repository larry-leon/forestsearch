# CC task — the primary vignette: the declaration at the recommended cap, and the release re-check

**File:** `dev/tasks/cc_task_vignette_cap_2026-10-31.md` · **Issued:** 2026-10-31 by chat, commissioned by Larry ("Yes, I would like to proceed", 10-31, after chat's comparison of `vignettes/forestsearch.qmd` with León and Anderson 2026b, arXiv v1).
**Repo:** `forestsearch` · **Base:** `master` at `748cec46` (part B's closing block, pushed). **New branch:** `feature/vignette-cap`.
**Transport:** `~/Downloads` → `dev/tasks/`: this document, and `vignette_cap_inputs_2026-10-31.zip` (md5 `0e4a71661e7b63f6209b397c1392343c`), unzipped to `dev/tasks/vignette_cap_2026-10-31/`. Filenames may arrive with hyphens or underscores altered; locate by pattern and restore the names.

**This document is the complete, authoritative task.** A gate failure stops the task at that gate: record it, leave uncommitted work uncommitted, write the report, and stop. Never ask, never work around, except under the one deviation clause of §3.3.

**What chat could not do.** Chat cannot run forestsearch. It parsed every chunk and inline expression of the new vignette and the gate script with R 4.3.3 (all parse), and checked the two cutoffs by arithmetic: the rounded-rate rule gives 1.6211 at p\* = 0.90 (the vignette prints the package's 1.621) and 2.4324 at p\* = 0.99 (the paper's 2.432). The build in §3 is the new chunk's first execution.

---

## ⚠ CATEGORY

**Documentation only. `R/` is not touched; no package behaviour changes.** One file of the tarball changes:

| file | in the tarball | what changes |
|---|---|---|
| `vignettes/forestsearch.qmd` | yes (rebuilt into `inst/doc` by `R CMD build`) | (1) the Inference paragraph names field-s as the studentized complement bound instead of listing it as a fourth output; (2) a new subsection, "The recommended cap", after the declaration table: the screen's cutoff at p\* = 0.99, the region's margin over it, and the diagnostic (8) at the capped cutoff, all read from the fit's existing draws (no refit), followed by the paper's S1.8.4 realized rates, quoted and labelled as quoted, and the reading of §6.2 of 2026b ("marginal on either reading"). The sentence that said the S1.8.4 simulation "is not run here" moves into that subsection. |

Diff against `master`: 63 lines added, 5 removed, all in `vignettes/forestsearch.qmd`. Table 1, the declaration table and every other chunk are unchanged; §3.2 gates that on the built HTML.

Also written: this document and the inputs under `dev/tasks/`; the record `dev/reports/vignette_cap_2026-10-31/` (`logs/`, `results/`); the report `dev/reports/REPORT_vignette_cap_2026-10-31.md`.

**Compute, its own line.** `R CMD build` once (it re-renders the vignette; part B measured 93 s; cap 15 min) and `R CMD check --as-cran --run-donttest` once on the built tarball (part B measured 294 s; cap 45 min). The new chunk adds no model fit. No simulation. Approved by Larry with the scope above. This is the release re-check that the check policy reserves for release work; run it once, at §4, not per edit.

---

## 0. Provenance, branch and first commit — GATE

1. No `git fetch`. `git rev-parse --abbrev-ref HEAD` must be `master`; `git rev-parse HEAD` must start `748cec46`; `git status --porcelain | grep -v '^?? '` must print nothing. Record `git status --porcelain | grep '^?? '` verbatim (**UNTRACKED_0**; untracked files are never staged). *GATE:* otherwise stop.
2. `md5sum vignettes/forestsearch.qmd` must be `037aad7ccfe9811bc59181356cc180a0` (the file chat edited). *GATE:* otherwise stop.
3. `git switch -c feature/vignette-cap`.
4. `md5sum` of the zip must equal the header value. Unzip into `dev/tasks/` so it creates `dev/tasks/vignette_cap_2026-10-31/`; `(cd dev/tasks/vignette_cap_2026-10-31 && md5sum -c MANIFEST.md5)` must report both files `OK`. *GATE:* otherwise stop.
5. Copy this document to `dev/tasks/cc_task_vignette_cap_2026-10-31.md`. `git add` it and `dev/tasks/vignette_cap_2026-10-31/` by path; commit them alone: `task: vignette — the declaration at the recommended cap — task document and inputs`.
6. Record `R --version | head -1`, `packageVersion()` of forestsearch, quarto (R), knitr, rmarkdown, xml2; the Quarto CLI path and version; `pandoc --version | head -1`. *GATE:* forestsearch is `0.4.1` and xml2 is installed. Otherwise stop; no installs.

Create `dev/reports/vignette_cap_2026-10-31/` with `logs/` and `results/`. Every log goes to `logs/`, every gate output to `results/`.

## 1. Place the vignette — GATE

1. `cp dev/tasks/vignette_cap_2026-10-31/forestsearch.qmd vignettes/forestsearch.qmd`; `md5sum` must be `f019e68cfa8d2538e889d834fc628d59`.
2. `git diff --numstat` must print exactly one line: `63	5	vignettes/forestsearch.qmd`. *GATE:* otherwise stop.
3. Parse check: `Rscript -e 'f <- knitr::purl("vignettes/forestsearch.qmd", output = tempfile(fileext = ".R"), quiet = TRUE); invisible(parse(f)); cat("parse OK\n")'` → `logs/parse.log`. *GATE:* `parse OK`.

## 2. Build the tarball from tracked files plus the placed vignette

1. `rm -rf /tmp/fs_cap && mkdir -p /tmp/fs_cap/src && git archive --format=tar HEAD | tar -x -C /tmp/fs_cap/src` (tracked files at commit 1, with no path quoting involved), then `cp vignettes/forestsearch.qmd /tmp/fs_cap/src/vignettes/forestsearch.qmd`.
2. From `/tmp/fs_cap`, run `R CMD build src` under `timeout 15m`, with the same environment part B §8.1 recorded (`RSTUDIO_PANDOC` exported; Quarto CLI on `PATH`, else `QUARTO_PATH`). Log to `logs/build.log`; record the wall-clock. *GATE:* exit 0 and `/tmp/fs_cap/forestsearch_0.4.1.tar.gz` exists.
3. Record its md5 (`results/tarball_md5.txt`). `tar -tzf` must list no path under `forestsearch/vignettes/articles/`, `forestsearch/dev/` or `forestsearch/docs/` (`results/tarball_paths.txt`). *GATE.*
4. The tarball's `forestsearch/vignettes/forestsearch.qmd` must have md5 `f019e68c…` (the placed file). *GATE.*

## 3. Gate the built vignette

1. Extract `forestsearch/inst/doc/forestsearch.html` to `/tmp/fs_cap/x/` (not into the check directory).
2. `Rscript dev/tasks/vignette_cap_2026-10-31/vignette_cap_gate.R /tmp/fs_cap/x/forestsearch/inst/doc/forestsearch.html dev/reports/vignette_cap_2026-10-31/results/gate_vignette.txt`. *GATE:* `RESULT: OK`. It checks:
   - Table 1, all ten rows × three identifiers, identical to the current vignette (which matches 2026b Table 1 to the digit);
   - the declaration table, all five rows, identical;
   - the cap table: cutoff `2.432`, region T `2.438`, a margin in (0, 0.010), uncapped cutoff `2.726`, and a diagnostic in (0, 0.667];
   - the new heading once; the new prose strings once each; the replaced wording absent; no `NaN`, `?@` or `#> Error`.
   It also **records, without gating**, the cap diagnostic against the paper's 0.186 (`MATCH` / `DIFFERS`). Put that line in the report's key values; a `DIFFERS` is a finding, not a stop.
3. **Deviation clause (the only one).** If the build in §2.2 fails inside the chunk `declaration-cap`, or §3.2 fails on the cap table alone, because of how the chunk reaches an object (for example `cal$c0$Mstar_c0` arriving as a vector or a data frame, or its columns not aligned with `tb`), you may make **up to two** edits to that chunk's code, then repeat §1.1's placement with the edited file, §1.2 with the new numstat recorded, and §2–§3. Each edit goes in the report verbatim (old line, new line, reason).
   - Not allowed: removing or loosening either `stopifnot()`; changing what is computed (the cutoff by the rounded-rate rule at p\* = 0.99; the diagnostic as the share of `Mstar_c0` draws at c0 = 0.75 above that cutoff); any edit to the prose, the quoted rates (0.432, 0.074, 0.022), any other chunk, or the gate script.
   - If Table 1 or the declaration table fails the gate, there is no deviation: stop.

## 4. Release re-check — GATE

1. From `/tmp/fs_cap`, `timeout 45m R CMD check --as-cran --run-donttest forestsearch_0.4.1.tar.gz`, same environment as §2.2. Log to `logs/check.log`; keep `forestsearch.Rcheck/00check.log`. *GATE:* `00check.log` ends `Status: OK`.
2. Record from `00check.log` the vignette re-build and examples timings beside `cran-comments.md`'s "about 75 s" / "about 40 s" (`results/gate4_check.txt`). If either is more than 50% above the stated figure, record an OPEN ITEM. It is not a stop, and `cran-comments.md` is not edited in this task.

## 5. Commits and post-conditions

1. Commit 2: `git add vignettes/forestsearch.qmd`; commit it alone: `docs: vignette — the declaration at the recommended cap (2026b §4.4, §6.2); field-s named as the studentized complement bound`.
2. Post-conditions, each recorded in `results/postconditions.txt`:
   - `git diff --name-only 748cec46 HEAD | grep -v '^dev/'` prints exactly `vignettes/forestsearch.qmd`;
   - `git show HEAD:vignettes/forestsearch.qmd | md5sum` equals the md5 of the tarball's `forestsearch/vignettes/forestsearch.qmd` (the tarball is built from the committed vignette);
   - `git ls-files 'vignettes/articles/*.html'` prints nothing.
3. Write the report (§6); commit 3: the report and `dev/reports/vignette_cap_2026-10-31/`, by path: `dev: report — vignette cap reading (2026-10-31)`.

No push, no fetch, no merge. `R/`, `man/`, `NEWS.md`, `README.md`, `cran-comments.md` and `_pkgdown.yml` are not touched.

## 6. Report

`dev/reports/REPORT_vignette_cap_2026-10-31.md`, in bullets:
- provenance and **UNTRACKED_0**;
- each gate and its result;
- the cap table's five values;
- the recorded diagnostic line (`MATCH` / `DIFFERS` against 0.186);
- build and check wall-clocks, and the check timings beside `cran-comments.md`;
- the tarball md5;
- the three commits;
- any §3.3 deviation, verbatim;
- **OPEN ITEMS**, findings only;
- a verdict block.

## 7. Closing block

Copy to `~/Downloads/vignette_cap_2026-10-31/`:
- the report;
- `forestsearch_0.4.1.tar.gz` (the tarball for win-builder);
- `/tmp/fs_cap/x/forestsearch/inst/doc/forestsearch.html` (the built vignette);
- `00check.log`;
- `dev/reports/vignette_cap_2026-10-31/results/`.

Then print the verdict block and stop.
