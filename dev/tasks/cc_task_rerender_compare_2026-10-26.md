# CC task — re-run the workstream's analyses under the lean-fit build and compare with the committed (legacy) records

**File:** `dev/tasks/cc_task_rerender_compare_2026-10-26.md` · **Issued:** 2026-10-26 by chat, commissioned by Larry ("re-produce the analyses described during this workstream and compare with the legacy results — don't worry if they may shift, just run the analyses and that will tell us")
**Repo:** `forestsearch`. Larry's `R CMD check` of the branch is clean; he merges `feature/vi-deprecation-lean-fits` into `master` in GitHub Desktop **before** this task starts, so §1 expects `master` (the feature branch is accepted if the merge has not happened — say which). CC never fetches, pulls, pushes or merges.
**Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone). Two legacy vignette HTML files also arrive in `~/Downloads` (§2). Filenames may lose hyphens; locate by pattern.

---

## CATEGORY — measurement only. No `R/` change. No edit to any document source. The only commits are this document and the comparison record under `dev/reports/`.

- **What runs.** The four documents the workstream's fixtures came from, rendered under the installed lean-fit build (0.4.0.9000, A+F): the two package vignettes and the two committed ACTG175 application documents. Each is compared, output line by output line and on its key reported quantities, with its legacy record: the CRAN 0.4.0 rendered vignettes, and the committed rendered HTML of the two application documents.
- **What is expected, so the reader can tell signal from noise.** Survival and continuous: identical numbers (the fits are bit-identical and the deprecation changed nothing for them). Binary: the selected subgroup's membership identical, but the candidate pool (145 → 160 passing the effect screen, 136 → 137 in the consistency pool) and the clause order of its label changed with the deprecation; anything downstream of the pool — bootstrap re-selections, p-hat, field bounds — may therefore move a little. The main vignette (`forestsearch.Rmd`) passed `vi.grf.min = -0.2` under 0.4.0 and no longer does, so it may move too. The two application documents still pass `vi.grf.min = -0.2` (left as committed — not edited here): they emit one deprecation warning per `forestsearch()` call and behave as `NULL`. The vignette legacy records were rendered on Larry's Mac (R 4.5.2, survival 3.8-9); this run is Linux R 4.6.1, survival 3.8.12, so a last-digit difference there is platform, not code — say so where it occurs rather than chasing it. **Nothing in the comparison is a gate.** The gates are provenance, a successful render, and a clean record.
- **Working tree at the end.** The two application documents are rendered **in place** (that is where their records live) and left **uncommitted**, listed by path in the report, for Larry's decision in GitHub Desktop: commit = accept as the new records; discard = keep the legacy ones. Vignette renders go to the record directory only (they are rebuilt at package build time). Untracked by-products of the vignette renders are removed.
- **Compute:** unknown until the documents are read — §2 prints each document's bootstrap size, multiplier draws, worker count and `eval` flags before anything renders. Hard caps: 20 min per vignette, 3 h per application document, sequential, nothing else running. Expect the binary document to be the long one.
- **Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

1. `git rev-parse --abbrev-ref HEAD` (expected `master`; `feature/vi-deprecation-lean-fits` accepted) and `git rev-parse HEAD`. **GATE:** `git merge-base --is-ancestor 2e2b7437 HEAD` succeeds (the lean-fits report commit is in HEAD's history), i.e. the tree carries both the deprecation and the lean fits. `git status --porcelain` verbatim (pre-existing untracked dirt tolerated, never staged).
2. `packageDescription("forestsearch")[c("Version", "Built")]`: **GATE** = 0.4.0.9000 built `2026-10-26 17:45:59 UTC` (the A+F build). If Larry has reinstalled since, any 0.4.0.9000 build from this tree is accepted — record it. If it is the CRAN 0.4.0 build, stop and report.
3. `quarto --version`, `R.version.string`, `packageVersion("survival")`, `requireNamespace("speff2trial")` — record.
4. Copy this document to `dev/tasks/cc_task_rerender_compare_2026-10-26.md`, `git add` by path, commit alone.
5. Record directory: `dev/reports/rerender_2026-10-26/`.

## 2. Inventory — print before rendering anything

The four documents:

| id | document | legacy record | provenance |
|---|---|---|---|
| V1 | `vignettes/survival-post-selection.qmd` | CRAN 0.4.0 `inst/doc/survival-post-selection.html`, MD5 `69effc528fdda3c2100074243f216878` | the survival fixture (Stage 0 "V") with `mr_inference = TRUE` |
| V2 | `vignettes/forestsearch.Rmd` | CRAN 0.4.0 `inst/doc/forestsearch.html`, MD5 `1fc55176f5506245a3c7e91220f06508` | not a workstream fixture; the package's bootstrap/CV demonstration (NB = 2, 2 cores), rendered because it is cheap and shows the `use_grf = TRUE` + former `vi.grf.min = -0.2` path |
| C | the continuous ACTG175 document named in the `source` field of the committed `dev/glm-continuous-sims/stage0b_glm_2026-10-26/fixture_continuous.R` | its committed sibling `.html` | the continuous fixture |
| B | `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd` (also the `source` of `fixture_binary.R`) | its committed sibling `.html` | the binary fixture |

Legacy files: the two CRAN HTMLs are in `~/Downloads` as `CRAN_0.4.0_survival-post-selection.html` and `CRAN_0.4.0_forestsearch.html`; verify the MD5s above (**GATE**; if a file is missing, `utils::download.packages("forestsearch", destdir = <tmp>, type = "source", repos = "https://cloud.r-project.org")`, untar, take `inst/doc/*.html`, verify the same MD5s — this is the one permitted download). The application legacy HTMLs are taken from git, not the working tree: `git show HEAD:<path>.html > <record dir>/legacy_<id>.html` (confirm each is tracked; if a document has no committed HTML, say so and compare against the working-tree file if one exists, else outputs only).

For each document read the source and print to `inventory.txt`: every chunk name with its `eval`/`cache` options; `nb_boots`, `mr_inference_args$draws`, `fs.splits`, `parallel_args` (plan, workers), `n_cores`-type settings; every `vi.grf.min` / `max_n_confounders` line (with line numbers — expected: V1 none, V2 none now, B at its two known lines, C to be found); every file the document writes (`saveRDS`, `write.csv`, `ggsave`, `qs::qsave`, …) with its path; and the rendered output path(s). Commit `inventory.txt` with the record at the end (not now).

## 3. Render — sequential, logged, capped

Order V1, V2, C, B. Each under `timeout` with its cap, `nohup … > <record dir>/render_<id>.log 2>&1`, wall-clock recorded (`date +%s` before/after). Nothing else runs meanwhile.

- V1, V2: render **copies**. Copy the source (and any sibling files it reads — check the inventory) into `<record dir>/vignette_build/`, render there (`quarto render` for the `.qmd`; `Rscript -e 'rmarkdown::render("forestsearch.Rmd")'` for the `.Rmd`), keep the HTML as `new_V1.html` / `new_V2.html`. The `vignettes/` directory itself is not touched.
- C, B: render **in place** with `quarto render <path>` from the repo root (the documents are self-contained by Larry's rule; if a render fails on a missing input, stop and report — do not patch the document). Before each, snapshot `git status --porcelain` so the by-products can be listed afterwards. Copy the new HTML to `<record dir>/new_<id>.html` as well, and to `~/Downloads/new_<id>.html`.

A render that errors or hits its cap is reported with the log's last 40 lines and the task continues with the next document (the comparison for that document is then "not rendered").

## 4. Compare

Script `<record dir>/compare.R` (xml2 or rvest; both are on the machine if quarto is — check `requireNamespace`).

1. **Output extraction.** From each legacy and new HTML, extract the text of every code-output block, in document order: quarto HTML → `div.cell-output pre` (stdout and stderr); rmarkdown HTML → `<pre>` blocks whose text begins with the comment prefix `#>`. Unescape entities, strip tags, split into lines, trim. Drop lines that are timing or environment noise: those matching `fit time|Computation time|elapsed|seconds|minutes|Using .* cores|R version|Platform:|Running under|BLAS|LAPACK|time zone|tzcode|attached base|other attached|loaded via|^\[ *[0-9]+\] [a-zA-Z0-9._]+_[0-9]` (the sessionInfo package lists) and ISO dates. Write `outputs_<id>_legacy.txt` and `outputs_<id>_new.txt`, then a unified diff `diff_<id>.txt` and the counts: lines compared, lines differing.
2. **Key quantities.** From the same blocks, pull the values the document reports, by the labels its own result chunks print (read the source to know them). For every document: selected subgroup label; N; effect estimate (HR / MD / OR); Pcons; candidates evaluated and passed; `family_status`. Where present: the post-selection block (naive, de-biased, field one-sided lower, field-s one-sided upper, Bonferroni pair, IJ interval, p-hat, `n_family`); the recovery diagnostics (sens_H, ppv_H, containment q10/q50/q90, share); bootstrap results (`SG_CIs` rows, H/Hc estimates, bias-corrected values); CV agreement metrics; for B, every method's selected subgroup (FS, GRF, DINA, whatever the document reports) and their effects. Table per document: quantity, legacy, new, identical / changed (numeric: equal to the printed precision or not; note the printed precision). Write `keyvalues_<id>.csv`.
3. **Binary membership check (B only).** The document reports the FS subgroup's label; if its label clause order differs from the legacy but the clause set is the same, say "same clauses, reordered". If the document prints N and HR/OR for the subgroup, identical N and effect with a reordered label means identical membership (the expectation). If it saves a membership payload (`saveRDS`), compare `setequal()` of the ids.
4. **Warnings.** Count the deprecation warnings visible in each new render (`grep -c "deprecated" new_<id>.html`) and in the legacy.

## 5. Report — `dev/reports/REPORT_rerender_compare_2026-10-26.md`

1. Provenance (§1 values; the installed build; branch/HEAD).
2. Inventory summary: per document, bootstrap size / draws / workers / eval flags, the `vi.grf.min` lines found, files written, render wall-clock, legacy-render platform where known (V1/V2: Mac R 4.5.2; C/B: read it from the legacy HTML's sessionInfo if printed).
3. Per document: the line-diff counts and the key-quantity table (legacy → new, identical/changed). Then one paragraph each: what changed and why, sorted into (a) platform/last-digit, (b) deprecation-driven family change (binary, and V2's dropped `vi.grf.min`), (c) anything else — anything in (c) is the finding of this task and goes first in the verdict.
4. Deprecation-warning counts per document.
5. Working-tree state: the in-place rendered files (modified tracked, new untracked) listed by path with sizes, and the one-line instruction for Larry: commit to accept as the new records, or discard.
6. Deviations verbatim. Ten-line verdict.
7. Closing block (§6).

## 6. Commit and post-conditions — GATE

Commit by explicit path: `dev/tasks/cc_task_rerender_compare_2026-10-26.md` (§1, alone) and, at the end, `dev/reports/REPORT_rerender_compare_2026-10-26.md` plus everything under `dev/reports/rerender_2026-10-26/` **except** the HTML files (`legacy_*.html`, `new_*.html` and `vignette_build/` stay untracked — add a `.gitignore` inside the record directory listing `*.html` and `vignette_build/`, and commit that `.gitignore`). Before the commit `git diff --cached --name-only` must be a subset of those paths. Then print:

```
git diff --name-only <start HEAD> HEAD | grep -v -E '^dev/(tasks/cc_task_rerender_compare_2026-10-26\.md|reports/)'   # must print nothing
git diff --name-only <start HEAD> HEAD | grep -E '^(R/|tests/|man/|vignettes/|NAMESPACE|DESCRIPTION|NEWS)'               # must print nothing
git status --porcelain | grep -v '^?? '     # may print ONLY the in-place rendered application records (the .html of C and B); list them
git status --porcelain | grep '^?? '        # untracked: pre-existing dirt + render by-products of C and B only; vignette_build/ and the record HTMLs are ignored
```

If the first two print anything, stop and report. Finally copy `REPORT_rerender_compare_2026-10-26.md`, every `keyvalues_<id>.csv`, every `diff_<id>.txt` and `new_C.html`, `new_B.html` to `~/Downloads` and print `ls -la` of the copies.

## 7. Out of scope

- No edit to any document source (the two `vi.grf.min = -0.2` lines of B stay; they are Larry's call after he sees this comparison).
- No simulation campaign re-runs; no re-render of any other application document.
- No commit of rendered HTML anywhere; no `R/` change; no install; no `R CMD check`; no push or merge.
