# TASK — forestsearch 0.3.5.9000: R/ change review and the fs-glms-interpretable revalidation questions Q1–Q3 (read-only)

**Date:** 2026-09-19 · **Spec:** chat · **Executor:** CC on pop-os, from the repo root · **Approver:** Larry — pasting the kickoff is the go
**Kind:** read-only review. **`R/` is not touched. No compute.**

---

## 0. Limits

- Read-only over `R/`, `DESCRIPTION`, `NAMESPACE`, `tests/`, `vignettes/`, `quarto/`, payloads and data.
- This task writes exactly three things:
  1. this file → `dev/tasks/TASK_code_review_revalidation_2026-09-19.md` (first commit);
  2. the report → `dev/reports/REPORT_code_review_revalidation_2026-09-19.md` (last commit);
  3. a copy of the report → `~/Downloads/REPORT_code_review_revalidation_2026-09-19.md`.

  Scratch goes only in `mktemp -d` directories, removed at the end.
- Never done, whatever happens: editing any tracked file; changing any default; installing or reinstalling any package; touching, moving or regenerating any committed result; calling any data generator; any simulation or re-run; `R CMD check` / `rcmdcheck`; building vignettes; running testthat (test files are read, never run); regenerating any `current_status.md` (this task writes no simulation directory); `git fetch` / `pull` / `push` / `stash` / `checkout` / `reset` / `rebase` / `merge`.
- A fact a committed record already establishes is quoted (path:line, commit), not re-derived.
- No stops for gaps. Anything not found goes in OPEN ITEMS with what was searched, and the run continues. The one hard rule: no action outside the write list above; if a step seems to need one, record it as an OPEN ITEM instead.

## 1. Concurrency — the binary stage-2 run may be live in another CC session in this working tree

- Pin first: `PIN=$(git rev-parse HEAD)`. Read committed content at the pin (`git show $PIN:<path>`, `git grep -n '<pattern>' $PIN -- <paths>`). Commits other sessions land during this task are listed in the report (SHA, subject), not reviewed.
- Read committed bundles from the pin (`git show $PIN:<path> > "$T/<file>"`), never from the working tree, so a file another session is writing is never opened.
- Stage and commit only this task's two paths, by name: `git add <path>` then `git commit -m "<msg>" -- <path>` (pathspec-only, so anything another session has staged is left alone). Pre-existing untracked files are never staged.
- Before every `git add` / `git commit`: wait while `.git/index.lock` exists — poll every 15 s, up to 10 min; if it persists, record an OPEN ITEM, continue, and retry at closeout.
- Never write under `quarto/simulations/actg175/binary*/`; never signal or inspect a running process.

## 2. Step 0 — commit this document

Copy this file from `~/Downloads/` to `dev/tasks/TASK_code_review_revalidation_2026-09-19.md` and commit it alone (§1 rules). If it is not in `~/Downloads/`, say so and stop.

## 3. Part 0 — pin

1. `PIN`, branch, `git log -1 --format='%H %ci %s' $PIN`; `DESCRIPTION` `Version` at `PIN`.
2. Installed forestsearch: `Version`, `Built`, `Packaged` (`packageDescription("forestsearch")`); `R.version.string`; platform.
3. Installed namespace against `PIN`'s `R/`. Skip and record why if `PIN` has a `src/` directory or `pkgload` is not installed (do not install it). Otherwise: `T=$(mktemp -d)`; `git archive $PIN DESCRIPTION NAMESPACE R | tar -x -C "$T"`; in one R session under `timeout 1200`: capture every closure in `asNamespace("forestsearch")` (formals, and body after `utils::removeSource()`); `pkgload::load_all("$T", export_all = FALSE, quiet = TRUE)`; capture again; compare by name. Report counts identical / differing / present on one side only, with the names that differ.

## 4. Part 1 — the `R/` changes (job 1)

1. **Anchors.** C0 = the commit that introduced `Version: 0.3.5` in `DESCRIPTION`; C1 = the commit that introduced `Version: 0.3.5.9000`. Quote each (SHA, date, subject, its `DESCRIPTION` diff line).
2. **Inventory table** — one row per commit in `C0..$PIN` touching `R/`, oldest first:
   SHA · date · `R/` files · functions touched · class · governing task document · before/after C1 · in the §A inventory (Y/N).
   - Governing task document: the `dev/tasks/` file that instructed the commit; state the basis (commit message, or the task document committed at the start of the same session).
   - Class — one of, taking the strongest if several apply: changes the method · changes behaviour (say what: a default, a new error or warning, an estimate or selection path) · add-only (new function or formal; existing defaults byte-identical) · moves existing code · dead-code removal · messages/display only · docs/roxygen only.
3. **Summary lines.** One bullet per governing task document: what it changed in `R/`, its class, its commit range.
4. **NEWS.md.** Quote the development section's header and every entry (path:line at `PIN`); map each entry to its commit(s); list the `R/` commits in `C0..$PIN` with no NEWS entry.
5. **§A claims.** For each claim in §A below: CONFIRMED (path:line at `PIN`) · CONTRADICTED (quote what the source says) · NOT FOUND.
6. **Check coverage.** From the committed records (e.g. `dev/reports/CHECK_ascran_2026-09-18.md`, `dev/reports/STATUS_thresholds_workstream_2026-09-18.md`): the commit at which the last `--as-cran` run and the last full-suite run were made, and the `R/` commits after each. Run nothing.
7. **Static CRAN/style scan** of the lines added or changed in `R/` over `C0..$PIN` (`git diff C0 $PIN -- R/`). Report path:line per hit; fix nothing:
   - `set.seed()`, `RNGkind()` or assignment to `.Random.seed` without restoring the caller's RNG state;
   - `options()`, `par()`, `Sys.setenv()`, `setwd()` without `on.exit()` restoration;
   - `print()` / `cat()` outside print/format/summary methods and not behind a verbosity flag;
   - `<<-`; bare `T` / `F`; non-ASCII characters (`tools::showNonASCIIfile()` on each changed file);
   - `library()` / `require()` inside functions; `forestsearch:::` self-calls;
   - each function newly exported in the range (from the `NAMESPACE` diff): `@param` for every formal, `@return`, `@examples` (present, and how wrapped), `@export` consistent with `NAMESPACE`.

### §A — the handoff's inventory, to check against source

From the 2026-09-19 code-review handoff, which is not on this machine. Line numbers are approximate.

- **A1.** `R/forestsearch_main.R` SECTION 1A3: `.fs_resolve_threshold_pair()` (~line 90) — `c2 > c1` is a `stop()`; an unset c2 derives 0.80·c1; worked values 1.25 → 1.00, 1.00 → 0.80, 0.90 → 0.72.
- **A2.** `R/forestsearch_main.R` SECTION 2B-ii: bootstrap and CV replicates resolve the parent fit's thresholds; parent-fit resolution unchanged.
- **A3.** Binary `effect_measure` default RD → OR at the live resolution site in `forestsearch_main.R`; a second, unreachable site deleted.
- **A4.** `R/glm_effect_estimators.R` and `R/consistency_resample.R`: binary default `"OR"`; `match.arg` order `c("OR","RD","RR","IRR","IRD")`.
- **A5.** `.dina_assert_ratio_estimand(family, effect_measure)` in `R/forestsearch_helpers.R`, called at the `use_dina` derivation site; DINA's floor is `if (family == "gaussian") hr.threshold else log(hr.threshold)`.
- **A6.** `.DINA_FRONTIER_KEYS` and a frontier-key warning in `.forestsearch_dina_select()`.
- **A7.** `R/dina_subgroup.R`: `max_per_covariate` / `max_subgroups` default `Inf`; `.dina_warn_cap_trim()` with condition class `dina_frontier_cap_trim`.
- **A8.** `.fs_c2_inert_note()` applied to both config banners; `threshold_config` `@return` documented.
- **A9.** Existence condition at the effect-estimator boundary: OR needs all four cells ≥ 1 (control/treated × events/non-events); RR, IRR, HR need ≥ 1 event per arm; RD, IRD, MD untouched. A non-estimable candidate returns `estimate = NA`, `se = NA`, `converged = FALSE` and a `reason` naming the empty cell; it is counted in `filter_counts`, shown in `fs_family_report()`'s stage map, and cannot rank.
- **A10.** `converged` was computed by every estimator and discarded at `R/subgroup_search.R:879` — say what that line does at `PIN`.
- **A11.** RD's tier-3 raw-proportions fallback returns `converged = FALSE` and is untouched.
- **A12.** `R/fs_dgm_feasibility.R`, exported, commit `f9b794f6`: `fs_dgm_feasibility(dgm, n, n.min, d0.min, d1.min, n_rep, tolerance)`; draws through the DGM's own generator; does not change the RNG kind; GLM path only.
- **A13.** `n.min` = 60 strict; `d0.min` / `d1.min` = 10 / 10 non-strict, skipped for continuous and count — unchanged over `C0..$PIN`.
- **A14.** The oracle helper keeps the pooled 5/5 criterion and adds the four-cell condition — say where it lives (`R/` or a template).
- **A15.** Tests `helper-threshold-sync.R`, `test-threshold-sync.R`, `test-threshold-pair-directive-a.R`, `test-directive-c.R`, `test-binary-default-or-entry-points.R` exist (read, don't run).
- **A16.** "Two workstreams (thresholds; admission floors) produced every `R/` change in 0.3.5.9000" — test against the Part 1 table.

## 5. Part 2 — Q1: the four MR default changes

**Question (fs-glms-interpretable revalidation, 2026-09-19):** were `ci_method` `"ij"` → `"field"`, and the defaults of `field_complement`, `field_scale_complement` and `return_reselection`, each Larry's decision? Where is each recorded?

1. **Sites.** Every place in `R/` at `PIN` that supplies a value for each of the four when the caller does not: formal defaults (including a `match.arg` first element), pass-through values in `forestsearch()` and other wrappers, internal fallbacks (the chat record points at `.fs_apply_mr()` falling back to `"ij"` for DINA/GRF). One row per site: function · path:line · value · which identifier paths reach it.
2. **Commits.** For each site: the commit that introduced the current value (`git log -L`, or `-S` / `-G` on the expression) and the value before it (`git show <sha>^:<path>`) — SHA · date · subject. If that commit changed any other MR default, list it too, one line.
3. **Authorization.** For each commit: the task document that instructed it and any decision record it cites. Quote verbatim, path:line, the lines instructing the default change and any line recording approval or a decision by Larry. Leads from the chat record — verify each (`git log --all --format='%h %ci' -- '<path>'`):
   - `dev/tasks/PROPOSAL_complement_field_scale_2026-09-08_v2.md`
   - `dev/notes/REVIEW_E1_fields_2026-09-08.md`
   - `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` — refers to "the recorded adoption" of `TRUE` / `"selected"` / `TRUE`
   - `REVIEW_certification_2026-09-09.md` and `SUMMARY_survival_properties_2026-09-10.md` — reported on 2026-09-13 as absent from every reachable commit; confirm with `--all`.

   Then search `dev/tasks/`, `dev/notes/`, `dev/reports/` and `claude/` for the four formal names around each commit date. If a lead is absent from forestsearch and `~/Documents/GitHub/fs-glms-interpretable` exists on this machine, search it read-only too (`git -C <path> log --all -- '<name>'`, `git -C <path> grep`); never write there.
4. **NEWS.md.** The entry line(s) for each default (path:line at `PIN`) and the commit that added each.
5. **Verdict per default**, one of: RECORDED DECISION (quote + path:line) · INSTRUCTED, NO DECISION ON RECORD (quote the instruction) · NO RECORD FOUND. Facts only — Larry confirms each.

## 6. Part 3 — Q2: the scope of threshold rule (A)

**Question:** was the HR/OR-only scope Larry's decision? **The chat record:** decided 2026-09-17 — rule (A) (the `c2 > c1` stop and the 0.80·c1 derivation) for hazard ratio and binary/OR only, nothing else; RD, MD and IRR get no derivation and no new error; the stop fires on the FS consistency path only, because c2 is inert under DINA/GRF. Check it in-repo:

1. `dev/tasks/TASK_directive_A_2026-09-18.md`, `dev/reports/STATUS_thresholds_workstream_2026-09-18.md` and the Directive A report in `dev/reports/`: quote, path:line, every line stating the scope, the exclusions, and the decision's date or source.
2. `.fs_resolve_threshold_pair()` and its call site(s), quoted with path:line. From the code, fill one row per estimand — HR, OR, RR, RD, IRR, IRD, MD — with columns FS consistency · dina · grf; each cell gives validation (the `c2 > c1` stop: yes/no) and derivation (0.80·c1: yes/no), citing the governing line(s).
3. The lines in `tests/testthat/test-threshold-pair-directive-a.R` that pin out-of-scope behaviour — path:line; read, don't run.
4. One line: record and implementation agree, or differ (how).

## 7. Part 4 — Q3: the GRF factor-covariate fix

**Question:** did any committed GRF run made before the fix, whose results feed the companion manuscript (fs-glms-interpretable), pass a factor-class covariate? The fs-glms chat pointed at the identifiers campaign in `quarto/simulations/gbsg_020/` (catalog `current_status.md`) and named "the dinamr stem".

1. **The fix.** The commit(s) that made GRF's membership evaluator handle factor covariates on the forest's coding. Lead from the chat record: `.grf_evaluate_subgroup()` in `R/grf_subgroup_labels.R`, landed 2026-09-17 from the md-field-rerun workstream, records `REPORT_grf_dina_fixes_2026-09-16.md` and `REPORT_md_grf_2026-09-16.md` — verify. Give SHA · date · subject · governing task document · report, and quote the before/after hunk.
2. **Reach.** Callers of the fixed function at the fix's parent (`git grep -n '<fn>' <sha>^ -- R/`): which identifiers can reach it. Say which `gbsg_020` stems are GRF, and whether DINA reaches the fixed code (the record identifies `dinamr` as DINA, `grfmr` as GRF).
3. **Run set.** The GRF result artifacts that feed the companion manuscript, established from the committed briefs for fs-glms-interpretable and the `current_status.md` catalogs — quote the lines naming each campaign. Expect at least the `gbsg_020` GRF stem(s), `actg175` continuous `mdgrf` and `actg175` binary `orgrf`; include `idsweep`'s GRF rows, and any GRF result under `quarto/applications/`, if a brief or catalog names them as manuscript input. Nothing outside the set is assessed.
4. **Vintage** per run: pre-fix · post-fix · indeterminate — from the build that produced it (bundle meta `pkg_version` / `Built` / install metadata, or the run's committed record), else from ancestry of the commit that added its bundles relative to the fix commit. Say which basis.
5. **Committed exposure records first.** Where a committed record already settles a run's exposure — the chat record describes a 2026-09-16 read-only check that found `grfmr` not exposed — quote it (path:line, commit), state which cells it covers and on what basis, and do not re-derive it.
6. **Two instruments**, only for pre-fix or indeterminate runs that no record covers:
   - *Classes as passed.* Read the campaign's generator and template at the producing commit (`git show <sha>:<path>`) and trace from data generation to the GRF call, including any coercion step; quote the lines that fix each covariate's class. Where a committed DGM or super-population object exists, corroborate with `vapply(<frame>[<covariates>], class, "")` on the stored object. Call no generator.
   - *Stored payloads.* Per cell: count of `NA` in the membership and re-selection columns; kept-family size (`n_family` or equivalent: min / median / max); the covariates appearing in stored candidate or selected labels, against the covariate list. If the pre-fix `mdgrf` Gate-1 smoke artifacts are committed, report the same statistics for them as the positive control; if not, say so. One bundle at a time, in R under `timeout 3600`.
7. **Answer:** either no run in the set passed a factor before the fix — with the classes as evidence — or the list of runs and the outputs that could be affected. No re-run.

## 8. Report and closeout

1. Write `dev/reports/REPORT_code_review_revalidation_2026-09-19.md`:
   - **§1 Answers:** Q1 (one line per default), Q2, Q3 — one short answer each, with its source as path:line, a commit, or a record entry; then the Part 1 summary lines.
   - §2 Part 0. §3 Part 1 (table, NEWS map, §A results, check coverage, scan hits). §4 Q1 evidence. §5 Q2 evidence. §6 Q3 evidence.
   - §7 OPEN ITEMS. §8 Commits landed by other sessions during this task (SHA, subject).
   - Facts only: no recommendations, proposals or tasks.
2. Commit it (§1 rules). Copy it to `~/Downloads/`.
3. Remove every `mktemp` directory this task made.
4. Post-conditions, each printed PASS / FAIL:
   - each of this task's two commits touches exactly its own path (`git show --name-only --format= <sha>`);
   - `git diff --quiet $PIN -- R/ DESCRIPTION NAMESPACE tests/`;
   - the installed forestsearch `Built` equals the Part 0 value;
   - no `mktemp` directory from this task remains;
   - `cmp` of `~/Downloads/REPORT_code_review_revalidation_2026-09-19.md` against the committed report.
5. Closing message: the §1 answers verbatim, the post-conditions, the two commit SHAs. Do not push.
