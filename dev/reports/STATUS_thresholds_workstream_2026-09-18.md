# STATUS — the thresholds workstream, 2026-09-18

Compiled from the **committed record only**: the five task documents in `dev/tasks/` and their reports in
`dev/reports/` and `quarto/simulations/actg175/binary_020/`, re-read for this document. No new analysis, no
fixes, no tasks proposed.

Branch `feature/glm-extension`, **nothing pushed**. Predecessors, for context, not part of the five:
`633e15b4`/`82b421e9` (minimum-events admission check) and `bb3492d2`/`ce5f1e84` (the criterion-and-defaults
audit that opened the workstream).

---

## 1. The five tasks, in order

| # | Task | Commit range | Class | Outcome |
|---|---|---|---|---|
| 1 | **Threshold naming docs** — `TASK_threshold_naming_docs_2026-09-18.md` | `216f3405` → `02491f4f` (edits `cc1714e0`, `0671525f`, `458cef48`) | documentation only | Named `c1`, `c2` and `p*` in every exported function carrying them; documented `threshold_config` and the two scales of "the screening threshold"; documented that the DINA frontier caps trim a report, not a search. `R/` carries roxygen changes and nothing else, asserted mechanically. |
| 2 | **Directive B — binary default estimand** — `TASK_binary_default_or_2026-09-18.md` | `6d0ce4df` → `9de4d0c5` (incl. `CLAUDE.md` policy `3427ca07`→`e68c5505`→`e500e5ff`) | **changes behaviour** + removes code | `forestsearch()` resolves an unset binary `effect_measure` to `"OR"`, not `"RD"`, at the live site (`5b1dac43`); the unreachable second resolution site deleted on proof (`1ba22ded`). Gate: 66 acceptance checks, confirmed failing pre-change (FAIL 30 / PASS 36 at `49893dec`). |
| 3 | **Threshold sync** — `TASK_threshold_sync_2026-09-18.md` | `d0708011` → `a389ebf8` | **fixes behaviour** + rider | Bootstrap replicates and CV folds now resolve the thresholds the parent fit resolved (`571b202c`, SECTION 2B-ii) — `missing()` on the legacy spellings was always `FALSE` inside a replay. Rider: the estimation entry points default binary to `"OR"` (`85847aa0`). Exposure gate was **not** empty: 11 hits, one unsafe class. 41 acceptance checks, 1.2 s. First `dev/reports/` report. |
| 4 | **Directive A — the threshold pair** — `TASK_directive_A_2026-09-18.md` | `ed34ce13` → `7cabca90` | **changes behaviour** ×2 | `c2 > c1` is an error on the forest-search consistency path (`822cd69a`); a silent `c2` derives `0.80 * c1` on the ratio scale (`0037ef6d`). Scope-gated to `subgroup_method = "consistency"` with survival `HR` or binary `"OR"`. The default pair (1.25, 1.0) is the rule's fixed point. 86 acceptance checks; 36 of 175 probe cells move, all in scope. |
| 5 | **Directive C — DINA guard and honest displays** — `TASK_directive_C_2026-09-18.md` | `bef180ce` → `4447a6c6` | behaviour ×3, messaging ×1, display ×1 | Five parts, five commits, each revertable alone (below). 182 acceptance checks, 5.9 s; all 175 probe cells byte-identical to baseline. |

A sixth, separately authorized item sits after the five: the one-off certification check,
`dev/reports/CHECK_ascran_2026-09-18.md` (`1daecaae`). §4 below.

---

## 2. The done-list, and where each item landed

| Item | Where it landed | Record |
|---|---|---|
| **Naming collision documented** — `c1` / `c2` / `p*` named in every exported function that carries them; `consistency.threshold` stated to be an *effect* threshold, not a proportion and not the consistency rate; `threshold_config` documented as a return component | `R/forestsearch_main.R` and the other exported blocks, `cc1714e0` · `0671525f` | `REPORT_threshold_docs_2026-09-18.md` §1-§3 |
| **Binary default `"OR"` at every resolution site** — the live `forestsearch()` site; the unreachable second site deleted; and the estimation entry points `make_effect_estimator()` and `.consistency_glm_pieces()` (the resolution `consistency_resample()` hands it) | `5b1dac43`, `1ba22ded` (B); `85847aa0` (sync rider); `787d138e` (A rider: `"OR"` first in `make_effect_estimator()`'s `match.arg` choices) | `REPORT_binary_default_or_2026-09-18.md` §2, §5; `REPORT_threshold_sync_2026-09-18.md` §6 |
| **Replicates resolve the parent fit's thresholds** — SECTION 2B-ii writes the resolved natural-scale values back into `effect.threshold` / `consistency.threshold` before `args_call_all` is captured, so bootstrap and CV replays cannot re-resolve | `R/forestsearch_main.R`, `571b202c` | `REPORT_threshold_sync_2026-09-18.md` §4, §5 |
| **`c2 > c1` stops, and a silent `c2` derives `0.80 * c1`** — HR and binary OR, consistency path only; disagreeing spellings of one threshold also stop; derivation on the ratio scale, announced once, carried into every replicate by the sync | `.fs_resolve_threshold_pair()` at `R/forestsearch_main.R:90`, SECTION 1A3; `822cd69a`, `0037ef6d`, roxygen `70aab091` | `REPORT_directive_A_2026-09-18.md` §3, §5 |
| **DINA refuses identity-scale estimands** — `RD` / `IRD` on a non-gaussian family stop at the `m_diff` derivation, on both routes (`.forestsearch_dina_select()` and the `use_dina` + `selected_only` screening branch); gaussian `MD` provably cannot fire | `.dina_assert_ratio_estimand()`, `R/forestsearch_helpers.R`; `c5bd3b5d` | `REPORT_directive_C_2026-09-18.md` Part 1 |
| **Frontier-key warning** — the seven `dina_frontier()` keys warn once per fit under `subgroup_method = "dina"`, naming every offending key | `R/forestsearch_helpers.R`, `.DINA_FRONTIER_KEYS`; `b8243432` | `REPORT_directive_C_2026-09-18.md` Part 2 |
| **`Inf` display caps** — `dina_frontier()`'s `max_per_covariate` / `max_subgroups` default to `Inf`; a finite cap that trims warns with kept-vs-available; muffled on the one internal call whose frontier is discarded | `R/dina_subgroup.R`, `R/forestsearch_main.R`; `31255426` | `REPORT_directive_C_2026-09-18.md` Part 3 |
| **Frontier print retitle** — the `details`-time table is titled as proposed single cuts, display only, beside the family counts | `R/forestsearch_helpers.R`; `c7e95160` | `REPORT_directive_C_2026-09-18.md` Part 4 |
| **c2 / p\* echo annotation under dina / grf** — seven echo sites annotated "not used on this path"; display only, no value changes | `forestsearch_main.R` (both banners), `interpret_search_config.R`, `forestsearch_methods.R`, `bootstrap_dofuture_main.R`, `forestsearch_cross_validation.R`, `fs_family_report.R`; `c8b475b4` | `REPORT_directive_C_2026-09-18.md` Part 5 |

---

## 3. Findings register — recorded, no task attached

Ids are this document's, assigned in report order; the docs-task and Directive C entries are unnumbered in
their own reports.

| Id | One sentence | Report |
|---|---|---|
| D1 | `R/subgroup_search.R:25` says `hr.threshold` is on the log scale for `HR`; on the survival path it is the natural scale, so a user passing `log(1.25)` would screen at `HR >= 1.13`. | docs §8.1 |
| D2 | Part 2's cited evidence, `quarto/gbsg/REPORT_dina_family_survey_2026-09-17.md` at `dba3073`, exists nowhere in the repository or its history. | docs §8.2 |
| D3 | `devtools::document()` errors on a pre-existing roxygen defect at `R/fs_bias_coverage.R:22` — a markdown code span beginning `` `r `` is parsed as inline R. | docs §8.3 |
| D4 | Two citations in the predecessor audit are off by one to two lines. | docs §8.4 |
| D5 | `dina_frontier()` is computed and discarded under the default screening mode; the wasted computation is behaviour, recorded without a fix. | docs §8.5 |
| D6 | `threshold_config` was an undocumented component of the `forestsearch()` return object. | docs §8.6 |
| F1 | Eight further `effect_measure` resolution sites carry two opposed binary defaults across layers; not reachable with `NULL` from `forestsearch()`. | B §9 |
| F2 | Exported entry points disagreed with `forestsearch()` on the binary default — dispositioned afterwards by the sync rider and A's `match.arg` rider. | B §9 |
| F3 | An explicit RD threshold becomes a log-ratio floor under `subgroup_method = "dina"` — **the defect Directive C Part 1 later refused.** | B §9 |
| F4 | Vignette build fails on an untouched tree ("Pandoc is required…") — **answered by the certification check, §4.** | B §9 |
| F5 | Pre-existing roxygen defect in `R/fs_bias_coverage.R` (same file as the non-ASCII warning). | B §9 |
| F6 | The truth layer `.fs_region_effect()` is already unconditionally OR for binary, ignoring `effect_measure`. | B §9 |
| S1 | The sync's exposure criterion is broader than the mechanism: only an *unsupplied* threshold on an identity-scale fit is exposed; 10 of 11 inventory hits were safe. | sync §8 |
| S2 | `consistency_resample_compare()` has no binary default to flip — it is survival-only, with neither an `outcome_type` nor an `effect_measure` formal. | sync §8 |
| S3 | `make_effect_estimator()`'s `match.arg` choice **order** was itself an `"RD"` default, unreachable at the time. | sync §8 |
| S4 | Five further `binary = "RD"` resolution sites remain: `frontier_cis.R:138`, `forestsearch_cross_validation.R:1449` and `:1498`, `plot_sg_glm_outcomes.R:147`. | sync §8 |
| S5 | The sync removes a per-replicate ratio-scale warning on RD / IRD; resolution was already equal there. | sync §8 |
| A1 | `mrct_region_sims()` already pairs the two thresholds (0.90 / 0.80) and passes both explicitly — the one place the pair was chosen rather than inherited. | A §6 |
| A2 | No committed caller is broken; the file-level scan behind this is weaker than the sync's parsed inventory and is reported as such. | A §6 |
| A3 | `quarto/simulations/gbsg_020/summary_grfmr.qmd:26` runs exactly the degenerate pair under `grf`, where c2 is inert — the motivating case Directive C Part 5 later annotated. | A §6 |
| A4 | The derivation makes an explicit default-valued `c1` newly chatty: one message per fit at the fixed point. | A §6 |
| A5 | Deviation from Directive A's Step 5 — the roxygen landed in its own commit rather than inside the behaviour commits. | A §6 |
| C1 | The `details = TRUE` frontier print on the dina path is now long: `Inf` caps show the whole non-dominated set (~200 rows in the micro-fit, formerly 10). | C, findings |
| C2 | The `use_dina` unoriented-floor residual is real and still present — `forestsearch_main.R:3054-3063` omits `tau_sign`, which `.forestsearch_dina_select()` computes and passes. | C, findings |
| C3 | `selection_rule` vocabularies disagree across the dina boundary: `forestsearch()` accepts `"hr"`/`"maxSG"`/…, `dina_subgroup()` accepts `"neighborhood"`/`"pareto"`/`"both"`. | C, findings |
| C4 | The 20 committed `dina_frontier()` call sites carry `# default 3L` / `# default 10L` comments that Part 3 made stale. | C, findings |

Two cross-references worth noting, both already closed by the record above: **F3** is what Directive C
Part 1 refuses, and **A3** is the case Directive C Part 5 annotates. **F2** was dispositioned by the sync
rider (`85847aa0`) and A's rider (`787d138e`). **F4** is answered in §4.

---

## 4. Standing policies, check reference sets, and the certification verdict

### Standing `CLAUDE.md` policies (unchanged; not edited by any of the five, nor by the check)

- **Verification per task is the task's own acceptance tests, and nothing else** — written into
  `tests/testthat/`, run with `devtools::load_all()` + `testthat::test_file()`, gated on those files being
  green. The full `devtools::test()` suite and `R CMD check` / `rcmdcheck` in any form are **not** run as
  part of a task. (Converged during Directive B: `3427ca07` → `e68c5505` → `e500e5ff`.)
- **The full suite and any CRAN check run only when Larry explicitly asks** — typically before an install
  or a release. The **certification surface** is `rcmdcheck::rcmdcheck(args = "--as-cran")`, which builds
  the PDF manual; `devtools::check()` is *not* certification, because its `manual = FALSE` default passes
  `--no-manual`. Either surface needs
  `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64` exported or on `PATH`.

### Check reference set 1 — the docs task (both runs, identical in text)

`REPORT_threshold_docs_2026-09-18.md` §6, at `216f3405` (pre-change) and `458cef48` (post-change):
**0 errors | 1 warning | 2 notes** — non-ASCII in `R/fs_bias_coverage.R`; 9 "no visible binding" in
`fs_plot_bias_coverage`; HTML manual validation skipped (no `tidy`).

### Check reference set 2 — Directive B

`REPORT_binary_default_or_2026-09-18.md` §7: **no `R CMD check` result exists for that task** — two
`--as-cran` runs and one reduced run were started and killed, producing no artifact; the docs reference set
stands untouched and B edited no file named in it. One completed data point is recorded there: a full
`devtools::test()` in **7.12 min** — 361 files, **FAIL 0 | ERROR 0 | PASS 5160 | SKIP 3 | WARN 32**.

### The one-off certification check (`1daecaae`, `dev/reports/CHECK_ascran_2026-09-18.md`)

Run at `4447a6c6` on the certification surface with the PDF manual and vignettes included; **10.3 min**
against a 90-minute timeout. Verdict: **1 error | 2 warnings | 2 notes**. The three reference findings
carry over identical in text; two findings are new, both naming files these tasks changed:

- **ERROR, `checking tests`** — `FAIL 14 | PASS 5097`. All 14 sit in the three probe-based test files added
  by the sync, A and C tasks, and are **test-harness artifacts, not package defects**: one is a `dev/` path
  read that does not exist in the built tarball (`test-directive-c.R:662`), and thirteen are the
  statement-lifting probe in `helper-threshold-sync.R` failing to match against the installed,
  byte-compiled package, which errors all 175 cells. Those three files therefore give no coverage under
  `R CMD check`. All five task-owned files are green under `devtools::load_all()` at the same tree.
- **WARNING, `checking Rd cross-references`** — `evaluate_combination_with_status.Rd` links
  `subgroup_search`; the alias carries a dot (`subgroup.search`). That Rd was regenerated by Directive B
  (`f391b816`).

**F4 answered: the vignettes built.** `creating vignettes ... OK`, `checking package vignettes ... OK`,
`checking re-building of vignette outputs ... OK`, and `checking PDF version of manual ... OK`. The cause
of F4 was the **`pandoc` PATH**, not stale or broken vignette content — confirming the `CLAUDE.md`
environment note and closing the B session's open question.
