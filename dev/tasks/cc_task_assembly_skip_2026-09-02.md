# CC task — the assembly skip, plus two read-only settlements: bootstrap outer parallelism and `stop_threshold`

**File:** `dev/tasks/cc_task_assembly_skip_2026-09-02.md` · **Issued:** 2026-09-02 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1). Hyphen-stripped stem `cc_task_assembly_skip_20260902.md` is the same file.
**Reads (in-repo):** `dev/glm-continuous-sims/REPORT_medians_on_survivors_2026-09-02.md` (the 0.3.4 state and residual-cost statement), `dev/glm-continuous-sims/medians_baseline_2026-09-02.R` (the battery to transplant), `dev/glm-continuous-sims/medians_postchange_2026-09-02.rds` (**the pre-change baseline for this task** — committed at `ef021f50`, produced by that battery at 0.3.4), `dev/glm-continuous-sims/medians_profile_F2_after_2026-09-02.Rprof` (the committed 0.3.4 profile that sizes this task's prize), `REPORT_bootstrap_reprofile_2026-09-01.md` (bucket definitions; the F3 anchor fixture and its `stop_threshold` finding).

**Why.** After 0.3.4 the gbsg per-candidate cost is coxph internals plus the per-candidate `df.x` assembly. The 09-01 profile put assembly at 25.9% of the fit bucket (≈1.03 s per gbsg call), and the medians move removed `df.x`'s only other consumer — the skip is now a pure substitution. **The 09-01 sub-split also misattributed `summary.coxph` by 8×**, so this task re-sizes the prize from the committed 0.3.4 Rprof **before** implementing, and the report states realized against both numbers.

Two recorded questions are settled here by measurement, no code change either way: whether `forestsearch_bootstrap_dofuture()` parallelizes replicates reproducibly (the difference between ~40–90 min sequential projections and a few minutes of wall on this machine), and whether the anchor's default first-passer `stop_threshold` selects identically to the drivers' full evaluation under `maxeffCons` (the parked semantics note from 09-01).

**The governing rule:** a performance change that alters a result is a defect regardless of its speed. Guard: `identical()` against the committed 0.3.4 baseline — no tolerance. Committed work is not re-run: **there is no Stage A in this task.**

---

## ⚠ CATEGORY — this task edits `R/`

- **Moves existing code:** the `df.x` assembly inside `fit_cox_for_subgroup()` relocates into the adjusted arm — the only arm that still needs a frame (covariate columns).
- **Adds:** in the unadjusted arm, three named local vectors (`Y <- yy[idx]`, `E <- dd[idx]`, `Treat <- tt[idx]`) and the same `coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)` call resolving them from the calling environment — same variable names, so the coefficient name, the `"Treat" %in% rownames` check, and the 0.3.3 CI arithmetic are untouched downstream.
- **Changes behaviour:** none on any output surface. Same values, same row order, same solver, same names.
- **Changes the method:** no.

Files expected to change: `R/subgroup_search.R` only, plus `DESCRIPTION` (→ 0.3.5) and `NEWS.md`. The medians site (`med_df` at the post-screen point) is **untouched**. If the restructure requires touching any other `R/` file, **STOP** and report why.

**§6 and §7 are read-only settlements** — measurement and source quotes, no edits of any kind, whatever they find.

**Compute:** Stage-C battery ≈5 min; one profiled F2 call ≈1 min; tests + `R CMD check` ≈10 min; §6 bootstrap pair ≈3 min; §7 anchor pair ≈20 s. **Estimate 30–50 minutes wall.** No Stage-A run, no renders.

**Unattended.** Gates stop the task; never ask, never work around. On any equality failure: case, component, first differing values verbatim, both sides — then STOP, tree committed and documented. Provenance gate dirt-tolerant; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block; `packageVersion("forestsearch")` expect **0.3.4**. *GATE:* branch `feature/glm-extension`; the medians-task commits `ff2f9ca2`, `4c923dfd`, `3921ffdd`, `ef021f50` in the log (note the fourth was amended — `ef021f50` is the hash that exists). Copy this document into `dev/tasks/`, commit **alone**, report filenames and hash. Record the `vi.grf.min` default in force.

## 2. Verify at HEAD, and size the prize — GATE

1. **`fit_cox_for_subgroup()` at HEAD, full body:** the assembly lines (quote them — constructor, columns, any coercion), both arms as they now stand (the 0.3.4 state: fast CI arithmetic unadjusted, `summary.coxph` adjusted, medians gone), and **every remaining consumer of `df.x`**. *GATE:* consumers are the two `coxph()` calls and nothing else. If anything else reads `df.x`, STOP and report it.
2. **Assembly semantics:** confirm from the quoted lines that assembly is pure column binding of the already-subset vectors — no NA handling, no sort, no filtering, no type coercion that changes values. If it does anything beyond binding, STOP and report — the substitution below assumes value-and-order identity.
3. **Upstream NA guarantee:** quote where `yy`, `dd`, `tt` are cleaned (or shown NA-free) before the candidate loop, so environment-resolved `coxph()` (default `na.action`) cannot diverge from the frame-based fit. If no guarantee exists in source, STOP and report.
4. **Size the prize from the committed record — before any edit.** Re-bucket `medians_profile_F2_after_2026-09-02.Rprof` (0.3.4, committed at `ef021f50`) with the established `BUCKETS` definitions, sub-classifying the fit bucket as the 09-01 report did (solve / assembly-and-subset / model.frame / vcov-summary / formula / wrapper). Report the assembly-and-subset seconds at 0.3.4. **This number is this task's predicted recovery**; the 09-01-derived ≈1.03 s is quoted beside it. Given the `summary.coxph` misattribution precedent, if the re-bucketed assembly cost is under ≈0.3 s, say so prominently in §9's report — the change still lands (it is harmless and simplifying), but the record must not claim a prize the profile does not show.
5. **Transplant the battery — author nothing.** Create `dev/glm-continuous-sims/assembly_battery_2026-09-02.R` as a **copy of the committed `medians_baseline_2026-09-02.R`**, changing only the named output lines: the `.rds`/log filenames (→ `assembly_postchange_2026-09-02.*`) and the env-var name if one is needed. Same for the profiler: `assembly_profile_2026-09-02.R` as a copy of `medians_profile_2026-09-02.R` with only the output-name lines changed. *GATE:* show the `diff` of each transplant against its committed original; anything beyond the named output lines is a stop. Commit the two transplanted scripts as the second commit.

**The pre-change baseline is not re-run:** `medians_postchange_2026-09-02.rds` is the 0.3.4 output of the identical fixture battery — provenance stated in the report (commit `ef021f50`, script, version, self-consistency 3/3 as recorded).

## 3. Implement

In `fit_cox_for_subgroup()`:

1. Compute `idx <- id.x == 1` once (if not already present in this form).
2. **Adjusted arm** (`length(adj) > 0L`): assemble `df.x` exactly as the current lines do — the assembly code moves here verbatim — and fit as now (`summary.coxph` path unchanged).
3. **Unadjusted arm:** no frame. `Y <- yy[idx]; E <- dd[idx]; Treat <- tt[idx]`, then the existing 0.3.3 fast block with the identical formula `survival::Surv(Y, E) ~ Treat` and **no `data =` argument** — `coxph()` resolves the three locals; the coefficient is named `Treat` exactly as before; the CI arithmetic, `try()` contract, and dimnames lines are byte-for-byte untouched.
4. Everything downstream of the two arms — the `trow` extraction, the `NA` medians return — unchanged.

**Housekeeping, same commit:** `DESCRIPTION` → **0.3.5**; `NEWS.md` entry under 0.3.5 (performance: unadjusted survival fits run on subset vectors, the per-candidate frame now built only for adjusted fits; a relocation plus substitution verified bit-identical); roxygen touched only if `fit_cox_for_subgroup()`'s docs mention the frame; `devtools::document()` clean; tidyverse style. Commit as the third commit; reinstall.

## 4. Equality gates — against the committed baseline

Run the transplanted battery at 0.3.5 → `assembly_postchange_2026-09-02.rds`. Compare against **`medians_postchange_2026-09-02.rds`** with `identical()` under the battery's own volatile-field exclusion (unchanged by the transplant):

- *GATE E-1 (F2, F5):* every retained component identical — F5 exercises the relocated adjusted-arm assembly.
- *GATE E-2 (F1, F3):* continuous guards identical — this change must be invisible there.
- *GATE E-3 (E-ties, E-named, E-finite, E-zero, E-degen):* identical per case — E-degen in particular re-proves the degenerate-fit warning/`try()` behaviour on the environment-resolved call.
- *GATE E-4:* the 20-replicate gbsg bootstrap payload identical.

The battery's within-stage self-consistency checks must pass as they did at both prior stages.

## 5. Substitution proof and realized recovery

- **Constructor count:** identify from §2.1 the assembly's constructor call; `trace()`-count it during one F2 call at 0.3.4 (installed from the log's prior state is unnecessary — cite the structural fact instead: one construction per fitted candidate, 1,410) and at 0.3.5 (expected: **0** — F2 is unadjusted; the 121 `med_df` `data.frame()` constructions at the medians site are a separate constructor and are counted and attributed separately). One F5 line: adjusted fits still construct.
- **Profile:** one F2 call at 0.3.5 via the transplanted profiler; bucket table beside the §2.4 re-bucketing of the committed 0.3.4 profile — assembly-and-subset seconds before → after, fit bucket, total wall (0.3.4 measured 4.82 s).
- **Bootstrap:** the battery's gbsg 20-replicate wall at 0.3.5 against the recorded 5.49 s/replicate; the survival B = 1000 line (91.5 min → measured projection).
- One line: realized against **both** predictions — §2.4's re-bucketed number and the 09-01-derived ≈1.03 s.

## 6. Read-only settlement A — can the bootstrap parallelize across replicates reproducibly?

No code change, whatever is found.

1. **Source:** quote from `forestsearch_bootstrap_dofuture()` (and `bootstrap_analysis_dofuture.R` as needed) how replicates are dispatched (the loop/backend construct), how the per-replicate RNG is established — specifically whether each replicate's seed is **pre-assigned deterministically from `seedit` and the replicate index** (or an equivalent parallel-safe design: `future.seed = TRUE` L'Ecuyer streams, `%dorng%`, etc.) — and where the caller's plan enters.
2. **Measurement:** the continuous F1 bootstrap configuration, `nb_boots = 40`, seed 8316952: once sequential, once under `plan(multisession, workers = 20)` (or the wrapper's own argument for the outer plan — use the mechanism the source shows; if the source shows the outer level cannot take a parallel plan at all, say so and skip the run). Compare payloads with `identical()` under the battery's exclusion; report both walls.
3. **Verdict for the record**, one of: **(i)** payloads identical — outer parallelism is reproducible; report the speedup and the implied B = 1000 wall at stated worker counts (this reframes every sequential projection in the record); **(ii)** payloads differ — quote the mechanism (which seeds/streams diverged) and connect it to the recorded seeding-unification item; **(iii)** structurally unavailable — quote why. No fix, no workaround, no plan-default change.

## 7. Read-only settlement B — is first-passer `stop_threshold` sound for `maxeffCons`?

No code change, whatever is found.

1. Run the F3 anchor fixture twice at 0.3.5: as-is (default `stop_threshold = pconsistency.threshold`, 1 candidate evaluated) and with `stop_threshold = NULL` (all evaluated). Compare the selected subgroup and `out.found` (the full-evaluation run will carry more consistency rows — compare the selection and the selected row's values, and say what extra the full run carries).
2. Quote from source the candidate ordering in force at the consistency loop under `sg_focus = "maxeffCons"` — the sort that would make the first passer the maximum-effect consistent candidate — and state whether the code establishes that ordering **before** the early-stop check.
3. Verdict for the record: soundness confirmed (matching winner + ordering quoted), or refuted (differing winner — report both selections verbatim; this would be a finding about the applied anchor's committed results and goes to chat before anything else happens), or ordering-unclear (quote what was found). The bootstrap-level early-stop question is **decided by Larry after this lands** — nothing here enables it.

## 8. Package health gates

- `devtools::test()` on 0.3.5 source: 0 failures; warning parity (31); `test-search-reproducibility.R` unmodified and passing. *GATE.*
- `devtools::check(document = FALSE, vignettes = FALSE, manual = FALSE)`: **0/0/0**, matching the recorded 0.3.4 result. *GATE.*
- The two standing guards (`fidelity_fs_oc_predict_2026-08-28.R`, `prerefactor_reference_2026-08-29.R` check) run on 0.3.5 source: both PASS bit-identical, as after every `R/` edit in this workstream.

## 9. Out of scope — explicit

- **Medians removal** (NA-out + display recompute): a decision pending with Larry; hygiene not speed (~0.22 s); separate task if commissioned.
- **Enabling `stop_threshold` in bootstrap replicates:** a semantics decision for Larry after §7; not touched here.
- **Changing bootstrap plan defaults or documentation** after §6: findings go in the record; Larry decides.
- **Deeper Cox surgery** (`coxph.fit`/`agreg.fit`), the **IPTW continuous defect**, batch-size defaults: recorded, untouched.
- **No renders, no push, no payload, application, or driver changes.**

## 10. Report

`dev/glm-continuous-sims/REPORT_assembly_skip_2026-09-02.md`:

1. Provenance; §2's HEAD quotes (assembly lines, consumers, NA guarantee), the transplant diffs, and **the §2.4 re-bucketed 0.3.4 prize beside the 09-01-derived figure**.
2. The restructure as implemented: lines moved, lines added, file:line.
3. The equality matrix against the committed baseline — every fixture × every retained component — with the baseline's provenance stated; self-consistency results.
4. §5's constructor counts, before/after bucket table, walls, bootstrap replicate and B = 1000 line, realized-vs-both-predictions.
5. §6's settlement: quotes, walls, payload verdict, and the implied production B = 1000 walls if reproducible.
6. §7's settlement: both selections, the ordering quote, the verdict.
7. §8's results against the recorded baselines, standing guards included.
8. Ten-line verdict.

Commits, in order: (1) this task doc alone; (2) the two transplanted scripts; (3) the implementation + version/NEWS; (4) all artefacts (`assembly_postchange_2026-09-02.*`, profile outputs, §6/§7 logs and `.rds`) + the report. Explicit paths; tree clean at close. **No push. No render. Nothing outside `R/subgroup_search.R`, `DESCRIPTION`, `NEWS.md`, and `dev/` artefacts.**
