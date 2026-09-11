# TASK — DINA MR campaign at the FS-analogous criterion, both prevalences (campaign `dinamr`)

Date: 2026-09-10. Author: chat (spec). Executor: Claude Code (Linux), unattended. Approver: Larry (kickoff paste = compute go). Reviewer: the Linux MR-field chat.
Predecessors: `REPORT_dinamr_stage0_2026-09-10.md` (the STOP), `REPORT_o1_forwarding_2026-09-10.md` (which cleared it), `REPORT_cert20_2026-09-08.md`, `REPORT_tier2_2026-09-08.md`, `REPORT_p12ext_2026-09-09.md` (the FS grid this mirrors), `dev/notes/NOTE_survival_products_2026-09-09.md`.

## Framing — state in the record's opening paragraph

- **This does not certify DINA.** DINA's candidates are read off a cross-fit surface a bootstrap would regenerate, so the manuscript's fixed-family condition (§2.1) does not hold. After the alignment repair the estimate and interval are first-order exact **conditional on the proposed family**; **every coverage number here is coverage of that conditional estimand** and each table must say so.
- **The criterion is FS-analogous by construction:** `effMaxSG` at ε = 0.20 reuses the inclusion-band logic shared with `forestsearch()` (`effect_neighborhood` honored, semantics matched), with DINA's effect floor log(0.90) = −0.105361 and **no consistency term** — the one structural difference from FS, because neither DINA nor GRF has a consistency floor.
- **Grid mirrors FS** (Larry, 2026-09-10): both prevalences, ε = 0.20, the same HR and n structure as `cert20` / `tier2` / `p12ext`.
- **GRF is not in scope.** Its `effMaxSG` uses the shared band **frontier-only, as a filter with no empty-band fallback**, where DINA uses it as a sort key, and `dmin.grf` is in DR-score units not comparable to the log-HR floor (`REPORT_o1_forwarding` V1/V5). Those are open decisions.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/` — **do not archive `HANDOFF_guohe_comparison_2026-09-09.md`**. Copy this file to `dev/tasks/` and commit. **Commit only; do not push.**
- **No `R/` change.** One document-level, add-only template change (Part T). If anything else appears to need `R/`, STOP and report.
- Gates stop per cell; a failed cell stops that cell and the next proceeds; the driver re-projects from realized walls of the same n and prevalence and defers, listing, any cell that would cross the ceiling. `.refuse_if_tracked()` live. `devtools::install(dependencies = FALSE)` if installed does not match HEAD.
- Standing conventions: bounds by location; Wilson intervals; marginal and error SDs side by side; NPV beside sens/spec/PPV; winner-only and winner-floor excluded. Leave the seven pre-existing untracked files alone. **No recommendation change; report and wait.**

## Part T — Engine knob (document-level, add-only, default-inert)

`subgroup_method` is hard-coded at `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:300`. Rather than duplicating the template, add `FS_S7_METHOD` (`.env_chr("FS_S7_METHOD", "consistency")`, `stopifnot` in `c("consistency","dina","grf")`), echo it in the knob audit line and the settings readout, and record it in batch and combined meta. **Nothing else changes.**

**Gate T** (5 replicates, the standing identity cell — effMaxSG ε 0.20, HR 1.50, n 500, `FS_S7_Z1Q=0.60`, seeds `8316951 + sim_id`, sim_id 1–5, tag `methknob`): with `FS_S7_METHOD` **unset**, every non-timing column and `truth` `identical()` to the committed `e1stud` rows 1–5. STOP and revert on failure.

## Part C — The campaign

All cells: `FS_S7_METHOD=dina FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=dinamr`, `dina_select_statistic = "effect"` (already the template default — confirm), seeds `8316951 + sim_id`, sim_id 1–2000, two seed-disjoint batches of 1,000 then combine, 100 workers. **`FS_S7_ER_JCUTS` is inert on DINA** (it feeds the consistency-only `method_args`) — do not set it, and say so in the record.

**Blocks, run in this order** (so a partial run still yields a complete picture at one prevalence):

| Block | Prevalence knob | Cells |
|---|---|---|
| **A** | `FS_S7_Z1Q` unset (12.4%) | HR 1.50 and 1.75, at n = 500, 1000, 1500 — six cells |
| **B** | `FS_S7_Z1Q=0.60` (31%) | HR 1.50 and 1.75, at n = 500, 1000, 1500 — six cells |
| **C** | both prevalences | HR 1.00 at n = 500, 1000, 1500 — six cells |

**Defer order if the ceiling threatens:** Block C first (all six, 31% before 12.4%), then Block B's n = 1500 cells, then Block B's n = 1000.

**Stage 1:** 5-replicate smokes at (12.4%, HR 1.50, n 500) and (31%, HR 1.50, n 1500): every construction finite; interval invariants on harm, complement, `_s`, joint and IJ; γ in range; realized prevalence stated for each block; p̂, ρᶜ and the nine recovery columns present (they are reachable now — **if any is absent, STOP**, since that would mean the O-1 fix did not take on this path); **the proposed-family size per replicate recorded**. Time per replicate at 100 workers. **Gate 1: proceed if the total projection is ≤ 10 h wall; hard timeout 12 h.** Project from the **family-size distribution**, not from a mean and never from the FS walls — the pilot measured 0.13–27.8 s per replicate with a ~100× right-skewed family (median 91, max 1745).

**Gate 2 per cell:** completeness (2,000 rows, sim_id 1–2000, no duplicates, no CONFIG-ERROR, meta knobs as set including `FS_S7_METHOD`, `n_workers` and `forestsearch_version` recorded); **detection rate recorded prominently** (DINA's differs from FS's and every summary conditions on it); proposed-family size distribution; realized prevalence; every harm / complement / `_s` / joint / IJ / β(Ĥ) / β(Ĥᶜ) / p̂ / ρᶜ / recovery quantity finite on detected replicates; interval invariants; γ ∈ [0.025, 0.05]; bound↔quantile identities ≤ 1e-12; **structurally-NA columns (`n_cons_qual`, `band_n`, `p_star`) reported as such, not as failures.**

## Stage 3 — Report

`summary_dinamr.qmd`, transplanted from the committed `summary_cert20.qmd` (globs, labels, comparator names only); `REPORT_dinamr_2026-09-10.md` beside the results with the Gate 2 record inside, every number verbatim:

1. **Standard tables** (per-cell constructions; across cells), both blocks, rows naive / field / field-s / IJ two-term: bias (log; marginal-SD and error-SD units), SDs, SE, r, SE/error-SD, one-sided on the exposed side [Wilson], two-sided [Wilson], Gaussian reference. **Every coverage column labelled as the conditional-on-proposed-family estimand.**
2. **Classification vs n and prevalence:** detection, sens, spec, PPV, NPV, mean |Ĥ|, |Ĥ|/|H| median and q90, β(Ĥ) and β(Ĥᶜ) against the planted values, naive optimism in both SD units.
3. **The proposed family:** size distribution per cell, and whether it stabilizes with n — the mechanism Supplementary S8.3 offers for conditional coverage recovering with n.
4. **Beside the FS grid** (`cert20` / `tier2` / `p12ext` read as they stand), **with the confound stated**: FS and DINA differ in identifier *and* in family construction, and each conditions on a different detection set. Descriptive, not a contest.
5. **The two-sided question:** does the harm-block two-sided decay with n at 12.4% appear on DINA as it does on FS (IJ 0.977 → 0.932 → 0.901 at HR 1.50; 0.961 → 0.917 → 0.913 at HR 1.75)? Report the same trajectory and the miss split by side.
6. **p̂ and the recovery diagnostics on a model-generated family:** distributions per cell; whether the low-p̂/high-p̂ bias structure found for FS appears; `sens_H` and containment quantiles. Describe; do not interpret beyond the record.
7. **No acceptance criteria.** Exploratory. FS certified ranges as reference lines only; no recommendation.

## Done means

Part T committed with Gate T PASS; Stage 1 and Gate 1 recorded; cells completed / deferred / dropped with walls; Gate 2 records; `summary_dinamr` rendered; the report committed; one-paragraph closing summary with the detection and family-size headline, the two-sided answer, cells and walls, and the commit range. **Out of scope:** GRF, any `R/` change, any recommendation change, the certification language.

---

## Appendix — kickoff amendments (added at receipt, 2026-09-11; not part of the document as authored)

The kickoff paste that authorized this campaign carried three amendments that **supersede the
body above where they differ**. They are recorded here so the task file and the executed
protocol do not diverge.

- **AMENDMENT 1 — Gate 1 ceiling.** The ceiling is **9 h wall, not 10 h**; the hard timeout
  stays 12 h. Source (quoted from the kickoff): FS projections have run +16% (`cert20`) to
  −29% (`p12ext`) against realized, and DINA's ~100× family-size skew is less predictable than
  FS's, so the wider margin is deliberate. Supersedes "Gate 1: proceed if the total projection
  is <= 10 h wall" under **Stage 1**.

- **AMENDMENT 2 — Block A checkpoint.** After Block A completes, run **one** checkpoint:
  re-project Blocks B and C from Block A's realized walls against the remaining budget to the
  9 h ceiling, record the revised projection, and defer per the stated defer order if it no
  longer fits. **Once, not per cell.** Adds to the Gate/defer machinery under **Protocol** and
  **Part C**.

- **AMENDMENT 3 — same-draws assertions at Gate 2.** Where a committed FS comparator shares the
  DGM draws — Block A against `tier2` / `p12ext` at 12.4%, Block B against `cert20` at 31% —
  assert `n_true` identical on all 2,000 rows and `truth` `identical()`, as `cert20` did.
  Report any cell where the draws do not match, and **do not treat a mismatch as a cell
  failure**: it is a finding about the DGM path, not about the campaign. Adds to **Gate 2 per
  cell**.

Also carried in the kickoff, restating the body: `FS_S7_ER_JCUTS` must **not** be set (it is
inert on DINA); GRF is out of scope; report and wait, with no acceptance criteria and no
recommendation.
