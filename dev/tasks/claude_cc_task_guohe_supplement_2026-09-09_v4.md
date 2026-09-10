# TASK AMENDMENT v4 — Guo & He supplement: probe resolved on the Mac; proceed to Gate 1a

Date: 2026-09-09. Governs **jointly with** `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` (v1, `179ae409`) and `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09_v2.md` (v2, `669c9ef5`). Where they differ, **this document wins**.

**Supersedes v3 in full.** A v3 amendment was drafted but **never delivered to CC and never committed**; it contained a contradiction between its M1 and M2 clauses. If a file named `claude_cc_task_guohe_supplement_2026-09-09_v3.md` reaches this repository by any route, **do not commit or execute it** — report its presence and stop.

Decision (Larry, 2026-09-09): **Phase A remains on the Mac Studio.** Cross-platform floating-point deviation at the 1e-16 level is expected and is not an obstacle. The A0 constraint stands unchanged: CC runs no `git fetch`, `pull`, or `push`; every add by explicit named path; out-of-scope untracked files never staged.

---

## N1. Probe and pairing standard (replaces A3's plain `identical()` requirement)

Precedent: `dev/tasks/TASK_tier2_mac_2026-09-08.md`. Applied to the Stage 1 probe and to every production replicate:

- **Integer, selection, flag and seed columns:** `identical()` **required**, every replicate; any failure is a STOP. Covers the selected cutpoint, engine selection index/label, harm flag, admission/reselection indicators, coverage flags, seeds, replicate index, and truth lookups keyed by the selection.
- **Floating-point columns:** `all.equal()` at tolerance `1e-8`. Record per cell the worst absolute and worst relative deviation against the stored bundle, and the fraction bit-identical.
- **Complement exoneration — established, not re-run:** the Mac isolation test (complement enabled vs disabled, same machine, same seeds, 29/29 shared columns `identical()`, both probe cells, all three replicates) is a stronger proof of non-perturbation than the cross-machine comparison A3 specified. Cite it with its commit in the T1 REPORT.
- The cross-machine comparison against Linux-built bundles is a **provenance measurement, not a gate**: stored bundles are x86_64 / R 4.6.1 / reference BLAS; this Mac is arm64 / R 4.5.2 / Accelerate.

## N2. Stage 1 probe — scored PASS on evidence already recorded; do not re-run

The probe outputs already satisfy N1 and re-running them would be a verification whose outcome is determined:

- 120/120 discrete comparisons `identical()` across the six probe replicates (selection, engine selection, admission/reselection, coverage flags, seeds, truth lookups).
- Worst absolute float deviation 4.44e-16 (`mr_upper_2s`, `t7_beta2_00`, m = 1); `all.equal` at 1e-8 TRUE throughout.

Action: append a short note to `quarto/GuoHe/REPORT_guohe_supp_stage0_2026-09-09_v2.md` formally scoring the probe PASS under N1, citing these measured values and the isolation-test commit. Then proceed to Gate 1a.

## N3. Selection stability — mismatch is a STOP; tally is reported

This supersedes any tolerance for selection mismatches. Rationale, to be recorded in the T1 REPORT:

- Top-1 minus top-2 oriented-score gap at selection (106 selections: six probe replicates plus 100 fresh `t7_beta2_00` draws at m = 2001–2100, disjoint from the stored 2000): minimum 1.438e-04, 1st percentile 3.961e-04, median 5.205e-03.
- Ratio of minimum gap to worst absolute float deviation: 3.2 × 10¹¹ (1st percentile 8.9 × 10¹¹). The floating-point deviations cannot reach a selection boundary.

Therefore: per production replicate, compare the selected cutpoint against the stored `mr_field_vs_guohe_<id>.rds` row; **tally per cell and report the tally in the T1 REPORT (expected zero)**. Any mismatch is a **STOP**, reported with the cell, replicate id, and that replicate's top-1/top-2 gap — on the recorded margin a mismatch indicates something structural, not arithmetic.

## N4. Gate 1a — unchanged

Pilot at reps = 20 on `t7_beta2_00`; proceed only if the 6 × 2000 projection is ≤ 40 core-h and ≤ 90 min Mac wall; otherwise STOP and report. Pilot scaffolding transplanted per A4.

## N5. Provenance for the new columns

The probe guards only columns with stored counterparts; the complement and joint columns have none. The T1 REPORT must quote, with file path and line numbers as read from this tree, the engine's own returned field names for the 0.95 one-sided complement upper bound and for `bonf_lower_H`, `bonf_upper_Hc`, `bonf_joint_prob` in `.fs_mr_field_joint` — demonstrating the emitted values are read from the engine, not reconstructed in the driver. No `R/` change is permitted (A2).

## N6. Two corrections to the record — append as a dated addendum; reopen nothing

- **Block ordering.** v1 §4 states the engine's complement block precedes the field block. It is inverted: the complement block **follows** the field block and re-reads the harm field's draws (`R/fs_mr_inference.R:820-826`, `:889-910`; re-verify when writing). Source of the error: `REPORT_mr_field_stage1_2026-09-05.md:15` describes the pre-`field_complement` layout and predates the 2026-09-06 addition. The executing code governs. **The 09-05 work is not reopened**; its numbers are unaffected.
- **Sync state.** No merge landed before the v2 reading — the last network operation is `a9c0d5e5`, already HEAD at Stage 0 v1, and `origin/feature/glm-extension` still points there. v2's "post-merge" framing is void. `REVIEW_certification_2026-09-09.md` remains absent, so the **A2 pending-citation marker path is live** for B6.

## N7. T3 dispositions (qmd already committed at `c1f7bc61`; edits are named-line only)

- **B6 range.** The computed range is kept over any transcribed figure. Change the wording so the sentence leads with the minimum: coverage **bottoms out at 0.933**, with the three lowest cells spanning the computed range. As drafted it reads as if 0.938 were a low point when it is mid-pack.
- **B5 classification metrics.** The closed forms are correct and are kept: with Ĥ = S(30) nested inside every S(ĉ), sensitivity ≡ 1 and NPV ≡ 1 **identically**, specificity = (80 − ĉ)/50, PPV = 30/ĉ. Report the closed forms, ĉ's distribution, and per-cell means of specificity and PPV. The standing NPV rule is satisfied by stating that NPV is degenerate by construction here — not by tabulating a column of ones. State plainly that sensitivity carries no information in this design and that the identification story is ĉ's overshoot past 30.
- **Provenance line.** One sentence stating that the T1 complement/joint bundles were produced on arm64/Accelerate while the B2 head-to-head bundles were produced on x86_64/reference BLAS, with the worst measured deviation and the N3 tally. If the tally is zero, no further caveat appears anywhere in the section.

## N8. Phase B — unchanged

T2 executes on the Linux box after Larry sequences the sync, at `B = 2000` per A5, with the Gate 1b envelope (≤ 2,500 core-h, ≤ 24 h wall) and its STOP-and-report-with-costed-options behavior. Its pairing proof runs against Linux-built bundles on Linux and therefore remains plain `identical()`; N1's float tolerance does not apply there.

## N9. STOP list (amended)

- **Replaced:** "any `identical()` failure" → the N1 standard (integer/selection/seed `identical()`; floats at `all.equal` 1e-8) and N3 (any selection mismatch).
- **Retained:** any missing bundle, truth cache or script; any Q-quote whose line number cannot be verified; any needed `R/` change; Gate 1a or 1b envelope exceedance; any interaction with out-of-scope untracked files; anything touching a closed line; any invocation of `git fetch`, `git pull`, or `git push`.
- **Added:** arrival of a v3 amendment file by any route (report and stop).

## N10. Order from here (each commit by explicit named path)

1. This amendment → `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09_v4.md`, committed alone, first.
2. N6 addendum + N2 probe-PASS note appended to `quarto/GuoHe/REPORT_guohe_supp_stage0_2026-09-09_v2.md`.
3. Gate 1a: pilot, record, then production — six T1 bundles.
4. `quarto/GuoHe/REPORT_mr_complement_vs_guohe_2026-09-09.md`, leading with the N3 tally and carrying the N5 quotes, the N1 deviation summary, and per-cell complement upper coverage (Wilson), mean bound location, joint coverage (Wilson), p̂(Ĥ), and marginal SD beside error SD.
5. N7 named-line edits to `guohe_supp_section.qmd`; re-render.
6. No action on the T2 driver (`ede97a2e`) or the T3 qmd creation (`c1f7bc61`) — already committed, not executed.
