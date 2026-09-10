# TASK AMENDMENT v2 — Guo & He supplement (Tier 1 + Tier 2(a)): post-merge Stage 0 and three resolutions

Date: 2026-09-09. Governs **jointly with** `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` (committed at `179ae409`). Where the two differ, **this document wins**. Everything in v1 not amended here stands unchanged — including the D1–D4 decisions, the section shape (B1–B6), the fixed sentences, the transplant-first rule, and the bullet-form reporting requirement.

Basis: `quarto/GuoHe/REPORT_guohe_supp_stage0_2026-09-09.md` (committed `33d0add3`), which stopped correctly at the absent certification record and surfaced three findings. Larry's dispositions on all three are below.

---

## A0. Constraint, restated and tightened

- **CC runs no network git at all**: no `fetch`, no `pull`, no `push`, in this phase or Phase B. Larry performs the Mac pull in GitHub Desktop before this task resumes, and sequences all pushes.
- Provenance first; every `git add` by explicit named path; `git add -A` / `git add .` forbidden; commit locally only.
- The seven out-of-scope untracked files (three `diag_*` bundles, three `diag_*.html`, the ACTG175 payload) are **expected to remain absent on the Mac** — untracked files do not travel through git. The rule stays in force and binds on the Linux box in Phase B. Record the untracked set as observed, whatever it is; never stage it.

## A1. Stage 0 v2 (Mac, post-merge) — new record, no overwrite

- Write `quarto/GuoHe/REPORT_guohe_supp_stage0_2026-09-09_v2.md`. Do **not** modify the committed v1 record; it correctly documents the pre-merge tree and the STOP.
- Open it with `git log -1 --oneline` and `git status -sb` verbatim, plus the observed untracked set.
- **All Q1–Q6 line numbers from the v1 record are void after the merge.** Re-verify and re-quote every one against the current tree, reporting old → new where they moved. A quote whose line number cannot be re-verified is a STOP.
- Re-assert the input inventory of v1 §2.2 (24 bundles, 6 truth caches, 11 scripts/`R/` sources) with byte sizes.
- Record whether `REVIEW_certification_2026-09-09.md` is now present. **Its absence is no longer a STOP** (see A2).

## A2. Q5 resolved from the engine; certification record demoted to a citation

- The joint-Bonferroni construction is transplanted from the engine's own `.fs_mr_field_joint` — `bonf_lower_H`, `bonf_upper_Hc`, `bonf_joint_prob` (v1 reading: `R/fs_mr_inference.R:1195–1220`; re-verify line numbers per A1). This supersedes v1's requirement to transplant from a certification harness. **No `R/` change is needed or permitted.**
- The certification record is now required only for the **B6 citations** in T3 (field-s upper 0.912–0.960; joint 0.939–0.963). If present: quote with file path and line numbers. If absent: render B6 with the literal marker `[certification citation pending sync]` at each of the two figures, note it in the T3 record, and continue. T1 and T2 are unblocked either way.

## A3. T1 driver calls `fs_mr_inference()` directly

- `mv_mr()` in `mr_vs_guohe_sim.R` has no `field_scale_complement` formal. **Do not modify `mr_vs_guohe_sim.R`** — it is a committed input to two completed campaigns.
- Instead: in `quarto/GuoHe/mr_field_complement_vs_guohe_run.R`, **copy `mv_mr()`'s argument-assembly block verbatim** into a local call to `fs_mr_inference()`, then add exactly `field_complement = TRUE`, `include_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`. Do not re-author the argument list from memory or from documentation; copy the executing lines and cite them in the record.
- The v1 Stage 1 identity probe is now **load-bearing**: it is the proof that the hand-assembled call reproduces the stored path. Run 3 replicates on each of `t7_beta2_00` and `t7_beta2_05`; assert `identical()` of recomputed naive, IJ, and field columns against the stored `mr_field_vs_guohe_<id>.rds` rows. Any mismatch: STOP, naming the differing column and the first differing replicate.

## A4. Pilot scaffolding transplant (Gate 1a)

- `mr_field_vs_guohe_run.R` has no `--pilot` flag. Transplant it from `guohe_sec52_run.R` (v1 reading: lines 70, 81–83, 144, 214–225; re-verify per A1) — named lines only, no re-authorship. Gate 1a's envelope is unchanged: proceed only if the 6 × 2000 projection is ≤ 40 core-h and ≤ 90 min Mac wall; otherwise STOP and report.

## A5. T2 adaptive bootstrap: **B = 2000**, superseding v1's B = 200

- v1 pinned `B = 200` for `guohe_adaptive_r()` on the stated grounds that it was "the validated reproduction setting." **That premise is false** — `--adaptive-B` is inert (A6), so the validated Tables 3–6 Adaptive columns executed at 2000. v1 §5 is amended accordingly.
- T2 calls `guohe_adaptive_r()` with `B = 2000` (one `B` serves both the inner CV fits and the final refit; that coupling is the function's own behavior and is **not** to be altered — splitting it would be an implementation change to their algorithm and is out of scope for this task). All other T2 settings stand: `orient = +1`, `r_grid = c(1/3, 1/12, 1/21, 1/30)`, `v = 5`, derived seed recorded, r̂ and per-candidate objectives stored.
- Rationale to record: B = 2000 is simultaneously the setting under which the function was validated, resolution parity with the fixed-r columns, and their method's strongest configuration — the honest-comparison requirement.
- **Cost consequence, to be stated in Gate 1b rather than assumed:** the ledger's 5.90 s/replicate adaptive figure may itself have executed at 2000 (flag inert) or at a genuine 200; the full-grid projection therefore spans ≈1,400–1,700 core-h (13–16 h wall) at best to roughly ten times that at worst. The Gate 1b pilot measures it. Envelope unchanged (≤ 2,500 core-h and ≤ 24 h wall): if exceeded, **STOP and report** the measured projection together with costed reduced options (500 reps/cell; three cells; B = 200 as a labelled sensitivity) — do not launch, and do not choose among them.

## A6. NOTE for the inert `--adaptive-B` flag — record only, no repair

- Write `dev/notes/NOTE_adaptive_B_inert_2026-09-09.md` (commit by named path). Content, bullet form, with verified line numbers:
  - `b_adapt` is parsed, printed, and stored as bundle metadata in `guohe_reproduction_run.R`, but never reaches `gh_one_rep`; line 109 passes `B = b_boot`, and one `B` serves both the fixed-r loop and `guohe_adaptive_r()`.
  - Consequence: committed Tables 3–6 Adaptive bundles record `adaptive_B = 200` while having executed at 2000. This is a **provenance-labelling defect, not a numerical error** in those columns; the published-comparison conclusions drawn from them are unaffected.
  - `guohe_reproduction_RUN.md:98–99` (184 vs 875 core-h) is therefore a projection ledger, not a record of what ran.
  - **No repair in this task. No re-run of committed work. Disposition is Larry's.**

## A7. Amended commit order (each by explicit named path)

1. This amendment → `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09_v2.md`, committed alone, first.
2. Stage 0 v2 record.
3. `NOTE_adaptive_B_inert_2026-09-09.md`.
4. T1 driver; then Stage 1 probe note appended to the Stage 0 v2 record.
5. Gate 1a record; then six T1 bundles + T1 REPORT.
6. T2 driver (authored, **not executed** on the Mac).
7. T3 qmd (B2 with the Adaptive column rendered "pending Phase B"; B6 per A2).
Phase B on Linux, after Larry's go: Gate 1b record, six T2 bundles, T2 REPORT.

## A8. STOP list (amended)

- Removed: absence of `REVIEW_certification_2026-09-09.md` (now A2's citation marker).
- Retained: any missing bundle, truth cache, or script; any Q-quote whose line number cannot be re-verified post-merge; any `identical()` failure in probe or production pairing; any needed `R/` change; Gate 1a or Gate 1b envelope exceedance; any interaction with out-of-scope untracked files; anything touching a closed line.
- Added: any invocation of `git fetch`, `git pull`, or `git push` — forbidden outright, in both phases.
