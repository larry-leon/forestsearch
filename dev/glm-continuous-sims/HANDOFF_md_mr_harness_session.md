# HANDOFF: GLM/continuous MR harness — session 2 start (v2, supersedes v1)

Next session: fetch from `github.com/larry-leon/forestsearch`, branch
`feature/mr-in-replicates`, path
`dev/glm-continuous-sims/HANDOFF_md_mr_harness_session.md`. Read fully.
Companion task document committed beside it: `CC_TASK_md_mr_harness.md`.
This v2 corrects v1, which repeated two chat-side fabrications (below).

## Workstream priority (standing)

Continuous/MD first — the tractable setting where closed forms confirm the
simulations. Sound theoretically AND computationally before binary; binary
before survival. Prior simulations will all be re-done; no main-thread time
on past results. Every task justified by its contribution to sound
GLM/continuous machinery.

## Decisions, closed — do not reopen

- `focus = "harm"`; DGM `calibrate_glm_interaction("continuous", "MD")`
  with `cal_target_md = -40`; no recalibration.
- Orientation `adverse_outcome = FALSE` per
  `dev/glm-continuous-sims/NOTE_threshold_sign_md.md` (traced for the SWEEP
  path — the MR path's orientation must be traced separately, see blockers).
- v1 identifier FS/consistency. Aim: unbiasedness + interval coverage of
  the MR-corrected estimate against EXACT β(Ĥ); θ† reference columns.
- NEW (this handoff): the harness adopts the MR-COMPATIBLE configuration
  in principle (the binary MR docs' pattern — no lasso/GRF screening, fixed
  candidate family). Rationale: MR is the aim; the −61.87 sweep was
  exploratory and config-specific, and the harness's own naive column
  re-measures selection under the config MR actually certifies. Adoption is
  CONDITIONAL on the blockers below being resolved.

## Ground truth (filesystem/git, reported by CC — verbatim)

- Remote tip `79bce8d`; local == remote; 0 unpushed.
- Committed and real: consolidation through T1–T10; `f8bdb80` (20 engines
  call `fs_attach_betaHhat()`; bundles carry status/nH/nHc); `ab53239`
  (attach-ITT); `bd11f87` (full-scale control); `79bce8d` (eval/theta SPEC).
- UNCOMMITTED: `R/betaHhat_truth.R` (the two functions
  `fs_build_eval_frame()` / `fs_betaHhat_theta_dagger_check()` implemented),
  `tests/testthat/test-betaHhat-contract.R` (T11–T14 drafted), `NAMESPACE`,
  untracked `man/` pages. **Test and check status: a clean result exists for
  exactly this file state, but re-verify.** Filesystem evidence:
  `/tmp/claude-1000/-home-larryleon-Documents-GitHub-forestsearch/c57cd555-9ef2-41ea-859a-e431519d800f/scratchpad/check5.log`,
  mtime `2026-08-07 19:41:05 -0700`, final line verbatim:
  `0 errors ✔ | 0 warnings ✔ | 0 notes ✔` (`Status: OK`,
  `Duration: 14m 46.8s`, with `checking tests ... Running 'testthat.R'
  [654s/484s] OK`). The log postdates all three uncommitted files —
  `R/betaHhat_truth.R` (19:24:19), `test-betaHhat-contract.R` (19:25:29),
  `NAMESPACE` (19:26:11) — so it certifies **this** tree state, not merely
  some earlier one. It is scratch, outside the repo, and will not survive;
  blocker 1 still runs the tests and check and commits the result. Treat a
  *disagreement* with this line as a signal worth investigating, not as an
  expected first run. The full engine re-point (dropping shim `source()` for
  frame/theta) NEVER HAPPENED; engines still source shims for those two
  pieces, which works.
- Harness `quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd`
  exists, untracked, 28,563 bytes. One pilot bundle
  (`mr_sweep_md_harm/md_harm_s50_pilot/fs_md_harm_n1000_res.rds`): pilot ran,
  **detected 0 of 50**, render halted in check-A. No grid, no HTML, no
  processes running.

## The incidents this handoff exists to prevent repeating

TWO chat-side stop-point approvals were fabricated on empty attachments:
the eval/theta "STOP 2 complete" review and the harness "STOP A approved /
proceed to full grid" (which CC never received — confirmed). The tell both
times: zero concrete numbers or hashes in the review. Standing rules:
1. No stop-point review without the report's verbatim numbers. Empty or
   number-free = ABSENT.
2. Transport: CC content as message-body text; long reports committed by CC
   and pushed, then fetched and read from the remote. Attachments retired.
3. Committing a report is transport, not approval.

## Blockers, in order (the next session's critical path)

1. **Eval/theta closeout.** Run T11–T14 and `R CMD check --as-cran`; report
   verbatim; commit implementation + tests + man (one commit). The engine
   re-point stays deferred/queued — not on the harness's critical path.
2. **Verify the MR incompatibility from source.** The exact guard or error
   (file:line) making `mr_inference = TRUE` incompatible with
   `use_lasso/use_grf = TRUE`. Establish whether it is structural (MR
   requires a fixed candidate family across resamples) or incidental. If
   structural, "MR + adaptive screening" is a different problem — record as
   a follow-on question, do not solve.
3. **Diagnose 0/50 on ONE replicate** under the MR-compatible config:
   candidate family size; the fixed cut grids vs the true boundaries
   (`age > 34`, `preanti <= 744.5` on raw scales — do the grids bracket
   them?); estimated effects on the truth-adjacent candidates; the threshold
   value, sign, and orientation actually applied on THIS executing path
   (the sign NOTE traced the other path); where detection failed. Then the
   minimal config fix, stated before applied.
4. **Re-pilot 50 replicates.** STOP A = the report with verbatim numbers
   (all readouts + the three closed-form assertions per the task document),
   committed to `dev/glm-continuous-sims/` and pushed; review against the
   fetched file. Then STOP B (grid n ∈ {500, 1000, 2000, 4000}, s = 1000)
   the same way. Then the binary port; survival last; campaign stays
   background-only.

## Key file map

- Task: `dev/glm-continuous-sims/CC_TASK_md_mr_harness.md`
- Estimand authority: `NOTE_target_is_collapsibility.md` (§2–3; the identity
  code lives in `verification/acceptance_betaHhat_md.qmd` criteria 3–4)
- Orientation: `NOTE_threshold_sign_md.md` (sweep path only — see blocker 3)
- Consolidation record: `dev/betaHhat-consolidation/`
- Larry-only, unblocked: Table S7 C_θ column mapping vs the supplement
  source, queued for the manuscript revision.
