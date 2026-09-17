# TASK — GRF on the MD design: campaign `mdgrf`, resumed after the membership fix

**File:** `dev/tasks/TASK_md_grf_resume_2026-09-16.md` · **Issued:** 2026-09-16 by chat, on Larry's yes to P1
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`
**Runs after:** `TASK_grf_dina_fixes_2026-09-16.md` has closed out
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (S1)
**Directory:** `quarto/simulations/actg175/continuous/`, written `<dir>`
**Governing task:** `dev/tasks/TASK_md_grf_2026-09-16.md`, written "the GRF task"
**Records:**
- `<dir>/REPORT_md_grf_stage1_2026-09-16.md` — the stopped Stage 1
- `<dir>/REPORT_grf_dina_fixes_2026-09-16.md` — the fix record

**What this is.** The GRF task stopped at its §1.6(c) on the membership defect. P1 has now landed and the package has been reinstalled. This task executes the GRF task from §1.6 through §3.5 as written, with the substitutions S1–S7 below.

Everything else in the GRF task stands:
- its dispositions, including `dmin.grf = 30` and `ci_method = "field"`;
- the governing constraint and the conditional-family labelling;
- its conventions;
- its gates;
- the Gate 1 advance go: projection under 8 hours → Stages 2 and 3 run.

**Compute authorized by this kickoff:**
- the smoke and calibration, with the GRF task's ceiling;
- Stages 2 and 3 under the advance go.

**No `R/` change. No install.**

## Substitutions

**S1 — preconditions (replace the GRF task's §1.1–§1.2 for this run).** *GATE* — all of the following hold, or stop:
- The host is `pop-os` and the branch is `feature/glm-extension`.
- HEAD contains the fix task's closeout commit.
- The fix record states that P1 landed.
- `git diff --quiet <the fix record's last R/ commit>..HEAD -- R/ DESCRIPTION NAMESPACE` succeeds.
- The installed `Built` equals the fix record's final `Built`, and two doFuture workers report the same.
- There are no tracked modifications, and no R, Rscript or quarto process is running.

Then copy this document from `~/Downloads` to `dev/tasks/TASK_md_grf_resume_2026-09-16.md` (exact name, else the single match for `*TASK_md_grf_resume_2026-09-16*.md`) and commit it alone.

**S2 — not repeated.**
- The GRF task's §1.3–§1.5 and the script creation in its §1.6 are already committed: the template at `894da993`, the scripts at `f0b9c844`. Quote both commits.
- *GATE:* both of these succeed, or stop:
  - `git diff --quiet 894da993..HEAD -- <dir>/sim_fs_maxeffCons_mr_field_md_template.qmd`
  - `git diff --quiet f0b9c844..HEAD -- <dir>/scripts_mdgrf/`

**S3 — clean the stopped smoke.** Delete the untracked smoke outputs the first Stage 1 record lists (tags `mdgrfsmokefs` and `mdgrfsmoke`: their result directories and logs), and nothing else.

**S4 — §1.6(a).** Re-run as written. The package was reinstalled, so FS must still be bit-identical to `mdsgnb20` within the GRF task's tolerance.

**S5 — §1.6(b) and (c).**
- Run (b) as written, with one change: report sim_id 1's selection beside Stage 0 S0.3's and beside the fix record's F4 after-P1 selection. These are facts, not gates.
- (c) is replaced by this *GATE:* every GRF smoke replicate records zero factor-comparison warnings in the captured-warnings column and zero NA-membership candidates. Any warning or NA stops the task.

**S6 — records.** The Gate 1 record is `<dir>/REPORT_md_grf_stage1_resume_2026-09-16.md` and cites the first Stage 1 record. Everything else in the GRF task's §1.8 stands, including applying the advance go.

**S7 — §3.4.** The DINA open-work line states the fix record's outcome for P2 as that record words it. Do not reword it.

## Closing

End with the GRF task's closing message for wherever this run stops, including the commit range to push. Then stop.

## Out of scope

- `R/`; installs.
- DINA.
- Any change to the GRF task's dispositions, cells, replicates, gates or knobs.
- Pushing.
