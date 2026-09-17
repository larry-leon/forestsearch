# TASK — ACTG175 binary Stage 0, resumed with the manuscript path restated

**File:** `dev/tasks/TASK_actg175_binary_stage0_resume_2026-09-17.md` · **Issued:** 2026-09-17 by chat, after v2 stopped at its §1 on the manuscript path
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `399e4b42`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (R2)
**Governing task:** `dev/tasks/TASK_actg175_binary_stage0_2026-09-17_v2.md`, written "v2" (committed at `d526f7ba`)
**Record:** `quarto/simulations/actg175/binary/REPORT_actg175_binary_stage0_2026-09-17.md` (committed at `399e4b42` with the §1 stop), completed in place

**What this is.** v2 stopped at its §1 because it stated the manuscript copy relative to forestsearch, where no such copy exists. The path Larry gave exists in the fs-glms-interpretable repository. This task executes v2 §2–§8 with the substitutions R1–R7 below. Everything else in v2 stands: its category (read-only, apart from this document and the record), its conventions, and §7's condition and 1.5 h ceiling.

## Substitutions

**R1 — the manuscript copy.**
- `<ms>` = `~/Documents/GitHub/fs-glms-interpretable/dev/reference/post_selection/manuscript_jrssb_initialsubmit/`, the path Larry gave.
- **Secondary**, consulted only for a file `<ms>` lacks: `~/Documents/GitHub/fs-post-selection/jrssb_submission/`.
- For every quoted file, state which of the two it came from.
- Neither repo is written. Run no `git` command that changes either one: no checkout, fetch, pull, or stash.

**R2 — preconditions** (replacing v2 §1).

*GATE:*
- the host is `pop-os` and the branch is `feature/glm-extension`;
- HEAD contains `399e4b42`;
- forestsearch has no tracked modifications;
- no R, Rscript or quarto process is running;
- `<ms>` exists.

Record `git -C <repo> rev-parse --short HEAD` and `git -C <repo> status --porcelain | wc -l` for fs-glms-interpretable and for fs-post-selection. Copy this document from `~/Downloads` to `dev/tasks/TASK_actg175_binary_stage0_resume_2026-09-17.md` (exact name, else the single match for `*TASK_actg175_binary_stage0_resume_2026-09-17*.md`) and commit it alone.

**R3 — record location** (v2 §2.5). The location is fixed at `quarto/simulations/actg175/binary/`, where v2's §1 found the sweep drivers and `build_actg175_glm_dgm.R`.

**R4 — the record.** Complete the committed record in place:
- Replace its stop section with R2's provenance, keeping how the DINA task ended.
- Add §2–§8's content.
- Commit it by explicit path.

The directory has no catalog generator, so there is no catalog step.

**R5 — post-conditions** (replacing v2 §8.3), printed in the closing message:
- `git diff --name-only <R2 HEAD>..HEAD` lists only this document and the record.
- For fs-glms-interpretable and fs-post-selection, HEAD and the `status --porcelain` count equal R2's.
- forestsearch has no tracked modifications.
- The installed `Built` is unchanged.
- The temporary directories are gone.

**R6 — bundle for the chat** (the last action, after the record's commit).
- Write one archive, `~/Downloads/bundle_md_and_binary_stage0_2026-09-17.zip`. If `zip` is unavailable, write `.tar.gz` instead.
- Copy each file from HEAD with `git show HEAD:<path>` into a temporary directory, then archive them:
  - the completed record from R4;
  - from `quarto/simulations/actg175/continuous/`: `md_field_metrics.csv`, `COLUMNS_md_field.md`, `md_grf_metrics.csv`, `COLUMNS_md_grf.md`, `md_dina_metrics.csv`, `COLUMNS_md_dina.md`, `REPORT_md_field_rerun_2026-09-15.md`, `REPORT_md_grf_2026-09-16.md`, `REPORT_md_dina_2026-09-17.md`, and `REPORT_grf_dina_fixes_2026-09-16.md`;
  - `quarto/applications/actg175/REPORT_actg175_continuous_intervals_2026-09-16.md`.
- A file that is not at HEAD is listed as missing. That is not a stop.
- Remove the temporary directory afterwards.

**R7 — closing message:** v2 §8.4's contents, plus the bundle's path and its listing (file names and sizes). Then stop.

## Out of scope

- Everything v2 puts out of scope.
- Any write to fs-glms-interpretable or fs-post-selection.
- Pushing.
