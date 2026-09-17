# REPORT — ACTG175 binary (OR) simulation study, Stage 0 (v2): **STOP at §1**

Date: 2026-09-17 (UTC). Machine: `pop-os`. Branch `feature/glm-extension`. Task: `dev/tasks/TASK_actg175_binary_stage0_2026-09-17_v2.md` (committed as received, `d526f7ba`). Read-only: no `R/`, template, script, document or payload edit; no install, render or computation. **Outcome: the §1 gate failed on "`<ms>` exists"; §2–§7 were not run.**

## 1. Provenance — GATE FAIL

```
pop-os
feature/glm-extension
1d122b21
1d122b21 actg175/continuous closeout: generate current_status.md at c767027d
c767027d actg175/continuous status_curated.md: campaign mddina (payloads, records, the conditional-family reading convention ...)
58be2fda mddina record (TASK_md_dina_campaign_2026-09-17 §3.3): DINA on the MD design, four cells x 2,000 replicates ...
34067b7d mddina extract (TASK_md_dina_campaign_2026-09-17 §3.2): md_dina_metrics.csv, 1,174 rows ...
1ffef474 mddina summary (TASK_md_dina_campaign_2026-09-17 §3.1): summary_continuous_field_mddina.qmd ...
[tracked modifications: none]
DINA task started
[ps | grep: two lines are the shell wrappers running this check; one is PID 2086198, an 11-day-old bash memory-sampling loop from another session whose command text contains "quarto"]
[by process name (pgrep -x R / Rscript / quarto): none]
ls: cannot access 'dev/reference/post_selection/manuscript_jrssb_initialsubmit': No such file or directory
bc82bc3 Create ujag122_supplemental_files.zip          (fs-post-selection, 2026-09-15)
R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix
```

| gate item | result |
|---|---|
| host `pop-os`, branch `feature/glm-extension` | holds |
| the DINA task document is at HEAD | holds |
| no tracked modifications | holds |
| no R, Rscript or quarto process (the DINA campaign finished or stopped) | holds (by process name; the grep pattern matches shell wrappers and the unrelated sampling loop) |
| `<ms>` = `dev/reference/post_selection/manuscript_jrssb_initialsubmit/` exists | **fails**: forestsearch has no `dev/reference/` directory at all, and `git log --all` shows no commit that ever touched `dev/reference/post_selection` |

**How the DINA task ended.** `TASK_md_dina_campaign_2026-09-17` closed out completely: Gate 1 passed and the advance go applied (`REPORT_md_dina_stage1_2026-09-17.md`, `37395393`). Gate 2 passed 65 of 65 checks in each of the four cells, with no halt file (`REPORT_md_dina_gate2_2026-09-17.md`; campaign complete `30b282eb`, 27,727 s). Stage 3 is recorded in `REPORT_md_dina_2026-09-17.md` (`58be2fda`); its catalog closeout is `1d122b21`, and `check_current_status.sh --commit` passes.

**Where copies of the manuscript exist on this machine** (located, not read; the choice is Larry's):
- `~/Documents/GitHub/fs-glms-interpretable/dev/reference/post_selection/manuscript_jrssb_initialsubmit`: the task's relative path, but in the repository `fs-glms-interpretable` (branch `main`; last commit touching it `47880a0`, 2026-09-11, "docs(theory): unlabel duplicate equations, complete bibliography").
- `~/Documents/GitHub/fs-post-selection/` (cloned; HEAD `bc82bc3`, 2026-09-15) holds `manuscript_jrssb_planning/`, `jrssb_submission/`, `manuscript/`, `submission/` and others; none is named `manuscript_jrssb_initialsubmit`.

The gate is explicit and more than one candidate source exists, so the task stops here rather than substituting a path.

## Record location

§2.5 was not reached. This record is at the task's default, `quarto/simulations/actg175/binary/`. That directory already holds the OR-0.75 sweep documents (`mr_coverage_sweep_or075*.qmd`, `effMaxSG_/maxcons_/maxeffCons_mr_coverage_sweep_or075.qmd`, `_sim_mr_coverage_or075.html`) and `build_actg175_glm_dgm.R`, so it is also where §2.2 would likely have pointed. It has no `current_status.md` generator, so no catalog was regenerated (§8.2).

## Findings

- **F1 (the stop).** `<ms>` is not in forestsearch. The same relative path exists in `fs-glms-interpretable`. To resume, the task needs `<ms>` restated as an absolute path to the intended copy, or the copy placed at the stated path.
- **F2.** The task's `ps | grep` pattern matches this session's own shell wrappers and PID 2086198, a bash `free -m` / `sleep 3` loop from another session, running since about 2026-09-06 (it runs no R). It was left alone.
- **F3.** Nothing in §2–§7 ran. No temporary directory was created, and no R process was started beyond reading `packageDescription()`.

## Commits

```
d526f7ba dev/tasks: TASK_actg175_binary_stage0_2026-09-17_v2.md as received (ACTG175 binary OR study, Stage 0, read-only)
<this record>
```
