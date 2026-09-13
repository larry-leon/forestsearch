# STATUS — p12x20 (Part A: FS effMaxSG eps 0.20, the nine 12.4% cells)

- Generated (UTC): 2026-09-13T13:07:46Z
- Branch: `campaign/p12x20`; HEAD at generation: `c30360f0` (informational; no pin claim on this file).
- Base SHA: `0ab5d1c560cfd0c536efd914f3a6266aa4511bab`.
- Task: `dev/tasks/TASK_p12x20_partA_2026-09-12_v2.md`; addendum `dev/tasks/ADDENDUM_p12x20_stage2_unattended_2026-09-12_v2.md`.
- Stage 0 record: `quarto/simulations/gbsg_020/REPORT_p12x20_stage0_2026-09-12.md`.
- Gate 1 record: `quarto/simulations/gbsg_020/REPORT_p12x20_gate1_2026-09-12.md`.
- Runner: `quarto/simulations/gbsg_020/scripts_p12x20/run_p12x20.sh`; payload directory `quarto/simulations/gbsg_020/p12x20_2026-09-12/`.
- HALT file: absent.
- Open items: none.

## Stage 0 install record (verbatim from the Stage 0 report, 0a-0f)

### 0a Install

- Command: `Rscript -e 'devtools::install(dependencies = FALSE, upgrade = FALSE)'` from the repo root at `f0a71ac3`.
- First attempt aborted before building on an invalid argument supplied by the executor (`upgrade = "never"`: "`upgrade` must be a single TRUE, FALSE, or NA"); nothing was installed; re-run with `upgrade = FALSE`.
- Built `forestsearch_0.3.5.tar.gz`; "installing to library '/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6'"; `* DONE (forestsearch)`; exit 0.
- Installed `DESCRIPTION`: `Version: 0.3.5`; `Packaged: 2026-09-13 06:52:20 UTC`; `Built: R 4.6.1; ; 2026-09-13 06:52:22 UTC; unix`.
- Install target (`.libPaths()[1]` in the installing session): `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6`.

### 0b O-1 forwarding commit

- SHA: `1d9401cb7765a6fdcf8239382a78598ec0c0e1db`.
- Subject: "O-1: .fs_apply_mr() forwards the full MR argument set (DINA + GRF operational)".
- `git merge-base --is-ancestor 1d9401cb 0ab5d1c5`: exit 0 — ancestor — PASS.

### 0c Named token

- Function: `forestsearch:::.fs_apply_mr` (`R/fs_mr_inference_methods.R`, the only R file in the O-1 diff).
- Token 1 (line added by O-1): `.mr_fml <- formals(fs_mr_inference)`.
- Token 2 (argument added by O-1): `field_recovery = .g(mr_inference_args$field_recovery, .d("field_recovery"))`.
- Match rule: `deparse(get(".fs_apply_mr", envir = asNamespace("forestsearch")), width.cutoff = 500)`, collapsed, all whitespace removed, `grepl(fixed = TRUE)` against the whitespace-stripped token.
- Also asserted: `environment(f)` is `asNamespace("forestsearch")` (the closure is the installed namespace's, not a sourced copy).
- Script: `stage0_tokens.R` (session scratchpad), run under `Rscript --vanilla`.

### 0d Fresh session, installed closure

- Before `library()`: `package:forestsearch` on search path FALSE; namespace loaded FALSE.
- `find.package("forestsearch")`: `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6/forestsearch`.
- `packageVersion("forestsearch")`: `0.3.5`.
- `.libPaths()[1]`: `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6`.
- token1_mr_fml: TRUE.
- token2_field_recovery: TRUE.
- fn_env_is_namespace: TRUE.
- Result: PASS.

### 0e Spawned parallel worker

- Mechanism: `future::plan(future::multisession, workers = 2L)`; one `future::future()` that calls `library(forestsearch)` and runs the same token check inside the worker.
- Main PID 111223; worker PID 111274; different: TRUE.
- Worker `packageVersion("forestsearch")`: `0.3.5`.
- Worker `.libPaths()[1]`: `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6`.
- Worker `find.package("forestsearch")`: `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6/forestsearch`.
- Worker token1_mr_fml: TRUE.
- Worker token2_field_recovery: TRUE.
- Worker fn_env_is_namespace: TRUE.
- Result: PASS.

### 0f Worker library matches install target

- Worker `.libPaths()[1]` `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6` vs install target `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6`: `identical()` TRUE — PASS.


## Cells

| cell | HR | n | status | cell wall (s) | replicates | gate counts | cell commit | gate record |
|---|---|---|---|---|---|---|---|---|
| A1 | 1.50 | 500 | done | 1575 | 2000 | run=66 passed=66 failed=0 | acbc02a5 | `GATE2_p12x20_A1_2026-09-12.txt` |
| A2 | 1.50 | 1000 | done | 2191 | 2000 | run=66 passed=66 failed=0 | 0fca1997 | `GATE2_p12x20_A2_2026-09-12.txt` |
| A3 | 1.50 | 1500 | done | 3352 | 2000 | run=66 passed=66 failed=0 | da252976 | `GATE2_p12x20_A3_2026-09-12.txt` |
| A4 | 1.75 | 500 | done | 1618 | 2000 | run=66 passed=66 failed=0 | 880f488c | `GATE2_p12x20_A4_2026-09-12.txt` |
| A5 | 1.75 | 1000 | done | 2213 | 2000 | run=66 passed=66 failed=0 | 6171aef4 | `GATE2_p12x20_A5_2026-09-12.txt` |
| A6 | 1.75 | 1500 | done | 3371 | 2000 | run=66 passed=66 failed=0 | 5b7cfe81 | `GATE2_p12x20_A6_2026-09-12.txt` |
| A7 | 1.00 | 500 | done | 1309 | 2000 | run=66 passed=66 failed=0 | bc5e6a69 | `GATE2_p12x20_A7_2026-09-12.txt` |
| A8 | 1.00 | 1000 | done | 1763 | 2000 | run=66 passed=66 failed=0 | bee11092 | `GATE2_p12x20_A8_2026-09-12.txt` |
| A9 | 1.00 | 1500 | done | 2396 | 2000 | run=66 passed=66 failed=0 | a2ec51c3 | `GATE2_p12x20_A9_2026-09-12.txt` |

- Gate counts over done cells: run=594 passed=594 failed=0.

## Payload inventory (from the directory)

- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_p12x20_combined_1_2000.rds`: 1448192 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_p12x20_res_1_1000.rds`: 731996 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_p12x20_res_1001_2000.rds`: 733460 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A1_2026-09-12.txt`: 7652 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A1_h150_n500_batch_1.html`: 4374858 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A1_h150_n500_batch_1001.html`: 4374284 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A1_h150_n500_combine_1.html`: 4381847 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A1_2026-09-12.md`: 9896 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n1000_nb20_p12x20_combined_1_2000.rds`: 1504071 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n1000_nb20_p12x20_res_1_1000.rds`: 760469 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n1000_nb20_p12x20_res_1001_2000.rds`: 760689 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A2_2026-09-12.txt`: 7660 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A2_h150_n1000_batch_1.html`: 4356879 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A2_h150_n1000_batch_1001.html`: 4366630 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A2_h150_n1000_combine_1.html`: 4377061 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A2_2026-09-12.md`: 9918 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n1500_nb20_p12x20_combined_1_2000.rds`: 1503618 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n1500_nb20_p12x20_res_1_1000.rds`: 759672 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n1500_nb20_p12x20_res_1001_2000.rds`: 760557 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A3_2026-09-12.txt`: 7660 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A3_h150_n1500_batch_1.html`: 4365331 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A3_h150_n1500_batch_1001.html`: 4364501 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A3_h150_n1500_combine_1.html`: 4366141 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A3_2026-09-12.md`: 9918 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n500_nb20_p12x20_combined_1_2000.rds`: 1484801 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n500_nb20_p12x20_res_1_1000.rds`: 751592 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n500_nb20_p12x20_res_1001_2000.rds`: 751701 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A4_2026-09-12.txt`: 7670 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A4_h175_n500_batch_1.html`: 4385713 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A4_h175_n500_batch_1001.html`: 4380982 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A4_h175_n500_combine_1.html`: 4387451 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A4_2026-09-12.md`: 9914 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n1000_nb20_p12x20_combined_1_2000.rds`: 1506351 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n1000_nb20_p12x20_res_1_1000.rds`: 761637 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n1000_nb20_p12x20_res_1001_2000.rds`: 762867 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A5_2026-09-12.txt`: 7677 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A5_h175_n1000_batch_1.html`: 4368468 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A5_h175_n1000_batch_1001.html`: 4357298 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A5_h175_n1000_combine_1.html`: 4374925 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A5_2026-09-12.md`: 9935 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n1500_nb20_p12x20_combined_1_2000.rds`: 1490480 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n1500_nb20_p12x20_res_1_1000.rds`: 753752 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n1500_nb20_p12x20_res_1001_2000.rds`: 752698 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A6_2026-09-12.txt`: 7675 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A6_h175_n1500_batch_1.html`: 4360691 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A6_h175_n1500_batch_1001.html`: 4366606 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A6_h175_n1500_combine_1.html`: 4353522 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A6_2026-09-12.md`: 9933 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_p12x20_combined_1_2000.rds`: 1163469 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_p12x20_res_1_1000.rds`: 583860 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_p12x20_res_1001_2000.rds`: 591748 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A7_2026-09-12.txt`: 7660 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A7_h100_n500_batch_1.html`: 4308551 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A7_h100_n500_batch_1001.html`: 4325876 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A7_h100_n500_combine_1.html`: 4322178 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A7_2026-09-12.md`: 9904 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_nb20_p12x20_combined_1_2000.rds`: 1135842 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_nb20_p12x20_res_1_1000.rds`: 576313 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_nb20_p12x20_res_1001_2000.rds`: 571051 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A8_2026-09-12.txt`: 7658 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A8_h100_n1000_batch_1.html`: 4348524 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A8_h100_n1000_batch_1001.html`: 4317353 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A8_h100_n1000_combine_1.html`: 4328128 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A8_2026-09-12.md`: 9914 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n1500_nb20_p12x20_combined_1_2000.rds`: 1085248 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n1500_nb20_p12x20_res_1_1000.rds`: 545226 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n1500_nb20_p12x20_res_1001_2000.rds`: 550170 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/GATE2_p12x20_A9_2026-09-12.txt`: 7656 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A9_h100_n1500_batch_1.html`: 4327475 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A9_h100_n1500_batch_1001.html`: 4328106 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/p12x20_A9_h100_n1500_combine_1.html`: 4333711 B, tracked
- `quarto/simulations/gbsg_020/p12x20_2026-09-12/REPORT_p12x20_A9_2026-09-12.md`: 9914 B, tracked

## Heartbeat (`LOG_p12x20_progress.txt`, verbatim)

```
2026-09-13T07:37:49Z	A1	start	elapsed_min=0
2026-09-13T08:04:05Z	A1	done	elapsed_min=26	commit=acbc02a5
2026-09-13T08:04:05Z	A2	start	elapsed_min=0
2026-09-13T08:40:37Z	A2	done	elapsed_min=36	commit=0fca1997
2026-09-13T08:40:37Z	A3	start	elapsed_min=0
2026-09-13T09:36:30Z	A3	done	elapsed_min=55	commit=da252976
2026-09-13T09:36:30Z	A4	start	elapsed_min=0
2026-09-13T10:03:29Z	A4	done	elapsed_min=26	commit=880f488c
2026-09-13T10:03:29Z	A5	start	elapsed_min=0
2026-09-13T10:40:23Z	A5	done	elapsed_min=36	commit=6171aef4
2026-09-13T10:40:23Z	A6	start	elapsed_min=0
2026-09-13T11:36:35Z	A6	done	elapsed_min=56	commit=5b7cfe81
2026-09-13T11:36:35Z	A7	start	elapsed_min=0
2026-09-13T11:58:25Z	A7	done	elapsed_min=21	commit=bc5e6a69
2026-09-13T11:58:25Z	A8	start	elapsed_min=0
2026-09-13T12:27:49Z	A8	done	elapsed_min=29	commit=bee11092
2026-09-13T12:27:49Z	A9	start	elapsed_min=0
2026-09-13T13:07:46Z	A9	done	elapsed_min=39	commit=a2ec51c3
```

