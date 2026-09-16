# REPORT — ACTG175 continuous (MD) re-run under the current field constructions: Stage 0, read-only

Date: 2026-09-15. Machine: `pop-os` (AMD Ryzen Threadripper PRO 5995WX, 64 physical cores / 128 logical, 251 GB). Branch `feature/glm-extension`; §1 HEAD `b15c8b7a`; task document committed as `70975bc4`. Task: `dev/tasks/TASK_md_field_rerun_stage0_2026-09-15.md`. Nothing was installed, rendered, simulated or timed; no `R/`, template, script, document, payload or bundle was edited. R was used to read the installed namespace and formals (§2.6, §3) and to summarise the committed `mdf1` bundles with base R (§7.3, §8.1; every such figure is labelled *computed here*).

Marking convention: **confirmed** / **differs** / **not determinable here**, each with `path:line` at HEAD `b15c8b7a` (the template, `R/` and the `mdf1` records are unchanged between `b15c8b7a` and `70975bc4`; the only change is the task document itself).

---

## 1. Provenance and first commit (output as run, verbatim except where marked)

```
pop-os
/home/larryleon/Documents/GitHub/forestsearch
feature/glm-extension
b15c8b7a
--- porcelain (tracked):            [empty: no tracked modifications]
--- log:
b15c8b7a gbsg_020 closeout: regenerate current_status.md at 6ce92a52 (merge of pop-os p12x20 with Mac idsweep)
6ce92a52 status_curated.md reconciled for idsweep and p12x20; current_status_regen.R covers idsweep
6b909809 Merge origin/feature/glm-extension (p12x20 Part A, status_curated restore) into local idsweep
872e06e5 gbsg_020 closeout: regenerate current_status.md at d62e1391
d62e1391 idsweep: 288 cell-runs (18 cells x 16 runs, 500 replicates, MR off), summary and report
behind ahead (vs upstream, as of last fetch):
0	0
-rw-r--r-- 1 larryleon larryleon 1817 Sep 14 14:17 .git/FETCH_HEAD
 21:50:48 up 37 days, 10:49,  1 user,  load average: 0.28, 0.39, 0.40
[ps: two lines are this stage's own §1 shell (pids 1207646, 1207685, etime 00:00) -- elided; the third line, in full:]
2086198 10-08:12:19 /bin/bash -c source /home/larryleon/.claude/shell-snapshots/snapshot-bash-1788626021452-65rwz3.sh 2>/dev/null || true && shopt -u extglob 2>/dev/null || true && { \builtin unalias -- 'unsetenv'; \builtin unset -f -- 'unsetenv'; } >/dev/null 2>&1 || true && eval 'cd ~/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020 && export PATH="/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH" && SCR=/tmp/claude-1000/-home-larryleon-Downloads/f57d1867-8de3-4a39-94eb-0ec33d02da97/scratchpad && \ ( while true; do free -m | awk "/^Mem:/{print \$3}"; sleep 3; done > "$SCR/mem_uburst.txt" ) & MP=$!; \ FS_S7_HR=1.75 FS_S7_NSIMS=100 FS_S7_START=1 FS_S7_FB=none FS_S7_CAMPAIGN=uburst FS_S7_UNIFORM=TRUE FS_S7_WORKERS=100 quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd --output u_burst_h175.html > "$SCR/u_burst_h175.log" 2>&1 && echo "UBURST H175 OK" || echo "UBURST H175 FAIL"; \ FS_S7_HR=1.0 FS_S7_NSIMS=100 FS_S7_START=1 FS_S7_FB=none FS_S7_CAMPAIGN=uburst FS_S7_UNIFORM=TRUE FS_S7_WORKERS=100 quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd --output u_burst_h100.html > "$SCR/u_burst_h100.log" 2>&1 && echo "UBURST H100 OK" || echo "UBURST H100 FAIL"; \ kill $MP 2>/dev/null; echo "peak mem MB: $(sort -n "$SCR/mem_uburst.txt" | tail -1) baseline: $(head -1 "$SCR/mem_uburst.txt")"' < /dev/null && pwd -P >| /tmp/claude-c9ff-cwd
128
               total        used        free      shared  buff/cache   available
Mem:             251          29         213           0          24         221
Swap:             19           1          18
0.3.5
R 4.6.1; ; 2026-09-13 06:52:22 UTC; unix
/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6/forestsearch
R version 4.6.1 (2026-06-24)
/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0
/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0
```

- **GATE:** neither `dev/tasks/TASK_md_field_rerun_stage0_2026-09-15.md` nor this record existed at HEAD `b15c8b7a`. Passed.
- **Running campaign:** pid 2086198 is a 10-day-old bash wrapper from the gbsg_020 `uburst` session; its only live child is the `sleep` of its memory-sampler loop. No `R`, `Rscript` or `quarto` process exists on the machine (`ps -eo comm` check). Reported and left alone.
- **First commit:** `70975bc4 Add TASK_md_field_rerun_stage0_2026-09-15 as received` (the task document alone, from `~/Downloads/TASK_md_field_rerun_stage0_2026-09-15.md`, the exact name, one match).

## 2. S0.1 — Branch, version, and what has changed since `mdf1`

**2.1 Carriers of `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd`** (every local and remote-tracking ref tested with `git cat-file -e <ref>:<path>`): `campaign/p12x20` (`c6d9c82a`), `feature/glm-extension` (`70975bc4`), `origin/campaign/p12x20` (`c6d9c82a`), `origin/feature/glm-extension` (`b15c8b7a`), `origin/feature/glm-extension-mac` (`540c16e8`). No other ref carries it.

**2.2 The Mac branch.** No local ref named `feature/glm-extension-mac`. One remote-tracking ref: `origin/feature/glm-extension-mac` at `540c16e8` (2026-09-08, "ACTG175 continuous intervals Stage 2: the record ..."). `git merge-base --is-ancestor`: ancestor of HEAD **yes**; of `feature/glm-extension` **yes**; of `origin/feature/glm-extension` **yes**. Everything the Mac produced for `mdf1` and the applied document is contained in HEAD.

**2.3 `mdf1`'s run commit and platform.**
- The record (`REPORT_continuous_field_2026-09-07.md:3`) states: "Machine: Mac Studio (M4 Max, 13 workers). Branch `feature/glm-extension-mac`; forestsearch 0.3.5 (with the add-only `scale` argument on `fs_sim_bias_coverage()`, commit `2f118042`)". `2f118042` (2026-09-07) is an ancestor of HEAD. It names the package state, not the template commit.
- The template was added at `c2c402e0` ("Mac continuous/field Stage 1: mr_field MD template ..."), which is the parent of the first bundle-adding commit `c13d0882` (cell 1); the other cells followed at `9114a435`, `5a169691`, `c80a1291` (all 2026-09-07). `git log 2f118042..HEAD -- R/` does not list `c2c402e0`, so `R/` at `c2c402e0` equals `R/` at `2f118042`. **Run commit for the re-run identity gate: `c2c402e0` (template) with the package at `2f118042` (same `R/`).** Source: bundle-adding commits' parent.
- Bundle meta (computed here from the four `*_combined_1_2000.rds`): `pkg_version 0.3.5`, `r_version 4.5.2`, `hostname Mac-Studio-3.local`, `built_at` 2026-09-07 13:27:23 / 14:14:02 / 14:53:58 / 15:41:22 (md40 n500 / md120 n500 / null n500 / md40 n700), `n_workers_by_batch 13/13`. Stage 0 record (`REPORT_continuous_field_stage0_2026-09-07.md:93`): "Apple M4 Max, 14 physical cores (10 performance + 4 efficiency), 14 logical, 36 GB RAM. R 4.5.2, Quarto 1.10.18".
- **BLAS: not recorded** in the `mdf1` record, its stage records or the bundle meta. The applied record for the same machine one day later (`quarto/applications/actg175/REPORT_actg175_continuous_intervals_2026-09-07.md:13`) states "R's BLAS here is `/System/Library/Frameworks/Accelerate.framework/.../libBLAS.dylib`" (Mac-Studio-3). This machine: reference BLAS/LAPACK 3.12.0 (§1).

**2.4 `R/` since `mdf1`** (`git log --oneline 2f118042..HEAD -- R/`; 4 files, +637/−44: `forestsearch_main.R`, `forestsearch_methods.R`, `fs_mr_inference.R`, `fs_mr_inference_methods.R`). Classification against what the template computes on the continuous FS path, given the arguments the template passes (§4.2):

| commit | what | can change a `mdf1` column? |
|---|---|---|
| `7245e898` | `field_decompose` (default FALSE), add-only complement-scale diagnostics | **cannot** (default-inert; not passed by the template) |
| `3880022e` | `field_scale_complement = c("none", "selected")` (default "none" at that commit): the field-s construction, `_s` companions, `joint_s` | **cannot at its own default**; the code path it added now runs under `fb62705c` (below). Its commit message: "G2b (decompose + scale) 131/131 identical" (survival rows). Its invariance argument is in code, `R/fs_mr_inference.R:1198–1207`: "on every path that reaches here Bc / fit_ok are exactly what the loop left, and Zo_c / Zi_c / lam_c below are unchanged (gates G2a-G2c)". The field-s root is computed beside, not instead: `:1223–1224` `lam_c <- rep(NA_real_, R_out); lam_cs <- rep(NA_real_, R_out)`; `:1234` `lam_c[r] <- Zo_c[G, r] - mean(...)`; `:1237–1239` `if (scale_on) lam_cs[r] <- (s[sel] / s[G]) * Zo_c[G, r] - mean(...)`. No draw is consumed by the scale step (`s <- sqrt(colSums(Bc * Bc))`, `:1219`). |
| `c0f48a7c` | merge of `origin/feature/glm-extension-mac` | no content of its own |
| `fb62705c` | **defaults flip**: `return_reselection` FALSE→TRUE, `field_complement` FALSE→TRUE, `field_scale_complement` "none"→"selected", in both `fs_mr_inference()` and `forestsearch()`'s fallbacks | **can change what is computed** (the field-s step now runs for the template, which does not pass `field_scale_complement`); **cannot change any recorded `mdf1` column** (the template passes `field_complement` and `return_reselection` explicitly at the same values `mdf1` used; the field-s step is add-beside per `3880022e`). Hunks: `R/forestsearch_main.R` `-        return_reselection = .g_mr(mr_inference_args$return_reselection, FALSE)` / `+        ... TRUE)`; `-        field_complement = .g_mr(mr_inference_args$field_complement, FALSE)` / `+        ... TRUE)`; `-        field_scale_complement = .g_mr(mr_inference_args$field_scale_complement, "none")` / `+        ... "selected")`; `R/fs_mr_inference.R` `-                           return_reselection = FALSE,` / `+ ... TRUE,`; `-                           field_complement = FALSE,` / `+ ... TRUE,`; `-                           field_scale_complement = c("none", "selected"),` / `+                           field_scale_complement = c("selected", "none"),`. |
| `beda000b` | `ci_method` default "ij"→"field" (`fs_mr_inference` formal order and the `forestsearch()` fallback) | **cannot**: the template passes `ci_method` explicitly (`FS_MD_CI`, default "field"). Hunk: `R/forestsearch_main.R` `-        ci_method     = .g_mr(mr_inference_args$ci_method,   "ij"),` / `+        ci_method     = .g_mr(mr_inference_args$ci_method,   "field"),`; `R/fs_mr_inference.R` `-                           ci_method = c("ij", "wald", "field"),` / `+                           ci_method = c("field", "ij", "wald"),`. |
| `df251a96`, `105186cd`, `0f1eb18c` | `print.forestsearch()` / `summary.forestsearch()` reporting | **cannot** (`R/forestsearch_methods.R` only; the template reads `fs.est$mr_inference` directly) |
| `e819a58a` | `field_recovery` (default FALSE): membership-agreement diagnostics | **cannot**: the only touched existing lines are the `G_out` allocation (`G_out <- if (rec_g) ...` with `rec_g <- fc || rec_on`; `fc` is TRUE for the template, so unchanged) and three assignments (`kept`, `Ncol`, `Nall`) hoisted out of a block with the same values; commit message "byte-identical defaults" |
| `1d9401cb` | `.fs_apply_mr()` forwards the full MR argument set (DINA + GRF) | **cannot** (`R/fs_mr_inference_methods.R`; the FS path calls `fs_mr_inference()` from `forestsearch()` directly, `R/forestsearch_main.R:3395`) |

Consequence for a re-run identity gate: with the template as committed, **every recorded `mdf1` column** (the 138 recorder columns, §5.1) is expected to reproduce; nothing on the FS continuous path changes at the template's passed arguments. Two caveats belong to the gate, not to `R/`: (i) the field-s output-invariance evidence (`3880022e` G2a–G2c) was gathered on survival rows, while the code is outcome-agnostic (§3); (ii) cross-platform (Mac `mdf1` → this machine) reproduction is to floating point, not bit-for-bit — the `mdf1` Gate 2 pairing used ≤ 1e-8 relative (`REPORT_continuous_field_gate2_2026-09-07.md:5`) and the applied record found max |d| 1e-13 Linux↔Mac (§9.3) — and R here is 4.6.1 against 4.5.2.

**2.5 Loading.** Template `sim_fs_maxeffCons_mr_field_md_template.qmd:96–99`:
```
# forestsearch must be INSTALLED (devtools::install()), not load_all(): the
# doFuture workers spawn separate R processes that only see the installed
# package.
library(forestsearch)
```
`scripts_mdf1/`: `library(forestsearch)` in `identity_postchange.R:1`, `fixture_check.R:3`, `fixture_baseline_save.R:4`, `reselection_check.R:3`, `stage3_aggregate.R:4`; the shell drivers call `quarto render` (`run_cell.sh:20`, `cal_renders.sh:7`, `id_renders.sh:9`). No `load_all()` and no `devtools::install()` appears in the template or in any script (the install is a precondition stated in the comment, not performed by them).

**2.6 Installed package against source** (`pkgload` is installed; `git archive HEAD DESCRIPTION NAMESPACE R` — the package has no `src/` — into a `tempfile()` directory, 111 files; one R session; `formals()`/`body()` after `utils::removeSource()`; extract removed afterwards): installed closures **666**, source closures **666**, common 666, **identical 666, differing 0**, installed-only 0, source-only 0. Installed build: 0.3.5, `Built: R 4.6.1; ; 2026-09-13 06:52:22 UTC; unix`. The tip of `feature/glm-extension` is HEAD (`70975bc4`, this session's task-document commit); `origin/feature/glm-extension` (`b15c8b7a`) differs from it only by the task document, with no `R/` difference, so no second comparison was needed. **Confirmed: the installed package is the source at HEAD.**

## 3. S0.2 — The loaded engine

Installed namespace (`formals(get("fs_mr_inference", asNamespace("forestsearch")))`; `fs_mr_inference` is not exported) and source `R/fs_mr_inference.R:549–570` agree:
- `ci_method = c("field", "ij", "wald")` (`:559`) — **confirmed**; `field_complement = TRUE` (`:566`) — **confirmed**; `field_scale_complement = c("selected", "none")` (`:568`) — **confirmed**; `return_reselection = TRUE` (`:561`) — **confirmed**; `field_R_out = 1000L` (`:562`), `field_R_in = 500L` (`:563`) — **confirmed**; `include_complement = FALSE` (`:558`) — **confirmed**.
- Multiplier-draw and seed formals: `draws = 2000L` (`:556`), `multiplier = c("poisson", "gaussian", "rademacher")` (`:557`), `seed = NULL` (`:560`); the field block re-seeds at `seed + 900000L` (`:864`, `seed_offset = 900000L` `:938`), the kappa sweep at `seed + 910000L` (`:981`, `:995`).
- **`forestsearch()` FS branch** (`R/forestsearch_main.R:3395–3453`, and the installed body read back): `include_complement = .g_mr(mr_inference_args$include_complement, TRUE)` (`:3412`), `ci_method = .g_mr(mr_inference_args$ci_method, "field")` (`:3421`), `field_complement = .g_mr(mr_inference_args$field_complement, TRUE)` (`:3436`), `field_scale_complement = .g_mr(mr_inference_args$field_scale_complement, "selected")` (`:3444`) — **confirmed** (these are the fallbacks when `mr_inference_args` omits the element). MR on a GLM outcome requires `consistency_method = "resample"`: `:3318` `.mr_glm_ok <- consistency_method == "resample" && !is.null(estimator_fn)`; `:3332–3336` the skip message "MR on a GLM outcome requires consistency_method = \"resample\"" — **confirmed**.
- **No outcome-type guard in the field code** — **confirmed**: the only `outcome_type`/`survival` test in `R/fs_mr_inference.R` is the pieces dispatch `:61–73` (`.fs_mr_pieces`: `if (identical(spec$outcome_type, "survival")) .consistency_cox_pieces(...) else .consistency_glm_pieces(df_sub, outcome_type = ..., effect_measure = ..., ..., adverse_outcome = if (is.null(spec$adverse_outcome)) TRUE else spec$adverse_outcome)`). Each candidate's coefficient and influence contributions on the GLM path come from `.consistency_glm_pieces()` (`R/consistency_resample.R:240`) — **confirmed**. MD bounds on the identity scale: `.consistency_glm_pieces()` returns `log_scale = effect_measure %in% c("OR", "RR", "IRR")` (`R/consistency_resample.R:323`), and `fs_mr_inference()` sets `to_eff <- function(x) if (log_scale) exp(x) else x` (`:607`) — **confirmed**. Under `adverse_outcome = FALSE` the refit negates Y: `R/consistency_resample.R:255–257` `if (!adverse_outcome) { if (outcome_type == "continuous") { df[[outcome.name]] <- -df[[outcome.name]] }` — **confirmed**; every bound is `to_eff(<oriented beta> - <quantile>)`, i.e. on the harm-oriented MD scale. (`forestsearch()` resolves `adverse_outcome` to FALSE for continuous when NULL, `R/forestsearch_main.R:1718–1720`; the template passes FALSE explicitly, §6.1.)
- **Returned** — **confirmed**: `field$lower_1s = to_eff(beta_deb - qs[5])` (`:932`); complement (`.fs_mr_field_complement`, `:1293–1316`): `upper_1s = to_eff(bdc - qs[1])` (`:1300`), and under `scale_on` the companions `lambda_mean_s`, `se_field_s = sd_cs`, `est2_s = to_eff(est2s_w)`, `upper_1s_s = to_eff(bdc - qss[1])`, `lower_1s_s`, `lower_2s_s = to_eff(bdc - qss[7])`, `upper_2s_s = to_eff(bdc - qss[6])`, `lower_se_s`, `upper_se_s` (`:1309–1315`). `upper_1s` and `upper_1s_s` are inverted about the same `bdc` from `lf <- lam_c[ok_c]` and `lfs <- lam_cs[ok_c]` over the same `ok_c` draws (`:1252`, `:1261`) — **confirmed, same draws**. `field$joint_s` (`:1291–1292`, attached at `:961`) is `.fs_mr_field_joint(lam_H[ok_c], lfs, ...)`, returning `gamma, joint_prob, alpha, lower_H, upper_Hc, bonf_gamma, bonf_lower_H, bonf_upper_Hc, bonf_joint_prob, corr, n_joint_draws, grid_gamma, grid_joint_prob` (`:1361–1367`) — **confirmed**. **The unstudentized `field$joint` is still returned** (`:1289–1290`, attached at `:960` `if (!is.null(fcres$joint)) field$joint <- fcres$joint`) — **confirmed**.

## 4. S0.3 — The template's MR call

**4.1** MR runs **inside `forestsearch()`** (`mr_inference = TRUE`, `mr_inference_args = mr_inference_args`); the template contains no direct `fs_mr_inference()` call. Identification + MR call, `sim_fs_maxeffCons_mr_field_md_template.qmd:590–613`:
```
    fs.est <- suppressWarnings(forestsearch(
      df.analysis = df, confounders.name = confs,
      outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
      outcome_type = "continuous", effect_measure = "MD",
      effect.threshold = md_threshold, consistency.threshold = md_consistency,
      pconsistency.threshold = pconsistency, fs.splits = fs_splits,
      n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
      vi.grf.min = vi_grf_min, sg_focus = sg_focus,
      selection_rule = selection_rule,
      effect_neighborhood = effect_neighborhood,
      stop_threshold = stop_threshold,
      consistency_method = consistency_method,
      conf.cont_jcuts = fs_conf.cont_jcuts,
      use_lasso = use_lasso, use_dina = use_dina, use_grf = use_grf,
      use_twostage = use_twostage, is.RCT = is_rct,
      adverse_outcome = adverse_outcome,
      details = FALSE, quiet = FALSE, seedit = sd_i,
      parallel_args = inner_parallel,
      mr_inference = TRUE,
      ...
      mr_inference_args = mr_inference_args)),
```
with `mr_inference_args` built at `:282–287`:
```
mr_inference_args <- c(
  list(ci_method = mr_ci_method, draws = mr_draws, include_complement = TRUE,
       confirm_rule = mr_confirm_rule, field_uniform = mr_field_uniform,
       field_complement = mr_field_complement, ij_residual = mr_ij_residual,
       return_reselection = mr_return_reselection),
  if (!is.null(mr_t_confirm)) list(t_confirm = mr_t_confirm))
```

**4.2** Passed or inherited:
- `ci_method`: **passed**, `mr_ci_method <- .env_chr("FS_MD_CI", "field")` (`:263`), so `"field"` unless `FS_MD_CI=ij`.
- `include_complement`: **passed**, literal `TRUE` (`:283`).
- `field_complement`: **passed**, `mr_field_complement <- identical(.env_chr("FS_MD_FIELD_COMPLEMENT", "TRUE"), "TRUE")` (`:271`), so `TRUE`.
- `field_scale_complement`: **not passed** (absent from `:282–287`) → **inherited from `forestsearch()`'s fallback** `"selected"` (`R/forestsearch_main.R:3444`), which `fs_mr_inference()` accepts by `match.arg` (`R/fs_mr_inference.R:573`).
- `return_reselection`: **passed**, `mr_return_reselection <- identical(.env_chr("FS_MD_RESELECTION", "TRUE"), "TRUE")` (`:278`), so `TRUE`.
- **Plainly:** the template as committed would run today with `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE` and `field_scale_complement = "selected"`: the harm field, **both** the unstudentized complement field and the studentized complement field (field-s), and **both** `field$joint` and `field$joint_s` are computed; the recorder captures only the harm field, the unstudentized complement field and the unstudentized joint pair (§5). `mdf1` ran the same call under the then-default `"none"` (`3880022e`), i.e. without field-s.

**4.3 `FS_MD_*` knobs** (line, default, effect): `FS_MD_WORKERS` (`:116`; physical cores − 1, capped there; worker count), `FS_MD_NSIMS` (`:119`; 1000; replicates in this batch), `FS_MD_KNOISE` (`:123`; 0; noise confounders), `FS_MD_MODE` (`:128`; "batch"; batch | combine), `FS_MD_START` (`:130`; 1; first `sim_id`), `FS_MD_MD` (`:142`; "40"; cell: 40 | 120 | null → `target_md_harm`, `null_cell`), `FS_MD_N` (`:146`; 500; `n_sample`), `FS_MD_BUILD` (`:151`; "direct" if MD 120 else "calibrate"; DGM route), `FS_MD_CAMPAIGN` (`:167`; "c1"; campaign tag in the stem and results dir), `FS_MD_QUICKRUN` (`:174`; FALSE; `_quickrun` stem), `FS_MD_SAVE_COMBINED` (`:188`; TRUE; write the pooled bundle), `FS_MD_CI` (`:263`; "field"; `ci_method`), `FS_MD_UNIFORM` (`:269`; FALSE; kappa sweep — excluded by convention 6), `FS_MD_FIELD_COMPLEMENT` (`:271`; TRUE; `field_complement`), `FS_MD_IJ_RESIDUAL` (`:274`; "two_term"; `ij_residual`), `FS_MD_RESELECTION` (`:278`; TRUE; `return_reselection`), `FS_MD_WINNER_ROWS` (`:281`; FALSE; show excluded rows), `FS_MD_FB` (`:294`; "none"; none | join), `FS_MD_FB_PATH` (`:296`; the committed FB bundle path), `FS_MD_JOIN_SKIP` (`:302`; ""; sim_ids to skip in the join). No knob sets `sg_focus`, `effect_neighborhood`, `selection_rule`, `field_scale_complement`, `mr_draws` (5000, `:122`) or the thresholds. Campaign-tag guard, `:167–171`:
```
campaign_tag <- .env_chr("FS_MD_CAMPAIGN", "c1")
stopifnot(grepl("^[A-Za-z0-9]+$", campaign_tag))
rds_stem <- sprintf("%s_%s_mr_field_md%02d_knoise%d_n%d_%s",
                    method_tag, focus_tag, abs(target_md_harm),
                    k_random_noise, n_sample, campaign_tag)
```

## 5. S0.4 — What it saves per replicate

**5.1** Per-replicate recorder: `.na_record()` (`:460–538`) builds **a fixed set of named columns** (`sim_id status detected mr_ok err_msg mr_msg n_sel n_harm n_true sg_def covs betaHhat_H betaHhat_Hc fb_secs fit_mr_secs fb_err fb_src1 fb_src2 fb_nres`; `or_H_{est,lo,hi,se}`, `or_Hc_*`; `nv_H_*`, `nv_Hc_*`; `mr_H_{est,lo,hi,se_ij}`, `mr_Hc_*`; then `label`, `mr_harm_flag`, `fld_H_*` (23), the eight uniform columns, `fld_Hc_*` (24), the twelve winner-variant columns, `fld_joint_*` (9), `p_hat_H`, `p_top1..3`, `p_lab1..3`), filled by `record_replicate()` (`:563–767`): naive `:646–649`, IJ `:650–653`, complement naive/IJ `:655–662`, field `:693–703`, complement field `:722–737`, joint `:740–751`, re-selection `:754–763`. Bundles have **138 columns** (computed here). Batch save `:890–929`:
```
out <- rds_path
.payload <- list(results = results, truth = truth,
             meta = list(n_sample = n_sample, n_sims = length(sim_ids), ... field_R = "R_out=1000/R_in=500 (package defaults; not forwarded by forestsearch())", ... pkg_version = ..., r_version = ..., hostname = ..., built_at = Sys.time()))
if (!is.null(dgm_scale)) .payload$scale <- dgm_scale
.payload$oc    <- fs_mr_oc_summary(.payload)
.refuse_if_tracked(out)
saveRDS(.payload[intersect(c("results", "truth", "scale", "oc", "meta"),
                           names(.payload))], out)
```
It saves the results data frame plus `truth`, `scale`, `oc`, `meta` — **not** the `forestsearch` result object. `meta` (`:892–921`) records `ci_method`, `field_uniform`, `field_complement`, `ij_residual`, `return_reselection`, `field_R` but **not `field_scale_complement`** (finding: the combine-mode poolability keys `:945–948` also omit it).

**5.2** `upper_1s_s`, `se_field_s`, `est2_s`, `lower_2s_s`, `upper_2s_s` and the `joint_s` fields are **not saved** (`fld_Hc_up1s_s` absent from all four bundles, computed here; recorder `:504–518`, `:526–530` carry no `_s` column). The template edit that would add them, **listed, not made**, transplanted from the survival template of record `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (last touched at `8fd89e1d`):
- (a) after `:518` (`fld_Hc_note = NA_character_,`) insert the source's `:931–936`:
  ```
    # same mr_Hc_est.  PRIMARY: fld_Hc_up1s_s, the one-sided 95% upper bound.
    fld_Hc_est2_s = NA_real_,
    fld_Hc_up1s_s = NA_real_, fld_Hc_lo1s_s = NA_real_,
    fld_Hc_lo2s_s = NA_real_, fld_Hc_hi2s_s = NA_real_,
    fld_Hc_lo_se_s = NA_real_, fld_Hc_hi_se_s = NA_real_,
    fld_Hc_se_s = NA_real_, fld_Hc_lam_mean_s = NA_real_,
  ```
- (b) after `:530` (`fld_joint_corr = NA_real_, fld_joint_n = NA_integer_,`) insert the source's `:973–978`:
  ```
    # joint_s mirrors (field-s; else all NA).
    fld_joint_s_gamma = NA_real_, fld_joint_s_prob = NA_real_,
    fld_joint_s_loH = NA_real_, fld_joint_s_upHc = NA_real_,
    fld_joint_s_bonf_loH = NA_real_, fld_joint_s_bonf_upHc = NA_real_,
    fld_joint_s_bonf_prob = NA_real_,
    fld_joint_s_corr = NA_real_, fld_joint_s_n = NA_integer_,
  ```
- (c) after `:733` (`rec$fld_Hc_secs <- fc$timing_seconds %||% NA_real_`, inside the `is.null(fc$note)` block) insert the source's `:1260–1268`:
  ```
          rec$fld_Hc_est2_s     <- fc$est2_s      %||% NA_real_
          rec$fld_Hc_up1s_s     <- fc$upper_1s_s  %||% NA_real_
          rec$fld_Hc_lo1s_s     <- fc$lower_1s_s  %||% NA_real_
          rec$fld_Hc_lo2s_s     <- fc$lower_2s_s  %||% NA_real_
          rec$fld_Hc_hi2s_s     <- fc$upper_2s_s  %||% NA_real_
          rec$fld_Hc_lo_se_s    <- fc$lower_se_s  %||% NA_real_
          rec$fld_Hc_hi_se_s    <- fc$upper_se_s  %||% NA_real_
          rec$fld_Hc_se_s       <- fc$se_field_s  %||% NA_real_
          rec$fld_Hc_lam_mean_s <- fc$lambda_mean_s %||% NA_real_
  ```
- (d) after `:751` (the closing `}` of the `f$joint` block) insert the source's `:1305–1316`:
  ```
        # ---- joint_s (field-s; absent -> all NA) ----
        if (is.list(f$joint_s) && is.null(f$joint_s$note)) {
          js <- f$joint_s
          rec$fld_joint_s_gamma     <- js$gamma
          rec$fld_joint_s_prob      <- js$joint_prob
          rec$fld_joint_s_loH       <- js$lower_H
          rec$fld_joint_s_upHc      <- js$upper_Hc
          rec$fld_joint_s_bonf_loH  <- js$bonf_lower_H
          rec$fld_joint_s_bonf_upHc <- js$bonf_upper_Hc
          rec$fld_joint_s_bonf_prob <- js$bonf_joint_prob
          rec$fld_joint_s_corr      <- js$corr
          rec$fld_joint_s_n         <- as.integer(js$n_joint_draws)
        }
  ```
- Companion lines a Stage 1 would consider with it (listed only): the interval invariant pairs at `:824–827` (`c("fld_Hc_lo2s_s","fld_Hc_hi2s_s")`, `c("fld_Hc_lo1s_s","fld_Hc_up1s_s")`), a `field_scale_complement` entry in `meta` (`:907–913`) and in the poolability keys (`:945–948`), and a `FS_MD_FIELD_SCALEC` knob mirroring the survival template's, so the construction is recorded rather than inherited. This is a template change, not an `R/` change.

**5.3** `.refuse_if_tracked()` (`:197–207`):
```
.refuse_if_tracked <- function(path) {
  tracked <- tryCatch(
    system2("git", c("ls-files", "--error-unmatch", shQuote(path)),
            stdout = FALSE, stderr = FALSE) == 0L,
    error = function(e) FALSE)
  if (isTRUE(tracked))
    stop("refusing to overwrite the git-tracked bundle '", path,
         "': committed bundles are read-only -- use a new campaign_tag.",
         call. = FALSE)
  invisible(path)
}
```
Live on the batch save (`:925 .refuse_if_tracked(out)` immediately before `saveRDS`, `:926`) and on the combine save (`:1006`). **Confirmed.** Every `mdf1` bundle path is tracked, so a re-run needs a new `FS_MD_CAMPAIGN` tag.

**5.4** `fs_sim_bias_coverage()` (`R/fs_bias_coverage.R:73–79`): `function(results, block = c("H", "Hc"), estimators = c("naive", "mr", "fld"), level = 0.95, target = "betaHhat", side = c("lower", "upper"), scale = c("log", "identity"))` — signature **confirmed**. Estimator-to-column mapping: `estimators <- match.arg(estimators, c("naive", "mr", "fld", "mr_w", "mr_wf"), several.ok = TRUE)` (`:85–86`); `fld_pre <- paste0("fld_", block, "_")` (`:90`); `cols("fld")` reads `f("est2"), f("lo2s"), f("hi2s"), f("se")` and `b1 = if (side == "lower") f("lo1s") else f("up1s")` (`:125–129`); `naive`/`mr` read `nv_|mr_<block>_est`, `_se`/`_se_ij`, `_lo`, `_hi` (`:133–148`). **It cannot read the field-s complement columns under the recorder's naming** (`fld_Hc_up1s_s`, `fld_Hc_est2_s`, ...) without an `R/` change: the estimator vocabulary has no field-s member and the prefix/suffix set is fixed. What is missing: an estimator key (e.g. a sixth vocabulary entry) mapping to the `_s` suffix. A document-level workaround (a copy of `results` with the `fld_Hc_*_s` columns renamed to `fld_Hc_*`) would reach it with no `R/` change; that is a summary-document choice, not drafted here. **Finding.**

**5.5** Transplant sources (paths only). Summary of `mdf1` (tables and the bias–coverage display): `quarto/simulations/actg175/continuous/summary_continuous_field_mdf1.qmd` (+ `.html`), the four `fs_maxeffCons_mr_field_<cell>_mdf1_combine_1_2000.html`, `scripts_mdf1/stage3_aggregate.R` (writes `stage3_tables.rds`, a list of wide tables, and the two display PNGs), `fig_mdf1_bias_coverage_display_H.png`, `fig_mdf1_bias_coverage_display_Hc.png`. **Long-format metrics file (one row per cell × block × estimator × metric with Monte Carlo SEs): none located** — a search of committed `quarto/simulations/**` `.qmd`/`.R` files for a writer (`fwrite`/`write.csv`/`saveRDS`) of a `metric`-keyed long table with MC SEs found no such document, and `git ls-files` lists no `*metrics*.csv|rds` or `*long*.csv|rds` under `quarto/simulations`. **Not determinable here** beyond that negative search; the re-run's extract would be new.

## 6. S0.5 — Settings, estimators, seeds, and the selection-rule question

**6.1 Identification, as the template sets it** (`sim_fs_maxeffCons_mr_field_md_template.qmd`):
- `sg_focus <- "maxeffCons"` (`:135`, hard-coded; `subgroup_method <- "consistency"` `:134`).
- `effect_neighborhood <- 0.10` (`:214`, hard-coded; comment "inert under maxcons-family foci; pinned"); `selection_rule <- "neighborhood"` (`:213`). The `maxeffCons` focus does not consult the band: `sort_subgroups()` applies `.compute_inclusion_band()` only `if (sg_focus %in% c("hrMaxSG", "hrMinSG"))` (`R/subgroup_consistency_helpers.R:587–606`), and MR's `.fs_mr_select()` calls `.inband()` only for `effMaxSG`/`effMinSG` (`R/fs_mr_inference.R:173–174`); `maxcons = passers[which.max(zcons[passers])]` (`:169`).
- Thresholds: `md_threshold <- 30` passed as `effect.threshold` (`:216`, `:594`); `md_consistency <- 10` passed as `consistency.threshold` (`:218`, `:594`); `pconsistency <- 0.90` passed as `pconsistency.threshold` (`:220`, `:595`).
- Cut grid: `fs_conf.cont_jcuts <- list(age = 10, preanti = 10)` passed as `conf.cont_jcuts` (`:227`, `:602`); `cut_type` is **not set** by the template → `forestsearch()`'s default `cut_type = "default"` (installed formals).
- `maxk <- 2L`, `n_min <- 60L`, arm minima `d0_min <- 12L; d1_min <- 12L` (`:221`, passed `:596`); `fs_splits <- 400L` (`:221`, `:595`).
- `consistency_method <- "resample"` (`:210`, `:601`).
- `adverse_outcome <- FALSE` (`:230`, `:605`).
- Effect screen: `effect.threshold = 30` on the oriented MD scale (comment `:216–217`: "true clearance <=> P(Q|H) >= 27.2%").
- `stop_threshold <- NULL` (`:215`, `:600`; "pinned so the resolved setting is visible").
- Also fixed: `use_lasso = use_dina = use_grf = FALSE`, `use_twostage = TRUE`, `is.RCT = TRUE`, `vi.grf.min = -0.2` (`:228–229`); `include_str2 <- TRUE` (`:248`); analysis covariates `age, preanti, wtkg, karnof, cd40, cd80` + `hemo, homo, drugs, race, gender, symptom, str2` (`:244–250`).

**6.2 The four cells.** DGM per cell (`:343–374`): the null cell `generate_glm_dgm(data = actg_df, factor_vars = z1..z12, outcome_var = "cd4_change", treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD", subgroup_vars = c("z1","z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L), model = "null", n_super = 5000L, seed = 8316951L)`; MD 120 (`dgm_build = "direct"`) `generate_glm_dgm(..., model = "alt", k_treat = 1, k_inter = stage2_forecast$k_inter, adverse_outcome = FALSE, ...)` with `k_inter` read from `dev/glm-continuous-sims/oc_breadth_ladder_2026-08-30_forecast120.rds` and asserted `abs(k_inter + 93.7447641240) < 1e-9` (`:351–356`); MD 40 (n 500 and n 700; `dgm_build = "calibrate"`) `calibrate_glm_interaction(..., target_effect = -40, k_inter_range = c(0, 120), grid_step = 2, n_super = 5000L, seed = 8316951L)`. ACTG175 arms 1/3, treat arm 1, age cut 34, preanti cut 744.5 (`:233–236`). The DGM is n-invariant; `n_sample` = 500 or 700. Recorded truths, harm-oriented (bundle `truth`, computed here; the record's line 7 agrees): planted region β(Q) = **+40** (md40 cells), **+120** (md120), undefined (null: `effect_Q NA`, `prevalence_Q 0`, `beta_inter 0`); complement β(Qᶜ) = **+26.255** in every cell; overall (ITT) = **+30.992** (md40), **+58.560** (md120), **+26.255** (null). Prevalence of Q 0.3446. Fitted (one super-population draw) `marg_H` +35.83 / +115.83 / NA, `marg_Hc` +22.54 / +22.54 / +22.41.

**6.3 Estimators, per replicate** (oriented, positive = harm; `Ĥ` = the selected subgroup):
- Unadjusted (naive): `rec$nv_H_est <- g$naive$est`, `nv_H_lo/hi <- g$naive$lower/upper`, `nv_H_se <- g$debiased$se_wald` (`:646–649`); complement `:655–658`. **Depends on Ĥ** (plug-in on the selected region).
- Oracle: `.oracle_md_on()` (`:541–555`) fits `lm(y_sim ~ treat_sim)` on `df[keep, ]` with `keep = df$flag_harm == 1L` (`.oracle_md`, `:556`) or `== 0L` (`:557`), `s <- if (isTRUE(adverse_outcome)) 1 else -1`, `c(est = s * est, lo = s * est - 1.96 * se, hi = s * est + 1.96 * se, se = se)`; recorded `:582–586`. **Does not depend on Ĥ** (refit on the true region Q / Qᶜ; NA under the null where Q is empty).
- IJ two-term: `rec$mr_H_est <- g$debiased$est`, `mr_H_lo/hi <- g$debiased$lower/upper`, `mr_H_se_ij <- g$debiased$se_ij` (`:650–653`); complement `:659–662`. **Depends on Ĥ.**
- Field: `rec$fld_H_est2 <- f$est2; fld_H_lo2s/hi2s <- f$lower_2s/upper_2s; fld_H_lo1s <- f$lower_1s; fld_H_se <- f$se_field; ...` (`:693–703`). **Depends on Ĥ.**
- Complement field: `rec$fld_Hc_est2 <- fc$est2; fld_Hc_up1s <- fc$upper_1s; fld_Hc_lo1s <- fc$lower_1s; fld_Hc_lo2s/hi2s <- fc$lower_2s/upper_2s; fld_Hc_se <- fc$se_field; ...` (`:722–737`); joint `:740–751`. **Depends on Ĥ** (through Ĥᶜ).

**6.4 Seeds.** `seed_base <- 8316951L` (`:131`); per replicate `sd_i <- seed_base + sim_id` with `RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)` (`:565–568`); data `simulate_from_glm_dgm(dgm, n = n_sample, seed = sd_i)` (`:570`); noise (if any) `set.seed(sd_i + 10^7)` (`:575`); search `seedit = sd_i` (`:606`); MR seed = `seedit` (`R/forestsearch_main.R:3453` `seed = .g_mr(mr_inference_args$seed, seedit)`; the template passes no `seed` in `mr_inference_args`); field offset `+900000L`, kappa `+910000L` (§3). `sim_id` range: `sim_id_start = FS_MD_START` to `sim_id_start + n_sims - 1` (`:182`); `mdf1` ran two batches per cell, `sim_id` 1–1000 and 1001–2000, then a combine (`scripts_mdf1/run_cell.sh:29–31`); bundle meta (computed here): `seed_base 8316951`, `sim_id 1–2000`, `n_batches 2`, sources `*_res_1_1000.rds`, `*_res_1001_2000.rds`. Multiplier draws `mr_draws <- 5000L` (`:122`; meta `mr_draws 5000`); field `R_out/R_in` 1000/500 (package defaults, not forwarded; meta `field_R`).

**6.5 `effMaxSG` at ε = 0.20 on this path** (facts only):
- (i) `sg_focus` and `effect_neighborhood` are **not** set by `FS_MD_*` knobs: `:135` `sg_focus <- "maxeffCons"` and `:214` `effect_neighborhood <- 0.10` are literals. The template edit that would make them so, in the file's own idiom (listed, not made): `sg_focus <- .env_chr("FS_MD_FOCUS", "maxeffCons")` at `:135` and `effect_neighborhood <- .env_num("FS_MD_NBHD", 0.10)` at `:214` (`.env_num` exists, `:107`); `focus_tag <- forestsearch::fs_focus_tag(subgroup_method, sg_focus)` (`:157`) then re-tags the output stem, and `meta$sg_focus`/`meta$selection_rule` (`:894–895`) already carry the values; `effect_neighborhood` is not in `meta` (a finding for that edit). The combine-mode poolability keys include `sg_focus` (`:946`).
- (ii) Identifier side, continuous outcome: `sort_subgroups()` (`R/subgroup_consistency_helpers.R:587–612`): `hr_vec <- as.numeric(result_new$hr); if (isTRUE(effect_log_scale)) hr_vec <- exp(hr_vec)`; `in_band <- .compute_inclusion_band(hr_vec = hr_vec, n_vec = N_vec, selection_rule = selection_rule, effect_neighborhood = effect_neighborhood)`; `ord <- if (sg_focus == "hrMaxSG") order(-in_band, -N_vec, -Pcons_vec, -hr_vec, K_vec)`. The band (`:778–799`): `hr_max <- max(hr_vec, na.rm = TRUE); hr_floor <- (1 - effect_neighborhood) * hr_max; as.integer(!is.na(hr_vec) & hr_vec >= hr_floor)`. `effect_log_scale` defaults FALSE (`R/subgroup_consistency_main.R:381`) and `forestsearch()` does not pass it (no `effect_log_scale` in `R/forestsearch_main.R`), and on the GLM path `hr = res$estimate` (`R/subgroup_search.R:887`, "effect estimate (RD, log-OR, etc.)") from `estimator_fn`, which for MD negates Y when `adverse_outcome = FALSE` (`R/glm_effect_estimators.R:815–820`). So the band acts on the **harm-oriented MD (identity) scale**: candidates with oriented MD ≥ 0.80 × max oriented MD (at ε = 0.20), largest N first. (`effMaxSG` is normalised to `hrMaxSG`, `:545–550`.)
- (iii) MR re-selection side: `.fs_mr_reselection_from_focus()` maps `hrMaxSG → "effMaxSG"` (`R/fs_mr_inference_methods.R:92`) and `forestsearch()` passes `effect_neighborhood = effect_neighborhood` and `selection_rule` through (`R/forestsearch_main.R:3405–3409`). `.fs_mr_select()` (`R/fs_mr_inference.R:132–176`): `eff <- if (log_scale) exp(beta[passers]) else beta[passers]  # natural effect`; `ib <- .compute_inclusion_band(hr_vec = eff, n_vec = sz, selection_rule = selection_rule, effect_neighborhood = nbhd) == 1L`; `if (!any(ib)) ib <- rep(TRUE, length(passers))`; `effMaxSG = { b <- .inband(); b[which.max(sizes[b])] }`. `log_scale` is FALSE for MD (§3), so `eff` is the harm-oriented MD; the same shared band helper, the same ε, natural scale, largest size in band. **Aligned** — same focus, same band expression, same scale — with two stated differences: (1) MR's band is over `passers` (the admitted set of the draw, floors resolved from the identification) rather than the identifier's full post-consistency table, and its empty-band fallback keeps all passers (`:165`); (2) the identifier's sort carries `-Pcons_vec` as a later key (`:609`), which MR's `effMaxSG` does not use. Which arm minima/`Pcons` floors reach `passers` is set by `admission_resolved` (`R/forestsearch_main.R:3402`), not re-derived here.

## 7. S0.6 — The `mdf1` record

**7.1 Locations** (all tracked at HEAD; all present on this machine):
- Record: `quarto/simulations/actg175/continuous/REPORT_continuous_field_2026-09-07.md` (208 lines); stage records `REPORT_continuous_field_stage0_2026-09-07.md` (105), `_stage1_` (76), `_gate2_` (202).
- Summary documents: `summary_continuous_field_mdf1.qmd` / `.html`; per cell `fs_maxeffCons_mr_field_<cell>_mdf1_combine_1_2000.html` (4); figures `fig_mdf1_bias_coverage_display_H.png`, `_Hc.png`.
- Per-replicate bundles (tracked): `mr_md_harm/fs_maxeffCons_mr_field_{md40_knoise0_n500,md120_knoise0_n500,mdnull_knoise0_n500,md40_knoise0_n700}_mdf1_d5000/` each with `*_res_1_1000.rds`, `*_res_1001_2000.rds`, `*_combined_1_2000.rds` and `gate2_flips.txt`; directory sizes 2.5 / 2.5 / 2.3 / 2.5 MB; combined bundles 1,271,333 / 1,237,302 / 1,149,263 / 1,254,403 bytes; **2,000 rows × 138 columns per cell** (computed here).
- Logs: the session's `campaign/<cell>_<label>.log`, `wall.txt`, `status.txt`, `peak_mem_mb.txt` under the Mac scratchpad path hard-coded in `scripts_mdf1/run_cell.sh:5–6` (`SP=/private/tmp/claude-501/.../scratchpad`, `C=$SP/campaign`). **Not committed and not on this machine.** The per-replicate timing columns (`fit_mr_secs`, `fld_H_secs`, `fld_Hc_secs`) in the bundles are the surviving timing record.

**7.2 Claims** (quotation, line numbers in `REPORT_continuous_field_2026-09-07.md`):
- Harm field one-sided lower coverage 0.947–0.950 in every cell — **confirmed**: `:188` "**0.950, 0.948, 0.950, 0.947** (md40 n500, md120 n500, null n500, md40 n700; Wilson half-widths ≈ 0.010)"; table `:20, :28, :36, :44` `cov1_w` 0.950 / 0.948 / 0.950 / 0.947.
- p̂(Ĥ) ≈ 0.16 on a family of about 1,842 candidates with many duplicate-membership labels — **confirmed**: `:196` "p̂(Ĥ) mean **0.155 / 0.189 / 0.155 / 0.163** ... The family has 1,842 candidates (36 J-quantile and default cuts × two directions, ≤ 2 conjunctions) with many duplicate-membership labels (`str2` ≡ `preanti > 0`; `karnof ≤ 90` ≡ `≤ 95`)"; regime table `:156–159` `p_hat_mean` 0.155 / 0.189 / 0.155 / 0.163 (md120 is 0.189, so "≈ 0.16" holds for three cells).
- Unstudentized complement field upper coverage 0.924–0.936 — **confirmed**: `:188` "field one-sided UPPER coverage of β(Ĥᶜ): **0.924, 0.934, 0.926, 0.936**"; `:24, :32, :40, :48`.
- IJ SE/SD 1.5–1.75 — **confirmed on Ĥ, differs on Ĥᶜ**: `:190` "SE/SD **1.74 / 1.49 / 1.75 / 1.66** on Ĥ ... and **1.95 / 1.83 / 1.97 / 1.90** on Ĥᶜ"; the harm block spans 1.485–1.751 (`:19, :27, :35, :43`), the complement 1.833–1.965 (`:23, :31, :39, :47`).
- Bonferroni joint coverage 0.936–0.943 — **confirmed**: `:198` "Bonferroni joint coverage **0.940 / 0.943 / 0.940 / 0.936**"; `:136, :139, :142, :145`.
- Retained bias in SD units: unadjusted +4, IJ +0.9, field +0.4 — **confirmed for the tie cells, differs at md120**: `:192` "naive → IJ → field = +4.1 → +0.92 → +0.39 (md40 n500), **+2.7 → +0.36 → +0.06** (md120), +4.2 → +0.95 → +0.41 (null), +4.4 → +0.96 → +0.41 (n700)".
- Field two-sided coverage 0.942–0.978 — **confirmed**: `:190` "The field's two-sided interval covers at **0.978 / 0.942 / 0.978 / 0.970**"; `:20, :28, :36, :44` `cov2_w`.
- Every display point within 0.008 of the Gaussian reference — **confirmed for the field's one-sided points; differs as an "every point" statement**: `:194` "Field points: |observed − reference| ≤ **0.008** for the one-sided coverage on both blocks in every cell ... and ≤ 0.015 for the two-sided. IJ points sit within 0.02 (the largest gap 0.942 vs 0.962, one-sided, n700)"; display table `:169–184`.
- The pairing proof enumerated 55 label-tie rows and no selection flips — **confirmed as 55 enumerated rows (label ties and MR-numerics), 0 flips**: Gate 2 record `:38` "14 enumerated and excluded (0 selection flips)", `:81` "8 enumerated", `:129` "16 enumerated", `:178` "17 enumerated"; 14 + 8 + 16 + 17 = 55, and the four committed `gate2_flips.txt` files hold 14 + 8 + 16 + 17 = 55 lines. Their classification (`:37, :80` and the cell-3/4 counts in the commit messages `5a169691`, `c80a1291`): 43 label ties (pure or with a super-population target move: 13 + 6 + 13 + 11) and 12 MR-numerics rows (1 + 2 + 3 + 6) (differences ≤ 4.5e-5 to 2.4e-4 relative), **0 selection flips** (`:202` "0 selection flips across 4000 anchored replicates"). So "55 label-tie rows" **differs**: 55 enumerated rows, of which 43 are label ties.
- 13 workers, 172.8 min, peak 17.5 GB, about 15 s per replicate with field and complement at n = 500 — **confirmed**: `:204` "14.9–17.6 s per replicate at 13 workers (field 11–13 s, complement 0.4 s); 8,000 replicates in 172.8 min cumulative wall; peak memory 17.5 GB"; `:3` "13 workers"; Gate 2 `:9` "Peak memory on the first render: 17.5 GB", `:146` "cumulative **10,371 s = 172.8 min**". At n = 500: 14.9 s (md40), 17.2 s (md120), 14.7 s (null) — md120 at n = 500 is 17 s, not 15.

**7.3 Per cell, both blocks.** The record gives declaration counts (`n_det`, `:156–159`), the mean oriented β(Ĥ) (`mean_pH`, `:156–159`; prose `:92`), median/q10/q90 of the field and IJ bounds (`:97–108`, `:116–127`) and the field/IJ/naive threshold shares **without** MC SEs and **without** oracle rows in the bound-location tables. Everything below not quoted from the record was **computed here** with base R from the four `*_combined_1_2000.rds` (detected replicates; oriented scale; the oracle's one-sided bound as the record defines it, `:8` "the normal-based rows use est ∓ 1.645 SE", from `or_H_est − 1.645·or_H_se` and `or_Hc_est + 1.645·or_Hc_se`; MC SE = √(p(1−p)/n)). The oracle refits on the **true region** Q / Qᶜ (`:8`, §6.3), so "the oracle lower bound on Ĥ" below is the oracle's bound, scored beside Ĥ's field bound, not a bound on Ĥ.

| cell | declared (rate) | mean β(Ĥ) | mean β(Ĥᶜ) |
|---|---|---|---|
| md40 n500 | 1998/2000 (0.999; record `:156` n_det 1998, `:202` "0.997–1.000") | 31.669 (record `:156` 72.078 is mean \|Ĥ\|; `:92` "+31.7") | 30.882 (`:111` "about +31") |
| md120 n500 | 2000/2000 (1.000) | 96.002 (`:92` "+96.0") | 51.578 (`:111` "+52") |
| null n500 | 1993/2000 (0.997) | 26.255 (`:92` "+26.3") | 26.255 (`:111` "+26") |
| md40 n700 | 1999/2000 (1.000) | 31.672 | 30.909 |

Quantiles (5 / 25 / 50 / 75 / 95 %), harm block, one-sided 95% LOWER bound (computed here; the record's medians `:99, :102, :105, :108` agree: −11.4 / 41.7 / −16.7 / −10.8):

| cell | field `fld_H_lo1s` | oracle lower |
|---|---|---|
| md40 n500 | −43.5 / −26.3 / −11.4 / 4.5 / 30.7 | −27.0 / −5.9 / 7.3 / 20.1 / 39.8 |
| md120 n500 | 1.5 / 23.7 / 41.7 / 61.8 / 93.5 | 53.0 / 74.1 / 87.3 / 100.1 / 119.8 |
| null n500 | −48.4 / −31.0 / −16.7 / −0.6 / 25.9 | undefined (Q empty; `or_H_*` NA on every row) |
| md40 n700 | −42.5 / −24.6 / −10.8 / 4.3 / 32.6 | −15.5 / 0.8 / 12.0 / 24.0 / 41.5 |

Complement block, one-sided 95% UPPER bound (computed here; record medians `:118, :121, :124, :127` agree: 48.8 / 71.5 / 44.2 / 46.2):

| cell | field `fld_Hc_up1s` | oracle upper |
|---|---|---|
| md40 n500 | 28.4 / 40.9 / 48.8 / 56.8 / 68.7 | 27.1 / 39.8 / 50.0 / 59.5 / 72.5 |
| md120 n500 | 49.5 / 62.0 / 71.5 / 80.0 / 92.7 | 26.8 / 39.8 / 49.9 / 59.5 / 72.5 |
| null n500 | 23.6 / 36.2 / 44.2 / 51.9 / 63.9 | 26.6 / 37.5 / 45.1 / 52.6 / 64.0 |
| md40 n700 | 29.5 / 39.4 / 46.2 / 53.3 / 64.2 | 27.2 / 38.4 / 46.6 / 54.6 / 66.5 |

Share of lower bounds at or above τ, field | oracle, with MC SE (computed here; the field shares equal the record's `:99, :102, :105, :108` to three decimals):

| cell | τ = 0 | 10 | 20 | 30 | 40 |
|---|---|---|---|---|---|
| md40 n500 | 0.307 (0.010) \| 0.639 (0.011) | 0.185 (0.009) \| 0.442 (0.011) | 0.107 (0.007) \| 0.252 (0.010) | 0.054 (0.005) \| 0.119 (0.007) | 0.025 (0.003) \| 0.048 (0.005) |
| md120 n500 | 0.956 (0.005) \| 1.000 (0.000) | 0.896 (0.007) \| 1.000 (0.000) | 0.799 (0.009) \| 0.999 (0.001) | 0.661 (0.011) \| 0.997 (0.001) | 0.524 (0.011) \| 0.986 (0.003) |
| null n500 | 0.240 (0.010) \| — | 0.138 (0.008) \| — | 0.076 (0.006) \| — | 0.037 (0.004) \| — | 0.017 (0.003) \| — |
| md40 n700 | 0.309 (0.010) \| 0.763 (0.010) | 0.185 (0.009) \| 0.552 (0.011) | 0.097 (0.007) \| 0.322 (0.010) | 0.059 (0.005) \| 0.165 (0.008) | 0.035 (0.004) \| 0.059 (0.005) |

Complement shares at or below τ, field | oracle (computed here; field shares equal the record's `:118, :121, :124, :127`): md40 n500 τ 0/10/20/30: 0.000 | 0.001; 0.002 | 0.006; 0.015 | 0.020; 0.068 | 0.080. md120: 0.000 | 0.001; 0.000 | 0.006; 0.001 | 0.021; 0.002 | 0.081. null: 0.000 | 0.000; 0.006 | 0.001; 0.029 | 0.015; 0.120 | 0.090. md40 n700: 0.000 | 0.000; 0.001 | 0.002; 0.006 | 0.012; 0.054 | 0.086 (MC SEs ≤ 0.007).

Complement thresholds the record used: `:111` "share at or below each threshold ("harm at most τ" supported)" with columns `P(U<=0) P(U<=10) P(U<=20) P(U<=30)` (`:114`); harm block `:90` "0 / 10 / 20 / 30 / 40 CD4 cells/mm³, oriented" (`:95`); `:9` "0 = no harm; 10 = the design's consistency threshold; 20; 30 = the search's effect threshold; 40 = the planted MD".

**7.4 The null cell.** Truth as recorded: `:7` "the null cell has no treatment-by-subgroup interaction (a homogeneous −26.26)"; `:92` "+26.3 under the null (no subgroup; homogeneous +26 on the harm-oriented scale, so the null cell's bound locations are read against +26, not as a false-claim rate)"; bundle `truth` (computed here): `effect_Q NA`, `effect_Qc −26.255`, `effect_ITT −26.255`, `beta_inter 0`, `prevalence_Q 0`. Declaration rate 1993/2000 = 0.997 (`:158` n_det 1993; `:202` "Detection 0.997–1.000 in every cell including the null"). Coverage of β(Ĥ) is conditioned on declaration and scored against the exact super-population target at the realized rule (`:8` "β(Ĥ), β(Ĥᶜ) are the exact super-population targets at the realized rule, oriented with −1"), which under the null equals +26.255 for every Ĥ (computed here: mean β(Ĥ) = 26.255 with zero spread); the oracle row is blank (`:8` "under the null the harm oracle is undefined (Q is empty) and its row is blank", `:34`). Bound-location figures for the null cell, with definitions: harm block, share of replicates whose one-sided 95% LOWER bound on β(Ĥ) is ≥ τ (`:105`): field 0.240 / 0.138 / 0.076 / 0.037 / 0.017 at τ = 0/10/20/30/40; IJ 0.211 / 0.120 / 0.065 / 0.030 / 0.015 (`:104`); naive 1.000 / 0.999 / 0.960 / 0.867 / 0.709 (`:103`); mean/median/q10/q90 of the field lower bound −14.844 / −16.708 / −42.663 / 15.641 (`:105`); `:200` "the field's lower bound sits above 0 on 24% of trials — a bound location against +26, not a false-claim rate ... above 10 on 14%, above 20 on 8%; the naive lower bound sits above 30 on 87%". Complement, share of replicates whose one-sided 95% UPPER bound on β(Ĥᶜ) is ≤ τ (`:124`): field 0.000 / 0.006 / 0.029 / 0.120 at τ = 0/10/20/30; mean/median 43.955 / 44.153. One-sided coverage: field lower on Ĥ 0.950 (`:36`), field upper on Ĥᶜ 0.926 (`:40`); Bonferroni joint 0.940 (`:142`).

## 8. S0.7 — Estimated cost

**8.1** From the Gate 2 record and the bundles (13 workers on the Mac; the render logs themselves are not committed, §7.1):

| cell | batch 1 / batch 2 / combine wall | replicates | peak memory |
|---|---|---|---|
| md40 n500 | 1,189 s / 1,209 s / 24 s (Gate 2 `:9`) | 2,000 | 17.5 GB (first render, `:9`; sampled every 5 s, summed RSS of R/Quarto) |
| md120 n500 | 1,371 s / 1,344 s / 20 s (`:58`) | 2,000 | not sampled (`run_cell.sh:28` samples only when `SAMPLE_MEM = 1`) |
| null n500 | 1,200 s / 1,163 s / 23 s (`:98`) | 2,000 | not sampled |
| md40 n700 | 1,391 s / 1,415 s / 22 s (`:146`) | 2,000 | not sampled |
| total | 10,371 s = 172.8 min (`:146`, record `:204`) | 8,000 | |

Per-replicate distributions (computed here from the bundles, all 2,000 rows per cell; seconds). **Nesting:** `fit_mr_secs` is the whole `forestsearch()` call including MR and both field blocks (template `:588`, `:619`; `:1765` "fit_mr_secs INCLUDES the field blocks"); `fld_H_secs` (`f$timing_seconds`, stamped at `R/fs_mr_inference.R:939` before the complement call at `:948`) and `fld_Hc_secs` (`fc$timing_seconds`, its own `t0c` `:1173`) are disjoint pieces nested inside `fit_mr_secs`. They are not summed here.

| cell | `fit_mr_secs` mean / median / p90 / max | `fld_H_secs` mean / median / p90 / max | `fld_Hc_secs` mean / median / p90 / max |
|---|---|---|---|
| md40 n500 | 15.0 / 15.0 / 17.3 / 20.0 | 11.0 / 11.0 / 12.5 / 14.2 | 0.4 / 0.4 / 0.5 / 0.7 |
| md120 n500 | 17.2 / 17.2 / 18.9 / 20.5 | 12.5 / 12.5 / 13.7 / 14.7 | 0.4 / 0.4 / 0.4 / 0.7 |
| null n500 | 14.7 / 14.7 / 17.2 / 20.9 | 10.9 / 10.9 / 12.5 / 15.1 | 0.4 / 0.4 / 0.5 / 0.7 |
| md40 n700 | 17.6 / 17.7 / 20.1 / 23.1 | 12.8 / 12.9 / 14.4 / 16.6 | 0.4 / 0.4 / 0.5 / 0.8 |

(Record `:156–159` `fit_secs` 14.953 / 17.151 / 14.744 / 17.648 and Gate 2 `:9, :58, :98, :146` means agree.) Effective throughput: 2,000 replicates in ≈ 2,400 s of batch wall ≈ 1.2 s/replicate at 13 workers, i.e. ≈ 15–17.6 s per replicate ÷ 13 with little loss.

**8.2** Per-replicate work today's defaults add that `mdf1` did not: **the field-s scale step** — for the template as committed (`field_scale_complement` inherited as `"selected"`, §4.2) `.fs_mr_field_complement()` additionally computes `s <- sqrt(colSums(Bc * Bc))` (`R/fs_mr_inference.R:1219`), the rescaled reading `lam_cs[r]` inside the existing per-draw loop (`:1237–1239`), the quantiles of `lfs` (`:1262–1263`) and a second `.fs_mr_field_joint()` on the same aligned draws (`:1291–1292`, a 26-point γ grid, `:1348`). It draws nothing and fits nothing new (the ensure-fit at `:1210–1218` can only run when `bdc` is finite, in which case `sel` was already fit, `:1198–1207`). No change to draw counts (`draws` 5000 from the template; `field_R_out/field_R_in` 1000/500 unchanged, `:562–563`) and no change to re-selection bookkeeping (`return_reselection` was already TRUE in `mdf1`, meta computed here). `field_recovery` (`e819a58a`) and `field_decompose` (`7245e898`) stay FALSE and are inert. The added arithmetic is O(Ncol × n + R_out × R_in) against a complement block that took 0.4 s of a 15 s replicate; **unmeasured here** (no timing run).

**8.3 Estimate** (labelled as such; `mdf1`'s machine, 13 workers): one rule (`mdf1`'s), 4 cells × 2,000 replicates = 8,000 replicates ≈ **173 min** (the measured `mdf1` wall, 172.8 min, plus the unmeasured field-s increment, which §8.2 argues is a small fraction of the 0.4 s complement block); both rules (`mdf1`'s and `effMaxSG` at ε = 0.20), 16,000 replicates ≈ **346 min ≈ 5.8 h**, assuming the second rule costs the same per replicate — unmeasured; the family, the admitted set and hence the field's lazy complement fits (`:1182–1196`; `mdf1` recorded ≈ 520–580 complement fits per replicate, record `:156–159` `nfit_mean`) could differ under a different winner. The `run_cell.sh` ceiling was 240 min (`CEIL=14400`, `:6`), so a two-rule campaign on the Mac exceeds the `mdf1` driver's ceiling as written. No estimate is made for this machine (no `mdf1`-configuration timing exists here; the template caps workers at physical cores − 1 = 63 by default, `:115–116`, and the survival memory notes 100-worker contention factors of 2.0–3.1 on gbsg_020, which this design has not measured).

**8.4** Render environment and timeout wrapper the `mdf1` scripts used, `scripts_mdf1/run_cell.sh:18–20`:
```
  FS_MD_MD=$MD FS_MD_N=$N FS_MD_MODE=$mode FS_MD_START=$start FS_MD_NSIMS=$nsims FS_MD_FB=$fb \
  FS_MD_CAMPAIGN=$TAG FS_MD_WORKERS=13 FS_MD_CI=field \
    bash /private/tmp/claude-501/-Users-larryleon-Downloads/476ca7fd-7bca-4851-bd29-048be32a2d48/scratchpad/tmo.sh $TMO quarto render sim_fs_maxeffCons_mr_field_md_template.qmd --output $out > $C/${CELL}_${label}.log 2>&1
```
with `CEIL=14400; TMO=5400` (`:6`), the batch sequence `:29–31` (`render batch1 batch 1 1000 $FB ...`, `render batch2 batch 1001 1000 none ...`, `render combine combine 1 1000 none ...`), the 1.5×-projection check `:25`, and the wrapper `scripts_mdf1/tmo.sh:3`:
```
perl -e 'my $t = shift @ARGV; my $pid = fork(); if ($pid == 0) { exec @ARGV or exit 127 } local $SIG{ALRM} = sub { kill "TERM", $pid; sleep 5; kill "KILL", $pid; exit 124 }; alarm $t; waitpid($pid, 0); exit($? >> 8);' "$@"
```
(README `:7`: "used because macOS has no `timeout`"; this machine has GNU `timeout`.) Memory sampler `scripts_mdf1/mem_sampler.sh` (5-second summed RSS of R/Rscript/quarto/deno). No `VECLIB_MAXIMUM_THREADS` or BLAS-thread setting appears in `run_cell.sh` (that variable appears only in the applied record for the OC document's `mclapply`, §9).

## 9. S0.8 — The applied document

`quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` (last change `f87518e9`, 2026-09-08, Mac; record `REPORT_actg175_continuous_intervals_2026-09-07.md`; payload `_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_payload.rds`, `forestsearch_version 0.3.5`, `built_at 2026-09-08 09:47:42`).

**9.1** Gate call's MR arguments, `:128–131`:
```
  mr_inference     = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L,
                           include_complement = TRUE, field_complement = TRUE,
                           return_reselection = TRUE)
```
(`field_scale_complement` not passed → inherited `"selected"` today, `"none"` when rendered.) Payload element: `iv` (`:207–239`), stored as `extras$intervals` (`:1285`), with fields `scale, thresholds, settings, H, Hc, joint, reselection, timing_seconds`; `Hc$field = list(est2, upper_1s, lower_1s, lower_2s, upper_2s, lambda_sd = fc$se_field, lambda_mean, q05, q95, n_out_used, n_complement_fits, timing_seconds)` (`:227–231`); `joint = list(gamma, lower_H, upper_Hc, joint_prob, bonf_lower_H, bonf_upper_Hc, bonf_joint_prob, corr, n_joint_draws)` (`:232–234`); `settings` records `ci_method = field, draws = 5000, multiplier = poisson, reselection = maxeff, ij_residual = two_term, seed = 8316951, field_R_out = 1000, field_R_in = 500, n_family = 4935, n_selected = 66` (computed here from the committed payload). **No field-s field is present** in the document's `iv` or in the committed payload (no `_s` leaf, no `joint_s`; computed here).

**9.2** Complement rows, `:277–290` (`tab2`): point `gc$naive$est`, `gc$debiased$est`, `fc$est2`; SE `gc$debiased$se_wald`, `gc$debiased$se_ij`, `fc$se_field`; two-sided `gc$naive$lower/upper`, `gc$debiased$lower/upper`, `fc$lower_2s/upper_2s`; one-sided upper `up1 <- c(naive = gc$naive$est + z95 * gc$debiased$se_wald, ij = gc$debiased$est + z95 * gc$debiased$se_ij, field = fc$upper_1s)` (`:203–204`). Joint rows, `:302–311` (`tab3`): `c(f$lower_1s, jt$bonf_lower_H, jt$lower_H)`, `c(fc$upper_1s, jt$bonf_upper_Hc, jt$upper_Hc)`, `jt$bonf_joint_prob`, `jt$joint_prob`, subtitle `jt$gamma, jt$corr, jt$n_joint_draws` — all read from the **unstudentized** `f$complement` and `f$joint` (`:182` `f <- g$field; fc <- f$complement; jt <- f$joint`). Re-rendering under today's defaults: the field-s step is add-beside and the unstudentized fields it displays are computed unchanged (§2.4 `3880022e`, §8.2), `ci_method` is passed, and the package version is the same 0.3.5, so **no displayed number changes by construction**; a re-render on this machine would differ from the Mac-built payload at floating-point level only (the record's Linux↔Mac finding, 9.3). Not verified by rendering.

**9.3** The intervals section **cannot be re-rendered without re-running the OC evaluation**: there is no `eval` flag, no `cache`/`freeze` option and no `Sys.getenv` switch in the document (grep), the header states `:30–33` "Everything below is computed in this document ... Nothing is read from disk", and a caching attempt was "attempted and reverted per gate" (`a0f732a0`). The intervals chunks (`:175–316`) only need `fs_anchor` (the gate call at `:92–132`, seconds), but the document runs the OC loop (`:566–568`) in the same render. Last recorded render (record `:39`): attempt 3 on the Mac at `n_workers = 1`, "wall 20 min 49 s; the document's own clock 'document compute wall-clock so far: 20.2 min'; 'evaluation loop: 1123.4 s (18.7 min) over 20 jobs, 1 workers'; peak summed R/quarto RSS 15,641 MB"; memory per worker `:116` "≈ 11 GB per worker" (`:14–17` measured 10.7 GB per `fs_oc_grid()` job; the document's own footer `:1229–1231` "the OC loop needs about 11 GB per worker"); the committed `params` are `draws: 20000`, `n_workers: 14` (`:5–6`), which the record calls "Linux's" (`:116`). Cross-platform tolerance for payload identity, `:115`: "the committed OC payload was built on Linux; a Mac render reproduces it to 1e-13 (last bit in T_obs), not bit-for-bit. Any future 'identity against the committed payload' gate on a document rendered on the other platform needs a floating-point tolerance stated in the task." The intervals-payload fields `settings$reselection = maxeff` show the applied document's gate re-selects under `maxeff`, i.e. the applied analysis's focus is not the simulation's `maxeffCons` (a fact, recorded for D5).

## Facts for Gate 0

**D1 — the re-run's selection rule** (`mdf1`'s `maxeffCons`, `effMaxSG` at ε = 0.20, or both)
- `mdf1` ran `sg_focus = "maxeffCons"`, `effect_neighborhood = 0.10` (inert for that focus), `selection_rule = "neighborhood"`, hard-coded (§6.1, §6.5(i)); bundle meta records `sg_focus maxeffCons`, `selection_rule neighborhood` (§6.4).
- Neither `sg_focus` nor `effect_neighborhood` is an `FS_MD_*` knob; the two-line template edit that would expose them is listed in §6.5(i); the stem re-tags via `fs_focus_tag()`; `meta` lacks `effect_neighborhood`.
- On the continuous path, `effMaxSG`'s band is `(1 − ε) × max` on the harm-oriented MD (identity) scale, largest N in band, on both the identifier and the MR re-selection, through the same `.compute_inclusion_band()`; aligned, with MR banding over the admitted `passers` with an all-passers fallback on an empty band (§6.5(ii)–(iii)).
- No `R/` commit since `mdf1` changes any recorded `mdf1` column at the template's passed arguments; the only computed addition is the field-s step, add-beside (§2.4, §8.2). A `mdf1`-rule re-run therefore has an identity anchor on every recorded column (to a floating-point tolerance across platforms and R 4.5.2→4.6.1); an `effMaxSG` re-run has none (different winner, different Ĥ).
- Cost on `mdf1`'s machine at 13 workers: ≈ 173 min for one rule, ≈ 346 min for both, the second rule's per-replicate cost unmeasured; the `mdf1` driver's ceiling was 240 min (§8.3).
- `mdf1` finding `:196`: every continuous cell is a tie regime under `maxeffCons` (p̂ 0.155–0.189; 97–99% of replicates below 0.5), with label ties of identical membership.

**D2 — whether the null cell is presented**
- Truth: no subgroup; `effect_Q` NA, `prevalence_Q` 0, `beta_inter` 0; homogeneous harm-oriented +26.255 everywhere (§7.4).
- Declaration 1993/2000 = 0.997 (the consistency screen is non-discriminating at thresholds 30/10, record `:202`, template `:218–219` "OPEN QUESTION D2").
- Coverage of β(Ĥ) is conditional on declaration against the exact target +26.255 for every Ĥ; the oracle row is blank; field one-sided lower coverage 0.950, complement upper 0.926, Bonferroni joint 0.940.
- Bound locations are read against +26, not as a false-claim rate: field lower bound ≥ 0 on 24.0% (MC SE 0.010), ≥ 10 on 13.8%, ≥ 20 on 7.6%, ≥ 30 on 3.7%, ≥ 40 on 1.7%; naive ≥ 30 on 86.7%; complement upper ≤ 30 on 12.0% (§7.3–7.4).
- The null cell's numbers coincide with the md40 n500 cell's on the harm block to within MC error (retained bias +0.41 vs +0.39 SD, coverage 0.950 both), record `:192`, `:196`.

**D3 — bound-location thresholds on the harm-oriented MD scale**
- `mdf1` used 0 / 10 / 20 / 30 / 40 for the harm block's lower bound and 0 / 10 / 20 / 30 for the complement's upper bound, defined as no harm / consistency threshold / 20 / effect threshold / planted MD (§7.3).
- md40 cells (β(Ĥ) ≈ +31.7): field lower-bound quantiles (5/25/50/75/95%) −43.5 / −26.3 / −11.4 / 4.5 / 30.7 (n500) and −42.5 / −24.6 / −10.8 / 4.3 / 32.6 (n700); shares ≥ 0/10/20/30/40: 0.307/0.185/0.107/0.054/0.025 and 0.309/0.185/0.097/0.059/0.035; the oracle's lower bound reaches ≥ 40 on 0.048 / 0.059.
- md120 (β(Ĥ) ≈ +96.0, q10–q90 64–120, record `:92`): field lower-bound quantiles 1.5 / 23.7 / 41.7 / 61.8 / 93.5; shares ≥ 0/10/20/30/40: 0.956/0.896/0.799/0.661/0.524 (MC SE ≤ 0.011); the top threshold used, 40, sits below the field bound's median (41.7) and far below the oracle's (median 87.3, 5% 53.0). Thresholds at 60, 80, 100 would fall within the field bound's inter-quantile range (25–95%: 23.7–93.5) and the oracle's (53.0–119.8); none is in the record.
- Complement (β(Ĥᶜ) ≈ +31 md40, +26 null, +52 md120): field upper-bound 5% quantiles 28.4 / 49.5 / 23.6 / 29.5, so shares ≤ 30 are 0.068 / 0.002 / 0.120 / 0.054 and ≤ 20 are ≤ 0.029; the oracle upper bound's shares are of the same size (≤ 0.09 at 30). A threshold at 40 or 50 would sit at the complement bound's 25% quantile (40.9 / 62.0 / 36.2 / 39.4); none is in the record.
- The field-s complement bound is not in any `mdf1` column; its location is unknown until a re-run records it (§5.2).

**D4 — the machine**
- This machine (`pop-os`): 64 physical / 128 logical cores, 251 GB, load 0.3, no R/quarto process; a stale 10-day bash wrapper with a `sleep` child (§1). forestsearch 0.3.5 built 2026-09-13 from source identical to HEAD (666/666 closures, §2.6); R 4.6.1 with reference BLAS/LAPACK 3.12.0.
- `mdf1`'s machine: Mac-Studio-3, M4 Max 14 cores / 36 GB, R 4.5.2, 13 workers, peak 17.5 GB, 172.8 min for 8,000 replicates; BLAS not recorded for `mdf1` (Accelerate per the applied record one day later) (§2.3, §7.1, §8.1).
- `origin/feature/glm-extension-mac` (`540c16e8`) is an ancestor of HEAD and of both `feature/glm-extension` refs; no local Mac branch exists (§2.2). The template and every `mdf1` artefact are in HEAD; the bundles are on this machine; the render logs are not (Mac scratchpad, §7.1).
- Cross-platform reproduction of `mdf1` rows is to floating point (Gate 2 pairing ≤ 1e-8 relative Mac vs the Linux-built twin bundles; the applied payload ≤ 1e-13), not bit-for-bit (§2.4, §9.3).
- The template caps `FS_MD_WORKERS` at physical cores − 1 (63 here) (§4.3); the survival memory of 100-worker contention on this machine is not a measurement of this design (§8.3). `mdf1`'s driver hard-codes a Mac scratchpad path and a Perl timeout because macOS lacks `timeout` (§8.4); GNU `timeout` exists here.

**D5 — the optional applied-document update**
- The applied gate call passes `ci_method = "field"`, `draws = 5000L`, `include_complement`, `field_complement`, `return_reselection` and not `field_scale_complement`; the committed payload's `extras$intervals` and the rendered tables carry only the unstudentized complement field and joint pair; no field-s field exists in the document or payload (§9.1–9.2).
- Re-rendering under today's defaults changes no displayed number by construction (add-beside field-s; same 0.3.5) but differs at floating point across platforms (≤ 1e-13, §9.3); displaying the field-s bound would need document edits to `iv` (`:227–239`) and `tab2`/`tab3` (`:277–316`), not listed here.
- The intervals section cannot be re-rendered without the OC loop: ≈ 21 min at `n_workers = 1` on the Mac (≈ 11 GB per worker; the committed `n_workers: 14` is the Linux setting) (§9.3).
- The applied gate re-selects under `maxeff` (`settings$reselection`), not the simulation's `maxeffCons`; its family has 4,935 candidates against the simulation's 1,842 (§9.1, §7.2).

## Findings (claims that differ, and gaps), in one place

1. "IJ SE/SD 1.5–1.75" holds on Ĥ only; on Ĥᶜ it is 1.83–1.97 (§7.2).
2. "Retained bias field +0.4 SD" holds in the three tie cells; md120 is +0.06 (§7.2).
3. "Every display point within 0.008" holds for the field's one-sided points; two-sided ≤ 0.015, IJ ≤ 0.02 (§7.2).
4. "55 label-tie rows": 55 enumerated rows, 43 label ties + 12 MR-numerics rows, 0 flips (§7.2).
5. "About 15 s per replicate at n = 500": 14.7–15.0 s in md40/null, 17.2 s in md120 (§7.2, §8.1).
6. `mdf1`'s record names `2f118042` as its commit; the template's commit is `c2c402e0` (same `R/`) (§2.3). BLAS for `mdf1` not recorded (§2.3). Render logs not committed (§7.1).
7. The template as committed inherits `field_scale_complement = "selected"` today: the field-s block runs but is not recorded, and `meta` and the poolability keys do not carry the setting (§4.2, §5.1). The recorder edit that would save it is listed in §5.2.
8. `fs_sim_bias_coverage()` cannot read field-s columns under the recorder's naming without an `R/` change (a document-level rename would) (§5.4).
9. No committed long-format metrics writer (cell × block × estimator × metric with MC SEs) was located (§5.5).
10. The `mdf1` record's bound-location tables have no MC SEs and no oracle rows; both are computed here (§7.3).
11. `effect_neighborhood` is not in the bundle `meta` (§6.5(i)).
12. `quarto/simulations/actg175/continuous/` has no `current_status.md` and no generator; nothing regenerated (§10.2 of the task).
