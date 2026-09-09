# REPORT — Stage 0, Guo & He supplement (Tier 1 + Tier 2(a)), Mac Studio, 2026-09-09

Task: `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` (committed at T0, `179ae409`).
Machine: Mac Studio. No `git fetch` and no `git pull` were run at any point in this task.

**OUTCOME: STOP at Stage 0 §2 (Records). One required input is absent from the tree:
`REVIEW_certification_2026-09-09.md`. Per §1 and the §7 STOP list, execution halts here;
no fetch was attempted. Everything else in Stage 0 verified clean and is recorded below,
so T1 can start on the same reading the moment the input arrives.**

> **2026-09-09 (post-merge `f221f75e`):** line numbers for `R/fs_mr_inference.R` refreshed (+8 throughout, from the expanded `ci_method` roxygen block); quoted content unchanged, re-verified by content; no result affected.

---

## 1. Provenance (verbatim)

`git log -1 --oneline` at the start of the task (before T0):

```
a9c0d5e5 Merge branch 'feature/glm-extension' of https://github.com/larry-leon/forestsearch into feature/glm-extension
```

`git log -1 --oneline` after T0 (the state this reading reflects):

```
179ae409 T0: commit the Guo & He supplement task document (TASK_guohe_supplement_2026-09-09)
```

`git status -sb` at the start of the task (before T0):

```
## feature/glm-extension...origin/feature/glm-extension
```

`git status -sb` after T0:

```
## feature/glm-extension...origin/feature/glm-extension [ahead 1]
```

### The seven untracked files — discrepancy, recorded not acted on

- The task's §1 states that seven pre-existing untracked files are out of scope: three `diag_*`
  bundles, three `diag_*.html`, and the ACTG175 payload.
- **This tree has zero untracked files.** `git status --porcelain --untracked-files=all` returns
  no `??` entries at all, both before and after T0. The working tree was clean at `a9c0d5e5`.
- Therefore the seven cannot be enumerated by name from this tree, and there was no opportunity
  to touch them. The out-of-scope guarantee holds vacuously here.
- Reading: those seven files live on the Linux box (the in-flight zero-compute session), not on
  this Mac. This is consistent with the missing certification record below — both point to the
  same un-synced Linux state.
- No `git add -A`, no `git add .`, no directory-level add was used. The single commit made
  (T0) named exactly one path.

---

## 2. Input inventory

### 2a. Bundles — ALL PRESENT (byte sizes as read)

| path | bytes |
|---|---|
| `quarto/GuoHe/guohe_repro_t7_beta2_00.rds` | 208906 |
| `quarto/GuoHe/guohe_repro_t7_beta2_01.rds` | 244231 |
| `quarto/GuoHe/guohe_repro_t7_beta2_02.rds` | 243147 |
| `quarto/GuoHe/guohe_repro_t7_beta2_03.rds` | 241623 |
| `quarto/GuoHe/guohe_repro_t7_beta2_04.rds` | 240678 |
| `quarto/GuoHe/guohe_repro_t7_beta2_05.rds` | 239156 |
| `quarto/GuoHe/mr_vs_guohe_t7_beta2_00.rds` | 491781 |
| `quarto/GuoHe/mr_vs_guohe_t7_beta2_01.rds` | 518900 |
| `quarto/GuoHe/mr_vs_guohe_t7_beta2_02.rds` | 517083 |
| `quarto/GuoHe/mr_vs_guohe_t7_beta2_03.rds` | 514834 |
| `quarto/GuoHe/mr_vs_guohe_t7_beta2_04.rds` | 512523 |
| `quarto/GuoHe/mr_vs_guohe_t7_beta2_05.rds` | 508011 |
| `quarto/GuoHe/mr_field_vs_guohe_t7_beta2_00.rds` | 717581 |
| `quarto/GuoHe/mr_field_vs_guohe_t7_beta2_01.rds` | 745629 |
| `quarto/GuoHe/mr_field_vs_guohe_t7_beta2_02.rds` | 744102 |
| `quarto/GuoHe/mr_field_vs_guohe_t7_beta2_03.rds` | 742337 |
| `quarto/GuoHe/mr_field_vs_guohe_t7_beta2_04.rds` | 740240 |
| `quarto/GuoHe/mr_field_vs_guohe_t7_beta2_05.rds` | 735427 |
| `quarto/GuoHe/guohe_sec52_truth_beta2_00.rds` | 3709 |
| `quarto/GuoHe/guohe_sec52_truth_beta2_01.rds` | 4055 |
| `quarto/GuoHe/guohe_sec52_truth_beta2_02.rds` | 4047 |
| `quarto/GuoHe/guohe_sec52_truth_beta2_03.rds` | 4041 |
| `quarto/GuoHe/guohe_sec52_truth_beta2_04.rds` | 4033 |
| `quarto/GuoHe/guohe_sec52_truth_beta2_05.rds` | 4036 |

### 2b. Scripts and package sources — ALL PRESENT

| path | bytes |
|---|---|
| `quarto/GuoHe/guohe_sec52_sim.R` | 19281 |
| `quarto/GuoHe/guohe_sec52_truth.R` | 18788 |
| `quarto/GuoHe/guohe_sec52_run.R` | 10042 |
| `quarto/GuoHe/guohe_reproduction_run.R` | 5648 |
| `quarto/GuoHe/guohe_reproduction_sim.R` | 16566 |
| `quarto/GuoHe/mr_vs_guohe_sim.R` | 24142 |
| `quarto/GuoHe/mr_field_vs_guohe_run.R` | 11268 |
| `quarto/GuoHe/mr_field_vs_guohe.qmd` | 14870 |
| `R/guohe_algorithm3.R` | 23642 |
| `R/guohe_adaptive_r.R` | 13961 |
| `R/fs_to_guohe.R` | 16961 |

### 2c. Records — ONE PRESENT, ONE ABSENT (**the STOP**)

| path | status |
|---|---|
| `dev/notes/NOTE_complement_product_2026-09-08.md` | present, 1402 bytes |
| `REVIEW_certification_2026-09-09.md` | **ABSENT — not in the working tree, not tracked, not in any local branch, not in any reachable commit** |

Searches run to establish the absence (all read-only, no fetch):

- `find . -name REVIEW_certification_2026-09-09.md -not -path "./.git/*"` — no hit.
- `git ls-files | grep -i -E "certification|REVIEW_"` — hits only
  `dev/notes/REVIEW_E1_fields_2026-09-08.md`, `dev/notes/REVIEW_partB_banddial_2026-09-08.md`,
  `dev/notes/REVIEW_partsAC_close_2026-09-08.md`, plus unrelated `dev/AI_review_updates/*`,
  `quarto/grf/grf_review_and_recommendations.qmd`, `quarto/guides/review_classification_calculations.qmd`.
- `git log --all --oneline -- "*REVIEW_certification*"` — empty.
- `git ls-tree -r` scan over reachable commits for any path matching `certification` — no hit.
- No `dev/notes` record of any kind is dated 2026-09-09; the newest are 2026-09-08.

Nearest-in-content documents present in this tree, for orientation only (they are **not**
substitutes, and no quote was taken from them for the task's purposes):

- `dev/notes/REVIEW_E1_fields_2026-09-08.md:78` — field-s Ĥᶜ upper `0.912 / 0.921 / 0.929`, and
  `:86` — `joint_s` `0.942 / 0.939 / 0.949 / 0.951–0.952` vs `joint` `0.932 / 0.933 / 0.935 / 0.935`.
- `dev/notes/HANDOFF_mr_field_linux_2026-09-08.md:21` — `joint_s` with field-s `0.939–0.952`.

The task's §2.4 asks for field-s one-sided upper **0.912–0.960** and joint **0.939–0.963**. The
upper end of both ranges (0.960 / 0.963) is not attained by anything in this tree; the certified
ranges evidently come from the cert20 campaign material that has not reached this Mac. This
corroborates that the document is on the Linux box, unpushed — the sync Larry sequences.

---

## 3. Source quotes (file path + line numbers as read from THIS tree)

### Q1 — disabled-complement lines in the campaign adapter — VERIFIED

`quarto/GuoHe/mr_vs_guohe_sim.R`, `mv_mr()` signature, lines 124-128:

```r
mv_mr <- function(df, cands, sel_label, spec, draws = MV_DRAWS,
                  multiplier = MV_MULTIPLIER, seed = NULL,
                  ci_method = "ij", field_R_out = 1000L, field_R_in = 500L,
                  field_uniform = FALSE, field_complement = FALSE,
                  include_complement = FALSE, ij_residual = "two_term") {
```

Where they are forwarded, lines 139 and 146-150:

```r
    include_complement = include_complement,
```

```r
  if (identical(ci_method, "field"))
    args <- c(args, list(field_R_out = field_R_out, field_R_in = field_R_in,
                         # kappa(Sigma-hat) sweep (TASK_mr_field_uniform_2026-09-05);
                         # FALSE keeps the 2026-09-05 field output byte-identical.
                         field_uniform = field_uniform,
                         # complement field (TASK_mr_field_complement_2026-09-06);
                         # FALSE keeps the field output byte-identical.
                         field_complement = field_complement))
```

Note for T1: the adapter's `FALSE / FALSE` are **adapter-level defaults that override the engine's
own `field_complement = TRUE`**. T1's named-line change is therefore at the `mv_mr()` call site in
the new driver (passing `field_complement = TRUE, include_complement = TRUE`), plus threading
`field_scale_complement` — which `mv_mr()` does **not** currently accept as a formal at all. That
is a change to `quarto/GuoHe/mr_vs_guohe_sim.R` (a shared campaign adapter), not to `R/`. Flagged
for Larry: T1 as specified needs either a new formal on `mv_mr()` or a bypass call to
`forestsearch:::fs_mr_inference` in the new driver. No change was made.

### Q2 — engine signature for `fs_mr_inference` — VERIFIED, NO CONTRADICTION

Defining file: `R/fs_mr_inference.R`, definition begins at line 522. Signature lines 522-542:

```r
fs_mr_inference <- function(df, candidates, spec, selected_members,
                           admission,
                           t_confirm = NULL, confirm_rule = c("point", "ci"),
                           reselection = c("maxcons", "maxeff", "maxSG",
                                           "minSG", "effMaxSG", "effMinSG"),
                           effect_neighborhood = 0.10,
                           selection_rule = c("neighborhood", "pareto", "both"),
                           draws = 2000L,
                           multiplier = c("poisson", "gaussian", "rademacher"),
                           include_complement = FALSE,
                           ci_method = c("ij", "wald", "field"),
                           seed = NULL,
                           return_reselection = TRUE,
                           field_R_out = 1000L,
                           field_R_in = 500L,
                           field_uniform = FALSE,
                           field_M_cap = NULL,
                           field_complement = TRUE,
                           field_decompose = FALSE,
                           field_scale_complement = c("selected", "none"),
                           ij_residual = c("two_term", "winner", "winner_floor")) {
```

Against the task's check:

| argument | required by task | read at this tree | verdict |
|---|---|---|---|
| `field_complement` | `TRUE` | `TRUE` (line 539) | match |
| `field_scale_complement` | `"selected"` | `c("selected", "none")`, `match.arg` at line 545 → `"selected"` | match |
| `return_reselection` | `TRUE` | `TRUE` (line 534) | match |

**No contradiction with the recorded adoption. Q2 does not trigger its STOP.**

Note (not a contradiction, but relevant to the adoption record): `NOTE_complement_product_2026-09-08.md`
states "The package default stays `"none"` so committed bundles remain byte-reproducible; flipping
the default is a separate decision." The installed/source default read here is `"selected"`, i.e.
the flip has since been made in `R/`. The task's expected value (`"selected"`) is what governs, and
it is what the tree has. Recorded because it bears on byte-reproducibility of older bundles.

**Interface for one-sided bounds at 0.95 and 0.975 — available with no `R/` change required.**

- Harm side, 0.95 one-sided lower: `field$lower_1s` (`R/fs_mr_inference.R:888`,
  `lower_1s = to_eff(beta_deb - qs[5])`).
- Complement side, 0.95 one-sided upper: `complement$upper_1s` (line 1160,
  `upper_1s = to_eff(bdc - qs[1])`); studentized companion `upper_1s_s` (line 1172).
- **Both 0.975 bounds and the joint indicator already exist as a matched pair** in the joint
  element, `R/fs_mr_inference.R:1219-1226`, `gamma = alpha/2 = 0.025`:

```r
  qh_b  <- stats::quantile(lh, 1 - alpha / 2, names = FALSE, type = 7)
  qc_b  <- stats::quantile(lc, alpha / 2, names = FALSE, type = 7)
  list(gamma = gamma, joint_prob = jp, alpha = alpha,
       lower_H = to_eff(beta_deb - qh_g), upper_Hc = to_eff(bdc - qc_g),
       bonf_gamma = alpha / 2,
       bonf_lower_H = to_eff(beta_deb - qh_b), bonf_upper_Hc = to_eff(bdc - qc_b),
       bonf_joint_prob = mean(lh <= qh_b & lc >= qc_b),
       corr = stats::cor(lh, lc), n_joint_draws = n,
       grid_gamma = grid, grid_joint_prob = probs)
```

  computed at lines 1146-1152 for both the unscaled (`joint`) and studentized (`joint_s`) fields:

```r
  # Joint (H lower, Hc upper) pair (method B): the harm field's lam and this
```
```r
    .fs_mr_field_joint(lam_H[ok_c], lf, beta_deb, bdc, to_eff, alpha) else NULL
    .fs_mr_field_joint(lam_H[ok_c], lfs, beta_deb, bdc, to_eff, alpha) else NULL
```

  So `bonf_lower_H` is the 0.975 one-sided lower on Ĥ, `bonf_upper_Hc` the 0.975 one-sided upper
  on Ĥᶜ, and the per-replicate both-correct indicator is formed from that pair. **T1 needs no
  `R/` edit to emit the 0.975 bounds** — the §4 "if emitting the 0.975 bounds requires any change
  under `R/`, STOP" condition is not triggered.

### Q3 — `guohe_adaptive_r()` signature — VERIFIED

`R/guohe_adaptive_r.R`, lines 178-195:

```r
guohe_adaptive_r <- function(data,
                             outcome = c("survival", "binary", "continuous"),
                             treatment,
                             candidates,
                             time = NULL,
                             event = NULL,
                             y = NULL,
                             orient = -1,
                             r_grid = c(0.03, 0.10, 0.20, 0.30, 0.40, 0.45),
                             v = 5L,
                             B = 200L,
                             level = 0.05,
                             seed = NULL,
                             min_events = 5L,
                             refit = TRUE,
                             adjust_covariates = NULL,
                             fast = NULL,
                             parallel = FALSE) {
```

- `r_grid` default `c(0.03, 0.10, 0.20, 0.30, 0.40, 0.45)` — **matches the value the task cites at
  `a9c0d5e`.**
- `v` default `5L`.
- **One-B behavior confirmed:** the same `B` is passed to the inner CV fits (line 255,
  `B = B, r = r_grid[l], level = level`) and to the final refit (line 285,
  `y = y, orient = orient, B = B, r = r_hat`). There is no separate refit-B formal.
- Constraint at lines 199-201: `all(r_grid > 0), all(r_grid < 0.5), v >= 2L, B >= 2L,
  orient %in% c(-1, 1)`. The task's pinned `r_grid = c(1/3, 1/12, 1/21, 1/30)` satisfies it.

### Q4 — the `--adaptive` path in `guohe_reproduction_run.R` — VERIFIED, WITH A MATERIAL FINDING

Flag parsing, lines 54-55:

```r
b_adapt <- as.integer(opt("adaptive-B", as.character(b_boot)))
adaptive <- flag("adaptive")
```

The replicate call, lines 108-111:

```r
      try(gh_one_rep(
        beta = s$beta, n = s$n, r_grid = GH_R_GRID, B = b_boot,
        v = 5L, adaptive = adaptive, seed = base + m
      ), silent = TRUE)
```

The recorded metadata, lines 127-128:

```r
      B = b_boot, adaptive_B = if (adaptive) b_adapt else NA_integer_,
      r_grid = GH_R_GRID, adaptive = adaptive, v = 5L,
```

`r_grid`: `quarto/GuoHe/guohe_reproduction_sim.R:52`

```r
GH_R_GRID <- c(1 / 3, 1 / 12, 1 / 21, 1 / 30)
```

The adaptive call itself, `quarto/GuoHe/guohe_reproduction_sim.R:184-190`, inside `gh_one_rep()`:

```r
    ar <- try(suppressWarnings(guohe_adaptive_r(
      data = df, outcome = "survival", treatment = "treat",
      candidates = cand, time = "time", event = "event",
      orient = -1, r_grid = r_grid, v = v, B = B, level = level,
      min_events = min_events, refit = TRUE
    )), silent = TRUE)
```

So the validated Tables 3–6 Adaptive columns used: `r_grid = c(1/3, 1/12, 1/21, 1/30)`, `v = 5L`,
`orient = -1`, `refit = TRUE`, `seed` **not** passed to `guohe_adaptive_r()`, and `B = b_boot`.

**Finding (footnoted, not reconciled): `--adaptive-B` is inert in this driver.** `b_adapt` is
parsed (line 54), printed (line 67) and stored in the bundle metadata as `adaptive_B` (line 127),
but it is **never passed to `gh_one_rep()`** — the call at line 109 passes `B = b_boot`, and
`gh_one_rep()` forwards that single `B` to both the fixed-r Algorithm-3 loop and the
`guohe_adaptive_r()` call. Consequence: the documented production commands in
`quarto/GuoHe/guohe_reproduction_RUN.md:57` and `:60`

```
Rscript guohe_reproduction_run.R --tables=35 --cores=120 --B=2000 --adaptive --adaptive-B=200
Rscript guohe_reproduction_run.R --tables=6  --cores=120 --B=2000 --adaptive --adaptive-B=200
```

would have recorded `adaptive_B = 200` in metadata while actually running the adaptive path at
`B = 2000`. Bearing on T2: the task pins `B = 200` as "the validated reproduction setting", citing
`--adaptive-B=200`. On this reading the *recorded* setting is 200 but the *executed* one was 2000,
which is also the ~10x cost line `guohe_reproduction_RUN.md:98-99` distinguishes
(`Adaptive at B = 200 | 184 core-h | ~1.7 h` vs `Adaptive at B = 2000 | 875 core-h | ~8.1 h`).
**This materially changes the T2 Gate 1b projection and the meaning of the "secondary" bound.**
Recorded here for Larry; nothing reconciled, nothing edited.

Two further Q4 differences against T2's t7 settings, footnoted not reconciled:

- `orient`: reproduction uses `-1`; T2 specifies `+1` (correct for t7 — the stored t7 bundles
  record `orient = 1`).
- `seed`: reproduction passes no `seed` to `guohe_adaptive_r()`; T2 requires a derived seed
  recorded.

### Q5 — certification harness's joint-Bonferroni computation — **BLOCKED**

The task directs: "locate the certification driver(s) behind `REVIEW_certification_2026-09-09.md`".
That document does not exist in this tree (§2c), so the driver set it names cannot be identified,
and the lines T1 is to transplant cannot be quoted from the authority the task specifies.

The engine-level pair T1 would ultimately be transplanting is quoted under Q2 above
(`R/fs_mr_inference.R:1203-1228`, `.fs_mr_field_joint()`), and the `joint` / `joint_s` wiring at
`:1146-1152`. That is recorded as orientation only. **It is not a substitute for Q5** — the task
is explicit that the certification harness's lines are transplanted and not re-derived, and
identifying them by inference from the engine would be exactly the re-derivation the task forbids.
**No transplant was attempted.**

### Q6 — transplant anchors in `mr_field_vs_guohe_run.R` — VERIFIED

- `--cells` flag block, lines 45-53:

```r
args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(args == paste0("--", nm))
opt <- function(nm, default) {
  hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
}
n_cores <- as.integer(opt("cores", "120"))
force <- flag("force")
cells_opt <- opt("cells", "")
```

- `MF_CELLS`, lines 55-58:

```r
MF_CELLS <- c(sprintf("t35_beta2_%02d", 0:5),
              sprintf("t6_k%02d", c(2L, 6L, 10L, 12L)),
              sprintf("t7_beta2_%02d", 0:5))
if (nzchar(cells_opt)) MF_CELLS <- strsplit(cells_opt, ",")[[1]]
```

  T1's named-line change is this block to `MF_CELLS <- sprintf("t7_beta2_%02d", 0:5)`.

- Seed structure. `mr_field_vs_guohe_run.R:157`:

```r
    id = id, m = m, seed_data = base + m, seed_mr = base + m + MV_SEED_MR,
```

  with `base <- if (sec52) mv_gh52_base(id) else mv_gh51_base(id)` at line 178, and from
  `quarto/GuoHe/mr_vs_guohe_sim.R`:

```r
MV_DRAWS      <- 5000L        # D3 default            (line 60)
MV_MULTIPLIER <- "poisson"    # D3 default (centred Poisson)   (line 61)
MV_SEED_MR    <- 700000L      # MR seed offset from the replicate data seed  (line 62)
```
```r
mv_gh52_base <- function(id) 1000000L + as.integer(sum(utf8ToInt(id))) * 100003L   (line 79)
```

  Field draws run under the derived seed `+ 900000L` inside the gate
  (`mr_field_vs_guohe_run.R:223`, `field_seed_offset = 900000L`).

- The `mf_rep_52` gate-call block, lines 151-155:

```r
  t0 <- proc.time()[["elapsed"]]
  mr <- mv_mr(cand$df, mv_cand_idx_52(cand), cand$names[sel], mv_spec52,
              seed = base + m + MV_SEED_MR, ci_method = "field",
              field_R_out = MF_R_OUT, field_R_in = MF_R_IN)
  t1 <- proc.time()[["elapsed"]]
```

  with `MF_R_OUT <- 1000L` / `MF_R_IN <- 500L` at lines 42-43.

- Row-wise `seed_data` identity assertions. The naive/selection pairing inside `mf_rep_52`,
  lines 142-149:

```r
  naive_ok <- identical(nv$point, row_r$naive_point) &&
    identical(nv$lower, row_r$naive_lower) &&
    identical(nv$c_hat, row_r$c_hat_naive) &&
    identical(nv$gamma_s, row_r$gamma_s_naive) &&
    identical(nv$cover, row_r$naive_cover)
  sel <- which.max(replace(fits$est, !is.finite(fits$est), -Inf))
  sel_ok <- identical(cand$cuts[sel], row_r$c_hat_gh) &&
    identical(gh52_truth_at(tru, cand$cuts[sel]), row_r$gamma_s)
```

  the MR-current pairing proof in `.mf_mr_cols`, lines 80-84:

```r
  cur_ok <- identical(est, row_c$mr_est) &&
    identical(mr$debiased$se_ij, row_c$mr_se_ij) &&
    identical(mr$selection_bias, row_c$mr_bias_sel) &&
    identical(mr$fixed_bias, row_c$mr_bias_fix) &&
    identical(unname(mr$reselection$p_hat[mr$selected_index]), row_c$p_hat_H)
```

  and the cell-level E5 join assertion, lines 206-210:

```r
  # E5 join: addendum-A columns from the 2026-09-04 bundle by (id, m);
  # row order is m for both (asserted above), seed equality asserted here.
  old <- cmp_bun$results[res$m, ]
  stopifnot(identical(old$seed_data, res$seed_data))
  res <- cbind(res, old[, MF_JOIN_COLS])
```

  plus the bundle-level base assertion at lines 179-181:

```r
  stopifnot(rep_bun$seed_base == base, cmp_bun$seed_base == base,
            nrow(cmp_bun$results) == nrow(rep_bun$results),
            identical(cmp_bun$results$m, seq_len(nrow(cmp_bun$results))))
```

- Gap noted for T1: `mr_field_vs_guohe_run.R` has **no `--pilot` flag**. Gate 1a's
  `--pilot` at reps = 20 requires transplanting that scaffolding from
  `quarto/GuoHe/guohe_sec52_run.R` (lines 70, 81-83: `pilot <- flag("pilot")`,
  `if (pilot) { ... n_rep <- as.integer(opt("reps", "20")) }`, output suffix `_pilot` at line 144,
  projection block at lines 214-225). Named-line addition, no re-authorship.

### §2.4 — certification quotes for B6 — **BLOCKED**

Requires `REVIEW_certification_2026-09-09.md` (§2c). Not quoted. No substitute used.

---

## 4. Supporting reads (recorded because T1/T2 depend on them)

Stored t7 reproduction bundle metadata, `quarto/GuoHe/guohe_repro_t7_beta2_00.rds`:

```
id t7_beta2_00 | target Table 7 | beta2 0 | n 400 | n_rep_requested 2000 | n_rep_used 2000
B 2000 | r_grid 0.333333333333333 0.0833333333333333 0.0476190476190476 0.0333333333333333
orient 1 | pilot FALSE | seed_base 92902757 | truth_identity 2e+06 0.25 20260721 TRUE exact
elapsed_sec 11569.456
results dim 2000 x 28
columns: beta2, n, n_cand, cens_rate, c_hat_gh, c_hat_naive, sel_agree, n_sel, gamma_s,
  naive_point, naive_lower, naive_cover, naive_dist, naive_bias, gamma_s_naive,
  r1_cover, r1_dist, r1_bias, r2_cover, r2_dist, r2_bias, r3_cover, r3_dist, r3_bias,
  r4_cover, r4_dist, r4_bias, rep_elapsed_sec
```

- The t7 bundles were produced by `guohe_sec52_run.R` (`orient = +1`, `pilot`, `truth_identity`,
  `c_hat_*`, `gamma_s`), **not** by `guohe_reproduction_run.R`. They carry no `adaptive`,
  `adaptive_B` or `v` field — consistent with the Adaptive column never having been run on t7,
  which is what T2 exists to supply.
- `quarto/GuoHe/guohe_sec52_sim.R:85`: `GH52_R_GRID <- c(1 / 3, 1 / 12, 1 / 21, 1 / 30)` — the same
  grid T2 pins, so T2's `r_grid` is exactly the stored fixed-r column set.
- **T2 lookup is feasible with the stored columns.** The Algorithm-3 bound itself is not stored,
  but `r{i}_dist = gamma_s - lower`, so the B = 2000 bound at any r is recovered exactly as
  `lower = gamma_s - r{i}_dist` — the same reconstruction `mr_field_vs_guohe_run.R:66-74` already
  performs in `.mf_gh_cols()`:

```r
    out[[sprintf("gh_r%d_low", i)]] <- theta - row[[sprintf("r%d_dist", i)]]
```

  So the "stored B = 2000 Algorithm-3 bound at r̂" is a `which(r_grid == r_hat)` index into those
  four columns. No re-run of the fixed-r columns is needed.

---

## 5. Stage 0 verdict

| item | verdict |
|---|---|
| §2.1 provenance | recorded; seven-untracked-file discrepancy documented (this tree has zero untracked files) |
| §2.2 bundles (24) | all present, sizes recorded |
| §2.2 scripts + `R/` sources (11) | all present, sizes recorded |
| §2.2 records | `NOTE_complement_product_2026-09-08.md` present; **`REVIEW_certification_2026-09-09.md` ABSENT** |
| §2.3 Q1 | verified |
| §2.3 Q2 | verified — defaults `TRUE` / `"selected"` / `TRUE`, no contradiction; 0.975 bounds available with no `R/` change |
| §2.3 Q3 | verified |
| §2.3 Q4 | verified, with the inert-`--adaptive-B` finding and the orient/seed differences footnoted |
| §2.3 Q5 | **BLOCKED** on the absent record |
| §2.3 Q6 | verified; `--pilot` scaffolding gap noted |
| §2.4 B6 quotes | **BLOCKED** on the absent record |

**STOP declared at §2c. No fetch, no pull, no push. T1, T2 and T3 not started.**

What unblocks: the Linux → Mac sync Larry sequences, bringing
`REVIEW_certification_2026-09-09.md` (and, apparently, the seven untracked files) into this tree.
On arrival, Q5 and §2.4 complete against it and T1 proceeds from the anchors recorded above.

## 6. Closing tree state

`git status -sb` at the close of Stage 0 (before this record is staged):

```
## feature/glm-extension...origin/feature/glm-extension [ahead 1]
```

- Untracked files: **zero**, identical to the opening state. The seven named in §1 were not
  present at any point and were therefore neither modified nor staged.
- Commits made in this session: one (`179ae409`, T0), by explicit named path.
- Nothing under `R/` was read-modified. No engine edit. No file outside
  `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` and this record was written.
