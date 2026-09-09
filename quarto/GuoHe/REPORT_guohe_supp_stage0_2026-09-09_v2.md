# REPORT — Stage 0 v2 (post-merge), Guo & He supplement (Tier 1 + Tier 2(a)), Mac Studio, 2026-09-09

Governing documents: `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` (v1, `179ae409`)
as amended by `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09_v2.md` (v2, `669c9ef5`).
Predecessor record: `quarto/GuoHe/REPORT_guohe_supp_stage0_2026-09-09.md` (`33d0add3`) —
**not modified**, per A1.

Machine: Mac Studio. **No `git fetch`, no `git pull`, no `git push` was run at any point** in
this session, per A0/A8.

**OUTCOME: Stage 0 v2 PASSES. No STOP condition is met.** All 24 bundles, 6 truth caches and 11
scripts/`R/` sources present; every Q1–Q6 line number re-verified; `REVIEW_certification_2026-09-09.md`
still absent, which under A2 is no longer a STOP and routes to the `[certification citation pending sync]`
marker in T3/B6.

> **2026-09-09 (post-merge `f221f75e`):** line numbers for `R/fs_mr_inference.R` refreshed (+8 throughout, from the expanded `ci_method` roxygen block) in this record and in both appendices below; quoted content unchanged, re-verified by content; no result affected. The Q1–Q6 comparison table in §3 is left as written: it is a dated snapshot of the v1 → v2 reading, and its `fs_mr_inference` entries are the pre-merge offsets. Add 8 to each for the current tree.

---

## 1. Provenance

`git log -1 --oneline` (at the opening of this record, after the v2 amendment commit):

```
669c9ef5 Amendment v2 to the Guo & He supplement task (TASK_guohe_supplement_2026-09-09_v2)
```

`git status -sb`:

```
## feature/glm-extension...origin/feature/glm-extension [ahead 3]
```

Observed untracked set (`git status --porcelain --untracked-files=all`, `??` entries):

```
(none)
```

- **Zero untracked files**, exactly as A0 anticipates ("expected to remain absent on the Mac —
  untracked files do not travel through git"). The three `diag_*` bundles, three `diag_*.html`
  and the ACTG175 payload are not present, were not staged, and could not have been touched.
  The rule stays in force and binds on the Linux box in Phase B.
- The three commits this session is ahead by are its own: `179ae409` (v1 task doc),
  `33d0add3` (Stage 0 v1 record), `669c9ef5` (v2 amendment). Nothing else is uncommitted.

### 1a. The merge referenced by A1 — factual observation for Larry

A1 is titled "Stage 0 v2 (Mac, post-merge)" and instructs that all v1 line numbers are "void
after the merge." **On this tree, no merge landed after the v1 reading.** Evidence, read-only:

- `git reflog` shows the most recent network operation as `HEAD@{2}`:

```
a9c0d5e5 HEAD@{2}: pull --ff --recurse-submodules --progress origin: Merge made by the 'ort' strategy.
```

  and the two entries above it are this session's own commits (`HEAD@{1}` = `179ae409`,
  `HEAD@{0}` = `33d0add3`, with `669c9ef5` added since).
- `.git/FETCH_HEAD` is dated `Sep 9 11:48`; the v1 T0 commit was made at 12:39.
- `git log -1 origin/feature/glm-extension` → `a9c0d5e5 … 2026-09-09 10:12:49 -0700`, i.e. the
  remote-tracking ref still points at the same merge the v1 record read from.

**Reading:** the Linux merge Larry pulled in GitHub Desktop is `a9c0d5e5` itself, which was
already the tree's HEAD when Stage 0 v1 ran. The v1 reading was therefore already post-merge.
This is recorded as an observation, not a STOP — A8 lists no such condition, and the required
re-verification was performed in full regardless (§3). Its practical consequence is benign:
every line number re-verifies at its v1 value, and the certification record is still absent
because it was never in the merge.

---

## 2. Input inventory (re-asserted, A1)

### 2a. Bundles and truth caches — 24 of 24 PRESENT, byte sizes unchanged from v1

| file (in `quarto/GuoHe/`) | bytes |
|---|---|
| `guohe_repro_t7_beta2_00.rds` | 208906 |
| `guohe_repro_t7_beta2_01.rds` | 244231 |
| `guohe_repro_t7_beta2_02.rds` | 243147 |
| `guohe_repro_t7_beta2_03.rds` | 241623 |
| `guohe_repro_t7_beta2_04.rds` | 240678 |
| `guohe_repro_t7_beta2_05.rds` | 239156 |
| `mr_vs_guohe_t7_beta2_00.rds` | 491781 |
| `mr_vs_guohe_t7_beta2_01.rds` | 518900 |
| `mr_vs_guohe_t7_beta2_02.rds` | 517083 |
| `mr_vs_guohe_t7_beta2_03.rds` | 514834 |
| `mr_vs_guohe_t7_beta2_04.rds` | 512523 |
| `mr_vs_guohe_t7_beta2_05.rds` | 508011 |
| `mr_field_vs_guohe_t7_beta2_00.rds` | 717581 |
| `mr_field_vs_guohe_t7_beta2_01.rds` | 745629 |
| `mr_field_vs_guohe_t7_beta2_02.rds` | 744102 |
| `mr_field_vs_guohe_t7_beta2_03.rds` | 742337 |
| `mr_field_vs_guohe_t7_beta2_04.rds` | 740240 |
| `mr_field_vs_guohe_t7_beta2_05.rds` | 735427 |
| `guohe_sec52_truth_beta2_00.rds` | 3709 |
| `guohe_sec52_truth_beta2_01.rds` | 4055 |
| `guohe_sec52_truth_beta2_02.rds` | 4047 |
| `guohe_sec52_truth_beta2_03.rds` | 4041 |
| `guohe_sec52_truth_beta2_04.rds` | 4033 |
| `guohe_sec52_truth_beta2_05.rds` | 4036 |

### 2b. Scripts and `R/` sources — 11 of 11 PRESENT, byte sizes unchanged from v1

| path | bytes | last commit touching it |
|---|---|---|
| `quarto/GuoHe/guohe_sec52_sim.R` | 19281 | `a9e099df` |
| `quarto/GuoHe/guohe_sec52_truth.R` | 18788 | — |
| `quarto/GuoHe/guohe_sec52_run.R` | 10042 | `a9e099df` |
| `quarto/GuoHe/guohe_reproduction_run.R` | 5648 | `8a98e05d` |
| `quarto/GuoHe/guohe_reproduction_sim.R` | 16566 | `d516a928` |
| `quarto/GuoHe/mr_vs_guohe_sim.R` | 24142 | `8fbb3bdc` |
| `quarto/GuoHe/mr_field_vs_guohe_run.R` | 11268 | `50d6b641` |
| `quarto/GuoHe/mr_field_vs_guohe.qmd` | 14870 | — |
| `R/guohe_algorithm3.R` | 23642 | — |
| `R/guohe_adaptive_r.R` | 13961 | `d516a928` |
| `R/fs_to_guohe.R` | 16961 | — |

Every byte size is identical to the v1 record — consistent with §1a (no post-v1 merge).

### 2c. `REVIEW_certification_2026-09-09.md` — STILL ABSENT (not a STOP under A2)

- `find . -name "REVIEW_certification_2026-09-09.md" -not -path "./.git/*"` — no hit.
- `git ls-files | grep -c "REVIEW_certification"` — `0`.
- Disposition per A2: T1 and T2 are unblocked; T3's B6 renders the literal marker
  `[certification citation pending sync]` at each of the two figures (field-s upper 0.912–0.960;
  joint 0.939–0.963), and the T3 record notes it.

---

## 3. Q1–Q6 re-verification (A1)

**Result: all six re-verify. No line number moved. One v1 citation range was imprecise and is
corrected below; the quoted content itself was correct.**

| quote | v1 line(s) | v2 line(s) | verdict |
|---|---|---|---|
| Q1 `mv_mr()` signature | 124–128 | 124–128 | unchanged |
| Q1 `include_complement` forward | 139 | 139 | unchanged |
| Q1 field-args block | "146–150" | **143–150** | **corrected** (content identical; v1 mislabelled the start of the 8-line block) |
| Q2 `fs_mr_inference` definition | 514 | 514 | unchanged |
| Q2 `include_complement = FALSE` | 523 | 523 | unchanged |
| Q2 `return_reselection = TRUE` | 526 | 526 | unchanged |
| Q2 `field_complement = TRUE` | 531 | 531 | unchanged |
| Q2 `field_scale_complement` | 533 | 533 | unchanged |
| Q2 its `match.arg` | 537 | 537 | unchanged |
| Q2/A2 `.fs_mr_field_joint` def | 1195 | 1195 | unchanged |
| Q2/A2 Bonferroni pair + indicator | 1211–1219 (fn to 1220) | 1211–1219 (fn to 1220) | unchanged |
| Q2/A2 `joint` / `joint_s` wiring | 1138–1144 | 1141–1144 (comment from 1138) | unchanged |
| Q2 harm `lower_1s` | 880 | 880 | unchanged |
| Q2 complement `upper_1s` | 1152 | 1152 | unchanged |
| Q2 complement `upper_1s_s` | 1164 | 1164 | unchanged |
| Q3 `guohe_adaptive_r` def | 178 | 178 | unchanged |
| Q3 defaults `orient`/`r_grid`/`v`/`B` | 185–188 | 185–188 | unchanged |
| Q3 `stopifnot` constraints | 199–201 | 199–201 | unchanged |
| Q3 inner-CV `B` | 255 | 255 | unchanged |
| Q3 refit `B` | 285 | 285 | unchanged |
| Q4 `b_adapt` parse / print / store | 54, 67, 127 | 54, 67, 127 | unchanged |
| Q4 `gh_one_rep` call (`B = b_boot`) | 109 | 109 | unchanged |
| Q4 `GH_R_GRID` | sim:52 | sim:52 | unchanged |
| Q4 adaptive call | sim:184–190 | sim:184–188 | unchanged (v1's 190 overran the closing paren at 188) |
| Q4 RUN.md commands | 57, 60 | 57, 60 | unchanged |
| Q4 RUN.md cost ledger | 98–99 | 98–99 | unchanged |
| Q6 `--cells` flag block | 45–53 | 45–53 | unchanged |
| Q6 `MF_CELLS` | 55–58 | 55–58 | unchanged |
| Q6 `MF_R_OUT`/`MF_R_IN` | 42–43 | 42–43 | unchanged |
| Q6 seed row | 157 | 157 | unchanged |
| Q6 `mf_rep_52` gate call | 151–155 | 151–155 | unchanged |
| Q6 `naive_ok`/`sel_ok` | 142–149 | 142–149 | unchanged |
| Q6 `cur_ok` | 80–84 | 80–84 | unchanged |
| Q6 bundle-level `stopifnot` | 179–181 | 179–181 | unchanged |
| Q6 E5 join assertion | 206–210 | 206–210 | unchanged |
| Q6 `field_seed_offset` | 223 | 223 | unchanged |
| Q6 `.mf_gh_cols` | 66–74 | 66–74 | unchanged |
| Q6 seed constants (`MV_*`, `mv_gh52_base`) | 60–62, 79 | 60–62, 79 | unchanged |
| A4 pilot scaffolding | 70, 81–83, 144, 214–225 | 70, **81–84**, 144, 214–225 | unchanged (v1 omitted the closing brace at 84) |
| `GH52_R_GRID` | sim:85 | sim:85 | unchanged |

Selected re-quotes, verbatim from the current tree:

**Q1** — `quarto/GuoHe/mr_vs_guohe_sim.R:124-128` and `:143-150`:

```r
mv_mr <- function(df, cands, sel_label, spec, draws = MV_DRAWS,
                  multiplier = MV_MULTIPLIER, seed = NULL,
                  ci_method = "ij", field_R_out = 1000L, field_R_in = 500L,
                  field_uniform = FALSE, field_complement = FALSE,
                  include_complement = FALSE, ij_residual = "two_term") {
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

**Q2** — `R/fs_mr_inference.R:531, 534, 539, 541, 545`; the three defaults A2 depends on are
unchanged and still match the recorded adoption (`field_complement = TRUE`,
`field_scale_complement` → `"selected"`, `return_reselection = TRUE`).

**A2's transplant source** — `R/fs_mr_inference.R:1219-1227`, inside `.fs_mr_field_joint`
(definition at `:1203`, closing brace `:1228`):

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
}
```

wired for both scalings at `R/fs_mr_inference.R:1149-1152`:

```r
  joint <- if (!is.null(lam_H) && is.finite(beta_deb))
    .fs_mr_field_joint(lam_H[ok_c], lf, beta_deb, bdc, to_eff, alpha) else NULL
  joint_s <- if (scale_on && !is.null(lam_H) && is.finite(beta_deb))
    .fs_mr_field_joint(lam_H[ok_c], lfs, beta_deb, bdc, to_eff, alpha) else NULL
```

So `bonf_lower_H` is the 0.975 one-sided lower on Ĥ, `bonf_upper_Hc` the 0.975 one-sided upper on
Ĥᶜ, and `bonf_joint_prob` the both-correct construction. **No `R/` change is needed**, consistent
with A2's prohibition on one.

**Q3** — `R/guohe_adaptive_r.R:185-188` (`orient = -1`, `r_grid = c(0.03, 0.10, 0.20, 0.30, 0.40,
0.45)`, `v = 5L`, `B = 200L`); one-B coupling confirmed again at `:255` (inner CV) and `:285`
(final refit). Constraints at `:199-201` admit T2's pinned grid and `orient = +1`.

**Q4** — `quarto/GuoHe/guohe_reproduction_run.R:54, 67, 109, 127`; `guohe_reproduction_sim.R:52`
(`GH_R_GRID <- c(1 / 3, 1 / 12, 1 / 21, 1 / 30)`) and `:184-188` (the adaptive call, `orient = -1`,
`B = B`, `refit = TRUE`, no `seed`). The inert-flag finding re-verifies exactly as in v1 and is
written up as `dev/notes/NOTE_adaptive_B_inert_2026-09-09.md` per A6.

**Q6 / A4** — all anchors re-verify; `quarto/GuoHe/guohe_sec52_run.R:70, 81-84, 144, 208-225`
supplies the `--pilot` transplant, including the projection block:

```r
if (pilot && !is.null(pilot_res)) {
  per_rep <- mean(pilot_res$rep_elapsed_sec, na.rm = TRUE)
  total_core_h <- per_rep * 6 * 2000 / 3600
```

---

## 4. Carried forward for T1 (A3) — the argument-assembly block to copy verbatim

A3 requires the T1 driver to copy `mv_mr()`'s executing argument-assembly lines rather than
re-author them. Those lines are `quarto/GuoHe/mr_vs_guohe_sim.R:132-151`, cited here so the T1
driver's provenance is on the record:

```r
  args <- list(
    df = df, candidates = cands, spec = spec,
    selected_members = cands[[sel_label]],
    admission = list(effect_floor = NULL, consistency = NULL),
    reselection = "maxeff",
    draws = draws, multiplier = multiplier,
    ci_method = ci_method, seed = seed, return_reselection = TRUE,
    include_complement = include_complement,
    # winner-only IJ variants (TASK_complement_refinements_2026-09-06);
    # "two_term" keeps the reported interval byte-identical.
    ij_residual = ij_residual)
  if (identical(ci_method, "field"))
    args <- c(args, list(field_R_out = field_R_out, field_R_in = field_R_in,
                         # kappa(Sigma-hat) sweep (TASK_mr_field_uniform_2026-09-05);
                         # FALSE keeps the 2026-09-05 field output byte-identical.
                         field_uniform = field_uniform,
                         # complement field (TASK_mr_field_complement_2026-09-06);
                         # FALSE keeps the field output byte-identical.
                         field_complement = field_complement))
  do.call(forestsearch:::fs_mr_inference, args)
```

`quarto/GuoHe/mr_vs_guohe_sim.R` is **not modified** by this task, per A3.

---

## 5. Carried forward for T2 (A5) — stored-bundle facts

`quarto/GuoHe/guohe_repro_t7_beta2_00.rds` metadata (read-only):

```
id t7_beta2_00 | target Table 7 | beta2 0 | n 400 | n_rep_requested 2000 | n_rep_used 2000
B 2000 | r_grid 0.333333333333333 0.0833333333333333 0.0476190476190476 0.0333333333333333
orient 1 | pilot FALSE | seed_base 92902757 | truth_identity 2e+06 0.25 20260721 TRUE exact
elapsed_sec 11569.456 | results 2000 x 28
```

- Produced by `guohe_sec52_run.R` (`orient = +1`), not `guohe_reproduction_run.R`; no `adaptive`,
  `adaptive_B` or `v` field — the Adaptive column has never been run on t7, which is what T2 supplies.
- `quarto/GuoHe/guohe_sec52_sim.R:85` — `GH52_R_GRID <- c(1 / 3, 1 / 12, 1 / 21, 1 / 30)` — is
  exactly the grid A5 pins, so T2's `r_grid` matches the stored fixed-r column set element for element.
- **The primary Adaptive bound is recoverable exactly from the stored columns.** The bound itself
  is not stored, but `r{i}_dist = gamma_s - lower`, so `lower = gamma_s - r{i}_dist` — the same
  reconstruction `mr_field_vs_guohe_run.R:70` already performs:

```r
    out[[sprintf("gh_r%d_low", i)]] <- theta - row[[sprintf("r%d_dist", i)]]
```

  The "stored B = 2000 Algorithm-3 bound at r̂" is therefore a `which(r_grid == r_hat)` index into
  those four columns; no re-run of the fixed-r columns is required.
- Under A5, T2's own `guohe_adaptive_r()` call uses `B = 2000`, matching the stored bundles' `B`
  and giving resolution parity with the fixed-r columns.

---

## 6. Stage 0 v2 verdict

| A1 requirement | verdict |
|---|---|
| new record, v1 not modified | done — v1 untouched at `33d0add3` |
| `git log -1 --oneline` + `git status -sb` + untracked set, verbatim | §1; untracked set empty |
| all Q1–Q6 line numbers re-verified post-merge | §3 — **all six verify**, none moved; two v1 range labels corrected, no content change |
| input inventory re-asserted with byte sizes | §2a/§2b — 24 bundles + 6 truth caches + 11 scripts/sources, all present, sizes identical to v1 |
| record whether the certification record is present | §2c — still absent; per A2 not a STOP |

**No A8 STOP condition is met.** Proceeding to A6 (the NOTE), then T1 per v1 as amended by A3/A4,
then T3. T2's driver is authored and committed but not executed on this Mac.

`git status -sb` at the close of this record (before it is staged):

```
## feature/glm-extension...origin/feature/glm-extension [ahead 3]
```

Untracked: none, unchanged from the opening. No network git was invoked.

---

# APPENDED — Stage 1 identity probe (A3), 2026-09-09

**OUTCOME: STOP. The probe's `identical()` assertion fails, which is a retained A8 STOP
condition. The cause is NOT the T1 driver and NOT the complement: the stored bundles were
built on Linux/x86_64/reference-BLAS and this Mac is arm64/Accelerate, so byte-identity to
them is unattainable here by any code path. Production T1 was not launched. Larry decides.**

## 1. What was run

- `devtools::install()` first, per A3 (workers see only the installed package): forestsearch
  0.3.5 installed from this tree, `* DONE (forestsearch)`.
- `Rscript quarto/GuoHe/mr_field_complement_vs_guohe_run.R --probe --cells=t7_beta2_00,t7_beta2_05 --reps=3`
- Columns compared by `identical()` against the stored `mr_field_vs_guohe_<id>.rds` rows:
  4 naive, 11 IJ, 19 field.

## 2. The literal result

| cell | naive | IJ | field | `naive_ok` | `cur_ok` |
|---|---|---|---|---|---|
| `t7_beta2_00` | IDENTICAL | MISMATCH | MISMATCH | 3/3 | 0/3 |
| `t7_beta2_05` | IDENTICAL | MISMATCH | MISMATCH | 3/3 | 0/3 |

- 46 column instances differ. Every difference is in the last one or two bits of the double.
  Examples as printed:
  - `IJ/mr_est`, `t7_beta2_05`, m = 1 — new `0.66115334364706191` vs stored `0.6611533436470618`.
  - `field/fld_lower_1s`, `t7_beta2_00`, m = 1 — new `-0.13077819374229155` vs stored `-0.13077819374229166`.
- **The naive columns are bit-identical in both cells.** They are computed by an independent
  `coxph` loop with no multiplier draws and no IJ algebra.

## 3. Diagnosis — three tests, in the order they were run

### 3a. Provenance of the stored bundles

`readRDS("mr_field_vs_guohe_t7_beta2_00.rds")$sessionInfo`:

```
R version 4.6.1 (2026-06-24)
Platform: x86_64-pc-linux-gnu
Running under: Pop!_OS 24.04 LTS
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
```

This Mac: arm64, R 4.5.2, Accelerate. Different architecture, different BLAS/LAPACK,
different R minor version.

### 3b. Control — the ORIGINAL complement-disabled path, on this Mac

The untouched `mv_mr()` (its own `field_complement = FALSE, include_complement = FALSE`
defaults), same seeds, same 3 replicates, same two cells — i.e. exactly what the committed
`mr_field_vs_guohe_run.R` does:

```
=== CONTROL (complement DISABLED, original mv_mr) t7_beta2_00 ===
  m=1 naive_est     identical=TRUE  reldiff=0
  m=1 mr_se_ij      identical=FALSE reldiff=3.39e-16
  m=1 fld_lower_1s  identical=FALSE reldiff=8.49e-16
  m=1 fld_lambda_sd identical=FALSE reldiff=6.4e-16
  ...
```

**The complement-disabled path fails `identical()` against the stored bundle too, at the same
~1e-16.** The probe's failure is therefore not attributable to enabling the complement, and
not to the hand-assembled call: no code path on this machine reproduces those bundles bit-for-bit.

### 3c. Isolation — complement enabled vs disabled, SAME machine, SAME seeds

The test that actually answers A3's question. 29 shared columns (naive, IJ, field, `p_hat_H`)
compared between the two calls:

```
=== ISOLATION t7_beta2_00 ===        === ISOLATION t7_beta2_05 ===
  m=1  29/29 shared columns identical    m=1  29/29 shared columns identical
  m=2  29/29 shared columns identical    m=2  29/29 shared columns identical
  m=3  29/29 shared columns identical    m=3  29/29 shared columns identical

=== ISOLATION VERDICT: enabling the complement perturbs NOTHING
    (all shared columns identical()) ===
```

**Enabling the complement perturbs no stored column.** This is the substantive proof A3 wanted,
and it matches the engine's own design claim at `R/fs_mr_inference.R:830-834` — the field block
re-seeds at `seed + 900000L` (`:828`) and draws `Xo`/`Xi_f` in the same order and sizes whether
or not the complement runs; the complement block re-reads them (`:897-918`) and consumes no RNG
of its own (`:710-733` is Cox fits and matrix algebra on already-drawn `Xi`).

**There is no RNG-stream finding.** The v1 §4 phrasing ("the engine's complement block precedes
the field block") is also inverted on this reading: the complement block *follows* the harm field
and re-reads its draws.

## 4. Magnitude of the cross-machine deviation

Over 180 compared values (6 replicates × 30 columns), new Mac run vs stored Linux bundle:

| statistic | value |
|---|---|
| compared values | 180 |
| bit-identical | 63 (35.0%) |
| max relative deviation | **9.51e-14** (`fld_q50`, `t7_beta2_05`, m = 1 — a quantile whose value is 7.3e-4, so the relative figure is inflated by a small denominator) |
| typical deviation | ~1e-16 relative |
| `all.equal` at 1e-8 | **TRUE** |

`mr_mean_r` is excluded from the relative figure and scored on absolute difference: it is ~0 by
construction (`-2.2e-17` vs `1.2e-18`, both numerical noise about zero; the engine documents it
as "exposure only").

## 5. Why this is Larry's call, not CC's

- A8 retains "any `identical()` failure in probe or production pairing" as a STOP. The literal
  condition is met, so CC stops.
- But the condition as written cannot be satisfied on this machine for any code, and D2 places
  Phase A on the Mac. That is a conflict between the spec and the hardware, not a defect in the
  work — and resolving it is a decision about evidentiary standard.
- Precedent exists in this repo: `TASK_tier2_mac_2026-09-08` hit the same wall and Larry gated
  truth at `~1e-8` with `all.equal` while keeping `identical()` for integer-valued quantities
  ("truth `identical()` FALSE but `all.equal` TRUE … the cross-machine BLAS difference").

### Options, costed — CC does not choose among them

1. **Gate the T1 pairing at a tolerance and run Phase A on the Mac as planned (recommended).**
   Adopt the tier2 precedent: `identical()` for the integer/selection columns (`naive_cover`,
   `sel`, `c_hat`, `n_sel`, `fld_n_out_used`), `all.equal` at 1e-8 for the floating-point ones.
   The measured worst case is 9.5e-14, six orders inside that. Requires a named-line change to
   the driver's probe/pairing check and Larry's word on the tolerance. Cost: minutes, then
   Gate 1a proceeds.
2. **Move Phase A to the Linux box.** Byte-identity to the stored bundles is then achievable
   and A3's assertion stands as literally written. Cost: T1 waits for the Linux box, which
   currently holds the Phase-B queue; no code change.
3. **Accept §3c as the pairing proof and record the cross-machine deviation as a quantified
   caveat.** No tolerance is introduced; the T1 record states that the isolation test, not
   byte-identity, is the evidence. Cost: none, but the T1 bundles are then not byte-comparable
   to the committed 16-cell record, which weakens D1's cross-reference sentence.

## 6. State at the STOP

- T1 driver `quarto/GuoHe/mr_field_complement_vs_guohe_run.R` authored and committed
  (`2e427309`); it parses, runs, and produces every required column. The probe's new columns
  computed cleanly, e.g. `t7_beta2_05`: `c_upper_95` 0.3357 / 0.6201 / −0.0704,
  `c_upper_975` 0.3875 / 0.6750 / −0.0299, `joint_cover` 1 / 1 / 0.
- **Gate 1a was not run. No production bundle was written. No T1 REPORT was written.**
- T2 driver and T3 qmd are downstream of this STOP in the A7 order and were not authored.
  Neither depends on the resolution above; both can proceed on Larry's word.
- Nothing under `R/` was changed. `mr_vs_guohe_sim.R` was not modified.
- No `git fetch`, `git pull` or `git push` at any point.

---

# APPENDED — v4 addendum, 2026-09-09: probe scored PASS under N1, and two corrections to the record

Governing amendment: `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09_v4.md` (`fb65b4a3`),
N10 item 2. Nothing above this line is modified.

## 1. Stage 1 probe — PASS under the N1 standard

The probe was **not re-run** (N2: its outcome is determined by evidence already recorded). Scored
against N1:

| N1 requirement | measured | verdict |
|---|---|---|
| integer / selection / flag / seed columns `identical()`, every replicate | **120/120** discrete comparisons across the six probe replicates | **PASS** |
| floating-point columns `all.equal()` at 1e-8 | worst absolute deviation **4.44e-16** (`mr_upper_2s`, `t7_beta2_00`, m = 1); worst relative **9.51e-14** (`fld_q50`, `t7_beta2_05`, m = 1, on a value of 7.3e-4); 63/180 (35.0%) bit-identical | **PASS** — six orders inside tolerance |
| complement exoneration | isolation test: complement enabled vs disabled, same machine, same seeds, **29/29 shared columns `identical()`**, both probe cells, all three replicates | **PASS** |

The 120 discrete comparisons cover: the selected cutpoint (`sel`, `c_hat` against both the repro
and stored-T1 bundles, `c_hat_naive`, `n_sel`), the engine's selection (`selected_index`,
`selected_label`, `sel_agree_mr`), admission/reselection indicators (`mr_ij_source`,
`mr_ij_draws`, `fld_n_out_used`), the coverage flags (`naive_cover`, `mr_cover`, `fld_cover_1s`,
`fld_cover_2s`, `gh_r1..r4_cover`), the seeds (`seed_data`, `seed_mr`), and the truth lookups
keyed by the selection (`gamma_s`, `gamma_s_naive`).

**Isolation-test citation for the T1 REPORT (N1):** recorded in the Stage 1 probe note appended
to this file, committed at **`7b2ed976`**.

**The cross-machine comparison is a provenance measurement, not a gate** (N1). Stored bundles:
`x86_64-pc-linux-gnu`, R 4.6.1, reference BLAS/LAPACK 3.12.0. This Mac: arm64, R 4.5.2,
Accelerate.

**Selection-stability margin (N3), for the record.** Top-1 minus top-2 oriented-score gap at
selection over 106 selections (the six probe replicates plus 100 fresh `t7_beta2_00` draws at
m = 2001–2100, disjoint from the stored 2000): minimum **1.438e-04**, 1st percentile
**3.961e-04**, median **5.205e-03**. Against the worst absolute float deviation of 4.44e-16 that
is a ratio of **3.2 × 10¹¹** (1st percentile **8.9 × 10¹¹**). The floating-point deviations
cannot reach a selection boundary, which is why every discrete column is bit-identical across the
two platforms while the continuous ones differ in the last bits.

**Verdict: Stage 1 probe PASSES under N1. Proceeding to Gate 1a.**

## 2. Correction 1 — block ordering (N6)

Re-verified against the current tree. **The correction stands, with one refinement: there are
two distinct complement blocks, and the v1 §4 error is a conflation of them, not a stale
description.**

- **The gate complement block** (`include_complement`), `R/fs_mr_inference.R:710-803`, **does
  precede** the field block. It consumes no RNG: Cox fits on complement row-sets and matrix
  algebra on the already-drawn `Xi`.
- **The complement *field* block** (`field_complement`), `R/fs_mr_inference.R:897-918`,
  **follows** the harm field's construction (`:882-895`) and re-reads its draws:

```r
      # -- Complement field (field_complement = TRUE) -- add-only and drawn
      # from NOTHING: it re-reads the harm field's xi (Xo / Xi_f) and winners
      # (G_out / W_in), so the harm field above and the uniform sweep below
      # are byte-identical whether or not it runs
```

- The field block re-seeds before drawing, `R/fs_mr_inference.R:828`, and its comment states the
  invariance directly (`:830-834`):

```r
    if (!is.null(seed)) set.seed(as.integer(seed) + 900000L)
```
```r
    # The raw N(0, I) multipliers are held in Xo / Xi_f so the complement
    # field (field_complement, TASK_mr_field_complement_2026-09-06) can
    # project the SAME xi through the complement's own influence.  Two rnorm
    # calls of the same sizes in the same order as before: the stream, Zo and
    # Zi are byte-identical whether or not the complement block runs.
```

**On the attributed source.** `quarto/GuoHe/REPORT_mr_field_stage1_2026-09-05.md:15` reads:

> The field block sits after the complement block, before the return assembly, and runs only
> under `ci_method == "field"`, drawing under the derived seed `seed + 900000L` **after** the
> main multiplier stream is fully consumed

That sentence is **accurate**, then and now — it refers to the *gate* complement block, which was
the only complement block in existence on 2026-09-05 (`field_complement` arrives with
`TASK_mr_field_complement_2026-09-06`). The defect is not in the 09-05 record but in v1 §4's
reading of it: applying a true statement about the gate complement block to the 09-06 complement
field block, and inferring a perturbation risk that the code forecloses. The executing code
governs; the isolation test measured the consequence and found none.

**The 09-05 work is not reopened. Its numbers are unaffected.**

## 3. Correction 2 — sync state (N6)

Re-verified:

- `git log -1 --format='%H %ci' origin/feature/glm-extension` →
  `a9c0d5e531d5a6d9d17c93340c24e830b9e2f4cc 2026-09-09 10:12:49 -0700` — still the same merge.
- The reflog's most recent network operation remains `HEAD@{...}: pull ... a9c0d5e5`; every entry
  since is this session's own commit.

**v2's "post-merge" framing is void**: no merge landed between the v1 and v2 readings, which is
why every line number and byte size re-verified unchanged.

`REVIEW_certification_2026-09-09.md` **remains absent** from the tree. The **A2 pending-citation
marker path is live** for B6, and `guohe_supp_section.qmd` already renders
`[certification citation pending sync]` at both figures (verified by executing its chunks:
`CERT_PRESENT` is `FALSE`).
