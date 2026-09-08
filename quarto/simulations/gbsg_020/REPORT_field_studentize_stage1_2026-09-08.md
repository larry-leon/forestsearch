# REPORT — Studentized complement field: Stage 0 (discovery) and Stage 1 (instrumentation, add-only `R/` change; identity gates)

**Task:** `dev/tasks/TASK_field_studentize_stage1_e0_2026-09-08.md` (1bb8bdf1), Stage 0 + Stage 1. Governing proposal: `dev/tasks/PROPOSAL_complement_field_scale_2026-09-08_v2.md` (committed alongside). Decisions D-1–D-4 at defaults.
**Date:** 2026-09-08. Executor: Claude Code (Linux), unattended. One `R/` change authorized, classified **adds code; byte-identical defaults**. Winner-only and winner-floor excluded from every table and line.

---

## GATE 1: PASS — G1a (decompose off): 131 / 131 pre-existing non-timing columns `identical()` to the committed nb20-A HR 1.50 rows 1–5, `truth` identical, the four new columns present and all NA; G1b (decompose on): 131 / 131 pre-existing columns `identical()` to committed, new columns finite on every detected row, `fld_Hc_scale_ratio` 0.97–1.14 (> 0 everywhere). Installed bodies verified by `deparse()` against source. Stage 0: every cited line re-quoted from HEAD (1bb8bdf1) matches the proposal's snapshot; no STOP.

## Stage 0 — Discovery (quotes from HEAD at 1bb8bdf1, `R/fs_mr_inference.R`, line numbers of that tree)

**Signature and `sel` (989–992):**
```r
.fs_mr_field_complement <- function(df, spec, kept, Bc, bh_c, winset, sel, bdc,
                                    G_out, W_in, Xo, Xi_f, to_eff, z975,
                                    lam_H = NULL, beta_deb = NA_real_,
                                    alpha = 0.05) {
```
**Lazy-fit loop (1004–1015):**
```r
  readings <- c(G_out[!is.na(G_out)], W_in[!is.na(W_in)])
  need     <- sort(unique(readings))
  new_fit  <- setdiff(need, winset)
  for (w in new_fit) {
    comp_idx <- setdiff(seq_len(Nall), kept[[w]])
    if (length(comp_idx) < 6L) next
    pcc <- tryCatch(.fs_mr_pieces(df[comp_idx, , drop = FALSE], spec),
                    error = function(e) NULL)
    if (is.null(pcc) || length(pcc$dfbeta) != length(comp_idx)) next
    Bc[comp_idx, w] <- pcc$dfbeta
    bh_c[w] <- pcc$beta_hat
  }
```
**`Zo_c` / `Zi_c` and the `lam_c` loop (1016–1031):**
```r
  fit_ok <- is.finite(bh_c)
  Zo_c <- crossprod(Bc, Xo)                  # Ncol x R_out : zeta^c (outer)
  Zi_c <- crossprod(Bc, Xi_f)                # Ncol x R_in  : zeta'^c (inner)
  lam_c  <- rep(NA_real_, R_out)
  n_in_c <- rep(NA_real_, R_out)
  n_drop_unfit <- 0L
  for (r in which(!is.na(G_out))) {
    G <- G_out[r]
    if (!fit_ok[G]) { n_drop_unfit <- n_drop_unfit + 1L; next }
    wi <- W_in[r, ]
    ok_in <- which(!is.na(wi))
    ok_in <- ok_in[fit_ok[wi[ok_in]]]
    if (!length(ok_in)) { n_drop_unfit <- n_drop_unfit + 1L; next }
    lam_c[r]  <- Zo_c[G, r] - mean(Zi_c[cbind(wi[ok_in], ok_in)])
    n_in_c[r] <- length(ok_in)
  }
```
**Bound assembly (1042–1045, 1057–1064):**
```r
  lf <- lam_c[ok_c]
  qs <- stats::quantile(lf, c(.05, .25, .50, .75, .95, .025, .975),
                        names = FALSE, type = 7)
  est2_w <- bdc - mean(lf)
  sd_c <- stats::sd(lf)
  ...
    est2 = to_eff(est2_w),
    # Primary: the one-sided 95% UPPER bound (benefit claim "at most U").
    upper_1s = to_eff(bdc - qs[1]),
    lower_1s = to_eff(bdc - qs[5]),
    lower_2s = to_eff(bdc - qs[7]), upper_2s = to_eff(bdc - qs[6]),
    se_field = sd_c,
    lower_se = to_eff(est2_w - z975 * sd_c),
    upper_se = to_eff(est2_w + z975 * sd_c)),
```
**Call site in the `ci_method == "field"` block (868–876), passing `sel`:**
```r
        fcres <- .fs_mr_field_complement(
          df = df, spec = spec, kept = kept, Bc = Bc, bh_c = bh_c,
          winset = winset, sel = sel, bdc = bdc_w,
          G_out = G_out, W_in = W_in, Xo = Xo, Xi_f = Xi_f,
          to_eff = to_eff, z975 = z975,
          # Joint (H lower, Hc upper) pair from the aligned outer draws
          # (TASK_complement_refinements_2026-09-06, method B); add-only,
          # attached as field$joint, no new draws.
          lam_H = lam, beta_deb = beta_deb)
```
**RNG discipline (794–800; raw ξ held in `Xo` / `Xi_f`):**
```r
    # The raw N(0, I) multipliers are held in Xo / Xi_f so the complement
    # field (field_complement, TASK_mr_field_complement_2026-09-06) can
    # project the SAME xi through the complement's own influence.  Two rnorm
    # calls of the same sizes in the same order as before: the stream, Zo and
    # Zi are byte-identical whether or not the complement block runs.
    Xo   <- matrix(stats::rnorm(Np * field_R_out), Np, field_R_out)
    Xi_f <- matrix(stats::rnorm(Np * field_R_in), Np, field_R_in)
```
All match the proposal's snapshot quotes (lines 989, 1004–1018, 1022–1031; 794–800) — no material difference; no STOP.

**Template recorder anchor:** `fld_Hc_se = NA_real_` at `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:829` (record: `rec$fld_Hc_se <- fc$se_field` at `:1073`). **Committed nb20-A bundles:** stems `fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_p30sgnb20` and `..._h175_..._nb20_p30sgnb20` (`_res_1_1000`, `_res_1001_2000`, `_combined_1_2000`); meta `seed_base = 8316951`, `campaign_tag = p30sgnb20`, `sg_focus = effMaxSG`, `effect_neighborhood = 0.20`, `er_jcuts = 10`, `harm_z1_quantile = 0.60`, `field_complement = TRUE`, `ij_residual = two_term`, `fb_mode = none`, `ci_method = field`, `forestsearch_version = 0.3.5`, `n_workers = 100`; per-replicate seed `8316951 + sim_id`. 136 columns, 1000 rows per batch; the HR 1.75 batch 1 has all of rows 1–200 detected.

**One forwarding finding (recorded, acted on):** the template reaches `fs_mr_inference()` through `forestsearch(mr_inference = TRUE, mr_inference_args = ...)`, and `forestsearch_main.R:3388–3421` forwards each MR knob **by name** (`field_complement = .g_mr(mr_inference_args$field_complement, FALSE)` at `:3416`); a knob absent from that list never reaches the gate. The task's S1a names only `fs_mr_inference.R`, but S1b's knob → `field_decompose` path is unreachable without the same one-line pass-through that every prior knob received. Added (`forestsearch_main.R`, 3 lines incl. comment, default `FALSE`, the identical add-only pattern) and covered by G1a. Flagged here as the one departure from S1a's file list; it is inside the task's classification ("adds code; byte-identical defaults").

## Stage 1 — The add-only edit (diff summary: 5 files, 87 insertions, 4 deletions; every deletion a line re-formed by an insertion)

| File | Change |
|---|---|
| `R/fs_mr_inference.R` (+40 / −1) | `fs_mr_inference()` gains `field_decompose = FALSE` (roxygen `@param`); forwarded at the complement call site (`field_decompose = field_decompose`); `.fs_mr_field_complement()` gains `field_decompose = FALSE`; the decomposition block inserted **after** `sd_c <- stats::sd(lf)` and **before** the joint / `complement <- c(list(...))` assembly, exactly as specified (selected-complement lazy fit only if `!fit_ok[sel]`, then `s <- sqrt(colSums(Bc * Bc))`, `s_win`, `zg`, `mi`, the seven scalars); `complement$decomp_fields` appended when TRUE, absent when FALSE. Comment records that the `sel` lazy fit runs after `lam_c` is final. |
| `R/forestsearch_main.R` (+3) | `field_decompose = .g_mr(mr_inference_args$field_decompose, FALSE)` pass-through (see the finding above). |
| `man/fs_mr_inference.Rd` (+10) | regenerated by `devtools::document()`. |
| `NEWS.md` (+8) | new "development version" header with the one bullet. |
| Template (+27 / −2) | knob `FS_S7_FIELD_DECOMP` (default FALSE) → `mr_field_decompose` → `mr_inference_args$field_decompose`; echoed in the "Template knobs:" audit line and the MR-settings-readout table; recorder gains `fld_Hc_scale_sel`, `fld_Hc_scale_win`, `fld_Hc_scale_cv`, `fld_Hc_scale_ratio` beside `fld_Hc_se`, filled `%||% NA_real_` from `fc$decomp_fields`. Two "deletions" are the re-formed `cat(sprintf(...))` line and the readout's `setting` vector. |

`devtools::install(dependencies = FALSE)` (never `load_all()`): installed 0.3.5; `deparse()` of the installed `fs_mr_inference` and `.fs_mr_field_complement` identical to the source tree's; installed `forestsearch()` carries the pass-through.

## Identity gates (template-driven, nb20-A HR 1.50 n500 config: `FS_S7_FOCUS=effMaxSG FS_S7_Z1Q=0.60 FS_S7_NBHD=0.20 FS_S7_N=500 FS_S7_HR=1.50 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_NSIMS=5 FS_S7_START=1`, 5 workers, seeds `8316951 + sim_id`, sim_id 1–5; `gate_identity.R`, session scratchpad; compared against the committed `..._nb20_p30sgnb20_res_1_1000.rds` rows 1–5, all 131 non-timing columns of the committed bundle — the five timing columns `fb_secs`, `fit_mr_secs`, `fld_H_secs`, `fld_Hc_secs`, `fld_H_uniform_secs` excluded as always)

**G1a — `FS_S7_FIELD_DECOMP` unset (tag `stud1inert`, 68 s render): PASS.** 131 / 131 pre-existing columns `identical()`; `truth` `identical()`; fresh bundle 140 columns = 136 + the four new, all NA.

**G1b — `FS_S7_FIELD_DECOMP=TRUE` (tag `stud1decomp`, 67 s render): PASS.** 131 / 131 pre-existing columns `identical()` (the block reads existing objects only); `truth` `identical()`; new columns finite on 5 / 5 detected rows:

| sim_id | \|Ĥ\| | p̂(Ĥ) | `nv_Hc_se` | `fld_Hc_se` | `fld_Hc_scale_sel` | `fld_Hc_scale_win` | `fld_Hc_scale_cv` | `fld_Hc_scale_ratio` (ρᶜ) |
|---|---|---|---|---|---|---|---|---|
| 1 | 161 | 0.013 | 0.1509 | 0.1340 | 0.1509 | 0.1320 | 0.050 | 1.143 |
| 2 | 69 | 0.186 | 0.1261 | 0.1264 | 0.1261 | 0.1297 | 0.057 | 0.973 |
| 3 | 113 | 0.064 | 0.1394 | 0.1335 | 0.1394 | 0.1324 | 0.047 | 1.053 |
| 4 | 180 | 0.037 | 0.1485 | 0.1400 | 0.1485 | 0.1343 | 0.062 | 1.105 |
| 5 | 87 | 0.159 | 0.1409 | 0.1347 | 0.1409 | 0.1384 | 0.033 | 1.018 |

Two reads, not checks: `fld_Hc_scale_sel` equals `nv_Hc_se` to four digits on every row (the influence-norm √Σ dfbeta² is the robust SE — the scale object is the naive SE itself, as §3 of the proposal says); and on these five rows ρᶜ is largest where p̂ is smallest and |Ĥ| largest (1.14 / 1.11 at p̂ 0.013 / 0.037) and ≈ 1 at p̂ 0.16–0.19 — the direction §3 predicts, to be tested on 200 at E0.

Gate bundles and renders committed beside the campaign's: `..._nb20_stud1inert_res_1_5.rds` / `_batch_1_5.html`, `..._nb20_stud1decomp_res_1_5.rds` / `_batch_1_5.html`.
