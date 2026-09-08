# TASK — Studentized complement field, Stage 1 (instrumentation, add-only `R/` change) and E0 (instrumented smoke; report-and-wait)

Date: 2026-09-08. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (go of 2026-09-08, defaults P-1–P-6 of the proposal). Reviewer: the Linux MR-field chat.
Governing proposal: `PROPOSAL_complement_field_scale_2026-09-08_v2.md` (committed alongside this file). Predecessors: `REPORT_complement_variance_2026-09-07.md` (A0–A4), `REPORT_banddial_2026-09-07.md` + Gate 2, nb20 records. Diagnosis and theory are settled there; this task implements §4 (Stage 1) and §6 (E0) of the proposal only.

## Protocol

- First action: archive any stale task/proposal variants in `~/Downloads` to `~/Downloads/cc_archive/`, then copy this file and `PROPOSAL_complement_field_scale_2026-09-08_v2.md` to `dev/tasks/` and commit. Do not push.
- **Scope: Stage 1 + E0 only.** E1 (the field-s campaign) is a separate later session; no R1/R2/R0 variant decision is made in this session — the record reports, the Linux chat decides.
- **This task authorizes one `R/` change**, classified **adds code; byte-identical defaults** — the add-only edit specified below and nothing else under `R/`. The identity gates are mandatory and stop on failure without asking.
- Standing conventions: winner-only and winner-floor excluded from every table and line; verify from source (Stage 0 quotes from HEAD, executing bodies by `deparse()` after install); records beside results; `.refuse_if_tracked()` respected; `devtools::install()` (never `load_all()`) before any parallel run; fail-fast; committed work never re-run (identity to committed rows is the pairing proof).

## Stage 0 — Discovery (quote from HEAD; STOP on material difference)

Quote in the record, with line numbers from the current tree: the `.fs_mr_field_complement()` signature and its `sel` argument; the lazy-fit loop; `Zo_c`/`Zi_c` construction and the `lam_c` loop (the proposal's cited lines 989, 1004–1018, 1022–1031 at the 2026-09-08 snapshot); the bound assembly (`est2`, `upper_1s`, `se_field`); the call site in the `ci_method == "field"` block that passes `sel`; and the RNG-discipline comment (raw ξ in `Xo`/`Xi_f`). STOP if any differs materially from the proposal's snapshot quotes. Confirm the template's recorder location for `fld_Hc_se` (anchor for the new columns) and the committed nb20-A bundles' stems and seeds (`8316951 + sim_id`).

## Stage 1 — The add-only edit, recorder, and identity gates

**S1a. `R/fs_mr_inference.R`:**

1. `fs_mr_inference()` gains `field_decompose = FALSE` (roxygen `@param`, one sentence), forwarded at the complement call site (`field_decompose = field_decompose`).
2. `.fs_mr_field_complement()` gains `field_decompose = FALSE`. When TRUE, insert a block **after** `lf <- lam_c[ok_c]` and the quantile/`sd_c` lines and **before** the `complement <- c(list(...))` assembly — it reads existing objects only and touches nothing upstream (so `lam_c`, `lf`, quantiles, `est2`, and every existing output are unchanged even when TRUE):

```r
if (isTRUE(field_decompose)) {
  if (!fit_ok[sel]) {                      # ensure the selected complement's column is populated
    comp_idx <- setdiff(seq_len(Nall), kept[[sel]])
    if (length(comp_idx) >= 6L) {
      pcc <- tryCatch(.fs_mr_pieces(df[comp_idx, , drop = FALSE], spec), error = function(e) NULL)
      if (!is.null(pcc) && length(pcc$dfbeta) == length(comp_idx)) {
        Bc[comp_idx, sel] <- pcc$dfbeta; bh_c[sel] <- pcc$beta_hat; fit_ok[sel] <- TRUE
      }
    }
  }
  s  <- sqrt(colSums(Bc * Bc))             # per-candidate complement noise scale; 0 for unfit
  gG <- G_out[ok_c]; s_win <- s[gG]
  zg <- Zo_c[cbind(gG, ok_c)]; mi <- zg - lf
  decomp_fields <- list(
    scale_sel      = if (fit_ok[sel]) s[sel] else NA_real_,
    scale_win_mean = mean(s_win),
    scale_win_cv   = stats::sd(s_win) / mean(s_win),
    scale_ratio_c  = if (fit_ok[sel]) s[sel] / mean(s_win) else NA_real_,
    var_zeta_G = stats::var(zg), var_m_in = stats::var(mi), cov_zeta_m = stats::cov(zg, mi))
}
```

3. Append `decomp_fields` to the `complement` list when TRUE (absent when FALSE). Note in a comment: the `sel` lazy fit runs after `lam_c` is final, so it cannot alter any existing output.
4. `NEWS.md`, one development-version bullet: `fs_mr_inference()` gains `field_decompose` (add-only complement-field scale diagnostics; PROPOSAL_complement_field_scale_2026-09-08_v2).

**S1b. Template (document-level, add-only):** env knob `FS_S7_FIELD_DECOMP` (default FALSE) → `field_decompose`; knob echoed in the settings readout; recorder gains four columns beside `fld_Hc_se`, `%||% NA_real_`: `fld_Hc_scale_sel`, `fld_Hc_scale_win`, `fld_Hc_scale_cv`, `fld_Hc_scale_ratio`.

**S1c. `devtools::install()`, then the identity gates** (template-driven, nb20-A HR 1.50 n500 config, seeds `8316951 + sim_id`, sim_id 1–5):

- **G1a (default off):** `FS_S7_FIELD_DECOMP` unset — every pre-existing column `identical()` to the committed nb20-A bundle rows 1–5; the four new columns present, all NA.
- **G1b (decompose on):** `FS_S7_FIELD_DECOMP=TRUE`, same replicates — every pre-existing column `identical()` to committed; new columns finite with `fld_Hc_scale_ratio > 0` on detected rows.

STOP on any mismatch. Commit Stage 1 (R/, template, NEWS) with `REPORT_field_studentize_stage1_2026-09-08.md` recording the diff summary and both gate results with concrete values.

## E0 — Instrumented smoke (compute approved for this stage only)

**Run:** nb20-A HR 1.75 configuration exactly (`FS_S7_FOCUS=effMaxSG`, `FS_S7_Z1Q=0.60`, `FS_S7_NBHD=0.20`, `FS_S7_N=500`, HR 1.75, J = 10, `FS_S7_FIELD_COMPLEMENT=TRUE`, `FS_S7_IJ_RESIDUAL=two_term`, `FS_S7_FB=none`, `return_reselection = TRUE`) plus `FS_S7_FIELD_DECOMP=TRUE`, `FS_S7_CAMPAIGN=e0stud`, seeds `8316951 + sim_id`, **sim_id 1–200**, one batch, 100 workers. **Ceiling 30 min wall; hard timeout 1 h.** No other compute in this session.

**Gate E0-a (pairing, before any analysis):** every pre-existing column of the 200 rows `identical()` to the committed nb20-A HR 1.75 bundle rows sim_id 1–200 — the proof that the instrumentation is inert on everything old under production settings.

**E0 analysis** (script beside the results; every number verbatim in the report; tertiles computed within the 200 detected replicates; n ≈ 67 per tertile, so report directional findings with SEs/Wilson where applicable, not precision claims):

1. corr(ρᶜ, `fld_Hc_se`/`nv_Hc_se`) — Pearson and Spearman — and mean ρᶜ by p̂(Ĥ) tertile and by |Ĥ|/|H| tertile. Prediction on record: ρᶜ > 1 concentrated in low-p̂ / large-Ĥ; ≈ 1 at the other end.
2. Tertile table of (ρᶜ · `fld_Hc_se`)/`nv_Hc_se` beside the uncorrected `fld_Hc_se`/`nv_Hc_se` (current pattern ≈ 1.00 / 0.92 / 0.79–0.89). Prediction: the corrected ratio ≈ 1 in every tertile.
3. `fld_Hc_scale_cv`: mean, q50, q90 — the R1-vs-R2 per-draw stability input.
4. Context row: mean and median ρᶜ, share ρᶜ > 1, against the cell's committed average deficit (A1: λ/naive SE = 0.945 at nb20-A HR 1.75, i.e. implied mean ρᶜ ≈ 1.06).

Output: `REPORT_field_studentize_e0_2026-09-08.md` beside the e0stud results (Stage 0 quotes, Stage 1 reference, Gate E0-a result, tables 1–4). **No R1/R2/R0 recommendation in the record — report and wait.**

## Decisions (defaults in brackets)

- D-1 E0 replicates: 200 [default].
- D-2 Recorder columns: the four listed [default].
- D-3 E0 cell: nb20-A HR 1.75 n500 [default].
- D-4 E0 ceiling 30 min / timeout 1 h [default].

## Done means

Stage 0 quotes in the record; Stage 1 committed with G1a/G1b PASS and concrete values; e0stud bundle, rendered document, and the E0 report committed; branch left unpushed for Larry; session ends with a one-paragraph summary listing gates passed and the commit range. Out of scope, deferred to the next documents: E1, the variant decision, any bound computation with scaling enabled.
