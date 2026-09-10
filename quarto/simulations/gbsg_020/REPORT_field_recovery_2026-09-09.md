# REPORT — Field re-selection recovery: membership agreement against the observed subgroup (add-only)

Date: 2026-09-09. Task: `dev/tasks/TASK_field_recovery_2026-09-09.md`, governing proposal `dev/tasks/PROPOSAL_field_recovery_2026-09-09.md` (both committed as received, b69fffab). Executor: Claude Code (Linux). **Report-and-wait.** Classification: **adds code; byte-identical defaults**. Compute: two 5-replicate verification renders, one 5-replicate scratch re-run for the arithmetic check, one GBSG vignette rebuild, one validation read of committed bundles. No campaign.

Decisions honoured: **R-1 yes** (membership agreement is the deliverable), **R-2 deferred** (the rule-name family is out of scope — no covariate-name comparison anywhere in this change), **R-3** `field_recovery` / `fld_recov_`, **R-4** now, **R-5** validation against the metrics available on a committed cell.

Predecessors: `REPORT_print_vignette_2026-09-09.md`, `REPORT_cimethod_flip_2026-09-09.md`, `dev/notes/NOTE_survival_products_2026-09-09.md`, `R/forestsearch_cross_validation.R` (the vocabulary this transplants).

**First action, as directed.** `~/Downloads/TASK_print_vignette_2026-09-09.md` (the completed predecessor task) was moved to `~/Downloads/cc_archive/`. `HANDOFF_guohe_comparison_2026-09-09.md` was **not** archived. `R_09Sep2026/` and `R_09Sep2026.zip` were left in place: they are the uploaded source snapshot the proposal cites for provenance, not a stale spec variant. The seven pre-existing untracked files were left alone.

---

## Stage 0 — Discovery (quoted from HEAD before the change)

### 1. `kept <- candidates[asm$keep]` and its enclosing branch

`R/fs_mr_inference.R:702–713` at HEAD (`105186cd`):

```r
  # ---------------------------------------------------------------------------
  # Complement subgroup (optional).  The complement is induced by the selection,
  # not chosen independently, so its bias is the perturbation of the complement
  # of the re-selected winner on each draw.  Fit complements only for candidates
  # that win across draws (plus the selected one) to keep the cost small.
  # ---------------------------------------------------------------------------
  complement <- NULL
  mean_r_c   <- NA_real_        # stays NA when no complement is fit
  bdc_w      <- NA_real_        # beta-tilde^c on the working scale, read by the
                                # complement field block; NA when unfit
  if (isTRUE(include_complement)) {
    kept   <- candidates[asm$keep]              # aligns with asm columns
    Ncol   <- length(asm$names)
    Nall   <- nrow(df)
```

**Both inputs are in scope well before that branch.** `candidates` is normalised and, when the observed pick is absent from the family, extended, at `R/fs_mr_inference.R:550–560`:

```r
  if (!length(candidates)) candidates <- list()
  if (is.null(names(candidates)))
    names(candidates) <- paste0("cand", seq_along(candidates))

  # Ensure the observed selected subgroup is in the family and is the target.
  H_lab <- ".selected_H"
  hit <- if (length(candidates))
    which(vapply(candidates, function(ix) setequal(ix, selected_members), logical(1)))
  else integer(0)
  if (length(hit)) sel_lab <- names(candidates)[hit[1]]
  else { candidates[[H_lab]] <- selected_members; sel_lab <- H_lab }
```

and `asm` (hence `asm$keep`) is built immediately after, with `sel` two lines later, at `R/fs_mr_inference.R:562–564`:

```r
  asm <- .fs_mr_assemble(df, candidates, spec)
  if (is.null(t_confirm)) t_confirm <- if (asm$log_scale) 1 else 0
  sel <- match(sel_lab, asm$names)
```

Neither is reassigned between line 564 and line 711 (`grep -n 'candidates\|asm <-'` over the file shows no further assignment to either), and `df` is fixed at line 549 (`df <- as.data.frame(df)`).

**Confirmation that `kept` is a pure subset.** `candidates[asm$keep]` is a logical subset of a list — no fitting, no drawing, no RNG consumption, no side effects. `asm$keep` is itself a plain logical vector produced by `.fs_mr_assemble()` (`R/fs_mr_inference.R:88–106`), which does all its fitting *there* and returns the keep flags:

```r
.fs_mr_assemble <- function(df, candidates, spec) {
  ...
  for (g in seq_len(S)) {
    idx <- candidates[[g]]
    if (length(idx) < 6L) next
    pc <- tryCatch(.fs_mr_pieces(df[idx, , drop = FALSE], spec), error = function(e) NULL)
    if (is.null(pc) || length(pc$dfbeta) != length(idx)) next
    B[idx, g] <- pc$dfbeta
    bh[g] <- pc$beta_hat; sdv[g] <- pc$sigma_D; sz[g] <- length(idx)
    keep[g] <- TRUE; log_scale <- isTRUE(pc$log_scale)
  }
  list(B = B[, keep, drop = FALSE], ..., keep = keep, ...)
}
```

`Ncol <- length(asm$names)` and `Nall <- nrow(df)` are likewise a length and a row count. **All three are pure. The hoist is byte-identical, and the task's STOP condition is not met.** Gate R2a proves it empirically.

### 2. `G_out`, `sel`, `ok_c`, and the `ci_method` field gate

The gate, `R/fs_mr_inference.R:818–822`:

```r
  field <- NULL
  if (ci_method == "field") {
    t0f <- proc.time()
    if (!is.null(seed)) set.seed(as.integer(seed) + 900000L)
    Np <- nrow(B)
```

`G_out` — allocated and written inside the outer loop, `R/fs_mr_inference.R:846–874`:

```r
    fc <- isTRUE(field_complement) && isTRUE(include_complement)
    G_out <- if (fc) rep(NA_integer_, field_R_out) else NULL
    W_in  <- if (fc) matrix(NA_integer_, field_R_out, field_R_in) else NULL
    for (r in seq_len(field_R_out)) {
      v <- w + Zo[, r]
      G <- if (fast) which.max(v) else sel_one(v)
      if (is.na(G)) next
      if (fast) {
        win <- max.col(t(v + Zi), ties.method = "first")
        lam[r] <- Zo[G, r] - mean(Zi[cbind(win, ii)])
        n_in_used[r] <- field_R_in
        if (fc) { G_out[r] <- as.integer(G); W_in[r, ] <- as.integer(win) }
      } else {
        wi <- vapply(ii, function(j) sel_one(v + Zi[, j]), integer(1))
        ok_in <- which(!is.na(wi))
        if (!length(ok_in)) next
        lam[r] <- Zo[G, r] - mean(Zi[cbind(wi[ok_in], ok_in)])
        n_in_used[r] <- length(ok_in)
        if (fc) { G_out[r] <- as.integer(G); W_in[r, ] <- wi }
      }
    }
```

`G_out[r]` is the index, **within the kept family**, of the candidate the gate's own re-selection map picks on outer draw *r* — the perturbed field `v = w + zeta*_r` pushed through `S`. `NA` where a draw produced no winner.

`sel` is the observed pick's index in the same family (line 564, quoted above).

**`ok_c` is not in scope in the harm field block.** It lives inside the complement helper, `R/fs_mr_inference.R:1102–1103`:

```r
  ok_c <- which(is.finite(lam_c))
  counts <- list(n_out_used = length(ok_c),
```

— it is the complement field's usable-draw set, computed from the *complement's* `lam_c`, and it exists only when `field_complement = TRUE` **and** `include_complement = TRUE`. **Reported deviation from the task's wording:** using `ok_c` literally would gate the recovery diagnostics on the complement path, which is not what R-1 asks for. The harm field's own analogue is `ok_f <- which(is.finite(lam))` (`R/fs_mr_inference.R:876`), and the two loop bodies above show that `G_out[r]` is written on exactly the same statement group that writes `lam[r]` — a draw that `next`s writes neither. So `{r : !is.na(G_out[r])}` **is** `ok_f`, and the recovery block uses that set. The invariant `n_used + n_skipped = |ok_c|` is therefore tested as `n_used + n_skipped = n_draws = |ok_f|`, and R2b additionally checks `n_used` against the field's own recorded `n_out_used` (`fld_H_nout`) row by row.

Two lines are needed to make `G_out` reach the recovery block when the complement field is off; both are side assignments and both collapse to the previous behaviour when `field_recovery = FALSE`. See §R1.

### 3. The CV membership cross-tab this transplants

`R/forestsearch_cross_validation.R:1409–1424`:

```r
    if (n_sgfound == 2) {
      tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
      sensH <- tabit[1, 1] / sum(tabit[, 1])
      sensHc <- tabit[2, 2] / sum(tabit[, 2])
      ppvH <- tabit[1, 1] / sum(tabit[1, ])
      ppvHc <- tabit[2, 2] / sum(tabit[2, ])
    } else {
      tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
      sensH <- if (nrow(tabit) > 0 && ncol(tabit) > 0) tabit[1, 1] / sum(tabit[, 1]) else NA
      sensHc <- NA
      ppvH <- if (nrow(tabit) > 0) tabit[1, 1] / sum(tabit[1, ]) else NA
      ppvHc <- NA
    }

    sens_metrics_original <- c(sensH, sensHc, ppvH, ppvHc)
    names(sens_metrics_original) <- c("sens_H", "sens_Hc", "ppv_H", "ppv_Hc")
```

Rows are the CV assignment, columns the **original full-data** assignment; row/column 1 is the harm cell. So, with the observed subgroup in the role of `treat.recommend.original` and each draw's re-selection in the role of `treat.recommend`, and writing `a = |G_r ∩ Ĥ|`, `b = |G_r|`, `c = |Ĥ|`, `d = n − b − c + a`:

| CV name | CV formula | field analogue |
|---|---|---|
| `sens_H` | `tabit[1,1] / sum(tabit[,1])` | `a / c` |
| `ppv_H` | `tabit[1,1] / sum(tabit[1,])` | `a / b` |
| `sens_Hc` | `tabit[2,2] / sum(tabit[,2])` | `d / (n − c)` |
| `ppv_Hc` | `tabit[2,2] / sum(tabit[2,])` | `d / (n − b)` |

**On `npv_Hc`.** Read with `Ĥ` as the positive class, CV's `ppv_Hc` *is* the negative predictive value — `d / (n − b)` is the share of the patients a re-selection excludes that the observed pick also excludes. It is returned under both names (`ppv_Hc` for CV parity, `npv_Hc` for the standing convention that sensitivity, specificity, PPV and NPV are reported together); they are the same number by construction, and the roxygen says so. The recorder column list in the task names `fld_recov_npv_Hc`, so that is the name the recorder writes.

The rule-name family at `R/forestsearch_cross_validation.R:1396–1407` (`Any`, `Exact`, `At least 1`, `Cov1`, `Cov2`, `Cov 1 & 2`, `Cov1 exact`, `Cov2 exact`) is R-2 and is **not** transplanted.

### 4. The template's `fld_Hc_scale_*` recorder rows, and the `field_decompose` precedent

The survival template is `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (the `gate_d2_cim_*.qmd` drivers are gate-local copies of it; `diff` against the template is the two `ci_method` hunks and nothing else).

Knob, template line 600 at HEAD:

```r
mr_field_decompose <- identical(.env_chr("FS_S7_FIELD_DECOMP", "FALSE"), "TRUE")
```

forwarded at line 642 inside `mr_inference_args`:

```r
       field_decompose = mr_field_decompose,
       field_scale_complement = mr_field_scalec,
```

NA initialisers, template lines 874–881:

```r
  # Scale decomposition (field_decompose = TRUE; else all NA;
  # TASK_field_studentize_stage1_e0_2026-09-08): s_sel = the selected
  # complement's influence-norm scale sqrt(sum(Bc[, sel]^2)); s_win = the
  # draw-weighted mean over the used outer winners' scales, with its CV;
  # scale_ratio = s_sel / s_win (rho^c).
  fld_Hc_scale_sel = NA_real_, fld_Hc_scale_win = NA_real_,
  fld_Hc_scale_cv = NA_real_, fld_Hc_scale_ratio = NA_real_,
```

fill site, template lines 1140–1144:

```r
          # ---- Scale decomposition (field_decompose = TRUE; absent -> NA) ----
          rec$fld_Hc_scale_sel   <- fc$decomp_fields$scale_sel      %||% NA_real_
          rec$fld_Hc_scale_win   <- fc$decomp_fields$scale_win_mean %||% NA_real_
          rec$fld_Hc_scale_cv    <- fc$decomp_fields$scale_win_cv   %||% NA_real_
          rec$fld_Hc_scale_ratio <- fc$decomp_fields$scale_ratio_c  %||% NA_real_
```

and the `.g_mr` pass-through the knob needs to reach the gate, `R/forestsearch_main.R:3439`:

```r
        field_decompose = .g_mr(mr_inference_args$field_decompose, FALSE),
```

**The pass-through is required.** The template reaches `fs_mr_inference()` only through `forestsearch(mr_inference_args = ...)`, and `forestsearch()` names every MR argument explicitly at that call site; an argument absent from the `.g_mr` block is never forwarded. So `field_recovery` gets the same one line.

### 5. R-5 validation cell

**No committed cell carries FB or CV *recovery* metrics on the same replicates as any field render.** Stated plainly, because the proposal's §4 assumed one might exist. What is available:

1. **The gate's own comparator cell**, `results/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_e1stud_{res_1_1000, res_1001_2000, combined_1_2000}.rds` (2,000 replicates; `meta`: `ci_method=field`, `field_decompose=TRUE`, `field_scale_complement=selected`, `ij_residual=two_term`, `nbhd=0.20`, `hr=1.50`, `n=500`, `z1q=0.60`, `mr_draws=5000`, `seed_base=8316951`, package 0.3.5). It carries the field's exact-match rate `p_hat_H` **and** per-replicate `sens` / `spec` / `ppv` / `npv` — but the latter score `Ĥ` against the **true planted harm flag**, not against a re-run's re-discovery. Its FB columns (`fb_H_*`) are estimates, not agreement rates: the cell ran `fb_mode = "none"` and joined nothing. This is the cell used below, because it is the only one on **exactly the same replicates** as Gate R's five rows.

2. **The closest CV comparison in the matching vocabulary** is `quarto/applications/gbsg/_payloads_2026-09-01_complete/analysis_gbsg_survival_frozen_family/` — a committed real-data payload whose `extras$loo$sens_metrics` are CV's own `sens_H`/`sens_Hc`/`ppv_H`/`ppv_Hc` against the **observed** full-data subgroup, i.e. the definitionally matched comparator. All four read **1.000**, and `find_metrics[["Exact"]]` reads **1.000** as well — but that analysis runs a *frozen* candidate family, so every leave-one-out fold re-selects the same rule by construction, and the payload's MR block predates both `return_reselection` and `ci_method = "field"`, so no `p̂` sits beside it. It anchors the vocabulary; it does not quantify a gap.

Both are reported in Part V2, with what each can and cannot say.

---

## Part R — The construction

### R1. The hoist

`R/fs_mr_inference.R:738–748` after the change:

```r
  # HOIST (TASK_field_recovery_2026-09-09, R1).  These three were built inside
  # the include_complement branch below; the field-recovery diagnostics need
  # `kept` and `Nall` whether or not that branch runs.  The hoist is
  # byte-identical: all three are pure functions of objects already fixed above
  # (`candidates` and `asm` at the top of the function, `df` unchanged since),
  # with no fitting, no drawing, no RNG consumption and no side effects --
  # a list subset, a length and a row count.  Nothing between here and the
  # branch reads or writes them.
  kept   <- candidates[asm$keep]              # aligns with asm columns
  Ncol   <- length(asm$names)
  Nall   <- nrow(df)
  if (isTRUE(include_complement)) {
    winset <- sort(unique(c(winner[!is.na(winner)], sel)))
```

The branch is otherwise untouched. `Ncol` and `Nall` were hoisted with `kept` because they are equally pure and both are read inside the branch immediately after.

The winner record needed two further lines. `R/fs_mr_inference.R:893–900`:

```r
    fc <- isTRUE(field_complement) && isTRUE(include_complement)
    # The recovery diagnostics (TASK_field_recovery_2026-09-09) read the same
    # outer-winner record, so G_out is allocated when EITHER consumer wants it.
    # W_in is the complement's alone.  With field_recovery = FALSE, rec_g is
    # fc and both allocations are exactly as before.
    rec_on <- isTRUE(field_recovery)
    rec_g  <- fc || rec_on
    G_out <- if (rec_g) rep(NA_integer_, field_R_out) else NULL
    W_in  <- if (fc) matrix(NA_integer_, field_R_out, field_R_in) else NULL
```

and inside both loop arms the single compound assignment splits into two guarded ones — `if (rec_g) G_out[r] <- as.integer(G)` and `if (fc) W_in[r, ] <- ...`. With `field_recovery = FALSE`, `rec_g` is `fc`, so both statements fire together exactly as the compound one did.

### R2. The new argument

`field_recovery = FALSE` is appended **last** in the signature of `fs_mr_inference()`, after `ij_residual`, so no positional call anywhere changes meaning. `R/forestsearch_main.R:3445–3448`:

```r
        # Add-only pass-through (TASK_field_recovery_2026-09-09): the field's
        # membership-agreement diagnostics; FALSE is the default and the block
        # is default-inert (nothing reads it, and it draws nothing).
        field_recovery = .g_mr(mr_inference_args$field_recovery, FALSE),
```

The roxygen `@param` states, in one paragraph, that these are descriptive diagnostics computed from draws already made — no new fits, no new randomness, no RNG consumption, no construction reading them — and that they answer a narrower question than FB/CV: re-selection within the fixed kept family under perturbation, not re-discovery from scratch, "related but not interchangeable, and a report that quotes one should say which."

### R3. The metrics

Attached as `field$recovery` when `field_recovery = TRUE`, absent when `FALSE`, by `.fs_mr_field_recovery()` (`R/fs_mr_inference.R:1094`). Called once, last, after every construction is complete, on **both** branches of the `length(ok_f) >= 2L` split, so a fit with fewer than two usable outer draws returns the all-`NA` shape rather than a missing element:

```r
    if (rec_on)
      field$recovery <- .fs_mr_field_recovery(kept, G_out, sel, Nall)
```

Returned elements: `sens_H`, `ppv_H`, `sens_Hc`, `ppv_Hc`, `npv_Hc`, `q10`, `q50`, `q90`, `share_equal_1`, `n_draws`, `n_used`, `n_skipped`, `n_selected`, `n_all`, `timing_seconds`. Guards: a draw is used when it recorded a winner that indexes the kept family and that candidate is non-empty; the rest count into `n_skipped`, and `n_used + n_skipped = n_draws` by construction. With no usable draw every metric is `NA_real_`. Sizes and intersections are computed once per **distinct** re-selected candidate and mapped back to draws.

### R4. Recorder

Template gains `FS_S7_FIELD_RECOV` (default `FALSE`, default-inert) beside the other MR knobs, forwards `field_recovery = mr_field_recovery` in `mr_inference_args`, echoes it in the "Template knobs" line and the MR-settings table, records it in `meta` (batch and combine), and carries the nine columns beside the `fld_Hc_scale_*` block: `fld_recov_sens_H`, `fld_recov_ppv_H`, `fld_recov_sens_Hc`, `fld_recov_npv_Hc`, `fld_recov_q10`, `fld_recov_q50`, `fld_recov_q90`, `fld_recov_share1`, `fld_recov_n_used`, each filled `%||% NA_real_`.

### R5. `NEWS.md`

One bullet under the development header.

---

## Gate R

Config: effMaxSG ε 0.20, HR 1.50, n = 500, `z1q` 0.60, `er_jcuts` 10, `seed_base = 8316951` (per-replicate seed `seed_base + sim_id`), `sim_id` 1–5, `FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FIELD_SCALEC=selected FS_S7_FB=none`, `mr_draws = 5000`, tags `recov_off` / `recov_on`. Timing columns (`fb_secs`, `fit_mr_secs`, `fld_H_secs`, `fld_Hc_secs`, `fld_H_uniform_secs`) excluded from every comparison.

**Driver.** Both arms rendered the **committed template itself**, through a byte-identical working copy (`diff` empty) whose only purpose was to keep the rendered HTML off the template's own output name. No gate-local driver is committed, because there is nothing in it to record: the knob is in the template.

### R2a — the hoist is inert

```
new columns: 9 fld_recov_sens_H, fld_recov_ppv_H, fld_recov_sens_Hc, fld_recov_npv_Hc,
             fld_recov_q10, fld_recov_q50, fld_recov_q90, fld_recov_share1, fld_recov_n_used
ncol off/on/e1stud: 167 167 158

pre-existing non-timing columns compared: 153
identical(): 153 / 153
truth identical(): TRUE
new columns present and all NA: 9 / 9
R2a: PASS
```

### R2b — on

```
pre-existing non-timing columns compared: 153
identical(): 153 / 153
truth identical(): TRUE
detected rows: 5 ; new columns finite on them: 9 / 9
invariant 0<=metrics<=1: TRUE
invariant q10<=q50<=q90: TRUE
invariant n_used <= field n_out_used: TRUE
  (n_used: 1000,1000,1000,1000,1000 ; fld_H_nout: 1000,1000,1000,1000,1000)
R2b: PASS
```

`n_skipped` is 0 on every row, so `n_used + n_skipped = n_draws = 1000 = |ok_f|` throughout.

### R2c — arithmetic check, one replicate, by hand

The gate's run loop is a `%dofuture%` with `seed = TRUE`, so each replicate draws from its own L'Ecuyer substream; the loop was reproduced verbatim under `plan("sequential")` in a scratch script so `trace()` could see into the iteration, and the capture taken on `s == 1` only. `sens_H` and `ppv_H` were then recomputed with a plain `for` loop over the draws using `intersect()` and `length()`, **never calling `.fs_mr_field_recovery()`**:

```r
Hhat <- kept[[sel]]
for (r in which(!is.na(G_out))) {
  Gr <- kept[[ G_out[r] ]]
  a  <- length(intersect(Gr, Hhat))
  sens_terms <- c(sens_terms, a / length(Hhat))
  ppv_terms  <- c(ppv_terms,  a / length(Gr))
}
```

```
scratch sim_id 1 : label 'q8.1 & q27.0'  |Hhat| 161
gate    sim_id 1 : label 'q8.1 & q27.0'  |Hhat| 161
same replicate as the gate row: TRUE

|Hhat| = 161, n = 500, kept family = 1233, draws with a winner = 1000
sens_H  by hand              : 0.31784472049689438
sens_H  recorded (gate row 1): 0.31784472049689444
        |diff|               : 5.55e-17
ppv_H   by hand              : 0.55241883811445247
ppv_H   recorded (gate row 1): 0.55241883811445236
        |diff|               : 1.11e-16
R2c: PASS
```

Both differences are one floating-point ulp — well inside the 1e-12 tolerance, and they are the expected residue of summing 1,000 terms in a different order.

### GATE R: PASS (R2a PASS, R2b PASS, R2c PASS)

### Cost

```
.fs_mr_field_recovery() over 1000 outer draws (family 1233, |Hhat| 161):
  median 0.0010 s, max 0.0020 s over 20 calls
package-reported field$recovery$timing_seconds on that replicate: 0.0010 s
```

Against the standing per-replicate reference — `fit_mr_secs` 32.1–40.3 s on these five rows — that is about **0.003%**. The measured `fit_mr_secs` are indistinguishable between the arms (mean 36.15 s off, 36.15 s on; per-row differences ±0.3 s, i.e. load noise, in both directions). No campaign is warranted by the cost.

---

## Part V2 — Validation and reporting (R-5)

### The field's numbers beside what the cell actually carries

Same five replicates as Gate R, from the `recov_on` bundle and the committed `e1stud` cell:

| sim_id | label | \|Ĥ\| | p̂ | field `sens_H` | field `ppv_H` | containment q50 | share = 1 | truth `sens` | truth `ppv` |
|---|---|---|---|---|---|---|---|---|---|
| 1 | q8.1 & q27.0 | 161 | 0.013 | 0.318 | 0.552 | 0.267 | 0.038 | 1.000 | 0.957 |
| 2 | q10.0 & q26.0 | 69 | 0.186 | 0.307 | 0.236 | 0.159 | 0.155 | 0.000 | 0.000 |
| 3 | q24.0 & q27.0 | 113 | 0.064 | 0.309 | 0.378 | 0.257 | 0.057 | 0.531 | 0.761 |
| 4 | q7.1 & q12.1 | 180 | 0.037 | 0.406 | 0.633 | 0.306 | 0.004 | 0.895 | 0.850 |
| 5 | q21.0 & q28.0 | 87 | 0.159 | 0.567 | 0.595 | 0.460 | 0.187 | 0.445 | 0.793 |
| **mean** | | 122.0 | **0.092** | **0.381** | 0.479 | 0.290 | 0.088 | **0.574** | 0.672 |

The wider committed cell (`e1stud`, 2,000 replicates, 1,999 detected), for context:

```
p_hat_H     : mean 0.103 | q10 0.014 | median 0.069 | q90 0.241
truth sens  : mean 0.590 | q10 0.243 | median 0.574 | q90 0.993
truth ppv   : mean 0.706 | q10 0.363 | median 0.734 | q90 1.000
truth spec  : mean 0.893 | q10 0.771 | median 0.906 | q90 1.000
truth npv   : mean 0.837 | q10 0.708 | median 0.831 | q90 0.997
share of detected replicates with p_hat_H < 0.20: 0.850
```

The five gate rows sit squarely inside that cell (`p̂` mean 0.092 against the cell's 0.103, truth `sens` 0.574 against 0.590), so they are not an unrepresentative corner.

### The gap, quantified — and what it is and is not

- **`sens_H` against `p̂`.** Mean `sens_H` 0.381 against mean `p̂` 0.092: **a factor of 4.1**. Every row has `sens_H` above its `p̂`, and the ordering is not preserved — the correlation across the five rows is 0.33. Row 1 makes the point sharpest: `p̂` = 0.013, so *exact* re-selection essentially never happens, yet a typical draw retains 32% of `Ĥ` and 3.8% of draws contain all of it. Row 2 is the mirror image: `p̂` = 0.186 is the *highest* in the set, but `ppv_H` = 0.236 is the lowest, because `|Ĥ|` = 69 is small and the draws that miss it land far away. **`p̂` and `sens_H` are not two readings of the same quantity, and neither can be recovered from the other.**

- **`sens_H` against the cell's truth-referenced `sens`.** These are the two metrics available on the same replicates, and the correlation across the five rows is **0.02** — no relationship. That is the expected answer, not a disappointment: the comparators differ. The field's `sens_H` scores each *re-selection* against the **observed** `Ĥ`; the cell's `sens` scores `Ĥ` against the **true planted subgroup**. Row 2 shows why they must be allowed to disagree — `Ĥ` misses the truth entirely (truth `sens` = 0.000) while the field's own re-selections still agree with `Ĥ` at 0.307. **Reproducibility of a pick is not correctness of a pick.** Nothing in this task claims otherwise, and `sens_H` must never be read as evidence that the identified subgroup is the right one.

- **`sens_H` against a genuine re-discovery rate.** The matching-vocabulary anchor is the committed GBSG frozen-family payload, where CV's own `sens_H` / `sens_Hc` / `ppv_H` / `ppv_Hc` against the **observed** subgroup all read **1.000** and `find_metrics[["Exact"]]` reads **1.000**. That is a frozen candidate family, where leave-one-out re-selects the same rule by construction, so it bounds the comparison rather than calibrating it — and no `p̂` sits beside it, since that payload's MR block predates both `return_reselection` and `ci_method = "field"`. **The honest statement is that a same-replicates FB/CV-versus-field comparison does not exist in the committed record, and this task did not create one** (that would be a campaign, which is out of scope). What can be said from the definitions, and is said in the vignette and the roxygen: FB and CV **re-run the search** on resampled data and rebuild the family; the field **re-selects within the fixed kept family** under multiplier perturbation. The field's version answers the narrower question and is free in every analysis; FB and CV answer the stronger one at a full re-run's cost. **Neither substitutes for the other.**

### `summary.forestsearch()`

One line added to the post-selection block, inside the `long` branch only (so `print()` is untouched, as directed), guarded on presence. On the vignette's GBSG fit:

```
Post-selection inference (certified products):
  Re-selection frequency  p-hat(H) = 0.006
    top-3 re-selection mass:  q1.1 & q18.1 0.182 | q1.1 & q27.1 0.156 | q10.0 & q30.0 0.056
    p_hat_sum = 1.000 over a family of 1744 candidates
    membership recovery: sens_H = 0.558 (mean share of Hhat retained over 998 draws)
```

**Absent-MR invariance re-run (Gate Pa from the print/vignette task): PASS.**

```
object: forestsearch(mr_inference = FALSE) on GBSG; is.null(mr_inference) = TRUE
print()  : before 16 lines, after 16; identical: TRUE
summary(): before 41 lines, after 41; identical: TRUE
differing lines: print 0, summary 0
```

(`before` captured from the installed package built at the Part R commit — i.e. without the `summary()` change — and `after` from the working tree, same fit, same seed.) Line counts also match the counts recorded in `REPORT_print_vignette_2026-09-09.md`. `.fs_mr_products()` returns `NULL` when `mr_inference` is absent, so the new branch is unreachable on that path.

### The vignette

`vignettes/survival-post-selection.qmd` now passes `field_recovery = TRUE` in the fit's `mr_inference_args`, and its p̂ section gains **"What p-hat alone cannot say"**, with the numbers computed live:

```
sens_H  = 0.558   mean share of Hhat retained by a re-selection
ppv_H   = 0.563   mean share of a re-selection that lies in Hhat
containment q10 / q50 / q90 = 0.000 / 0.547 / 1.000
share of draws containing all of Hhat = 0.243  (over 998 draws)
```

So the worked example that motivated the whole proposal now reads: **p̂ = 0.006 over a family of 1,744, but `sens_H` = 0.558 and a quarter of the draws contain all of `Ĥ`.** The section says that the pick is unstable *in its exact boundary*, not in the region it points at; that the low-p̂ caveat on the bounds still stands, because the bias correction is driven by the exact re-selection map; and that the field's family-conditional agreement and the FB/CV re-discovery rates answer different questions.

**Rebuild: 34.4 s elapsed, 518 MB peak RSS** (`/usr/bin/time -v` on `quarto render`), against 34.55 s / 488 MB recorded for the previous build — unchanged within noise, and the new chunk reads three elements off an object the fit already carries.

---

## Verification tallies

- **Full test suite:** `FAIL 0 | WARN 32 | SKIP 3 | PASS 5051` — identical to the tally recorded at `105186cd`. No test was added, edited, or skipped for this change.
- **`R CMD check --as-cran`** (`rcmdcheck::rcmdcheck(args = "--as-cran")`, the certification surface, which builds the PDF manual): **0 errors | 1 warning | 2 notes** — **unchanged** from `REPORT_print_vignette_2026-09-09.md`. All three are pre-existing and none is in a file this task touched:
  - WARNING `checking code files for non-ASCII characters` — `R/fs_bias_coverage.R` (the same file as side issue 1);
  - NOTE `checking R code for possible problems` — `fs_plot_bias_coverage`'s nine NSE bindings (`b`, `cell`, `cov`, `cov1`, `cov2`, `estimator`, `obs`, `r`, `ref`);
  - NOTE `checking HTML version of manual` — `no command 'tidy' found`, an environment note.

  The vignette re-build step inside the check reports `[53s/56s] OK` (against `[54s/57s]` last time), so the new chunk costs nothing measurable there either.
- `tools::showNonASCIIfile()` reports 0 non-ASCII characters in both touched `R/` files.

---

## Files touched

| file | change |
|---|---|
| `R/fs_mr_inference.R` | R1 hoist; `field_recovery` formal, roxygen `@param` and `@return`; `rec_on`/`rec_g` winner record; the `field$recovery` attachment; `.fs_mr_field_recovery()` |
| `R/forestsearch_main.R` | one `.g_mr` pass-through line (plus its three-line comment) |
| `R/forestsearch_methods.R` | `recov_sens_H`/`recov_n_used` in `.fs_mr_products()`; one guarded `summary()` line; roxygen paragraph |
| `man/fs_mr_inference.Rd`, `man/summary.forestsearch.Rd` | regenerated by `devtools::document()`; NAMESPACE unchanged |
| `NEWS.md` | one bullet |
| `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` | `FS_S7_FIELD_RECOV`, forwarding, knob echo, MR-settings row, batch and combine `meta`, nine `fld_recov_*` columns and their fill |
| `vignettes/survival-post-selection.qmd` | `field_recovery = TRUE` in the fit; the "What p-hat alone cannot say" section |
| `quarto/simulations/gbsg_020/results/…_recov_off_res_1_5.rds`, `…_recov_on_res_1_5.rds`, `…_recov_on_batch_1_5.html` | the two gate bundles and the on-arm render |

Out of scope and untouched, as directed: the rule-name family (R-2), `.fs_apply_mr()`'s `ci_method` fallback, `print.forestsearch()`, every construction, bound and default, and the seven pre-existing untracked files.

---

## Side issues (not fixed here)

1. **`devtools::document()` emits a pre-existing roxygen error** on `R/fs_bias_coverage.R:17` — `@description failed to evaluate inline markdown code … Failed to parse the inline R code: 'r = se_mean / sd_emp'`. The backtick-`r` in that description is being read as an inline code chunk. It predates this task (`2f118042`), does not block `document()` (every other `.Rd` is written, including both regenerated here), and is not touched.

2. **`.fs_apply_mr()` still has no `field_recovery` fallback.** Same asymmetry already recorded for `ci_method` in `REPORT_cimethod_flip_2026-09-09.md` §7: the DINA and GRF hooks fall back independently in `R/fs_mr_inference_methods.R`. The consistency engine forwards `field_recovery`; the DINA/GRF hooks do not, so the diagnostics are unreachable from those branches. Out of scope by the task's own wording; flagged for a later decision.

3. **The recorder writes `fld_recov_n_used` as `NA_real_`-initialised**, matching the task's column list, so an integer count is stored in a double column. Harmless, and consistent with `fld_Hc_nin_mean` beside it; noted only because `fld_Hc_nout` next to it is `NA_integer_`.
