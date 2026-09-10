# guohe_sec52_adaptive_run.R
#
# T2 driver for Tier 2(a) of the fs-post-selection supplement: the Guo & He
# ADAPTIVE-r column on t7 (dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md
# section 5, as amended by ..._v2.md A5).
#
# AUTHORED ONLY -- NOT EXECUTED ON THE MAC. Phase B runs this on the Linux box
# after Larry declares the sync complete and gives the Phase-B pointer line
# (Gate 1b).
#
# TRANSPLANT: the pilot/flag/skip scaffolding of guohe_sec52_run.R (a9e099df,
# lines 63-84, 95-138, 141-145, 208-225) plus the --adaptive named lines of
# guohe_reproduction_run.R (8a98e05d, lines 54-55, 67, 109, 127-128) and its
# adaptive call in guohe_reproduction_sim.R:184-188. Named lines only.
#
# Per replicate m of cell id:
#   1. Regenerate the data from the STORED seed (base + m, base read from the
#      committed bundle and asserted against the guohe_sec52_run.R formula),
#      recompute the naive argmax and the naive one-sided lower bound, and
#      assert identical() to the stored guohe_repro_t7_<id>.rds row -- the
#      pairing proof. Any mismatch is a STOP.
#   2. Run guohe_adaptive_r() with orient = +1, r_grid = GH52_R_GRID
#      (= c(1/3, 1/12, 1/21, 1/30), the published Table-7 grid, pinned), v = 5,
#      B = 2000, under a derived seed (base + m + GHA_SEED_OFFSET) that is
#      recorded. Record r-hat and the per-candidate objective values.
#   3. PRIMARY Adaptive bound = the STORED B = 2000 Algorithm-3 bound at r-hat,
#      looked up by (id, m) from guohe_repro_t7_<id>.rds -- exact resolution
#      parity with the fixed-r columns. The bound itself is not stored, but
#      r{i}_dist = gamma_s - lower, so lower = gamma_s - r{i}_dist recovers it
#      exactly; this is the same reconstruction mr_field_vs_guohe_run.R:70
#      already performs.
#   4. SECONDARY = the adaptive function's OWN final refit bound. Under A5 both
#      are at B = 2000, so the primary/secondary distinction is no longer one of
#      resolution: it is the stored repro run's bootstrap stream (boot_seed =
#      seed + 500000L) against this run's own (GHA_SEED_OFFSET). Both are
#      recorded and both are scored.
#   5. Coverage of both against gamma_s from the committed truth caches.
#
# Output guohe_adaptive_t7_<id>.rds per cell, bundle format as guohe_sec52_run.R.
#
# ---- A5: B = 2000, superseding v1's B = 200 --------------------------------
# v1 pinned B = 200 as "the validated reproduction setting". That premise is
# false: --adaptive-B is inert in guohe_reproduction_run.R (b_adapt is parsed at
# :54, printed at :67 and stored as metadata at :127 but never reaches
# gh_one_rep; :109 passes B = b_boot), so the validated Tables 3-6 Adaptive
# columns executed at 2000 -- see dev/notes/NOTE_adaptive_B_inert_2026-09-09.md.
# B = 2000 is therefore simultaneously the setting under which the function was
# validated, resolution parity with the fixed-r columns, and their method's
# strongest configuration (the honest-comparison requirement).
#
# ONE B SERVES BOTH the inner CV fits (R/guohe_adaptive_r.R:255) and the final
# refit (:285). That coupling is the function's own behavior and is NOT altered
# here: splitting it would be an implementation change to their algorithm and is
# out of scope (A5). Their method is run, never revised.
#
# ---- Gate 1b (Phase B, Linux) ----------------------------------------------
# Precondition: Larry states the Linux->Mac->push sequencing is complete and
# gives the Phase-B pointer line. Then devtools::install(); pilot at reps = 20
# on t7_beta2_00; print the 6 x 2000 projection.
#
# The cost must be MEASURED, not inherited (A5): guohe_reproduction_RUN.md:79
# records 5.90 s per replicate labelled "B = 200", but with the flag inert that
# figure may have been measured at either setting, so the full-grid projection
# spans ~1,400-1,700 core-h (13-16 h wall) at best to roughly ten times that at
# worst. Envelope: proceed only if <= 2,500 core-h AND <= 24 h wall at the
# available cores. If exceeded: STOP and report the measured projection together
# with costed reduced options (500 reps/cell; three cells; B = 200 as a labelled
# sensitivity) -- do not launch, and do not choose among them.
#
# Two caveats to carry into the T2 record:
#   (i)  Guo & He's own Table 6 adaptive caution.
#   (ii) guohe_adaptive_r() uses INDEPENDENT draws across r (no common random
#        numbers across the grid), which inflates Var(r-hat) -- and the Adaptive
#        column is precisely a measurement of r-hat
#        (guohe_reproduction_RUN.md:131-136).
#
# Usage (Phase B, Linux):
#   Rscript quarto/GuoHe/guohe_sec52_adaptive_run.R --pilot --cores=120
#   Rscript quarto/GuoHe/guohe_sec52_adaptive_run.R --cores=120
#   Rscript quarto/GuoHe/guohe_sec52_adaptive_run.R --cells=t7_beta2_00 --force

suppressMessages(library(survival))
suppressMessages(library(forestsearch))

.gha_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
})
source(file.path(.gha_dir, "guohe_sec52_sim.R"))

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# Derived seed offset for the adaptive call. Distinct from gh52_one_rep()'s
# boot_seed (+500000L) so the adaptive run's stream cannot collide with the
# stored repro run's Algorithm-3 stream.
GHA_SEED_OFFSET <- 600000L
GHA_V <- 5L          # their tables' value (Q4: guohe_reproduction_run.R:128)
GHA_ORIENT <- +1     # t7 orientation (stored bundles record orient = 1)

# ---- arguments (transplant: guohe_sec52_run.R:63-84) -----------------------
args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(args == paste0("--", nm))
opt <- function(nm, default) {
  hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
}

pilot <- flag("pilot")
n_rep_opt <- as.integer(opt("reps", NA_character_))
b_boot <- as.integer(opt("B", "2000"))          # A5: 2000, not 200
force <- flag("force")
out_dir <- opt("out", .gha_dir)
truth_dir <- opt("truth-dir", out_dir)
cells_opt <- opt("cells", "")
# Gate 1b measurement flags (2026-09-09). BOTH default to the committed
# behaviour: the full published grid and the CV path. Nothing is changed
# permanently -- omit them and the driver runs exactly as committed.
.parse_r <- function(txt) vapply(strsplit(txt, ",")[[1]], function(z) {
  z <- trimws(z)
  if (grepl("/", z)) { ab <- as.numeric(strsplit(z, "/")[[1]]); ab[1] / ab[2] }
  else as.numeric(z)
}, numeric(1), USE.NAMES = FALSE)
rgrid_opt <- opt("r-grid", "")
GHA_RUN_GRID <- if (nzchar(rgrid_opt)) .parse_r(rgrid_opt) else GH52_R_GRID
fixed_r_opt <- opt("fixed-r", "")
GHA_FIXED_R <- if (nzchar(fixed_r_opt)) .parse_r(fixed_r_opt)[1] else NA_real_
n_cores <- as.integer(opt(
  "cores",
  as.character(max(1L, floor(0.80 * parallel::detectCores(logical = FALSE))))
))

GHA_CELLS <- sprintf("t7_beta2_%02d", 0:5)
if (nzchar(cells_opt)) GHA_CELLS <- strsplit(cells_opt, ",")[[1]]
if (pilot) {
  GHA_CELLS <- "t7_beta2_00"                    # the maximal-bias cell
  n_rep_opt <- as.integer(opt("reps", "20"))
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Guo & He Section 5.2 / Table 7 -- ADAPTIVE-r column (T2)\n")
cat(sprintf("  mode      : %s\n", if (pilot) "PILOT (projection only)" else "production"))
cat(sprintf("  cells     : %s\n", paste(GHA_CELLS, collapse = ", ")))
cat(sprintf("  B         : %d   (A5: one B serves inner CV and final refit)\n", b_boot))
cat(sprintf("  r_grid    : %s%s\n", paste(sprintf("%.6f", GHA_RUN_GRID), collapse = ", "),
            if (!identical(GHA_RUN_GRID, GH52_R_GRID)) "   [--r-grid override]" else ""))
if (is.finite(GHA_FIXED_R))
  cat(sprintf("  MODE      : FIXED r = %.6f (no CV; Algorithm 3 only)\n", GHA_FIXED_R))
cat(sprintf("  v         : %d   orient : %+d\n", GHA_V, GHA_ORIENT))
cat(sprintf("  cores     : %d\n", n_cores))
cat(sprintf("  out       : %s\n", out_dir))
cat(sprintf("  truth dir : %s\n\n", truth_dir))

# ---- truth caches: refuse to proceed without them (guohe_sec52_run.R:95-138)
truth_file <- function(b2) {
  file.path(truth_dir, sprintf("guohe_sec52_truth_beta2_%02d.rds", round(b2 * 10)))
}
gha_truth <- function(id) {
  b2 <- as.integer(sub("^t7_beta2_", "", id)) / 10
  tf <- truth_file(b2)
  if (!file.exists(tf))
    stop("Truth cache missing for ", id, " (", tf, ").\nRun production step 0 ",
         "first:  Rscript quarto/GuoHe/guohe_sec52_truth.R\nA coverage run ",
         "must never start without gate-cleared truth curves.", call. = FALSE)
  tr <- readRDS(tf)
  if (!inherits(tr, "gh52_truth") || !all(tr$gates$pass))
    stop("Truth cache for ", id, " is not a gate-cleared gh52_truth object; ",
         "recompute with guohe_sec52_truth.R.", call. = FALSE)
  # exact-null scoring basis: validate, never mutate (guohe_sec52_run.R:118-137)
  if (tr$beta2 == 0 && is.null(tr$beta_exact))
    stop("The beta2 = 0 truth cache predates the exact-null scoring basis ",
         "(no `beta_exact`); regenerate the caches:\n",
         "  Rscript quarto/GuoHe/guohe_sec52_truth.R --force", call. = FALSE)
  tr
}

# ---- one replicate ---------------------------------------------------------
gha_one_rep <- function(id, m, row_r, truth, base, B) {
  b2 <- as.integer(sub("^t7_beta2_", "", id)) / 10
  seed <- base + m
  # (1) regenerate from the stored seed, exactly as gh52_one_rep() does
  set.seed(seed)
  df0 <- gh52_sim_data(b2, 400L)
  cand <- gh52_candidates(df0)
  fits <- gh52_subgroup_fits(cand$df, cand)
  nv <- gh52_naive(fits, cand, truth, level = 0.05)

  # PAIRING PROOF: the naive argmax and naive bound must reproduce the stored
  # row exactly. Any mismatch is a STOP (v1 section 5, A8).
  pair_ok <- identical(nv$c_hat, row_r$c_hat_naive) &&
    identical(nv$point, row_r$naive_point) &&
    identical(nv$lower, row_r$naive_lower) &&
    identical(nv$cover, row_r$naive_cover) &&
    identical(nv$dist, row_r$naive_dist) &&
    identical(nv$bias, row_r$naive_bias) &&
    identical(nv$gamma_s, row_r$gamma_s_naive) &&
    identical(mean(df0$event == 0L), row_r$cens_rate) &&
    identical(length(cand$cuts), as.integer(row_r$n_cand))

  gamma_s <- row_r$gamma_s          # truth at the ENGINE's selected cutpoint
  # (2) the adaptive fit, under its own recorded derived seed
  seed_ad <- seed + GHA_SEED_OFFSET
  t0 <- proc.time()[["elapsed"]]
  if (is.finite(GHA_FIXED_R)) {
    # (a) FIXED r -- Algorithm 3 at one r, no cross-validation. A single r
    # leaves the CV objective nothing to select over, so the CV path is skipped
    # entirely rather than run degenerately.
    .fit <- try(suppressWarnings(guohe_algorithm3(
      data = cand$df, outcome = "survival", treatment = "treat",
      candidates = cand$names, time = "time", event = "event",
      orient = GHA_ORIENT, B = B, r = GHA_FIXED_R, level = 0.05,
      seed = seed_ad, min_events = 5L, diagnostics = FALSE
    )), silent = TRUE)
    ar <- if (inherits(.fit, "try-error")) .fit else
      list(r_hat = GHA_FIXED_R, r_grid = GHA_FIXED_R,
           objective = NA_real_, per_candidate = NA_real_, fit = .fit)
  } else {
    ar <- try(suppressWarnings(guohe_adaptive_r(
      data = cand$df, outcome = "survival", treatment = "treat",
      candidates = cand$names, time = "time", event = "event",
      orient = GHA_ORIENT, r_grid = GHA_RUN_GRID, v = GHA_V, B = B,
      level = 0.05, seed = seed_ad, min_events = 5L, refit = TRUE
    )), silent = TRUE)
  }
  t_ad <- proc.time()[["elapsed"]] - t0

  out <- data.frame(
    id = id, m = m, beta2 = b2, n = 400L,
    seed_data = seed, seed_adaptive = seed_ad,
    n_cand = length(cand$cuts), cens_rate = mean(df0$event == 0L),
    pair_ok = as.integer(pair_ok),
    c_hat_gh = row_r$c_hat_gh, c_hat_naive = nv$c_hat,
    n_sel = row_r$n_sel, gamma_s = gamma_s,
    naive_point = nv$point, naive_lower = nv$lower, naive_cover = nv$cover,
    stringsAsFactors = FALSE)

  if (inherits(ar, "try-error")) {
    out$ad_error <- 1L
    out$r_hat <- NA_real_; out$r_hat_index <- NA_integer_
    out$ad_sel_ok <- NA_integer_
    out$ad_lower_primary <- NA_real_; out$ad_cover_primary <- NA_integer_
    out$ad_dist_primary <- NA_real_
    out$ad_lower_secondary <- NA_real_; out$ad_cover_secondary <- NA_integer_
    out$ad_dist_secondary <- NA_real_; out$ad_bias_secondary <- NA_real_
    for (i in seq_along(GHA_RUN_GRID)) out[[sprintf("obj_r%d", i)]] <- NA_real_
    out$obj_min <- NA_real_
    out$t_adaptive_sec <- t_ad
    out$ad_note <- as.character(ar)
    return(out)
  }

  out$ad_error <- 0L
  r_hat <- ar$r_hat
  # index into the stored fixed-r columns; exact match on the pinned grid
  ri <- which(vapply(GH52_R_GRID, function(z) isTRUE(all.equal(z, r_hat)), logical(1)))
  ri <- if (length(ri)) ri[1] else NA_integer_
  out$r_hat <- r_hat
  out$r_hat_index <- ri

  # the adaptive refit must land on the same subgroup the stored run selected
  sel_i <- match(ar$fit$selected, cand$names)
  out$ad_sel_ok <- as.integer(!is.na(sel_i) &&
                                isTRUE(cand$cuts[sel_i] == row_r$c_hat_gh))

  # (3) PRIMARY: the stored B = 2000 Algorithm-3 bound at r-hat, recovered as
  #     lower = gamma_s - r{ri}_dist (mr_field_vs_guohe_run.R:70 pattern)
  if (!is.na(ri)) {
    dist_stored <- row_r[[sprintf("r%d_dist", ri)]]
    lower_p <- gamma_s - dist_stored
    out$ad_lower_primary <- lower_p
    out$ad_dist_primary <- dist_stored
    out$ad_cover_primary <- as.integer(lower_p <= gamma_s)
    out$ad_bias_primary <- row_r[[sprintf("r%d_bias", ri)]]
    # cross-check: the stored per-r cover flag at r-hat must agree
    out$ad_cover_primary_stored <- as.integer(row_r[[sprintf("r%d_cover", ri)]])
  } else {
    out$ad_lower_primary <- NA_real_; out$ad_dist_primary <- NA_real_
    out$ad_cover_primary <- NA_integer_; out$ad_bias_primary <- NA_real_
    out$ad_cover_primary_stored <- NA_integer_
  }

  # (4) SECONDARY: the function's own final refit bound at B = 2000
  lower_s <- gh52_to_score(ar$fit$bound_one_sided)
  deb_s <- gh52_to_score(ar$fit$debiased)
  out$ad_lower_secondary <- lower_s
  out$ad_cover_secondary <- as.integer(lower_s <= gamma_s)
  out$ad_dist_secondary <- gamma_s - lower_s
  out$ad_bias_secondary <- deb_s - gamma_s

  # per-candidate / per-r objective values
  for (i in seq_along(GHA_RUN_GRID))
    out[[sprintf("obj_r%d", i)]] <- if (length(ar$objective) >= i)
      ar$objective[i] else NA_real_
  out$obj_min <- suppressWarnings(min(ar$objective, na.rm = TRUE))
  out$per_cand_min <- suppressWarnings(min(ar$per_candidate, na.rm = TRUE))
  out$per_cand_at_sel <- if (!is.na(sel_i) && length(ar$per_candidate) >= sel_i)
    unname(ar$per_candidate[sel_i]) else NA_real_
  out$t_adaptive_sec <- t_ad
  out$ad_note <- NA_character_
  out
}

# ---- one cell (transplant: guohe_sec52_run.R:141-207) ----------------------
gha_run_cell <- function(id) {
  f_out <- file.path(out_dir, paste0("guohe_adaptive_", id,
                                     if (is.finite(GHA_FIXED_R))
                                       sprintf("_fixedr%s", sub("[.]", "", sprintf("%.4f", GHA_FIXED_R)))
                                     else if (length(GHA_RUN_GRID) != length(GH52_R_GRID))
                                       sprintf("_grid%d", length(GHA_RUN_GRID)) else "",
                                     if (pilot) "_pilot" else "", ".rds"))
  if (file.exists(f_out) && !force) {
    cat(sprintf("[skip] %s (exists)\n", id))
    return(invisible(NULL))
  }
  f_rep <- file.path(.gha_dir, paste0("guohe_repro_", id, ".rds"))
  if (!file.exists(f_rep))
    stop("Stored reproduction bundle absent: ", f_rep, call. = FALSE)
  rep_bun <- readRDS(f_rep)
  truth <- gha_truth(id)

  # the base must be the STORED one, and must match guohe_sec52_run.R's formula
  base <- 1000000L + as.integer(sum(utf8ToInt(id)) * 100003L)
  stopifnot(identical(rep_bun$seed_base, base),
            identical(rep_bun$orient, +1),
            isTRUE(all.equal(rep_bun$r_grid, GH52_R_GRID)))
  # the stored fixed-r columns must be at the same B this run uses (A5 parity)
  if (!identical(as.integer(rep_bun$B), as.integer(b_boot)))
    warning("stored bundle B = ", rep_bun$B, " but this run uses B = ", b_boot,
            "; primary/secondary resolution parity is broken for ", id)

  n_rep <- nrow(rep_bun$results)
  if (!is.na(n_rep_opt)) n_rep <- min(n_rep, n_rep_opt)

  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("[run ] %s  %d reps, B=%d, v=%d ...\n", id, n_rep, b_boot, GHA_V))
  utils::flush.console()
  rows <- parallel::mclapply(
    seq_len(n_rep),
    function(m) try(gha_one_rep(id, m, rep_bun$results[m, ], truth, base, b_boot),
                    silent = TRUE),
    mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- vapply(rows, function(z) inherits(z, "try-error"), logical(1))
  if (any(bad))
    warning(sum(bad), " replicate(s) errored in ", id, "; dropped. First: ",
            as.character(rows[[which(bad)[1]]]))
  res <- do.call(rbind, rows[!bad])
  el <- proc.time()[["elapsed"]] - t0

  gate <- list(
    n_rep_expected = n_rep, n_rep_done = nrow(res), n_errored = sum(bad),
    pair_mismatch = sum(res$pair_ok == 0L),
    ad_errors = sum(res$ad_error, na.rm = TRUE),
    ad_sel_mismatch = sum(res$ad_sel_ok == 0L, na.rm = TRUE),
    r_hat_offgrid = sum(is.na(res$r_hat_index) & res$ad_error == 0L),
    primary_flag_disagree = sum(res$ad_cover_primary !=
                                  res$ad_cover_primary_stored, na.rm = TRUE))

  saveRDS(list(
    id = id, target = "Table 7 -- Adaptive column", beta2 = rep_bun$beta2,
    n = 400L, n_rep_requested = n_rep, n_rep_used = nrow(res),
    B = b_boot, r_grid = GHA_RUN_GRID, r_grid_stored = GH52_R_GRID,
    fixed_r = GHA_FIXED_R, mode = if (is.finite(GHA_FIXED_R)) "fixed-r" else "adaptive-CV",
    v = GHA_V, orient = GHA_ORIENT,
    one_B_serves_cv_and_refit = TRUE,
    adaptive_seed_offset = GHA_SEED_OFFSET,
    source_bundle = paste0("guohe_repro_", id, ".rds"),
    primary_bound = "stored Algorithm-3 bound at r_hat (gamma_s - r{i}_dist)",
    secondary_bound = "guohe_adaptive_r() own final refit bound",
    pilot = pilot, seed_base = base,
    truth_identity = list(n_big = truth$n_big, c_step = truth$c_step,
                          seed = truth$seed,
                          beta_exact = !is.null(truth$beta_exact),
                          scoring_basis = if (is.null(truth$beta_exact))
                            "smooth" else "exact"),
    gate = gate, elapsed_sec = el,
    sessionInfo = utils::capture.output(utils::sessionInfo()),
    results = res), f_out)
  cat(sprintf("[done] %s  %d/%d reps in %.1f min  pair_mm %d  ad_err %d  sel_mm %d  -> %s\n",
              id, gate$n_rep_done, n_rep, el / 60, gate$pair_mismatch,
              gate$ad_errors, gate$ad_sel_mismatch, basename(f_out)))
  utils::flush.console()
  invisible(res)
}

all_res <- list()
for (id in GHA_CELLS) all_res[[id]] <- gha_run_cell(id)

# ---- Gate 1b projection (transplant: guohe_sec52_run.R:214-225) ------------
if (pilot) {
  r <- all_res[[GHA_CELLS[1]]]
  if (!is.null(r)) {
    per_rep <- mean(r$t_adaptive_sec, na.rm = TRUE)
    total_core_h <- per_rep * 6 * 2000 / 3600
    wall_h <- total_core_h / n_cores
    cat("\n==== GATE 1b PILOT PROJECTION ====\n")
    cat(sprintf("  mean per-replicate (serial, adaptive only): %8.2f s  at B = %d\n",
                per_rep, b_boot))
    cat(sprintf("  full study 6 x 2000 reps, single core     : %8.1f core-h\n",
                total_core_h))
    cat(sprintf("  projected wall-clock at %3d cores         : %8.1f h\n",
                n_cores, wall_h))
    cat(sprintf("  envelope <= 2500 core-h AND <= 24 h wall  : %s\n",
                if (total_core_h <= 2500 && wall_h <= 24) "WITHIN" else "EXCEEDED"))
    cat(sprintf("  r_hat distribution in the pilot           : %s\n",
                paste(sprintf("%.4f", sort(unique(r$r_hat))), collapse = ", ")))
    cat("\nREPORT THIS PROJECTION TO THE OWNER BEFORE LAUNCHING PRODUCTION.\n")
    cat("If EXCEEDED: do NOT launch. Report the measured projection with costed\n")
    cat("reduced options (500 reps/cell; three cells; B = 200 as a labelled\n")
    cat("sensitivity) and let Larry choose.\n")
  }
} else {
  cat("\n=== T2 GATE TALLY ===\n")
  for (id in GHA_CELLS) {
    f <- file.path(out_dir, paste0("guohe_adaptive_", id, ".rds"))
    if (!file.exists(f)) { cat(sprintf("%-16s MISSING\n", id)); next }
    g <- readRDS(f)$gate
    cat(sprintf("%-16s reps %d/%d  errored %d  pair_mm %d  ad_err %d  sel_mm %d  offgrid %d  flag_dis %d\n",
                id, g$n_rep_done, g$n_rep_expected, g$n_errored,
                g$pair_mismatch, g$ad_errors, g$ad_sel_mismatch,
                g$r_hat_offgrid, g$primary_flag_disagree))
  }
}
