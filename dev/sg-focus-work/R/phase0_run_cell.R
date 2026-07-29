# =============================================================================
# Phase 0 -- search-stage re-run with selection-rule instrumentation
# =============================================================================
# READ-ONLY with respect to the package: nothing here modifies R/, and every
# quantity is read off the object subgroup.consistency() already returns
# (early_stop_triggered, early_stop_candidate, n_candidates_evaluated,
# n_candidates_total, n_passed, out_sg$result).  No instrumentation of the
# package source is required or performed.
#
# Usage (from the repository root):
#   Rscript dev/sg-focus-work/R/phase0_run_cell.R <cell> <sim_from> <sim_to> \
#           <workers> <out.rds> [stop_threshold]
#
# The optional 6th argument overrides stop_threshold.  Passing "none" disables
# early stopping, forcing every candidate to be evaluated.  That is NOT the
# as-run configuration -- it is a cross-check: it produces the COMPLETE
# qualifier set, against which the prefix-only qualifier set recorded by an
# early-stopping run can be validated (see PHASE0_FINDINGS.md).  Leave it off
# to reproduce the simulations.
#
# What is re-run: forestsearch() only.  The Tier-1 bootstrap is simply never
# called (the .qmd gates it on nb_boots >= 1) and debias_gate = FALSE, per the
# brief.  Neither touches the search stage, so the selected subgroup is
# unaffected -- which the seed-fidelity check against the saved sg_def verifies
# empirically rather than by assertion.

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(data.table)
  library(future)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) stop("usage: <cell> <sim_from> <sim_to> <workers> <out.rds>")
cell_name <- args[[1]]
sim_from  <- as.integer(args[[2]])
sim_to    <- as.integer(args[[3]])
n_workers <- as.integer(args[[4]])
out_path  <- args[[5]]
stop_thr_override <- if (length(args) >= 6L) {
  if (identical(args[[6]], "none")) NA else as.numeric(args[[6]])
} else NULL   # NULL = do not pass; inherit the forestsearch() default (0.95)

source("dev/sg-focus-work/R/phase0_cells.R")
K <- PHASE0_COMMON
cell <- PHASE0_CELLS[[cell_name]]
if (is.null(cell)) stop("unknown cell: ", cell_name)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# -----------------------------------------------------------------------------
# DGM -- built exactly as the .qmd build-dgm chunk does
# -----------------------------------------------------------------------------
k_inter <- calibrate_k_inter(target_hr_harm = cell$target_hr_harm,
                             model = K$dgm_model, use_ahr = FALSE)
dgm <- setup_gbsg_dgm(model = K$dgm_model, k_inter = k_inter,
                      n_super = K$n_super, seed = K$seed_base)
dgm <- compute_dgm_cde(dgm)

inner_parallel <- list(plan = "multisession", workers = n_workers)

# -----------------------------------------------------------------------------
# Per-replicate record
# -----------------------------------------------------------------------------
# Columns beyond the brief's four are kept because they cost nothing and make
# the branch classification auditable after the fact (n_qual, sel_m, min_m).

na_rec <- function(sim_id) data.frame(
  sim_id = sim_id, ok = 0L, detected = 0L,
  n_cand_total = NA_integer_, n_cand_evaluated = NA_integer_, n_passed = NA_integer_,
  early_stop_triggered = NA, early_stop_candidate = NA_integer_,
  stop_threshold = NA_real_,
  sel_Pcons = NA_real_, sel_hr = NA_real_, sel_N = NA_real_, sel_m = NA_integer_,
  n_qual = NA_integer_, max_Pcons_qual = NA_real_,
  first_qual_m = NA_integer_, first_qual_Pcons = NA_real_, first_qual_hr = NA_real_,
  hr_rank_selected = NA_integer_,
  n_qual_before_sel = NA_integer_,
  sg_def = NA_character_, n_sel = NA_integer_,
  secs = NA_real_, err = NA_character_,
  stringsAsFactors = FALSE)

record_replicate <- function(sim_id) {
  rec <- na_rec(sim_id)

  df <- tryCatch(
    simulate_from_dgm(dgm, n = cell$n_sample, analysis_time = K$analysis_time,
                      cens_adjust = K$cens_adjust, seed = K$seed_base + sim_id),
    error = function(e) NULL)
  if (is.null(df)) { rec$err <- "simulate_from_dgm failed"; return(rec) }
  df[[K$id_name]] <- seq_len(nrow(df))

  confs <- intersect(K$confounders_base, names(df))
  if (cell$k_random_noise > 0L) {
    noise_names <- paste0("noise", seq_len(cell$k_random_noise))
    set.seed(K$seed_base + sim_id + 1000000L)
    for (nm in noise_names) df[[nm]] <- rnorm(n = nrow(df))
    confs <- c(confs, noise_names)
  }

  fs_args <- list(
      df.analysis     = df,
      outcome.name    = K$outcome_name, event.name = K$event_name,
      treat.name      = K$treat_name,   id.name    = K$id_name,
      flag_harm.name  = K$harm_col,     confounders.name = confs,
      is.RCT          = TRUE,           seedit = K$seed_base + sim_id,
      quiet           = TRUE,
      sg_focus        = K$sg_focus,
      subgroup_method = K$subgroup_method,
      hr.threshold    = K$hr_threshold, hr.consistency = K$hr_consistency,
      pconsistency.threshold = K$pconsistency, n.min = K$n_min,
      selection_rule  = K$selection_rule,
      parallel_args   = inner_parallel,
      debias_gate     = FALSE,
      consistency_method = K$consistency_method,
      use_lasso = K$use_lasso, use_grf = K$use_grf,
      use_twostage = K$use_twostage, use_dina = K$use_dina,
      conf_force = K$fs_conf_force, conf.cont_jcuts = K$fs_conf.cont_jcuts,
      fs.splits = K$fs_splits, maxk = K$maxk,
      d0.min = K$d0_min, d1.min = K$d1_min)
  # `fs_args$x <- NULL` would DELETE the element; `fs_args["x"] <- list(NULL)`
  # is what actually passes stop_threshold = NULL through do.call().
  if (!is.null(stop_thr_override)) {
    if (is.na(stop_thr_override)) fs_args["stop_threshold"] <- list(NULL)
    else                          fs_args$stop_threshold <- stop_thr_override
  }

  t0 <- proc.time()[3]
  fs.est <- tryCatch(do.call(forestsearch, fs_args),
                     error = function(e) { rec$err <<- conditionMessage(e); NULL })
  rec$secs <- proc.time()[3] - t0
  if (is.null(fs.est)) return(rec)
  rec$ok <- 1L

  gc_ <- fs.est$grp.consistency
  if (!is.null(gc_)) {
    rec$n_cand_total         <- as.integer(gc_$n_candidates_total %||% NA_integer_)
    rec$n_cand_evaluated     <- as.integer(gc_$n_candidates_evaluated %||% NA_integer_)
    rec$n_passed             <- as.integer(gc_$n_passed %||% NA_integer_)
    rec$early_stop_triggered <- isTRUE(gc_$early_stop_triggered)
    rec$early_stop_candidate <- as.integer(gc_$early_stop_candidate %||% NA_integer_)
    rec$stop_threshold       <- as.numeric(gc_$stop_threshold %||% NA_real_)
  }

  found <- !is.null(fs.est$sg.harm)
  rec$detected <- as.integer(found)
  if (found) rec$sg_def <- paste(fs.est$sg.harm, collapse = " & ")
  if (!is.null(gc_$sg.harm.id))
    rec$n_sel <- sum(gc_$sg.harm.id == 1L, na.rm = TRUE)

  # --------------------------------------------------------------------------
  # Selection-rule instrumentation
  # --------------------------------------------------------------------------
  # out_sg$result is the post-consistency table AFTER sort_subgroups(), i.e.
  # setorder(-Pcons, -hr, K) for sg_focus = "hr".  Row 1 is therefore the
  # selected candidate.  Every row in it cleared pconsistency.threshold
  # (evaluate_subgroup_consistency() returns NULL otherwise), so the table IS
  # the qualifier set among the candidates that were actually evaluated.
  #
  # `m` is the candidate's GLOBAL index in the preview order, which for
  # sg_focus = "hr" is setorder(-HR, K) -- so m ascending IS effect-descending.
  # The intended rule ("argmax effect among Pcons >= threshold") therefore
  # picks the qualifier with the smallest m.  Candidates after an early stop
  # were never evaluated, but they all have lower HR than the stopping
  # candidate, so they cannot be the intended rule's winner whenever at least
  # one qualifier exists in the evaluated prefix -- which early stopping
  # guarantees.  The comparison below is thus exact, not an approximation.
  res <- gc_$out_sg$result
  if (!is.null(res) && nrow(res) > 0L) {
    res <- as.data.frame(res)
    Pc <- as.numeric(res$Pcons); hr <- as.numeric(res$hr)
    mm <- as.integer(as.numeric(res$m)); KK <- as.numeric(res$K)

    rec$sel_Pcons <- Pc[1]; rec$sel_hr <- hr[1]
    rec$sel_N <- as.numeric(res$N)[1]; rec$sel_m <- mm[1]

    rec$n_qual         <- length(mm)
    rec$max_Pcons_qual <- max(Pc, na.rm = TRUE)

    ord <- order(-hr, KK, mm)          # HR-descending, ties by K then index
    rec$hr_rank_selected <- which(ord == 1L)[1]
    j <- ord[1]
    rec$first_qual_m     <- mm[j]
    rec$first_qual_Pcons <- Pc[j]
    rec$first_qual_hr    <- hr[j]
    rec$n_qual_before_sel <- sum(mm < mm[1])
  }

  attr(rec, "qual_table") <- if (!is.null(res) && nrow(res) > 0L)
    data.frame(m = as.integer(as.numeric(res$m)), Pcons = as.numeric(res$Pcons),
               hr = as.numeric(res$hr), N = as.numeric(res$N),
               K = as.numeric(res$K), stringsAsFactors = FALSE) else NULL
  rec
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
sim_ids <- sim_from:sim_to
plan("sequential"); gc()
t_all <- proc.time()[3]
recs <- vector("list", length(sim_ids))
quals <- vector("list", length(sim_ids))
for (i in seq_along(sim_ids)) {
  r <- tryCatch(record_replicate(sim_ids[i]), error = function(e) {
    z <- na_rec(sim_ids[i]); z$err <- conditionMessage(e); z })
  # `quals[[i]] <- NULL` would DELETE the slot and misalign the list against
  # sim_ids; `quals[i] <- list(x)` stores a NULL without shrinking.
  quals[i] <- list(attr(r, "qual_table"))
  attr(r, "qual_table") <- NULL
  recs[[i]] <- r
  if (i %% 25L == 0L)
    cat(sprintf("[%s] %d/%d  (%.0f s elapsed)\n", cell_name, i, length(sim_ids),
                proc.time()[3] - t_all))
}
elapsed <- proc.time()[3] - t_all
plan("sequential"); gc()

out <- do.call(rbind, recs)
names(quals) <- as.character(sim_ids)

saveRDS(list(cell = cell_name, cell_spec = cell, results = out, qual = quals,
             meta = list(sim_from = sim_from, sim_to = sim_to,
                         n_workers = n_workers, elapsed = elapsed,
                         stop_threshold_override =
                           if (is.null(stop_thr_override)) "inherited(0.95)"
                           else if (is.na(stop_thr_override)) "NULL (no early stop)"
                           else as.character(stop_thr_override),
                         fs_version = as.character(utils::packageVersion("forestsearch")),
                         timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
        out_path)
cat(sprintf("[%s] done: %d reps, %.0f s (%.1f s/rep) -> %s\n",
            cell_name, nrow(out), elapsed, elapsed / nrow(out), out_path))
