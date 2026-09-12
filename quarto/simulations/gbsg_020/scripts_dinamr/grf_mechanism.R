# Part B mechanism diagnostic: WHY a GRF replicate records no subgroup, and how
# often the frontier BAND (as opposed to the admission floor) is what emptied.
# TASK_dinamr_blockC_grfprobe_2026-09-11, Part B.  Cost/mechanism only -- no
# coverage, no comparison, no acceptance criterion.
#
# There are three distinct no-selection mechanisms on this path, and the
# recorder collapses all three into the same all-NA NO-DETECTION row (the
# template returns at `if (!found) { rec$status <- "NO-DETECTION"; return(rec) }`
# before n_family is written).  They are separable from the forestsearch object,
# which carries grf_res even on a non-detection
# (.forestsearch_grf_select()'s not-found branch returns grf_res; forestsearch()
# puts it in out$grf_res):
#
#   (1) NO DR CANDIDATES      grf_res$candidates NULL or 0 rows -- no threshold
#       or pair candidate met n.min on the DR scores.
#   (2) DR FLOOR / BAND EMPTY native .grf_frontier_select() returned NULL inside
#       grf.subg.harm.survival(), so the object has no $candidates at all.
#       With dmin.grf = 0.0 the eligible set is {DR effect >= 0}; the band
#       (1 - 0.20) * max cannot then exclude the maximum, so this is the FLOOR
#       emptying, not the band.
#   (3) EFFECT ADMISSION EMPTY .grf_reselect_on_effect() found no candidate with
#       HR >= 0.90 (the resolved admission floor, hr.threshold, NOT dmin.grf)
#       and set grf_res$admitted_n <- 0L, sg_def <- NULL.
#
# The BAND is checked explicitly at both stages by recomputing
# (1 - nbhd) * max(effect) over the eligible set, which is what
# .compute_inclusion_band() applies (subgroup_consistency_helpers.R:783-785).
suppressPackageStartupMessages(library(forestsearch))
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES     <- file.path(QMD_DIR, "results")

# Template constants, quoted by line (sim_fs_maxeffCons_fb_mr_field_m1_template.qmd).
SEED_BASE      <- 8316951L      # seed_base; seeds are seed_base + sim_id
ANALYSIS_TIME  <- 84            # line 489
CENS_ADJUST    <- log(1.5)      # line 490
N_SUPER        <- 100000L       # line 491
DGM_MODEL      <- "alt"         # line 488
GRF_SELECTION  <- "frontier"    # line 503
GRF_SELECT_ST  <- "effect"      # line 504
GRF_DEPTH      <- 2L            # line 505
DMIN_GRF       <- 0.0           # line 506  (Larry, 2026-09-11)
HR_THRESHOLD   <- 0.90          # line 531
HR_CONSISTENCY <- 0.80          # line 532
NBHD           <- 0.20          # FS_S7_NBHD
SG_FOCUS       <- "effMaxSG"    # FS_S7_FOCUS
SEL_RULE       <- "neighborhood"

probe_file <- function(hr, n, z1q)
  file.path(RES, sprintf("grf_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_grfprobe_res_1_36.rds",
                         round(100*hr), n, if (z1q) "_z1q60" else ""))

band_card <- function(eff, nbhd) {           # |band| over an eligible effect vector
  if (!length(eff)) return(0L)
  sum(eff >= (1 - nbhd) * max(eff, na.rm = TRUE), na.rm = TRUE)
}

one_corner <- function(hr, n, z1q, sim_ids = NULL, label = "") {
  cat(sprintf("\n===== MECHANISM: %s (HR %.2f, n %d, %s) =====\n", label, hr, n,
              if (z1q) "31%" else "12.4%"))
  pf <- probe_file(hr, n, z1q)
  if (!file.exists(pf)) { cat("  probe bundle absent:", basename(pf), "\n"); return(NULL) }
  pr <- readRDS(pf)$results
  nd <- pr$sim_id[!(pr$detected %in% 1L)]
  # ALL replicates, not only the non-detections: "how often the frontier band
  # comes back empty" is a rate over every replicate, and on a DETECTED
  # replicate the band cardinality is the direct evidence that it did not empty.
  if (is.null(sim_ids)) sim_ids <- pr$sim_id
  cat(sprintf("  probe: %d replicates, %d detected, %d NOT detected -> re-running %d\n",
              nrow(pr), sum(pr$detected %in% 1L), length(nd), length(sim_ids)))
  if (!length(sim_ids)) { cat("  no replicates in the bundle.\n"); return(NULL) }
  detmap <- setNames(pr$detected %in% 1L, pr$sim_id)

  z1qv <- if (z1q) 0.60 else 0.25
  k    <- calibrate_k_inter(target_hr_harm = hr, model = DGM_MODEL,
                            use_ahr = FALSE, z1_quantile = z1qv)
  dgm  <- setup_gbsg_dgm(model = DGM_MODEL, k_inter = k, z1_quantile = z1qv,
                         n_super = N_SUPER, seed = SEED_BASE)
  confs <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")  # template line 667

  rows <- lapply(sim_ids, function(sid) {
    df <- simulate_from_dgm(dgm, n = n, analysis_time = ANALYSIS_TIME,
                            cens_adjust = CENS_ADJUST, seed = SEED_BASE + sid)
    df$id <- seq_len(nrow(df))
    cf <- intersect(confs, names(df))
    t0 <- proc.time()[3]
    fs <- tryCatch(forestsearch(
      df.analysis = df, outcome.name = "y_sim", event.name = "event_sim",
      treat.name = "treat_sim", id.name = "id", flag_harm.name = "flag_harm",
      confounders.name = cf, is.RCT = TRUE, seedit = SEED_BASE + sid,
      quiet = TRUE, sg_focus = SG_FOCUS, subgroup_method = "grf",
      hr.threshold = HR_THRESHOLD, hr.consistency = HR_CONSISTENCY,
      pconsistency.threshold = 0.90,  # template line 533 n.min = NULL,
      selection_rule = SEL_RULE, effect_neighborhood = NBHD,
      grf_selection = GRF_SELECTION, grf_depth = GRF_DEPTH,
      dmin.grf = DMIN_GRF, grf_select_statistic = GRF_SELECT_ST,
      mr_inference = FALSE),
      error = function(e) structure(list(err = conditionMessage(e)), class = "fserr"))
    el <- proc.time()[3] - t0
    if (inherits(fs, "fserr"))
      return(data.frame(sim_id = sid, mech = "ERROR", n_cand = NA_integer_,
                        n_elig_eff = NA_integer_, n_band_eff = NA_integer_,
                        admitted_n = NA_integer_, max_hr = NA_real_,
                        secs = el, note = substr(fs$err, 1, 90)))
    g  <- fs$grf_res
    cd <- if (is.null(g)) NULL else g$candidates
    adm <- if (is.null(g$admitted_n)) NA_integer_ else as.integer(g$admitted_n)
    if (is.null(cd) || !nrow(cd)) {
      # No $candidates on the object: either the DR candidate pool was empty or
      # the native frontier returned before candidates were attached.
      return(data.frame(sim_id = sid, mech = "NO_DR_CANDIDATES_OR_DR_FLOOR",
                        n_cand = 0L, n_elig_eff = NA_integer_, n_band_eff = NA_integer_,
                        admitted_n = adm, max_hr = NA_real_, secs = el, note = ""))
    }
    # Effect-stage cardinalities, recomputed exactly as .grf_reselect_on_effect
    # -> .grf_frontier_select does: eligible = {HR >= exp(log(hr.threshold))},
    # band = {HR >= (1 - nbhd) * max(eligible HR)}.
    hrv  <- if ("sel_effect" %in% names(cd)) exp(cd$sel_effect) else rep(NA_real_, nrow(cd))
    hrv  <- hrv[is.finite(hrv)]
    elig <- hrv[hrv >= HR_THRESHOLD]
    data.frame(sim_id = sid,
               mech = if (!length(hrv)) "NO_SCORABLE_CANDIDATE"
                      else if (!length(elig)) "EFFECT_ADMISSION_EMPTY"
                      else if (band_card(elig, NBHD) == 0L) "BAND_EMPTY"
                      else "ADMITTED_BAND_NONEMPTY",
               n_cand = nrow(cd), n_elig_eff = length(elig),
               n_band_eff = band_card(elig, NBHD), admitted_n = adm,
               max_hr = if (length(hrv)) max(hrv) else NA_real_,
               secs = el, note = "")
  })
  out <- do.call(rbind, rows)
  out$detected_in_probe <- unname(detmap[as.character(out$sim_id)])
  cat("\n  mechanism split over ALL replicates (probe detection beside it):\n")
  print(as.data.frame(table(mechanism = out$mech, detected = out$detected_in_probe)),
        row.names = FALSE)
  cat(sprintf("\n  >> FRONTIER BAND EMPTY on %d of %d replicates (%.4f)\n",
              sum(out$mech == "BAND_EMPTY"), nrow(out), mean(out$mech == "BAND_EMPTY")))
  cat(sprintf("  >> effect-admission empty on %d (%.4f); admitted_n == 0 recorded on %d\n",
              sum(out$mech == "EFFECT_ADMISSION_EMPTY"),
              mean(out$mech == "EFFECT_ADMISSION_EMPTY"),
              sum(out$admitted_n %in% 0L)))
  ok <- out$mech == "ADMITTED_BAND_NONEMPTY"
  if (any(ok)) {
    cat(sprintf("  >> band cardinality where the admission set is non-empty: min %d q10 %g med %g q90 %g max %d\n",
        min(out$n_band_eff[ok]), quantile(out$n_band_eff[ok], .1, names = FALSE),
        median(out$n_band_eff[ok]), quantile(out$n_band_eff[ok], .9, names = FALSE),
        max(out$n_band_eff[ok])))
    cat(sprintf("  >> DR-candidate pool size : min %d med %g max %d ; admitted (HR >= %.2f): min %d med %g max %d\n",
        min(out$n_cand[ok]), median(out$n_cand[ok]), max(out$n_cand[ok]), HR_THRESHOLD,
        min(out$n_elig_eff[ok]), median(out$n_elig_eff[ok]), max(out$n_elig_eff[ok])))
    cat(sprintf("  >> admitted_n agrees with the recomputed eligible count on %d of %d rows\n",
        sum(out$admitted_n[ok] == out$n_elig_eff[ok], na.rm = TRUE), sum(ok)))
  }
  cat(sprintf("  >> identification-only re-run wall: median %.2f s, p90 %.2f s, max %.2f s\n",
              median(out$secs), quantile(out$secs, .9, names = FALSE), max(out$secs)))
  if (any(out$mech == "ERROR")) {
    cat("  >> ERRORS (reported, NOT fixed):\n"); print(out[out$mech == "ERROR", ], row.names = FALSE)
  }
  out$corner <- label
  out
}

corners <- list(
  list(1.50,  500L, FALSE, "g_p124_h150_n500"),
  list(1.50, 1500L, FALSE, "g_p124_h150_n1500"),
  list(1.50,  500L, TRUE,  "g_p31_h150_n500"),
  list(1.50, 1500L, TRUE,  "g_p31_h150_n1500"),
  list(1.00,  500L, FALSE, "g_p124_h100_n500"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) corners <- Filter(function(cc) cc[[4]] %in% args, corners)
ALL <- do.call(rbind, lapply(corners, function(cc)
  one_corner(cc[[1]], cc[[2]], cc[[3]], label = cc[[4]])))
if (!is.null(ALL)) {
  cat("\n\n===== ALL CORNERS, MECHANISM SPLIT =====\n")
  print(table(ALL$corner, ALL$mech))
  cat(sprintf("\nFRONTIER BAND EMPTY, ALL CORNERS: %d of %d replicates (%.4f)\n",
              sum(ALL$mech == "BAND_EMPTY"), nrow(ALL), mean(ALL$mech == "BAND_EMPTY")))
  cat("\nWhat a replicate records when nothing is selected: the template returns at\n")
  cat("`if (!found) { rec$status <- \"NO-DETECTION\"; return(rec) }` BEFORE n_family is\n")
  cat("written, so every no-selection replicate is the all-NA NO-DETECTION row and the\n")
  cat("three mechanisms are indistinguishable from the committed columns.  That is why\n")
  cat("this diagnostic re-runs the identifier rather than reading the bundle.\n")
  saveRDS(ALL, file.path(SCRATCH, "grf_mechanism.rds"))
}
