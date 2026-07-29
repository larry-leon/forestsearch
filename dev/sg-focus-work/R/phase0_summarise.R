# =============================================================================
# Phase 0 -- seed fidelity check + per-cell branch/divergence summary
# =============================================================================
# Usage (from the repository root):
#   Rscript dev/sg-focus-work/R/phase0_summarise.R <glob-of-run-rds> [tag]

source("dev/sg-focus-work/R/phase0_cells.R")

args <- commandArgs(trailingOnly = TRUE)
pattern <- if (length(args) >= 1L) args[[1]] else "dev/sg-focus-work/out/run_*.rds"
tag     <- if (length(args) >= 2L) args[[2]] else "run"

files <- Sys.glob(pattern)
if (!length(files)) stop("no files matching ", pattern)

# Cut labels render the same membership two ways depending on whether the cut
# came in as a forced cut ("{meno == 0}") or as an auto-dichotomised binary
# ("!{meno}").  Normalise before comparing so a pure relabelling is not counted
# as a fidelity failure; membership is checked separately via n_sel.
normalise_sgdef <- function(x) {
  vapply(x, function(s) {
    if (is.na(s)) return("<none>")
    parts <- trimws(strsplit(s, "&", fixed = TRUE)[[1]])
    parts <- sub("^\\{ *([^=<> ]+) *== *0 *\\}$",  "!{\\1}", parts)
    parts <- sub("^\\{ *([^=<> ]+) *== *1 *\\}$",   "{\\1}", parts)
    parts <- sub("^!\\{ *([^=<> ]+) *== *0 *\\}$",  "{\\1}", parts)
    parts <- sub("^!\\{ *([^=<> ]+) *== *1 *\\}$", "!{\\1}", parts)
    parts <- gsub(" +", " ", parts)
    paste(sort(parts), collapse = " & ")
  }, character(1), USE.NAMES = FALSE)
}

saved_for <- function(cell) {
  p <- file.path(SIM_DIR, sprintf("fs_t1_t2_m1_%s_combined_1_500.rds", cell))
  if (!file.exists(p)) return(NULL)
  readRDS(p)$results
}

# Pool shards belonging to the same cell.
bundles <- lapply(files, readRDS)
cells <- vapply(bundles, `[[`, character(1), "cell")

fid_rows <- list(); sum_rows <- list(); pcons_rows <- list(); div_detail <- list()

for (cl in unique(cells)) {
  r <- do.call(rbind, lapply(bundles[cells == cl], `[[`, "results"))
  r <- r[order(r$sim_id), ]
  r <- r[!duplicated(r$sim_id), ]

  # ---- seed fidelity -------------------------------------------------------
  s <- saved_for(cl)
  fid <- data.frame(cell = cl, n_compared = NA_integer_, n_match_sgdef = NA_integer_,
                    n_match_sgdef_norm = NA_integer_,
                    n_match_det = NA_integer_, n_match_nsel = NA_integer_,
                    stringsAsFactors = FALSE)
  if (!is.null(s)) {
    ss <- s[match(r$sim_id, s$sim_id), ]
    keep <- !is.na(ss$sim_id)
    a <- ifelse(is.na(r$sg_def[keep]),  "<none>", r$sg_def[keep])
    b <- ifelse(is.na(ss$sg_def[keep]), "<none>", ss$sg_def[keep])
    fid$n_compared     <- sum(keep)
    fid$n_match_sgdef  <- sum(a == b)
    fid$n_match_sgdef_norm <- sum(normalise_sgdef(a) == normalise_sgdef(b))
    fid$n_match_det    <- sum(r$detected[keep] == ss$detected[keep])
    # n_sel is recorded only on gate-successful replicates in the saved bundle;
    # compare on the rows where both are present.
    both <- keep & !is.na(r$n_sel) & !is.na(ss$n_sel)
    fid$n_match_nsel   <- sum(r$n_sel[both] == ss$n_sel[both])
    fid$n_nsel_compared <- sum(both)
    mism <- which(normalise_sgdef(a) != normalise_sgdef(b))
    if (length(mism))
      fid$example_mismatch <- paste0("sim ", r$sim_id[keep][mism[1]], ": new='",
                                     a[mism[1]], "' old='", b[mism[1]], "'")
    else fid$example_mismatch <- ""
  }
  fid_rows[[cl]] <- fid

  # ---- branch / divergence -------------------------------------------------
  # A replicate contributes to the branch mix only if the search actually ran to
  # a selection (detected == 1).  Replicates with no qualifying candidate at all
  # (n_passed == 0) select nothing, so no selection rule was exercised.
  d <- r[r$ok == 1L & r$detected == 1L & !is.na(r$sel_m), ]

  b1 <- d$early_stop_triggered %in% TRUE
  b2 <- !b1
  div <- d$sel_m != d$first_qual_m          # selected is NOT the HR-argmax qualifier

  sum_rows[[cl]] <- data.frame(
    cell = cl,
    n_reps = nrow(r),
    n_run_ok = sum(r$ok == 1L),
    n_detected = sum(r$detected == 1L, na.rm = TRUE),
    n_selected = nrow(d),
    branch1_n = sum(b1), branch1_rate = mean(b1),
    branch2_n = sum(b2), branch2_rate = mean(b2),
    diverge_n = sum(div), diverge_rate = mean(div),
    diverge_b1_n = sum(div & b1), diverge_b1_rate = if (sum(b1)) mean(div[b1]) else NA_real_,
    diverge_b2_n = sum(div & b2), diverge_b2_rate = if (sum(b2)) mean(div[b2]) else NA_real_,
    med_cand_total = stats::median(d$n_cand_total),
    med_cand_eval  = stats::median(d$n_cand_evaluated),
    med_n_qual     = stats::median(d$n_qual),
    stringsAsFactors = FALSE)

  q <- stats::quantile(d$sel_Pcons, c(0, .10, .25, .50, .75, .90, 1), na.rm = TRUE)
  pcons_rows[[cl]] <- data.frame(
    cell = cl, n = sum(is.finite(d$sel_Pcons)),
    min = q[1], p10 = q[2], p25 = q[3], median = q[4],
    p75 = q[5], p90 = q[6], max = q[7],
    mean = mean(d$sel_Pcons, na.rm = TRUE),
    pct_ge_095 = mean(d$sel_Pcons >= 0.95, na.rm = TRUE),
    pct_in_090_095 = mean(d$sel_Pcons >= 0.90 & d$sel_Pcons < 0.95, na.rm = TRUE),
    stringsAsFactors = FALSE)

  if (any(div))
    div_detail[[cl]] <- data.frame(
      cell = cl, sim_id = d$sim_id[div], branch = ifelse(b1[div], 1L, 2L),
      sel_m = d$sel_m[div], sel_Pcons = d$sel_Pcons[div], sel_hr = d$sel_hr[div],
      passed_over_m = d$first_qual_m[div],
      passed_over_Pcons = d$first_qual_Pcons[div],
      passed_over_hr = d$first_qual_hr[div],
      hr_rank_selected = d$hr_rank_selected[div],
      n_qual = d$n_qual[div], stringsAsFactors = FALSE)
}

fid <- do.call(rbind, fid_rows)
sm  <- do.call(rbind, sum_rows)
pc  <- do.call(rbind, pcons_rows)
dd  <- if (length(div_detail)) do.call(rbind, div_detail) else NULL

ord <- match(names(PHASE0_CELLS), sm$cell); ord <- ord[!is.na(ord)]
sm <- sm[ord, ]; pc <- pc[match(sm$cell, pc$cell), ]; fid <- fid[match(sm$cell, fid$cell), ]

cat("\n================ SEED FIDELITY (vs saved combined bundles) ================\n")
print(fid, row.names = FALSE)
cat("\n================ BRANCH MIX / DIVERGENCE ================\n")
print(sm, row.names = FALSE, digits = 4)
cat("\n================ DISTRIBUTION OF SELECTED Pcons ================\n")
print(pc, row.names = FALSE, digits = 4)
if (!is.null(dd)) {
  cat("\n================ DIVERGENT REPLICATES (detail) ================\n")
  print(utils::head(dd, 60), row.names = FALSE, digits = 4)
  cat(sprintf("\n(%d divergent replicates total)\n", nrow(dd)))
}

saveRDS(list(fidelity = fid, summary = sm, pcons = pc, divergences = dd,
             files = files),
        sprintf("dev/sg-focus-work/out/summary_%s.rds", tag))
cat(sprintf("\nsaved -> dev/sg-focus-work/out/summary_%s.rds\n", tag))
