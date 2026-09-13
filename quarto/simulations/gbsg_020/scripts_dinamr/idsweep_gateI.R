# Gate I -- integrity, PER CELL, stop-on-failure (TASK_idsweep_2026-09-12).
#
# For the cell's sixteen bundles (consistency x 6, DINA x 5, GRF x 5):
#   FAIL-able:
#     all sixteen present; REPS rows with sim_id 1..REPS;
#     status only DETECTED / NO-DETECTION, and no replicate error (err_msg is
#       written only by the template's .safe_record() error path);
#     detected == (status == "DETECTED");
#     the identification and classification columns populated on EVERY
#       detected replicate: detected, n_sel, label, sg_def, sens, spec, ppv,
#       npv, and admitted_n on GRF; the four rates within [0, 1].
#   REPORTED, never a failure:
#     non-detections, counted;
#     realized trial prevalence mean(n_true / n) against the super-population
#       value in meta;
#     n_family (NA with MR off: MR's own fitted family) and n_cons_qual /
#       band_n on DINA and GRF (no consistency table) -- structural NA;
#     mr_ok (0 with MR off);
#     n_true identical across the sixteen runs (same seeds within the cell);
#     SAME DRAWS against every committed (git-tracked), MR-on bundle of this
#       template at the same coordinate that covers sim_id 1..500: n_true
#       identical() on 1..500.  A mismatch is a FINDING about the DGM path, not
#       a cell failure.  Only n_true is read from those bundles; none is used
#       as data in this sweep.
#
# usage: Rscript idsweep_gateI.R <cell-tag> <z1q|-> <n> <hr>
#   IDSWEEP_TAG / IDSWEEP_REPS / IDSWEEP_OUT_PREFIX / IDSWEEP_GATEI_RDS exist only so
#   the gate can be exercised on the committed pBoc smoke before the sweep.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4L)
CELL <- args[1]; Z <- args[2]; N <- as.integer(args[3]); H <- as.numeric(args[4])
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES  <- file.path(QMD_DIR, "results")
TAG  <- Sys.getenv("IDSWEEP_TAG", unset = "idsweep")
REPS <- as.integer(Sys.getenv("IDSWEEP_REPS", unset = "500"))
GI_RDS <- Sys.getenv("IDSWEEP_GATEI_RDS", unset = file.path(SCRATCH, "idsweep_gateI.rds"))
SAME_RANGE <- min(500L, REPS)   # 500 in the sweep; the replicate count when exercised on the smoke
z1q <- if (Z == "-") 0.25 else as.numeric(Z)

fail <- 0L; pass <- 0L
say <- function(ok, msg) { ok <- isTRUE(ok)
  if (ok) pass <<- pass + 1L else fail <<- fail + 1L
  if (!ok) cat(sprintf("[FAIL] %s\n", msg)); invisible(ok) }

runs <- rbind(data.frame(engine = "consistency", focus = c("effMaxSG", "effMinSG", "maxeffCons", "maxeff", "maxSG", "minSG")),
              data.frame(engine = "dina", focus = c("effMaxSG", "effMinSG", "maxSG", "minSG", "maxeffCons")),
              data.frame(engine = "grf",  focus = c("effMaxSG", "effMinSG", "maxSG", "minSG", "maxeffCons")))
stem <- function(e, f) sprintf("%s_%s_fb_mr_field_m1_h%03d_knoise0_n%d%s%s_nomr_%s_res_1_%d.rds",
  if (e == "consistency") "fs" else e, forestsearch::fs_focus_tag(e, f), round(100 * H), N,
  if (abs(z1q - 0.25) < 1e-12) "" else sprintf("_z1q%02d", round(100 * z1q)),
  if (f %in% c("effMaxSG", "effMinSG")) "_nb20" else "", TAG, REPS)
cat(sprintf("=== GATE I: cell %s (z1q %s, n %d, HR %.2f), tag %s, %d replicates ===\n", CELL, Z, N, H, TAG, REPS))

B <- list(); rows <- list()
for (j in seq_len(nrow(runs))) {
  e <- runs$engine[j]; f <- runs$focus[j]; k <- paste(e, f)
  p <- file.path(RES, stem(e, f))
  if (!say(file.exists(p), sprintf("%s: bundle missing: %s", k, basename(p)))) next
  b <- readRDS(p); r <- b$results; B[[k]] <- b
  say(nrow(r) == REPS && !anyNA(r$sim_id) && identical(as.integer(r$sim_id), seq_len(REPS)),
      sprintf("%s: %d rows, sim_id 1..%d", k, nrow(r), REPS))
  say(identical(b$meta$mr_inference, FALSE), sprintf("%s: meta$mr_inference FALSE", k))
  st <- r$status
  say(all(st %in% c("DETECTED", "NO-DETECTION")),
      sprintf("%s: status outside DETECTED/NO-DETECTION: %s", k,
              paste(unique(st[!st %in% c("DETECTED", "NO-DETECTION")]), collapse = ",")))
  nerr <- sum(!is.na(r$err_msg))
  say(nerr == 0L, sprintf("%s: %d replicate error(s): %s", k, nerr, paste(unique(na.omit(r$err_msg)), collapse = " | ")))
  d <- st %in% "DETECTED"
  say(identical(as.integer(r$detected %in% 1L), as.integer(d)), sprintf("%s: detected != (status == DETECTED)", k))
  idc <- c("detected", "n_sel", "label", "sg_def", "sens", "spec", "ppv", "npv", if (e == "grf") "admitted_n")
  # AMENDMENT 1 (2026-09-13, after Gate I stopped the sweep at p124_h150_n500).
  # The template's .classify() (template :848) returns NA for a rate whose
  # denominator is zero.  On a detected replicate the recorded counts say
  # exactly when: NPV when n_sel == n (the selection is the whole trial, no
  # predicted negatives), PPV when n_sel == 0, sensitivity when n_true == 0,
  # specificity when n_true == n.  Such an NA is UNDEFINED, not unrecorded; it is
  # counted and reported as structural.  Any NA not explained this way still
  # fails.  (First seen: DINA maxSG, 8 of 436 detected replicates, each with
  # n_sel = 500 and NPV NA, and no other NA on any detected row of the cell.)
  zero_den <- list(sens = r$n_true == 0L, spec = r$n_true == N, ppv = r$n_sel == 0L, npv = r$n_sel == N)
  und <- setNames(integer(4), names(zero_den))
  for (cc in idc) {
    present <- cc %in% names(r)
    na_d <- if (present) is.na(r[[cc]]) & d else d
    if (present && cc %in% names(zero_den)) {
      allowed <- na_d & (zero_den[[cc]] %in% TRUE)
      und[[cc]] <- sum(allowed); na_d <- na_d & !allowed
    }
    ok <- present && !any(na_d) &&
          (!is.character(r[[cc]]) || all(nzchar(r[[cc]][d & !is.na(r[[cc]])])))
    say(ok, sprintf("%s: `%s` NA on %d detected replicate(s) not explained by a zero denominator (%d detected)",
                    k, cc, sum(na_d), sum(d)))
  }
  for (cc in c("sens", "spec", "ppv", "npv")) {
    x <- r[[cc]][d]; x <- x[!is.na(x)]
    say(all(x >= 0 & x <= 1), sprintf("%s: `%s` outside [0,1] on a detected replicate", k, cc))
  }
  struct <- c(n_family = all(is.na(r$n_family)),
              n_cons_qual = all(is.na(r$n_cons_qual)), band_n = all(is.na(r$band_n)))
  rows[[length(rows) + 1L]] <- data.frame(
    engine = e, sg_focus = f, detected = sum(d), no_detection = sum(st %in% "NO-DETECTION"),
    errors = nerr, prev_realized = mean(r$n_true / N), prev_super = b$meta$harm_prevalence_super,
    n_sel_NA_on_nondetected = sum(is.na(r$n_sel[!d])),
    n_family_allNA = struct[["n_family"]],
    n_cons_qual = if (struct[["n_cons_qual"]]) "all NA" else sprintf("populated %d/%d det", sum(!is.na(r$n_cons_qual[d])), sum(d)),
    band_n = if (struct[["band_n"]]) "all NA" else sprintf("populated %d/%d det", sum(!is.na(r$band_n[d])), sum(d)),
    admitted_n = if (e == "grf") sprintf("finite %d/%d", sum(is.finite(r$admitted_n)), nrow(r)) else "--",
    mr_ok_max = suppressWarnings(max(r$mr_ok, na.rm = TRUE)),
    undefined_rates = if (any(und > 0)) paste(sprintf("%s:%d", names(und)[und > 0], und[und > 0]), collapse = " ") else "none",
    stringsAsFactors = FALSE)
}
TB <- do.call(rbind, rows)
cat("\n--- per run ---\n")
op <- options(width = 250); print(transform(TB, prev_realized = round(prev_realized, 4), prev_super = round(prev_super, 4)),
                                  row.names = FALSE); options(op)

cat("\n--- structural NA (by design, not failures) ---\n")
for (e in unique(TB$engine)) {
  s <- TB[TB$engine == e, ]
  cat(sprintf("  %-11s n_family all NA on %d/%d runs (MR-only quantity, MR off)", e, sum(s$n_family_allNA), nrow(s)))
  if (e != "consistency") cat(sprintf("; n_cons_qual all NA %d/%d, band_n all NA %d/%d (no consistency table)",
                                      sum(s$n_cons_qual == "all NA"), nrow(s), sum(s$band_n == "all NA"), nrow(s)))
  cat("\n")
}
if (any(!TB$n_family_allNA)) cat("  NOTE: n_family populated on a run with MR off -- unexpected; reported, not gated\n")
cat(sprintf("  mr_ok max over all runs: %s\n", paste(unique(TB$mr_ok_max), collapse = ",")))
ud <- TB[TB$undefined_rates != "none", ]
cat(sprintf("  undefined classification rates on detected replicates (zero denominator; Amendment 1): %s\n",
            if (nrow(ud)) paste(sprintf("%s %s [%s]", ud$engine, ud$sg_focus, ud$undefined_rates), collapse = "; ") else "none"))

cat("\n--- counts ---\n")
cat(sprintf("  detections %d, non-detections %d, errors %d over %d replicates x %d runs\n",
            sum(TB$detected), sum(TB$no_detection), sum(TB$errors), REPS, nrow(TB)))
cat(sprintf("  realized prevalence (all runs): %s ; super-population: %s\n",
            paste(unique(round(TB$prev_realized, 4)), collapse = ","), paste(unique(round(TB$prev_super, 4)), collapse = ",")))

cat("\n--- same draws within the cell ---\n")
nt <- vapply(B, function(b) paste(b$results$n_true, collapse = ","), character(1))
tr <- vapply(B, function(b) paste(deparse(b$truth), collapse = ""), character(1))
within_nt <- length(unique(nt)) == 1L; within_tr <- length(unique(tr)) == 1L
cat(sprintf("  n_true identical across the %d runs: %s ; truth identical: %s%s\n", length(B), within_nt, within_tr,
            if (!within_nt) "  <-- FINDING (reported, not gated)" else ""))

cat("\n--- SAME DRAWS against committed bundles at this coordinate (sim_id 1-500; a finding, not a failure) ---\n")
trk <- basename(system2("git", c("-C", shQuote(normalizePath(QMD_DIR)), "ls-files", "results"), stdout = TRUE))
rx <- sprintf("^[a-z]+_[A-Za-z]+_fb_mr_field_m1_h%03d_knoise0_n%d%s(_nb[0-9]+)?(_j[0-9]+)?_([A-Za-z0-9]+)_(res|combined)_([0-9]+)_([0-9]+)\\.rds$",
              round(100 * H), N, if (abs(z1q - 0.25) < 1e-12) "" else sprintf("_z1q%02d", round(100 * z1q)))
cand <- trk[grepl(rx, trk)]
SD <- NULL
if (length(B)) {
  ref <- B[[1]]$results; ref_nt <- ref$n_true[match(seq_len(SAME_RANGE), ref$sim_id)]
  if (length(cand)) {
    cc <- data.frame(file = cand, ctag = sub(rx, "\\3", cand), kind = sub(rx, "\\4", cand),
                     lo = as.integer(sub(rx, "\\5", cand)), hi = as.integer(sub(rx, "\\6", cand)), stringsAsFactors = FALSE)
    cc$covers <- cc$lo <= 1L & cc$hi >= SAME_RANGE
    for (ct in unique(cc$ctag)) {
      s <- cc[cc$ctag == ct, ]
      pick <- s[s$covers & s$kind == "combined", ]; if (!nrow(pick)) pick <- s[s$covers, ]
      if (!nrow(pick)) { cat(sprintf("  %-14s skipped: no tracked bundle covers sim_id 1-%d (%s)\n", ct, SAME_RANGE,
                                     paste(unique(sprintf("%d-%d", s$lo, s$hi)), collapse = ", "))); next }
      f <- pick$file[1]; cb <- readRDS(file.path(RES, f)); r2 <- cb$results
      nt2 <- r2$n_true[match(seq_len(SAME_RANGE), r2$sim_id)]
      idn <- identical(ref_nt, nt2)
      eqv <- isTRUE(all.equal(as.numeric(ref_nt), as.numeric(nt2)))
      sb <- cb$meta$seed_base %||% NA
      cat(sprintf("  %-14s %s  identical()=%s%s  seed_base=%s  engine=%s focus=%s\n", ct, f, idn,
                  if (!idn) sprintf("  equal-as-numbers=%s  <-- FINDING", eqv) else "", format(sb),
                  format(cb$meta$subgroup_method %||% NA), format(cb$meta$sg_focus %||% NA)))
      SD <- rbind(SD, data.frame(cell = CELL, campaign = ct, file = f, identical = idn, equal_numeric = eqv,
                                 seed_base = as.character(sb), stringsAsFactors = FALSE))
    }
  } else cat("  no committed MR-on bundle of this template at this coordinate\n")
}

if (nzchar(GI_RDS)) {
  acc <- if (file.exists(GI_RDS)) readRDS(GI_RDS) else list()
  acc[[CELL]] <- list(per_run = TB, within_n_true = within_nt, within_truth = within_tr, same_draws = SD,
                      pass = pass, fail = fail, at = Sys.time())
  saveRDS(acc, GI_RDS)
}
cat(sprintf("\nGATE I %s: %d checks passed, %d failed -- %s\n", CELL, pass, fail, if (fail) "FAIL: STOP" else "PASS"))
if (fail) quit(status = 1L)
