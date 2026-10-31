# ladder.R -- ACTG 175 article stage 0 (2026-10-30), Part 2 §2.4.
# Usage (from the forestsearch repo root):
#   timeout 1200 /usr/bin/time -v Rscript dev/reports/actg175_article_stage0_2026-10-30/ladder.R <rung>
# Rungs: M-bin M-cont D-bin D-cont M-S3-FS M-S3-DINA M-S3-GRF B-S3-FS B-S3-DINA B-S3-GRF
# Runs the installed forestsearch, single process. The transplanted chunks under calls_2026b/ and
# calls_2026a/ are evaluated verbatim; the only changes are the overrides below, applied by
# rewriting the call (symbols kept), never by editing the transplants.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
rung <- args[[1L]]
REC  <- "dev/reports/actg175_article_stage0_2026-10-30"
stopifnot(dir.exists(REC))
out_dir <- file.path(REC, "rungs")
dir.create(out_dir, showWarnings = FALSE)

# ---- chunk runner: evaluate one named chunk of a transplant file in the global env ----
read_chunks <- function(file) {
  ln <- readLines(file, encoding = "UTF-8")
  hdr <- grep("^## ==== chunk '", ln)
  labs <- sub("^## ==== chunk '(.*)' ====$", "\\1", ln[hdr])
  ends <- c(hdr[-1] - 1L, length(ln))
  setNames(lapply(seq_along(hdr), function(i) ln[(hdr[i] + 1L):ends[i]]), labs)
}
run_chunk <- function(file, label) {
  ch <- read_chunks(file)
  stopifnot(label %in% names(ch))
  exprs <- parse(text = ch[[label]], keep.source = FALSE)
  for (e in exprs) eval(e, envir = globalenv())
  invisible(NULL)
}
run_file <- function(file) for (lab in names(read_chunks(file))) run_chunk(file, lab)

# ---- the overrides (§2.4): quiet, no details, no plots, sequential ----
fs_overrides <- list(quiet = TRUE, details = FALSE, plot.sg = FALSE, plot.grf = FALSE,
                     parallel_args = list(plan = "sequential"))
boot_overrides <- list(nb_boots = 5L, details = FALSE, parallel_args = list(plan = "sequential"))
override_call <- function(cl, ov, fn) {
  cl[[1L]] <- fn
  for (nm in names(ov)) cl[[nm]] <- ov[[nm]]
  cl
}
CALLS <- list()
forestsearch <- function(...) {
  cl <- override_call(sys.call(), fs_overrides, quote(forestsearch::forestsearch))
  CALLS[[length(CALLS) + 1L]] <<- cl
  eval(cl, parent.frame())
}
forestsearch_bootstrap_dofuture <- function(...) {
  cl <- override_call(sys.call(), boot_overrides, quote(forestsearch::forestsearch_bootstrap_dofuture))
  CALLS[[length(CALLS) + 1L]] <<- cl
  eval(cl, parent.frame())
}

suppressPackageStartupMessages(library(forestsearch))
cat("RUNG", rung, "| forestsearch", as.character(packageVersion("forestsearch")),
    "| built", packageDescription("forestsearch")$Built, "\n")
B2 <- file.path(REC, "calls_2026b"); A3 <- file.path(REC, "calls_2026a")
t_rung <- proc.time()[["elapsed"]]
res <- list(rung = rung)

# ---- the documents' selection chunk carries their own anchor checks (stop()/stopifnot()); a failure
# there is recorded, and the membership the chunk would have built is rebuilt the chunk's own way ----
run_selection <- function(file) {
  st <- tryCatch({ run_chunk(file, "selection"); "ok" },
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  if (!identical(st, "ok")) {
    label_fs <<- paste(fit$sg.harm, collapse = " & ")
    est <- fit$df.est
    in_H <<- est$treat.recommend[match(df$id, est$id)] == 0
    n_H <<- sum(in_H); n_Hc <<- sum(!in_H)
  }
  cat("selection chunk:", st, "\n")
  st
}

# ---- 2026b cells, by the paths of the documents' extract chunk ----
cells_2026b <- function(binary) {
  g <- fit$mr_inference; cc <- g$complement; fb <- g$field; fc <- fb$complement; js <- fb$joint_s
  whole <- if (binary) reference_or[1, ] else reference_md[1, ]
  w_est <- if (binary) whole$or else whole$md
  out <- list(
    region = paste(fit$sg.harm, collapse = " & "), N_H = n_H, N_Hc = n_Hc,
    whole_n = whole$n, whole_est = w_est, whole_lo = whole$lo, whole_hi = whole$hi,
    naive_est = g$naive$est, naive_lo = g$naive$lower, naive_hi = g$naive$upper,
    mr_est = g$debiased$est, mr_lo = g$debiased$lower, mr_hi = g$debiased$upper,
    field_lower_H = fb$lower_1s,
    c_naive_est = cc$naive$est, c_naive_lo = cc$naive$lower, c_naive_hi = cc$naive$upper,
    c_mr_est = cc$debiased$est, c_mr_lo = cc$debiased$lower, c_mr_hi = cc$debiased$upper,
    field_s_upper_Hc = fc$upper_1s_s,
    bonf_lower_H = js$bonf_lower_H, bonf_upper_Hc = js$bonf_upper_Hc,
    selection_bias = g$selection_bias, n_family = g$n_family)
  if (binary) {
    out$whole_events <- whole$ev_ctrl + whole$ev_trt
    out$events_H <- sum(df$y_binary_neg[in_H]); out$events_Hc <- sum(df$y_binary_neg[!in_H])
  }
  out
}

if (rung %in% c("M-bin", "M-cont")) {
  binary <- rung == "M-bin"
  stem <- if (binary) "binary" else "continuous"
  run_file(file.path(B2, paste0(stem, "_data_prep.R")))
  t0 <- proc.time()[["elapsed"]]
  run_file(file.path(B2, paste0(stem, "_fit.R")))
  res$fit_sec <- proc.time()[["elapsed"]] - t0
  res$selection_status <- run_selection(file.path(B2, paste0(stem, "_assembly.R")))
  res$cells <- cells_2026b(binary)
  res$membership <- in_H
  saveRDS(list(args = CALLS, fit = fit), file.path(out_dir, paste0(rung, "_fit.rds")))
}

if (rung %in% c("D-bin", "D-cont")) {
  binary <- rung == "D-bin"
  stem <- if (binary) "binary" else "continuous"
  run_file(file.path(B2, paste0(stem, "_data_prep.R")))
  fit <- readRDS(file.path(out_dir, paste0(if (binary) "M-bin" else "M-cont", "_fit.rds")))$fit
  g <- fit$mr_inference
  res$selection_status <- run_selection(file.path(B2, paste0(stem, "_assembly.R")))
  t0 <- proc.time()[["elapsed"]]
  dc <- tryCatch({
    for (lab in c("declcal", "declcal-table", "declcal-observed", "declcal-decision"))
      run_chunk(file.path(B2, paste0(stem, "_declcal.R")), lab)
    "ok"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  res$declcal_sec <- proc.time()[["elapsed"]] - t0
  res$declcal_status <- dc
  if (exists("cal") && exists("cal_c2")) {
    tb <- cal$c0$table; ip <- which(tb$c0 == c0_primary)
    adm <- c(cal$c0$admitted_calibrated, list(c2 = cal_c2$admitted_calibrated))
    res$cells <- list(
      T_H = if (exists("T_H")) T_H else NA_real_,
      max_T = if (exists("max_T")) max_T$T else max(cal$T_hat),
      c0_primary = c0_primary, kappa_primary = tb$kappa_hat[ip],
      implied_level = 2 * pnorm(tb$kappa_hat[ip]) - 1,
      fw_primary = tb$fw_size[ip], kappa_claim = cal_c2$kappa_hat, fw_claim = cal_c2$fw_size,
      admitted_any = any(lengths(lapply(adm, intersect, cal$screened)) > 0),
      admitted_counts = lengths(adm), table = tb, z_pstar = cal$z_pstar)
  }
  saveRDS(list(cal = if (exists("cal")) cal else NULL, cal_c2 = if (exists("cal_c2")) cal_c2 else NULL),
          file.path(out_dir, paste0(rung, "_cal.rds")))
}

s3_fitname <- c(FS = "fs", DINA = "fs_dina", GRF = "fs_grf")
s3_idfile  <- c(FS = "identifier_fs.R", DINA = "identifier_dina.R", GRF = "identifier_grf.R")
s3_bootchunk <- c(FS = "bootstrap", DINA = "fs-dina-bootstrap", GRF = "fs-grf-bootstrap")
s3_bootname  <- c(FS = "fs_bc", DINA = "fs_dina_bc", GRF = "fs_grf_bc")

if (grepl("^M-S3-", rung)) {
  m <- sub("^M-S3-", "", rung)
  run_file(file.path(A3, "data_prep.R"))
  t0 <- proc.time()[["elapsed"]]
  run_file(file.path(A3, s3_idfile[[m]]))
  res$fit_sec <- proc.time()[["elapsed"]] - t0
  f <- get(s3_fitname[[m]])
  g <- f$mr_inference; cc <- g$complement
  tr <- f$df.est$treat.recommend
  res$cells <- list(
    region = paste(f$sg.harm, collapse = " & "), n_H = sum(tr == 0), n_Hc = sum(tr != 0),
    naive_est = g$naive$est, naive_lo = g$naive$lower, naive_hi = g$naive$upper,
    mr_est = g$debiased$est, mr_lo = g$debiased$lower, mr_hi = g$debiased$upper,
    c_naive_est = cc$naive$est, c_naive_lo = cc$naive$lower, c_naive_hi = cc$naive$upper,
    c_mr_est = cc$debiased$est, c_mr_lo = cc$debiased$lower, c_mr_hi = cc$debiased$upper,
    n_family = g$n_family, selection_bias = g$selection_bias, mr_sec = g$timing_seconds)
  ids <- f$df.est[[id.name]][tr == 0]
  res$membership_ids <- sort(ids)
  saveRDS(list(args = CALLS, fit = f), file.path(out_dir, paste0(rung, "_fit.rds")))
}

if (grepl("^B-S3-", rung)) {
  m <- sub("^B-S3-", "", rung)
  run_file(file.path(A3, "data_prep.R"))
  assign(s3_fitname[[m]], readRDS(file.path(out_dir, paste0("M-S3-", m, "_fit.rds")))$fit)
  if (m == "DINA") run_chunk(file.path(A3, s3_idfile[[m]]), "dina-config")   # dina-config is setup only
  t0 <- proc.time()[["elapsed"]]
  bs <- tryCatch({ run_chunk(file.path(A3, "bootstrap_calls.R"), s3_bootchunk[[m]]); "ok" },
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  res$boot_sec <- proc.time()[["elapsed"]] - t0
  res$boot_status <- bs
  if (exists(s3_bootname[[m]])) {
    b <- get(s3_bootname[[m]])
    r <- as.data.frame(b$results)
    res$per_rep_min <- r$tmins_iteration
    res$timing_attr <- attr(b$results, "timing")
    res$H_estimates <- b$H_estimates; res$Hc_estimates <- b$Hc_estimates
    saveRDS(list(args = CALLS, boot = b), file.path(out_dir, paste0(rung, "_boot.rds")))
  }
}

res$rung_sec <- proc.time()[["elapsed"]] - t_rung
res$args <- CALLS
saveRDS(res, file.path(out_dir, paste0(rung, "_result.rds")))
str(res[setdiff(names(res), c("args", "membership", "membership_ids"))], max.level = 2, give.attr = FALSE)
cat("RUNG_DONE", rung, sprintf("%.1f", res$rung_sec), "\n")
