# Timing ladder, post-selection vignette stage 0 (2026-10-28) -- one rung per invocation.
#   Rscript ladder.R <rung>      (run from this directory; see run_ladder.sh)
# Rungs: M-FS M-DINA M-GRF M-FS-ij (forestsearch() fits with MR), then
#        B-FS B-DINA B-GRF (the document's bootstrap call on the M fit, nb_boots = 5).
# table_s2_calls/data_prep.R is evaluated verbatim (parallel-setup, load-packages,
# data-setup, dina-config chunks). For an M rung the transplanted call file is
# parsed, every expression before `<obj> <- forestsearch(...)` is evaluated
# verbatim (the chunk's `t0 <- proc.time()` line), and the forestsearch() call
# is evaluated with exactly the task's overrides applied to the call object.
# For a B rung the matching forestsearch_bootstrap_dofuture() call is taken from
# table_s2_calls/bootstrap_calls.R and evaluated with nb_boots = 5 and the
# overrides parallel_args = list(plan = "sequential"), details = FALSE, on the
# fit saved by the M rung. Nothing after the call in either chunk is evaluated.

rung <- commandArgs(trailingOnly = TRUE)[1]
spec <- list(
  "M-FS"    = list(id = "fs",   kind = "M", obj = "fs",      ci = NULL),
  "M-DINA"  = list(id = "dina", kind = "M", obj = "fs_dina", ci = NULL),
  "M-GRF"   = list(id = "grf",  kind = "M", obj = "fs_grf",  ci = NULL),
  "M-FS-ij" = list(id = "fs",   kind = "M", obj = "fs",      ci = "ij"),
  "B-FS"    = list(id = "fs",   kind = "B", obj = "fs_bc",      fit = "M-FS",   fit_obj = "fs"),
  "B-DINA"  = list(id = "dina", kind = "B", obj = "fs_dina_bc", fit = "M-DINA", fit_obj = "fs_dina"),
  "B-GRF"   = list(id = "grf",  kind = "B", obj = "fs_grf_bc",  fit = "M-GRF",  fit_obj = "fs_grf"))[[rung]]
if (is.null(spec)) stop("unknown rung: ", rung)
dir.create("rungs", showWarnings = FALSE)
dir.create("fits", showWarnings = FALSE)

suppressPackageStartupMessages(library(forestsearch))
cat(sprintf("rung %s | forestsearch %s from %s | built %s\n", rung,
            packageVersion("forestsearch"), find.package("forestsearch"),
            packageDescription("forestsearch")$Built))
stopifnot(packageVersion("forestsearch") == "0.4.0.9000")

# ---- data preparation, verbatim ----------------------------------------------
invisible(capture.output(suppressPackageStartupMessages(
  for (e in parse("table_s2_calls/data_prep.R", keep.source = FALSE)) eval(e, envir = globalenv()))))
stopifnot(nrow(df.analysis) == 686L)

assign_target <- function(e, obj, fn) is.call(e) && identical(e[[1]], as.name("<-")) &&
  identical(e[[2]], as.name(obj)) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name(fn))

warns <- character(0)
run_timed <- function(cl) {
  gc(reset = TRUE)
  t0 <- proc.time()
  val <- withCallingHandlers(eval(cl, envir = globalenv()),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  tt <- proc.time() - t0
  mem <- gc()
  cat(sprintf("---- call elapsed %.2f s (user %.2f, sys %.2f); R max used %.1f Mb\n",
              tt[["elapsed"]], tt[["user.self"]], tt[["sys.self"]], sum(mem[, ncol(mem)])))
  cat(sprintf("---- warnings (%d):\n", length(warns))); if (length(warns)) cat(paste0("  W: ", unique(warns)), sep = "\n")
  list(val = val, sec = unname(tt[["elapsed"]]))
}
f6 <- function(x) if (is.null(x) || !length(x) || is.na(x[1])) NA_character_ else sprintf("%.6f", x[1])

if (spec$kind == "M") {
  # ---- the transplanted forestsearch() call --------------------------------
  src <- file.path("table_s2_calls", paste0(spec$id, "_call.R"))
  exprs <- parse(src, keep.source = FALSE)
  is_fit <- vapply(exprs, assign_target, logical(1), obj = spec$obj, fn = "forestsearch")
  stopifnot(sum(is_fit) == 1L)
  idx <- which(is_fit)
  for (i in seq_len(idx - 1L)) eval(exprs[[i]], envir = globalenv())
  cl <- exprs[[idx]][[3]]
  call_verbatim <- cl
  cl$parallel_args <- quote(list(plan = "sequential"))
  cl$quiet <- TRUE
  cl$details <- FALSE
  cl$plot.sg <- FALSE
  if ("plot.grf" %in% names(formals(forestsearch::forestsearch))) cl$plot.grf <- FALSE
  if (!is.null(spec$ci)) cl$mr_inference_args$ci_method <- spec$ci
  cat("---- call as evaluated ----\n"); cat(deparse(cl, width.cutoff = 100L), sep = "\n")
  arg_vals <- lapply(as.list(cl)[-1], eval, envir = globalenv())
  arg_vals[[which(names(arg_vals) == "")[1]]] <- NULL   # positional df.analysis
  saveRDS(list(rung = rung, source = src, call_verbatim = call_verbatim, call_evaluated = cl,
               args_resolved = arg_vals), file.path("rungs", paste0("args_", rung, ".rds")))

  r <- run_timed(cl); fit <- r$val
  saveRDS(fit, file.path("fits", paste0("fit_", rung, ".rds")))

  label <- paste(fit$sg.harm, collapse = " & ")
  est <- fit$df.est
  in_H <- if (!is.null(est) && length(fit$sg.harm)) est$treat.recommend[match(df.analysis$id, est$id)] == 0 else rep(FALSE, nrow(df.analysis))
  nH <- sum(in_H); eH <- sum(df.analysis$status[in_H]); nHc <- sum(!in_H); eHc <- sum(df.analysis$status[!in_H])
  g <- fit$mr_inference; cc <- g$complement; fb <- g$field; fc <- fb$complement; js <- fb$joint_s
  q <- list(naive_H = g$naive$est, naive_lo_H = g$naive$lower, naive_hi_H = g$naive$upper,
            mr_H = g$debiased$est, mr_lo_H = g$debiased$lower, mr_hi_H = g$debiased$upper,
            naive_Hc = cc$naive$est, naive_lo_Hc = cc$naive$lower, naive_hi_Hc = cc$naive$upper,
            mr_Hc = cc$debiased$est, mr_lo_Hc = cc$debiased$lower, mr_hi_Hc = cc$debiased$upper,
            se_ij_H = g$debiased$se_ij, se_ij_Hc = cc$debiased$se_ij,
            field_lower_1s_H = fb$lower_1s, field_lower_2s_H = fb$lower_2s, field_upper_2s_H = fb$upper_2s,
            field_est2_H = fb$est2, field_s_upper_1s_Hc = fc$upper_1s_s,
            bonf_lower_H = js$bonf_lower_H, bonf_upper_Hc = js$bonf_upper_Hc)
  row <- data.frame(rung = rung, identifier = spec$id, fit_sec = round(r$sec, 2),
                    ci_method = if (is.null(g$ci_method)) NA else g$ci_method,
                    region = label, N_H = nH, E_H = eH, N_Hc = nHc, E_Hc = eHc,
                    n_family = if (is.null(g$n_family)) NA else g$n_family,
                    mr_draws = if (is.null(g$settings$draws)) NA else g$settings$draws,
                    reselection = if (is.null(g$settings$reselection)) NA else g$settings$reselection,
                    mr_sec = f6(g$timing_seconds), field_sec_H = f6(fb$timing_seconds), field_sec_Hc = f6(fc$timing_seconds),
                    as.data.frame(lapply(q, f6), stringsAsFactors = FALSE),
                    n_warnings = length(warns), stringsAsFactors = FALSE)
  saveRDS(list(rung = rung, row = row, sg.harm = fit$sg.harm, in_H = in_H, mr = q,
               settings = g$settings, warnings = warns,
               filter_counts = fit$find.grps$filter_counts),
          file.path("rungs", paste0("rung_", rung, ".rds")))
  qf <- "ladder_quantities_M.csv"
  write.table(row, qf, sep = ",", row.names = FALSE, col.names = !file.exists(qf), append = file.exists(qf))
  cat("---- record ----\n"); print(t(row))
  cat("sg.harm:\n"); print(fit$sg.harm)
  cat("MR settings:\n"); str(g$settings)

  # ---- GATE A (M-FS only) -------------------------------------------------
  if (rung == "M-FS") {
    r2 <- function(x) sprintf("%.2f", x)
    ok_region <- identical(label, "{er <= 0} & {size <= 35}") && nH == 61L
    ok_naive  <- identical(r2(c(q$naive_H, q$naive_lo_H, q$naive_hi_H)), c("2.54", "1.25", "5.17"))
    ok_mr     <- identical(r2(c(q$mr_H, q$mr_lo_H, q$mr_hi_H)), c("1.44", "0.62", "3.36"))
    cat(sprintf("GATE A: %s (region %s, N %d [%s]; naive %s [%s]; MR %s [%s])\n",
                if (ok_region && ok_naive && ok_mr) "PASS" else "FAIL", label, nH, ok_region,
                paste(r2(c(q$naive_H, q$naive_lo_H, q$naive_hi_H)), collapse = "/"), ok_naive,
                paste(r2(c(q$mr_H, q$mr_lo_H, q$mr_hi_H)), collapse = "/"), ok_mr))
    if (!(ok_region && ok_naive && ok_mr)) {
      cat("---- GATE A FAIL detail ----\nfull sg.harm:\n"); print(fit$sg.harm)
      cat("filter_counts:\n"); str(fit$find.grps$filter_counts)
      cat("args_call_all (selection/MR):\n")
      str(fit$args_call_all[intersect(names(fit$args_call_all),
          c("sg_focus", "selection_rule", "effect_neighborhood", "hr.threshold", "hr.consistency",
            "pconsistency.threshold", "conf_force", "conf.cont_jcuts", "n.min", "d0.min", "d1.min",
            "maxk", "mr_inference_args", "consistency_method", "use_twostage", "fs.splits"))])
      quit(status = 3L)
    }
  }
} else {
  # ---- the transplanted bootstrap call on the M fit -------------------------
  assign(spec$fit_obj, readRDS(file.path("fits", paste0("fit_", spec$fit, ".rds"))), envir = globalenv())
  exprs <- parse("table_s2_calls/bootstrap_calls.R", keep.source = FALSE)
  is_b <- vapply(exprs, assign_target, logical(1), obj = spec$obj, fn = "forestsearch_bootstrap_dofuture")
  stopifnot(sum(is_b) == 1L)
  cl <- exprs[[which(is_b)]][[3]]
  call_verbatim <- cl
  cl$nb_boots <- 5L
  cl$parallel_args <- quote(list(plan = "sequential"))
  cl$details <- FALSE
  cat("---- call as evaluated ----\n"); cat(deparse(cl, width.cutoff = 100L), sep = "\n")
  saveRDS(list(rung = rung, call_verbatim = call_verbatim, call_evaluated = cl, fit = spec$fit),
          file.path("rungs", paste0("args_", rung, ".rds")))
  r <- tryCatch(run_timed(cl), error = function(e) { cat("BOOTSTRAP ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) quit(status = 4L)
  bc <- r$val
  saveRDS(bc, file.path("fits", paste0("boot_", rung, ".rds")))
  res <- as.data.frame(bc$results)
  it <- res$tmins_iteration * 60
  H <- bc$H_estimates; Hc <- bc$Hc_estimates
  row <- data.frame(rung = rung, identifier = spec$id, call_sec = round(r$sec, 2),
                    per_rep_sec_total_over_5 = round(r$sec / 5, 2),
                    rep1_sec = round(it[res$boot_id == 1][1], 2),
                    rep_sec_iter = paste(round(it[order(res$boot_id)], 2), collapse = "/"),
                    n_found = sum(!is.na(res$H_biasadj_2)),
                    H2 = f6(H$H2), H2_lower = f6(H$H2_lower), H2_upper = f6(H$H2_upper),
                    H0 = f6(H$H0), Hc2 = f6(Hc$H2), Hc2_lower = f6(Hc$H2_lower), Hc2_upper = f6(Hc$H2_upper),
                    n_warnings = length(warns), stringsAsFactors = FALSE)
  saveRDS(list(rung = rung, row = row, results = res, H_estimates = H, Hc_estimates = Hc,
               SG_CIs = bc$SG_CIs, warnings = warns), file.path("rungs", paste0("rung_", rung, ".rds")))
  qf <- "ladder_quantities_B.csv"
  write.table(row, qf, sep = ",", row.names = FALSE, col.names = !file.exists(qf), append = file.exists(qf))
  cat("---- record ----\n"); print(t(row))
  cat("per-replicate results:\n"); print(res[, intersect(names(res), c("boot_id", "H_biasadj_1", "H_biasadj_2",
        "Hc_biasadj_2", "tmins_search", "tmins_iteration", "Pcons", "any_found"))], row.names = FALSE)
  cat("H_estimates:\n"); print(H); cat("Hc_estimates:\n"); print(Hc)
}
