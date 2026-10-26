# s1_profile.R -- the single sequential forestsearch() call on fixture V:
# wall-clock (3 timed calls after a warm-up) and an Rprof split (3 profiles)
# attributing the call to the Cox fit machinery vs the fit proper.
#
# Nothing is modified; the installed package is called as a user would.

source("stage0_common.R")
fx <- stage0_fixture()
kv(fixture = fx$label, fixture_n = nrow(fx$df))

# warm-up (package code loading, lazy-load DB, JIT) -- recorded, not used
t_warm <- system.time(fs <- run_fixture_call(fx))[["elapsed"]]
kv(call_warmup_s = t_warm)
cc <- candidate_counts(fs)
kv(search_n_evaluated     = cc$n_evaluated_search,
   search_n_reached_cox   = cc$n_reached_cox_fit,
   search_n_cox_ok        = cc$n_cox_fit_ok,
   search_n_passed_effect = cc$n_passed_effect,
   consistency_n_total    = cc$n_consistency_total,
   consistency_n_evaluated = cc$n_consistency_evaluated,
   consistency_n_passed   = cc$n_consistency_passed,
   consistency_early_stop = cc$early_stop_triggered,
   identified_subgroup    = cc$sg_harm)

# GATE: the candidate flow is fixed by the data and the 0.4.0 code, not by the
# machine.  These are the counts CRAN 0.4.0 produced on fixture V in the chat's
# dry run (and the 1,410 / 121 the 0.3.5 record carries).  A mismatch means the
# installed package is not computing what 0.4.0 computes here -- stop.
expected <- list(n_reached_cox_fit = 1410L, n_passed_effect = 121L,
                 n_consistency_total = 120L, sg_harm = "{er <= 0} & {pgr <= 26}")
bad <- names(expected)[!mapply(function(k, v) identical(as.character(cc[[k]]), as.character(v)),
                               names(expected), expected)]
if (length(bad)) {
  cat("GATE FAILED: candidate flow differs from CRAN 0.4.0 on fixture V for:",
      paste(bad, collapse = ", "), "\n")
  quit(status = 5)
}

# --- plain timings -----------------------------------------------------------
t_call <- vapply(1:3, function(i) system.time(run_fixture_call(fx))[["elapsed"]], numeric(1))
kv(call_s_rep = t_call, call_s_median = median(t_call))

# --- Rprof split -------------------------------------------------------------
# Functions whose total time is reported.  "coxph" is the formula interface;
# "coxph.fit" is the fit proper it delegates to (the lean path of candidate A
# calls it directly).  Their difference is the formula-interface machinery:
# model.frame / model.matrix / Surv / aeqSurv / concordancefit / object assembly.
watch <- c("forestsearch", "subgroup.search", "search_combinations_parallel",
           "evaluate_combination_with_status", "fit_cox_for_subgroup",
           "coxph", "coxph.fit", "agreg.fit", "model.frame.default",
           "model.matrix.default", "Surv", "aeqSurv", "concordancefit",
           "concordance", "survfit", "survfit.formula", "summary.survfit",
           "subgroup.consistency", "evaluate_consistency_twostage",
           "evaluate_subgroup_consistency", "consistency_resample",
           ".consistency_cox_pieces", "residuals.coxph",
           "remove_near_duplicate_subgroups", "sort_subgroups_preview",
           "get_FSdata", "dummy", "get_dfpred", "merge", "future_lapply")

prof_one <- function(i) {
  f <- out_path(sprintf("Rprof_call_%d.out", i))
  # the platform may raise a too-short interval (e.g. to 0.01 s); harmless
  suppressWarnings(Rprof(f, interval = 0.002))
  run_fixture_call(fx)
  Rprof(NULL)
  sp <- summaryRprof(f)
  bt <- sp$by.total
  # Rprof labels namespaced calls "pkg::fun"; strip the prefix so "survival::coxph"
  # and "coxph" are the same row.  Functions reached through do.call() show as
  # "<Anonymous>" and are not attributable by name (forestsearch(), get_FSdata(),
  # subgroup.consistency()) -- their children are.
  bt$fun <- sub("^[A-Za-z0-9._]+::", "", gsub('"', "", rownames(bt)))
  bt <- aggregate(cbind(total.time, total.pct, self.time, self.pct) ~ fun, data = bt, FUN = sum)
  bt$rep <- i
  bt$sampling_time_s <- sp$sampling.time
  bt
}
prof <- do.call(rbind, lapply(1:3, prof_one))
utils::write.csv(prof, out_path("rprof_by_total_all.csv"), row.names = FALSE)

# median over reps for the watched functions
w <- prof[prof$fun %in% watch, c("fun", "total.time", "total.pct", "rep")]
med <- aggregate(cbind(total.time, total.pct) ~ fun, data = w, FUN = median)
med <- med[order(-med$total.time), ]
med$fun <- factor(med$fun, levels = watch)
med <- med[order(med$fun), ]; med$fun <- as.character(med$fun)
utils::write.csv(med, out_path("rprof_watched_median.csv"), row.names = FALSE)
cat("\nRprof by.total, median of 3 profiles (seconds, % of samples):\n")
print(med, row.names = FALSE)

get_t <- function(f) { v <- med$total.time[med$fun == f]; if (length(v)) v else 0 }
cox_total <- get_t("coxph")
fit_total <- get_t("coxph.fit") + get_t("agreg.fit")
kv(rprof_sampling_s_median   = median(prof$sampling_time_s[!duplicated(prof$rep)]),
   rprof_coxph_total_s       = cox_total,
   rprof_coxph_fit_total_s   = fit_total,
   rprof_formula_machinery_s = cox_total - fit_total,
   rprof_machinery_share_of_coxph = if (cox_total > 0) (cox_total - fit_total) / cox_total else NA,
   rprof_fit_cox_for_subgroup_s = get_t("fit_cox_for_subgroup"),
   rprof_search_stage_s      = get_t("search_combinations_parallel"),
   rprof_survfit_s           = get_t("survfit"),
   rprof_concordancefit_s    = get_t("concordancefit"),
   rprof_consistency_stage_s = get_t("evaluate_consistency_twostage") + get_t("evaluate_subgroup_consistency"),
   rprof_consistency_cox_pieces_s = get_t(".consistency_cox_pieces"),
   rprof_residuals_coxph_s   = get_t("residuals.coxph"))
save_result(list(t_warm = t_warm, t_call = t_call, counts = cc, prof = prof, med = med), "s1_profile")
cat("\ns1 complete\n")
