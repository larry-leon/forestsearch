# b1_profiles.R -- per fixture (survival V, continuous, binary): one warm-up
# call, candidate-flow counts, 3 timed sequential calls, 3 Rprof profiles
# bucketed by stage, and a sequential nb = 5 bootstrap for the per-replicate
# denominator.  Nothing is modified; the installed package is called as a
# user would.

source("stage0b_common.R")

watch <- c(
  # search stage
  "search_combinations_parallel", "evaluate_combination_with_status",
  "fit_cox_for_subgroup", "fit_glm_for_subgroup", "fit_glm_for_subgroup_fast",
  ".fs_lm_two_group", ".estimate_or", ".estimate_rd", ".estimate_rr",
  "coxph", "coxph.fit", "glm", "glm.fit", "lm", "lm.fit",
  "model.frame.default", "model.matrix.default", "summary.glm", "vcov", "concordancefit",
  "survfit", "survfit.formula",
  # consistency stage
  "evaluate_consistency_twostage", "evaluate_subgroup_consistency",
  "consistency_resample", ".consistency_glm_pieces", ".consistency_cox_pieces",
  ".dfbeta_glm", "lm.influence", "residuals.coxph",
  "as.data.table", "[.data.table", "FS_labels",
  # selection / dedup
  "remove_near_duplicate_subgroups", ".maxeff_membership_dedup",
  "sort_subgroups_preview", "sort_subgroups",
  # screens / cuts
  "get_FSdata", "lasso_selection", "glmnet", "cv.glmnet",
  "causal_forest", "causal_survival_forest", "fit_causal_forest_glm",
  "variable_importance", "policy_tree", "dummy", "get_dfpred", "merge")

prof_one <- function(fx, tag, i) {
  f <- out_path(sprintf("Rprof_%s_%d.out", tag, i))
  suppressWarnings(Rprof(f, interval = 0.002))
  run_fixture_call(fx)
  Rprof(NULL)
  sp <- summaryRprof(f)
  bt <- sp$by.total
  bt$fun <- sub("^[A-Za-z0-9._]+::", "", gsub('"', "", rownames(bt)))
  bt <- aggregate(cbind(total.time, total.pct, self.time, self.pct) ~ fun, data = bt, FUN = sum)
  bt$rep <- i; bt$fixture <- tag; bt$sampling_time_s <- sp$sampling.time
  bt
}

all_med <- list(); all_prof <- list()
for (tag in names(FIXTURES)) {
  cat("\n==========", tag, "==========\n")
  fx <- FIXTURES[[tag]]()
  t_warm <- system.time(fs <- run_fixture_call(fx))[["elapsed"]]
  cc <- candidate_counts(fs)
  do.call(kv, setNames(c(list(t_warm), cc), paste0(tag, "_", c("call_warmup_s", names(cc)))))
  if (tag == "survival") {
    # GATE (as in Stage 0): the survival candidate flow is fixed by the data and
    # the 0.4.0 code, not by the machine.
    expected <- list(n_reached_fit = 1410L, n_passed_effect = 121L,
                     n_consistency_total = 120L, sg_harm = "{er <= 0} & {pgr <= 26}")
    bad <- names(expected)[!mapply(function(k, v) identical(as.character(cc[[k]]), as.character(v)),
                                   names(expected), expected)]
    if (length(bad)) {
      cat("GATE FAILED: survival candidate flow differs from CRAN 0.4.0 on fixture V for:",
          paste(bad, collapse = ", "), "\n")
      quit(status = 5)
    }
  }

  t_call <- vapply(1:3, function(i) system.time(run_fixture_call(fx))[["elapsed"]], numeric(1))
  do.call(kv, setNames(list(t_call, median(t_call)), paste0(tag, c("_call_s_rep", "_call_s_median"))))

  prof <- do.call(rbind, lapply(1:3, function(i) prof_one(fx, tag, i)))
  all_prof[[tag]] <- prof
  w <- prof[prof$fun %in% watch, c("fun", "total.time", "total.pct", "rep")]
  med <- aggregate(cbind(total.time, total.pct) ~ fun, data = w, FUN = median)
  med <- med[order(-med$total.time), ]
  all_med[[tag]] <- cbind(fixture = tag, med)
  cat("\nRprof by.total, median of 3 (s):\n"); print(med, row.names = FALSE)
  get_t <- function(f) { v <- med$total.time[med$fun == f]; if (length(v)) v else 0 }
  do.call(kv, setNames(list(
    median(prof$sampling_time_s[!duplicated(prof$rep)]),
    get_t("search_combinations_parallel"),
    get_t("fit_cox_for_subgroup") + get_t("fit_glm_for_subgroup") + get_t("fit_glm_for_subgroup_fast"),
    get_t("evaluate_consistency_twostage") + get_t("evaluate_subgroup_consistency"),
    get_t("consistency_resample"),
    get_t(".consistency_glm_pieces") + get_t(".consistency_cox_pieces"),
    get_t("remove_near_duplicate_subgroups") + get_t(".maxeff_membership_dedup") + get_t("sort_subgroups_preview"),
    get_t("causal_forest") + get_t("causal_survival_forest") + get_t("fit_causal_forest_glm") + get_t("policy_tree"),
    get_t("lasso_selection"), get_t("get_FSdata"), get_t("survfit"),
    get_t("coxph"), get_t("coxph.fit"), get_t("glm"), get_t("glm.fit"), get_t("lm"), get_t("lm.fit"),
    get_t("lm.influence"), get_t(".dfbeta_glm"), get_t("model.frame.default"), get_t("model.matrix.default"),
    get_t("as.data.table")),
    paste0(tag, "_rprof_", c("sampling_s", "search_stage_s", "search_fit_s", "consistency_stage_s",
                             "consistency_resample_s", "consistency_pieces_s", "dedup_sort_s", "grf_s",
                             "lasso_s", "get_FSdata_s", "survfit_s", "coxph_s", "coxph_fit_s", "glm_s", "glm_fit_s",
                             "lm_s", "lm_fit_s", "lm_influence_s", "dfbeta_glm_s", "model_frame_s", "model_matrix_s",
                             "as_data_table_s"))))

  # sequential replicate denominator: nb = 5 bootstrap, plan sequential
  if (!is.null(fs$sg.harm)) {
    t0 <- proc.time()[["elapsed"]]
    b <- suppressWarnings(suppressMessages(forestsearch::forestsearch_bootstrap_dofuture(
      fs.est = fs, nb_boots = 5L, seed = 8316951L,
      parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
    t_outer <- proc.time()[["elapsed"]] - t0
    future::plan("sequential")
    rep_s <- b$results$tmins_iteration * 60
    do.call(kv, setNames(list(t_outer, mean(rep_s), sum(b$results$any_found == 1L)),
                         paste0(tag, c("_boot5_outer_s", "_rep_s_mean", "_boot5_n_found"))))
  } else {
    do.call(kv, setNames(list(NA, NA, NA), paste0(tag, c("_boot5_outer_s", "_rep_s_mean", "_boot5_n_found"))))
    cat("NOTE: no subgroup identified on this fixture; bootstrap denominator unavailable\n")
  }
}
utils::write.csv(do.call(rbind, all_prof), out_path("rprof_by_total_all.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, all_med), out_path("rprof_watched_median.csv"), row.names = FALSE)
cat("\nb1 complete\n")
