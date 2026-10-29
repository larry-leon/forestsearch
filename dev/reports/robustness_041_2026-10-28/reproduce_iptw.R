library(forestsearch)
# (a) closure level: an explicit sw column, the slow continuous closure built as
#     forestsearch() builds it for ps_adjust_method = "iptw"
set.seed(11)
n <- 200L
d <- data.frame(treat = rep(0:1, each = n / 2))
d$y  <- 1 + 0.5 * d$treat + stats::rnorm(n)
ps   <- stats::runif(n, 0.2, 0.8)
d$sw <- ifelse(d$treat == 1L, mean(d$treat) / ps, (1 - mean(d$treat)) / (1 - ps))
est <- make_effect_estimator(outcome_type = "continuous", treat.name = "treat",
                             outcome.name = "y", ps_adjust_method = "iptw")
r <- est(d)
cat(sprintf("closure: estimate = %s, se = %s, converged = %s, method_used = %s\n",
            format(r$estimate), format(r$se), r$converged, r$method_used))
cat(sprintf("%s: closure-level weighted fit returns NA (the defect is present)\n",
            if (is.na(r$estimate)) "PASS" else "FAIL"))
# (b) end to end: the same call shape as tests/testthat/test-vi-deprecation.R's
#     .run(), with is.RCT = FALSE, ps_method = "logistic", ps_adjust_method = "iptw";
#     numeric confounders for the propensity model.
source("tests/testthat/helper-synthetic-dgm.R")
df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))
df$sex <- as.numeric(as.character(df$sex))
rec <- new.env(); rec$names <- NULL
suppressMessages(trace("fit_glm_for_subgroup",
  tracer = bquote(if (is.null(.(rec)$names)) assign("names", names(df_clean), envir = .(rec))),
  where = asNamespace("forestsearch"), print = FALSE))
fit <- suppressWarnings(suppressMessages(forestsearch(
  df.analysis = df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
  outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
  outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
  effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
  fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
  sg_focus = "eff", is.RCT = FALSE, ps_method = "logistic", ps_adjust_method = "iptw",
  details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"))))
suppressMessages(untrace("fit_glm_for_subgroup", where = asNamespace("forestsearch")))
tab <- fit$find.grps$out.found$hr.subgroups
cat(sprintf("end to end: sw in df_clean = %s; fitted candidates = %d; sg.harm = %s\n",
            "sw" %in% rec$names, if (is.null(tab)) 0L else nrow(tab),
            paste(fit$sg.harm, collapse = " & ")))
