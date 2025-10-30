# Main function becomes orchestrator
draw_sim_stratified <- function(dgm, ss = 1, Ndraw = nrow(dgm$df_super),
                                strata_rand = "stratum", wname = "ecog",
                                bw = 0, checking = FALSE, details = FALSE,
                                return_df = TRUE, time_eos = Inf,
                                keep_rand = FALSE) {

  # 1. Prepare data
  sim_data <- prepare_simulation_data(dgm, Ndraw, ss)

  # 2. Generate potential outcomes
  po_data <- generate_potential_outcomes(
    df = sim_data,
    gamma_true = dgm$gamma.true,
    mu = dgm$mu,
    tau_strata = sim_data$tau.strataO,
    gamma_w = calculate_gamma_w(bw, dgm$tau.approx),
    wname = wname
  )

  # 3. Apply randomization
  obs_data <- apply_randomization_treatment(
    po_data = po_data,
    strata_rand = strata_rand,
    keep_rand = keep_rand
  )

  # 4. Apply censoring
  obs_data <- apply_censoring(
    obs_data = obs_data,
    muC = dgm$muC,
    tauC = dgm$tauC,
    time_eos = time_eos,
    seed = ss
  )

  # 5. Validation (optional)
  if (checking) {
    validate_simulation_parameters(obs_data, dgm)
  }

  # 6. Reporting (optional)
  if (details && ss <= 10) {
    report_simulation_metrics(obs_data)
  }

  # 7. Format output
  return(format_simulation_output(obs_data, return_df))
}

# New modular functions:

prepare_simulation_data <- function(dgm, Ndraw, seed) {
  set.seed(8316951 + seed * 1000)

  df_super <- dgm$df_super

  if (Ndraw != nrow(df_super)) {
    id_sample <- sample(seq_len(nrow(df_super)), size = Ndraw, replace = TRUE)
    df_super <- df_super[id_sample, ]
  }

  return(df_super)
}

generate_potential_outcomes <- function(df, gamma_true, mu, tau_strata, gamma_w, wname) {
  N <- nrow(df)

  # Build covariate matrices
  zmat_1 <- build_covariate_matrix(df, treat = 1)
  zmat_0 <- build_covariate_matrix(df, treat = 0)

  w <- df[[wname]]
  epsilon <- log(rexp(N))

  # Calculate potential outcomes
  df$log.Y1 <- mu + c(zmat_1 %*% gamma_true) + w * gamma_w + tau_strata * epsilon
  df$log.Y0 <- mu + c(zmat_0 %*% gamma_true) + w * gamma_w + tau_strata * epsilon

  # Calculate treatment effects
  df$loghr.po <- (-1) * c((zmat_1 - zmat_0) %*% gamma_true) / tau_strata
  df$theta1.po <- -c(zmat_1 %*% gamma_true + w * gamma_w) / tau_strata
  df$theta0.po <- -c(zmat_0 %*% gamma_true + w * gamma_w) / tau_strata

  return(df)
}

build_covariate_matrix <- function(df, treat) {
  zmat <- as.matrix(df[, c("treat", "z", "z.treat", "z.k", "z.k.treat")])
  zmat[, "treat"] <- treat
  zmat[, "z.treat"] <- treat * df$z
  zmat[, "z.k.treat"] <- treat * zmat[, "z.k"]
  return(zmat)
}

apply_randomization_treatment <- function(po_data, strata_rand, keep_rand) {
  if (!keep_rand) {
    blocks <- po_data[[strata_rand]]
    po_data$treat.sim <- randomizr::block_ra(blocks = blocks)
  } else {
    po_data$treat.sim <- ifelse(po_data$treat == 1, 1, 0)
  }

  # Observed outcomes based on treatment assignment
  po_data$y.sim <- exp(
    po_data$treat.sim * po_data$log.Y1 +
      (1 - po_data$treat.sim) * po_data$log.Y0
  )

  return(po_data)
}

apply_censoring <- function(obs_data, muC, tauC, time_eos, seed) {
  set.seed(8316951 + seed * 1000 + 100)

  N <- nrow(obs_data)
  epsilonC <- log(rexp(N))
  censor_time <- exp(muC + tauC * epsilonC)
  censor_time <- pmin(censor_time, time_eos)

  obs_data$event.sim <- ifelse(obs_data$y.sim <= censor_time, 1, 0)
  obs_data$y.sim <- pmin(obs_data$y.sim, censor_time)

  return(obs_data)
}

format_simulation_output <- function(obs_data, return_df) {
  obs_data <- data.table::setorder(obs_data, z)

  if (!return_df) {
    return(calculate_population_summaries(obs_data))
  }

  return(obs_data)
}

calculate_population_summaries <- function(df) {
  list(
    AHR = exp(mean(df$loghr.po)),
    AHR_W1 = with(subset(df, w == 1), exp(mean(loghr.po))),
    AHR_W0 = with(subset(df, w == 0), exp(mean(loghr.po))),
    CDE = mean(exp(df$theta1.po)) / mean(exp(df$theta0.po)),
    # ... other summaries
    zpoints = seq(min(df$z), max(df$z), by = 1),
    HR.zpoints = calculate_ahr_by_threshold(df, ">="),
    HRminus.zpoints = calculate_ahr_by_threshold(df, "<=")
  )
}


get_SGanalyses <- function(dgm, subgroups_name, subgroups_id, sims,
                           wname, bw, Ndraw = nrow(dgm$df_super),
                           bmcut = log(1), outfile = NULL) {

  # Validate inputs
  validate_subgroup_inputs(subgroups_name, subgroups_id)

  # Define analysis specifications
  analyses <- define_cox_analyses(wname)

  # Initialize storage
  results <- initialize_results_storage(sims, subgroups_name, analyses)

  # Run simulations
  for (ss in seq_len(sims)) {
    df_sim <- generate_simulation_dataset(dgm, ss, Ndraw, wname, bw, bmcut)

    # Fit all analyses for all subgroups
    for (gg in seq_along(subgroups_id)) {
      df_sg <- subset_by_expression(df_sim, subgroups_id[gg])
      results <- fit_all_analyses(df_sg, analyses, results, ss, gg)
    }

    # Optional: Plot first 10
    if (ss <= 10) {
      plot_km_curve(df_sim)
    }
  }

  # Save and return
  if (!is.null(outfile)) {
    save_results(results, dgm, subgroups_name, outfile)
  }

  return(results)
}

# New helper functions:

define_cox_analyses <- function(wname) {
  list(
    list(
      name = "sR",
      formula = "Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR)"
    ),
    list(
      name = "none",
      formula = "Surv(y.sim, event.sim) ~ treat.sim"
    ),
    list(
      name = "sW",
      formula = sprintf("Surv(y.sim, event.sim) ~ treat.sim + strata(%s)", wname)
    ),
    list(
      name = "sBM",
      formula = "Surv(y.sim, event.sim) ~ treat.sim + strata(cut_opt)",
      min_strata_size = 5
    ),
    list(
      name = "sW+sR",
      formula = sprintf("Surv(y.sim, event.sim) ~ treat.sim + strata(%s) + strata(strata.simR)", wname)
    ),
    list(
      name = "W+age+sR",
      formula = sprintf("Surv(y.sim, event.sim) ~ treat.sim + %s + age + strata(strata.simR)", wname)
    )
  )
}

fit_all_analyses <- function(df_sg, analyses, results, sim_idx, sg_idx) {
  results$ns[sim_idx, sg_idx] <- nrow(df_sg)

  for (i in seq_along(analyses)) {
    analysis <- analyses[[i]]

    # Check if analysis is feasible
    if (!is_analysis_feasible(df_sg, analysis)) {
      next
    }

    # Fit Cox model
    fit_result <- fit_cox_safely(analysis$formula, df_sg)

    if (!is.null(fit_result)) {
      hr_col <- paste0("hrs", i)
      ub_col <- paste0("ubs", i)

      results[[hr_col]][sim_idx, sg_idx] <- fit_result$hr
      results[[ub_col]][sim_idx, sg_idx] <- fit_result$upper_bound
    }
  }

  return(results)
}

fit_cox_safely <- function(formula_str, data) {
  tryCatch({
    fit <- coxph(as.formula(formula_str), data = data)
    summary_fit <- summary(fit)

    list(
      hr = summary_fit$conf.int[1, 1],
      upper_bound = summary_fit$conf.int[1, 4],
      lower_bound = summary_fit$conf.int[1, 3]
    )
  }, error = function(e) {
    NULL
  })
}

is_analysis_feasible <- function(df, analysis) {
  if (!is.null(analysis$min_strata_size)) {
    # Check stratification variable has enough observations
    strata_var <- extract_strata_var(analysis$formula)
    if (!is.null(strata_var)) {
      strata_counts <- table(df[[strata_var]])
      if (any(strata_counts < analysis$min_strata_size)) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}


# Separate computation from visualization
cox_cs_fit <- function(df, tte.name = "os_time", event.name = "os_event",
                       treat.name = "treat", strata.name = NULL,
                       z.name = "bm", alpha = 0.20, max_z = Inf,
                       z_grid = seq(min_z, max_z, by = 1)) {

  # 1. Prepare data
  model_data <- prepare_spline_data(df, tte.name, event.name, treat.name,
                                    strata.name, z.name)

  # 2. Fit model
  fit <- fit_spline_cox(model_data, strata.name)

  # 3. Generate predictions
  predictions <- predict_treatment_effects(fit, z_grid, model_data, alpha)

  # 4. Calculate counts
  counts <- calculate_z_counts(model_data$z, z_grid, window = 0.25)

  structure(
    list(
      fit = fit,
      predictions = predictions,
      counts = counts,
      data = model_data
    ),
    class = "cox_spline_fit"
  )
}

# Separate plotting function
plot.cox_spline_fit <- function(x, show_true_beta = FALSE,
                                true_beta = NULL, ...) {

  pred <- x$predictions
  counts <- x$counts

  # Set up plot region
  ylim <- calculate_plot_limits(pred, true_beta)

  # Main plot
  plot(pred$z_grid, pred$log_hr,
       type = "l", lwd = 3, col = "black",
       ylim = ylim, xlab = "Biomarker", ylab = "log(HR)",
       axes = FALSE, ...)

  # Add confidence bands
  lines(pred$z_grid, pred$lower, lty = 2, col = "black")
  lines(pred$z_grid, pred$upper, lty = 2, col = "black")

  # Add reference lines
  add_reference_lines(pred$standard_cox_hr)

  # Add true beta if provided
  if (show_true_beta && !is.null(true_beta)) {
    lines(x$data$z, true_beta, col = "orange", lwd = 2)
  }

  # Add counts
  add_count_labels(pred$z_grid, counts, ylim)

  # Finalize
  add_axes_and_legend(pred$z_grid, ylim)
}

# Core computation functions
fit_spline_cox <- function(model_data, strata.name) {
  # Create spline basis
  z_basis <- splines::ns(model_data$z, df = 3)

  # Build model matrix
  X <- build_interaction_matrix(z_basis, model_data$treat)

  # Fit Cox model
  formula <- build_cox_formula(strata.name)
  fit <- survival::coxph(formula, data = model_data)

  list(
    coefficients = coef(fit),
    vcov = vcov(fit),
    basis = z_basis
  )
}

predict_treatment_effects <- function(fit, z_grid, model_data, alpha) {
  # Generate basis for prediction points
  z_basis_pred <- predict(fit$basis, z_grid)

  # Calculate log(HR) and standard errors
  X_pred <- build_prediction_matrix(z_basis_pred)
  log_hr <- c(X_pred %*% fit$coefficients[treatment_indices(fit)])

  # Calculate standard errors efficiently
  se_log_hr <- calculate_treatment_effect_se(X_pred, fit$vcov)

  # Confidence intervals
  c_alpha <- qnorm(1 - alpha / 2)

  data.frame(
    z_grid = z_grid,
    log_hr = log_hr,
    se = se_log_hr,
    lower = log_hr - c_alpha * se_log_hr,
    upper = log_hr + c_alpha * se_log_hr,
    standard_cox_hr = calculate_standard_cox_hr(model_data)
  )
}

calculate_treatment_effect_se <- function(X_pred, vcov_full) {
  # Extract only treatment-related covariance
  treatment_idx <- c(1, 5:7)  # Indices for treatment effects
  vcov_treatment <- vcov_full[treatment_idx, treatment_idx]

  # Efficient variance calculation
  sqrt(rowSums((X_pred %*% vcov_treatment) * X_pred))
}
