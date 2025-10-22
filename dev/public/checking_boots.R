check_boot <- subset(fs_bc$results, is.na(H_biasadj_1))

check <- subset(check_boot, boot_id == 515)

#check15


# Return boot = 515 dataset

bootstrap_reproduce_aboot <- function(this_boot, fs.est, cox.formula.boot) {

  parallel_args = list(NULL)

  args_forestsearch_call <- fs.est$args_call_all
  parallel_args <- resolve_bootstrap_parallel_args(parallel_args, args_forestsearch_call)
  cox.formula.boot <- do.call(build_cox_formula,filter_call_args(args_forestsearch_call, build_cox_formula))


  df_boot_analysis <- fs.est$df.est
  NN <- nrow(df_boot_analysis)
  id0 <- seq_len(NN)
  seq_boots <- seq_len(this_boot)

  set.seed(8316951)

  # Do not modify seed it needs to align with ystar
  foreach::foreach(
    boot = seq_boots,
    .options.future = list(
      seed = TRUE,
      add = get_bootstrap_exports()
    ),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {

    if(boot == this_boot){

      in_boot <- sample.int(NN, size = NN, replace = TRUE)
      df_boot <- df_boot_analysis[in_boot, ]
      df_boot$id_boot <- seq_len(nrow(df_boot))

      # =================================================================
      # Extract variable names from formula
      # =================================================================
      outcome_var <- all.vars(cox.formula.boot[[2]])[1]
      event_var <- all.vars(cox.formula.boot[[2]])[2]
      treat_var <- all.vars(cox.formula.boot[[3]])[1]

      # =================================================================
      # Bootstrap data evaluated at ORIGINAL subgroup H
      # =================================================================

      # Check events in subgroup H (treat.recommend == 0)
      df_H <- subset(df_boot, treat.recommend == 0)

      events_H_0 <- sum(df_H[df_H[[treat_var]] == 0, event_var], na.rm = TRUE)
      events_H_1 <- sum(df_H[df_H[[treat_var]] == 1, event_var], na.rm = TRUE)

      cat("H_0 and H_1 events",c(events_H_0,events_H_1),"\n")


      # Note that any identified subgroups are of size at least n.min with minimum
      # number of events in the control and experimental arms of d0.min and d1.min
      # In addition, for any bootstrap forestsearch analysis the idenfitied subgroups
      # satisfy the same criteria;
      # However, for the bootstrap sample, the observed df_H and df_Hc analyses
      # do not have any constraints

      system.time({fitH_star <- get_Cox_sg(
        df_sg = df_H,
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      })


      # Check events in subgroup Hc (treat.recommend == 1)
      df_Hc <- subset(df_boot, treat.recommend == 1)
      events_Hc_0 <- sum(df_Hc[df_Hc[[treat_var]] == 0, event_var], na.rm = TRUE)
      events_Hc_1 <- sum(df_Hc[df_Hc[[treat_var]] == 1, event_var], na.rm = TRUE)

      cat("Hc_0 and Hc_1 events",c(events_Hc_0,events_Hc_1),"\n")

      fitHc_star <- get_Cox_sg(
        df_sg = df_Hc,
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )

      # =================================================================
      # Prepare bootstrap dataframes - drop confounders and treat.recommend
      # =================================================================
      drop.vars <- c(fs.est$confounders.candidate, "treat.recommend")
      dfnew <- df_boot_analysis[, !(names(df_boot_analysis) %in% drop.vars)]
      dfnew_boot <- df_boot[, !(names(df_boot) %in% drop.vars)]

      return(dfnew_boot)

      # =================================================================
      # Configure forestsearch arguments for bootstrap
      # =================================================================
      args_FS_boot <- fs.est$args_call_all
      args_FS_boot$df.analysis <- dfnew_boot
      args_FS_boot$df.predict <- dfnew

      # CATEGORY 1: OUTPUT SUPPRESSION
      args_FS_boot$details <- TRUE
      args_FS_boot$showten_subgroups <- FALSE
      args_FS_boot$plot.sg <- FALSE
      args_FS_boot$plot.grf <- FALSE

      # CATEGORY 2: VARIABLE RE-SELECTION
      args_FS_boot$grf_res <- NULL
      args_FS_boot$grf_cuts <- NULL

      # CATEGORY 3: SEQUENTIAL EXECUTION
      args_FS_boot$parallel_args$plan <- "sequential"
      args_FS_boot$parallel_args$workers <- 1L
      args_FS_boot$parallel_args$show_message <- FALSE

      # =================================================================
      # Run forestsearch on bootstrap sample
      # =================================================================
      run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)

      if (inherits(run_bootstrap, "try-error")) {
        warning("Bootstrap ", boot, " failed: ", as.character(run_bootstrap))
      }

      # Extract prediction datasets from bootstrap ForestSearch run
      df_PredBoot <- run_bootstrap$df.predict
      dfboot_PredBoot <- run_bootstrap$df.est

      # ==============================================================
      # Check events in NEW subgroups found by bootstrap
      # ==============================================================

      # NEW subgroup H* on ORIGINAL data
      df_Hstar <- subset(df_PredBoot, treat.recommend == 0)

      events_Hstar_0 <- sum(df_Hstar[df_Hstar[[treat_var]] == 0, event_var], na.rm = TRUE)
      events_Hstar_1 <- sum(df_Hstar[df_Hstar[[treat_var]] == 1, event_var], na.rm = TRUE)

      # NEW subgroup Hc* on ORIGINAL data
      df_Hcstar <- subset(df_PredBoot, treat.recommend == 1)

      events_Hcstar_0 <- sum(df_Hcstar[df_Hcstar[[treat_var]] == 0, event_var], na.rm = TRUE)
      events_Hcstar_1 <- sum(df_Hcstar[df_Hcstar[[treat_var]] == 1, event_var], na.rm = TRUE)

      out <- data.table::data.table(events_Hcstar_0, events_Hstar_0,fitH_star$est_obs)

    }
  }
}

df_badboot <- bootstrap_reproduce_aboot(this_boot = 4, fs.est = fs, cox.formula.boot = cox.formula.boot)


# WithOUT parallel processing
system.time({fs_bad <- forestsearch(df_badboot, confounders.name=confounders.name,
                                outcome.name = "tte", treat.name = "treat", event.name = "event", id.name = "id",
                                hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                                sg_focus = "hrMaxSG",
                                showten_subgroups = FALSE, details=TRUE,
                                conf_force = c("age <= 65", "er <= 0", "er <= 1", "er <= 2"," er <= 5"),
                                cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                                maxk = 2, n.min = 60, d0.min = 10, d1.min = 10,
                                plot.sg = TRUE, by.risk = 12,
                       parallel_args = list()
)
})

# With parallel
system.time({fs_bad <- forestsearch(df_badboot, confounders.name=confounders.name,
                       outcome.name = "tte", treat.name = "treat", event.name = "event", id.name = "id",
                       hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                       sg_focus = "hrMaxSG", fs.splits = 1000,
                       showten_subgroups = FALSE, details=TRUE,
                       conf_force = c("age <= 65", "er <= 0", "er <= 1", "er <= 2"," er <= 5"),
                       cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                       maxk = 2, n.min = 60, d0.min = 10, d1.min = 10,
                       plot.sg = TRUE, by.risk = 12,
                       parallel_args = list(plan = "multisession", show_message = TRUE, workers = 15)
)
})


