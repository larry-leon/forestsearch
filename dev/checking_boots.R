check_boot <- subset(fs_bc$results, is.na(H_biasadj_1))

check15 <- subset(check_boot, boot_id == 15)

check15


# Return boot = 515 dataset

df_badboot <- bootstrap_reproduce_aboot(this_boot = 1962, fs.est = fs, cox.formula.boot = cox.formula.boot)


# WithOUT parallel processing
fs_bad <- forestsearch(df_badboot, confounders.name=confounders.name,
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


