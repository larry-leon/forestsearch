
# Testing
dgm <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "pgr"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),   # er <= 25th percentile
    pgr = list(type = "quantile", value = 0.50)   # pgr <= median
  ),
  model = "alt"
)



# Your original request - using quantiles
dgm <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),   # er <= 25th percentile
    meno = 0
  ),
  model = "alt",
  k_inter = 0.0
)


result <- find_k_inter_for_target_hr(
  target_hr_harm = 2.0,
  data = gbsg,
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),
    meno = 0
  ),
  k_treat = 1.0
)

# Result: k_inter = 1.3009 achieves HR_harm = 2.0



base_params <- list(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),
    meno = 0
  ),
  k_treat = 1.0,
  n_super = 5000  # Using smaller for faster demonstration
)


sensitivity_results <- do.call(sensitivity_analysis_k_inter, c(
  list(
    k_inter_range = c(-1.5, 1.5),
    n_points = 11,
    model = "alt"
  ),
  base_params
))

cat("\nSensitivity results:\n")
print(round(sensitivity_results, 3))



df_gbsg <- gbsg
df_gbsg$tte <- with(gbsg, rfstime/30.4375)

dgm_null <- generate_aft_dgm_flex(
  data = df_gbsg,
  n_super = 5000,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "tte",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),   # er <= 25th percentile
    meno = 0
  ),
  model = "alt",
  k_inter = 0.0
)

# Testing simulation module

# Drawing treatment 1:1
draw_null <- simulate_from_dgm(dgm = dgm_null, n = 700, rand_ratio = 1, draw_treatment = TRUE, max_follow = Inf, seed = 123)
with(draw_null, mean(treat_sim))

# Retaining original study treatment assignment
draw_null <- simulate_from_dgm(dgm = dgm_null, n = 700, draw_treatment = FALSE, max_follow = Inf, seed = 123)
with(draw_null, mean(treat_sim))


library(weightedsurv)

dfcount <- df_counting(df = draw_null, tte.name = "y_sim", event.name = "event_sim", treat.name = "treat_sim")

par(mfrow=c(1,1))

plot_weighted_km(dfcount, conf.int=TRUE, show.logrank = TRUE, ymax = 1.05)
title(main="Simulated data")

library(forestsearch)

names(draw_null)

confounders.name <- c("z_age","z_er","z_pgr", "z_meno", "z_grade_1", "z_grade_2", "size", "nodes")

library(doFuture)
library(doRNG)

registerDoFuture()
registerDoRNG()

# conf_force 'forces' a specific covariate cut, eg. 'age <= 65'

system.time({fs <- forestsearch(draw_null,  confounders.name=confounders.name,
                                outcome.name = "y_sim", treat.name = "treat_sim", event.name = "event_sim", id.name = "id",
                                hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                                sg_focus = "hrMaxSG",
                                showten_subgroups = FALSE, details=TRUE,
                                conf_force = c("z_age <= 65", "z_er <= 0", "z_er <= 1", "z_er <= 2","z_er <= 5"),
                                cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                                maxk = 2, n.min = 60, d0.min = 12, d1.min = 12,
                                plot.sg = TRUE, by.risk = 12,
                                parallel_args = list(plan="multisession", workers = 12, show_message = TRUE)
)
})




dgm_alt <- generate_aft_dgm_flex(
  data = df_gbsg,
  n_super = 5000,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "tte",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),   # er <= 25th percentile
    meno = 0
  ),
  model = "alt",
  k_inter = 1.3009
)

# Testing simulation module

draw_alt <- simulate_from_dgm(dgm = dgm_alt, n = 700, rand_ratio = 1, max_follow = Inf, seed = 123)

dfcount <- df_counting(df = draw_alt, tte.name = "y_sim", event.name = "event_sim", treat.name = "treat_sim")

par(mfrow=c(1,1))

plot_weighted_km(dfcount, conf.int=TRUE, show.logrank = TRUE, ymax = 1.05)
title(main="Simulated data")


system.time({fs <- forestsearch(draw_alt,  confounders.name=confounders.name,
                                outcome.name = "y_sim", treat.name = "treat_sim", event.name = "event_sim", id.name = "id",
                                potentialOutcome.name = "loghr_po",
                                hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                                sg_focus = "hrMaxSG",
                                showten_subgroups = FALSE, details=TRUE,
                                conf_force = c("z_age <= 65", "z_er <= 0", "z_er <= 1", "z_er <= 2","z_er <= 5"),
                                cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                                maxk = 2, n.min = 60, d0.min = 12, d1.min = 12,
                                plot.sg = TRUE, by.risk = 12,
                                parallel_args = list(plan="callr", workers = 12, show_message = TRUE)
)
})


res_tabs <- sg_tables(fs, ndecimals = 3)

res_tabs$sg10_out
res_tabs$tab_estimates

plan("sequential")


NB <- 30

t.start <- proc.time()[3]

fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, show_three = FALSE, details = TRUE)

t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60
cat("Minutes (total) for bootstrap (boots,mins)",c(NB,t.min),"\n")
cat("Projected minutes for 1000",c(t.min*(1000/NB)),"\n")


sg_tab <- fs_bc$summary$table
sg_tab

plan("sequential")

