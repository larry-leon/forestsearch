# run_simulations.R
# Run this ONCE from the same directory as the .qmd
# Usage:  Rscript run_simulations.R
#    or:  source("run_simulations.R")  from the R console

library(forestsearch)
library(survival)
source("summarize_extreme_sims.R")   # for save_sim_results()

# ── Settings ────────────────────────────────────────────────────
n_sims_null <- 100L
hr_true     <- 0.70
output_file <- "results/extreme_sims_5000.rds"

# ── Data prep ───────────────────────────────────────────────────
data(gbsg, package = "survival")
gbsg$time_months <- gbsg$rfstime / 30.4375
gbsg$treat       <- gbsg$hormon

# ── Build null DGM ──────────────────────────────────────────────
set.seed(99)
dgm_null <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age","size","nodes","pgr","er"),
  factor_vars = c("meno","grade"),
  outcome_var = "time_months", event_var = "status",
  treatment_var = "hormon", subgroup_vars = NULL,
  model = "null", n_super = 5000, seed = 99, verbose = TRUE)

# ── Calibrate censoring ────────────────────────────────────────
cal_null <- calibrate_cens_adjust(
  dgm = dgm_null, target = "rate", n = 800,
  analysis_time = 84, max_entry = 24,
  seed = 42, n_eval = 1500, verbose = TRUE)

# ── Define subgroups (same list as in the .qmd) ────────────────
gbsg$flag_itt <- 1L
age_med  <- median(gbsg$age);    size_med <- median(gbsg$size)
nodes_q3 <- quantile(gbsg$nodes, 0.75)
pgr_med  <- median(gbsg$pgr);   er_cut   <- 20

subgroups <- list(
  # ... identical list from the .qmd ...
  list(id="flag_itt == 1", name="All Patients", grp="ITT"),
  list(id="meno == 1", name="Post-menopausal", grp="Clinical"),
  # ... all 56 subgroups ...
  list(id="random15==1", name="random15", grp="Random (N~15)")
)

# ── Cox helper ──────────────────────────────────────────────────
cox_sR <- function(data) {
  fit <- tryCatch(
    coxph(Surv(y_sim, event_sim) ~ treat_sim + strata(grade), data = data),
    error = function(e) NULL)
  if (is.null(fit) || nrow(summary(fit)$conf.int) == 0)
    return(c(NA_real_, NA_real_))
  ci <- summary(fit)$conf.int[1, ]
  c(ci["exp(coef)"], ci["upper .95"])
}

# ── Simulation loop ────────────────────────────────────────────
sg_names <- sapply(subgroups, `[[`, "name")
sg_ids   <- sapply(subgroups, `[[`, "id")
n_sg     <- length(subgroups)

sim_hrs <- matrix(NA_real_, n_sims_null, n_sg, dimnames=list(NULL, sg_names))
sim_ubs <- matrix(NA_real_, n_sims_null, n_sg, dimnames=list(NULL, sg_names))
sim_ns  <- matrix(NA_real_, n_sims_null, n_sg, dimnames=list(NULL, sg_names))

t0 <- Sys.time()
for (ss in seq_len(n_sims_null)) {
  if (ss %% 100 == 0) message("  sim ", ss, " / ", n_sims_null,
                              "  (", round(difftime(Sys.time(), t0, units="mins"), 1), " min)")

  df_s <- simulate_from_dgm(dgm=dgm_null, n=nrow(gbsg),
                            analysis_time=84, max_entry=24,
                            cens_adjust=cal_null$cens_adjust, seed=ss)
  df_s$flag_itt <- 1L
  df_s$age_med <- age_med; df_s$size_med <- size_med
  df_s$nodes_q3 <- nodes_q3; df_s$pgr_med <- pgr_med; df_s$er_cut <- er_cut

  set.seed(ss + 1e6L)
  n_s <- nrow(df_s)
  r_idx <- sample.int(n_s, min(60L, n_s), replace=FALSE)
  df_s$random60 <- as.integer(seq_len(n_s) %in% r_idx)
  df_s$random40 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(40,length(r_idx)))])
  df_s$random20 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(20,length(r_idx)))])
  df_s$random15 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(15,length(r_idx)))])

  for (gg in seq_len(n_sg)) {
    df_sg <- tryCatch(subset(df_s, eval(parse(text=sg_ids[[gg]]))),
                      error=function(e) NULL)
    if (is.null(df_sg) || nrow(df_sg) < 5) next
    sim_ns[ss,gg] <- nrow(df_sg)
    r <- cox_sR(df_sg)
    sim_hrs[ss,gg] <- r[1]; sim_ubs[ss,gg] <- r[2]
  }
}
message("Done: ", n_sims_null, " sims in ",
        round(difftime(Sys.time(), t0, units="mins"), 1), " minutes")

# ── Save ────────────────────────────────────────────────────────
save_sim_results(sim_hrs, sim_ubs, sim_ns, subgroups,
                 file = output_file, hr_true = hr_true,
                 metadata = list(dgm_seed = 99,
                                 cens_adjust = cal_null$cens_adjust,
                                 n_per_trial = nrow(gbsg)))
