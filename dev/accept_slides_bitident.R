#!/usr/bin/env Rscript
# =========================================================================
# SLIDES-DECK ACCEPTANCE: run_subgroup_sims() must reproduce the Beamer
# deck's former inline sequential loop BIT-IDENTICALLY (its own DGM: null
# model without k_treat calibration; its own cens_adjust calibration).
# All "verbatim" blocks are machine-extracted from the pre-port
# extreme_subgroups_slides.qmd -- regenerate rather than edit.
#
# Usage: Rscript dev/accept_slides_bitident.R [n_sims]     # default 50
# =========================================================================
suppressMessages({ library(forestsearch); library(survival) })
args  <- commandArgs(trailingOnly = TRUE)
n_acc <- if (length(args) >= 1) as.integer(args[1]) else 50L
.fails <- 0L
ok <- function(cond, msg) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
  if (!isTRUE(cond)) .fails <<- .fails + 1L
}
cat("forestsearch:", as.character(packageVersion("forestsearch")), "\n")
hr_true <- 0.70
n_sims_null <- n_acc

# ---- deck prelude: data, null DGM, censoring calibration, subgroups ----
# (verbatim from the pre-port run-or-load chunk)
data(gbsg, package = "survival")
gbsg$time_months <- gbsg$rfstime / 30.4375
gbsg$treat <- gbsg$hormon

set.seed(99)
dgm_null <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age","size","nodes","pgr","er"),
  factor_vars = c("meno","grade"),
  outcome_var = "time_months", event_var = "status",
  treatment_var = "hormon", subgroup_vars = NULL,
  model = "null", n_super = 5000, seed = 99, verbose = FALSE)

cal_null <- calibrate_cens_adjust(
  dgm = dgm_null, target = "rate", n = 800,
  analysis_time = 84, max_entry = 24,
  seed = 42, n_eval = 1500, verbose = FALSE)

gbsg$flag_itt <- 1L
age_med <- median(gbsg$age); size_med <- median(gbsg$size)
nodes_q3 <- quantile(gbsg$nodes, 0.75)
pgr_med <- median(gbsg$pgr); er_cut <- 20

subgroups <- list(
  list(id="flag_itt == 1", name="All Patients", grp="ITT"),
  list(id="meno == 1", name="Post-menopausal", grp="Clinical"),
  list(id="meno == 0", name="Pre-menopausal", grp="Clinical"),
  list(id="grade == 3", name="Grade 3", grp="Clinical"),
  list(id="grade != 3", name="Grade 1/2", grp="Clinical"),
  list(id="age <= age_med", name="Age (young)", grp="Continuous"),
  list(id="age > age_med", name="Age (older)", grp="Continuous"),
  list(id="age <= 50", name="Age <= 50", grp="Continuous"),
  list(id="age > 50", name="Age > 50", grp="Continuous"),
  list(id="size <= size_med", name="Tumour size (small)", grp="Continuous"),
  list(id="size > size_med", name="Tumour size (large)", grp="Continuous"),
  list(id="nodes == 0", name="Node-negative", grp="Continuous"),
  list(id="nodes > 0 & nodes <= 3", name="Nodes 1-3", grp="Continuous"),
  list(id="nodes > nodes_q3", name="High nodes (>Q3)", grp="Continuous"),
  list(id="pgr <= pgr_med", name="PGR low", grp="Continuous"),
  list(id="pgr > pgr_med", name="PGR high", grp="Continuous"),
  list(id="er <= er_cut", name="ER-low (<20)", grp="Continuous"),
  list(id="er > er_cut", name="ER-high (>=20)", grp="Continuous"),
  list(id="meno==0 & grade==3", name="Pre-meno / Grade 3", grp="Interaction: meno x grade"),
  list(id="meno==1 & grade==3", name="Post-meno / Grade 3", grp="Interaction: meno x grade"),
  list(id="meno==0 & grade!=3", name="Pre-meno / Grade 1-2", grp="Interaction: meno x grade"),
  list(id="meno==1 & grade!=3", name="Post-meno / Grade 1-2", grp="Interaction: meno x grade"),
  list(id="meno==0 & age<=50", name="Pre-meno / Age<=50", grp="Interaction: meno x age"),
  list(id="meno==0 & age>50", name="Pre-meno / Age>50", grp="Interaction: meno x age"),
  list(id="meno==1 & age<=age_med", name="Post-meno / Age(young)", grp="Interaction: meno x age"),
  list(id="meno==1 & age>age_med", name="Post-meno / Age(older)", grp="Interaction: meno x age"),
  list(id="meno==0 & er<=er_cut", name="Pre-meno / ER-low", grp="Interaction: meno x ER"),
  list(id="meno==0 & er>er_cut", name="Pre-meno / ER-high", grp="Interaction: meno x ER"),
  list(id="meno==1 & er<=er_cut", name="Post-meno / ER-low", grp="Interaction: meno x ER"),
  list(id="meno==1 & er>er_cut", name="Post-meno / ER-high", grp="Interaction: meno x ER"),
  list(id="grade==3 & nodes==0", name="Grade 3 / Node-neg", grp="Interaction: grade x nodes"),
  list(id="grade==3 & nodes>0", name="Grade 3 / Node-pos", grp="Interaction: grade x nodes"),
  list(id="grade!=3 & nodes==0", name="Grade 1-2 / Node-neg", grp="Interaction: grade x nodes"),
  list(id="grade!=3 & nodes>0", name="Grade 1-2 / Node-pos", grp="Interaction: grade x nodes"),
  list(id="grade==3 & nodes>nodes_q3", name="Grade 3 / High nodes", grp="Interaction: grade x nodes"),
  list(id="grade==3 & er<=er_cut", name="Grade 3 / ER-low", grp="Interaction: grade x ER"),
  list(id="grade==3 & er>er_cut", name="Grade 3 / ER-high", grp="Interaction: grade x ER"),
  list(id="grade!=3 & er<=er_cut", name="Grade 1-2 / ER-low", grp="Interaction: grade x ER"),
  list(id="grade!=3 & er>er_cut", name="Grade 1-2 / ER-high", grp="Interaction: grade x ER"),
  list(id="grade==3 & pgr<=pgr_med", name="Grade 3 / PGR low", grp="Interaction: grade x PGR"),
  list(id="grade==3 & pgr>pgr_med", name="Grade 3 / PGR high", grp="Interaction: grade x PGR"),
  list(id="grade!=3 & pgr<=pgr_med", name="Grade 1-2 / PGR low", grp="Interaction: grade x PGR"),
  list(id="grade!=3 & pgr>pgr_med", name="Grade 1-2 / PGR high", grp="Interaction: grade x PGR"),
  list(id="meno==0 & age<=50 & er<=er_cut", name="Pre-meno/Yng/ER-low", grp="3-way"),
  list(id="meno==0 & age<=50 & er>er_cut", name="Pre-meno/Yng/ER-high", grp="3-way"),
  list(id="meno==1 & grade==3 & er<=er_cut", name="Post-meno/G3/ER-low", grp="3-way"),
  list(id="meno==1 & grade==3 & nodes>0", name="Post-meno/G3/Node+", grp="3-way"),
  list(id="grade==3 & nodes>0 & er<=er_cut", name="G3/Node+/ER-low", grp="3-way"),
  list(id="size>size_med & nodes>0", name="Large/Node-pos", grp="Interaction: size x nodes"),
  list(id="size>size_med & nodes==0", name="Large/Node-neg", grp="Interaction: size x nodes"),
  list(id="size<=size_med & nodes>0", name="Small/Node-pos", grp="Interaction: size x nodes"),
  list(id="size<=size_med & nodes==0", name="Small/Node-neg", grp="Interaction: size x nodes"),
  list(id="random60==1", name="random60", grp="Random (N~60)"),
  list(id="random40==1", name="random40", grp="Random (N~40)"),
  list(id="random20==1", name="random20", grp="Random (N~20)"),
  list(id="random15==1", name="random15", grp="Random (N~15)"))

# ---- transcribed former inline loop (verbatim) --------------------------
cox_sR <- function(data) {
  fit <- tryCatch(
    coxph(Surv(y_sim, event_sim) ~ treat_sim + strata(grade), data = data),
    error = function(e) NULL)
  if (is.null(fit) || nrow(summary(fit)$conf.int) == 0)
    return(c(NA_real_, NA_real_))
  ci <- summary(fit)$conf.int[1, ]
  c(ci["exp(coef)"], ci["upper .95"])
}

sg_names <- sapply(subgroups, `[[`, "name")
sg_ids <- sapply(subgroups, `[[`, "id")
n_sg <- length(subgroups)
sim_hrs <- matrix(NA_real_, n_sims_null, n_sg, dimnames=list(NULL, sg_names))
sim_ubs <- matrix(NA_real_, n_sims_null, n_sg, dimnames=list(NULL, sg_names))
sim_ns  <- matrix(NA_real_, n_sims_null, n_sg, dimnames=list(NULL, sg_names))

for (ss in seq_len(n_sims_null)) {
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

V <- list(hr = sim_hrs, ub = sim_ubs, n = sim_ns)

# ---- package runner -----------------------------------------------------
cut_points <- list(age_med = age_med, size_med = size_med,
                    nodes_q3 = nodes_q3, pgr_med = pgr_med, er_cut = er_cut)
fit_v <- subgroup_cox(Surv(y_sim, event_sim) ~ treat_sim + strata(grade))
R <- run_subgroup_sims(
  dgm = dgm_null, subgroups = subgroups, n_sims = n_acc, fit = fit_v,
  baseline = "resample", n = nrow(gbsg),
  analysis_time = 84, max_entry = 24, cens_adjust = cal_null$cens_adjust,
  cutpoints = cut_points, hr_true = hr_true)

ok(identical(V$hr, R$sim_hrs), "sim_hrs bit-identical to the deck's loop")
ok(identical(V$ub, R$sim_ubs), "sim_ubs bit-identical")
ok(identical(V$n,  R$sim_ns),  "sim_ns  bit-identical")
S <- summary(R, hr_true = hr_true)
ok(inherits(S, "subgroup_sims_summary") && S$n_single + S$n_combo == sum(S$ok),
   "summary() runs; panels partition the valid subgroups")
cat("\n")
if (.fails == 0L) {
  cat("==== SLIDES ACCEPTANCE: ALL GATES PASSED ====\n")
} else {
  cat(sprintf("==== SLIDES ACCEPTANCE: %d GATE(S) FAILED ====\n", .fails))
  quit(status = 1L)
}
