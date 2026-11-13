# Part 1 unaffected case: OR, IRR, MD and survival values from the two
# extractors (hand data and simulated data) and run_simulation_analysis()
# replicates, saved for an identical() comparison between two trees.
args <- commandArgs(TRUE); tree <- args[1]; out_rds <- args[2]
which_part <- if (length(args) >= 3) args[3] else "all"
suppressMessages(pkgload::load_all(tree, quiet = TRUE))
future::plan("sequential")
source(file.path(tree, "tests/testthat/helper-synthetic-dgm.R"))
res <- list()
t0 <- proc.time()[["elapsed"]]

# ---- (a) function values on hand data and simulated data ------------------
set.seed(20261112)
N <- 240
hd <- data.frame(id = seq_len(N), treat = rep(0:1, N / 2),
                 flag_harm = rep(c(1L, 0L, 0L), N / 3),
                 Hhat = rep(c(1L, 1L, 0L, 0L), N / 4),
                 time = round(runif(N, 0.5, 3), 2))
hd$ybin  <- rbinom(N, 1, plogis(-0.5 + 0.6 * hd$treat * hd$flag_harm))
hd$ycnt  <- rpois(N, hd$time * exp(-0.7 + 0.5 * hd$treat * hd$flag_harm))
hd$ycont <- rnorm(N, 1 + 0.8 * hd$treat * hd$flag_harm)
hd$ysurv <- rexp(N, exp(-0.2 + 0.4 * hd$treat * hd$flag_harm))
hd$ev    <- rbinom(N, 1, 0.8)

ext <- function(eng, y, ev, em, ot, off) {
  if (eng == "FS") {
    .extract_fs_estimates_gen(
      df = hd, fs_res = NULL,
      fs_full = list(df.est = data.frame(treat.recommend = 1L - hd$Hhat),
                     sg.harm = "Hhat == 1"),
      dgm = NULL, cox_formula = NULL, cox_formula_adj = NULL,
      outcome_name = y, event_name = ev, treat_name = "treat",
      harm_col = "flag_harm", analysis = "FS", verbose = FALSE,
      outcome_type = ot, effect_measure = em, offset_name = off)
  } else {
    .extract_grf_estimates_gen(
      df = hd, grf_est = list(sg.harm.id = "Hhat == 1",
                              data = cbind(hd, treat.recommend = 1L - hd$Hhat)),
      dgm = NULL, cox_formula = NULL, cox_formula_adj = NULL,
      outcome_name = y, event_name = ev, treat_name = "treat",
      harm_col = "flag_harm", analysis = "GRF", verbose = FALSE, debug = FALSE,
      outcome_type = ot, effect_measure = em, offset_name = off, id_name = "id")
  }
}
cases <- list(
  OR          = list(y = "ybin",  ev = "ybin", em = "OR",  ot = "binary",     off = NULL),
  IRR_offset  = list(y = "ycnt",  ev = "ev",   em = "IRR", ot = "count",      off = "time"),
  IRR_nooff   = list(y = "ycnt",  ev = "ev",   em = "IRR", ot = "count",      off = NULL),
  MD          = list(y = "ycont", ev = "ev",   em = "MD",  ot = "continuous", off = NULL),
  NULL_em_bin = list(y = "ybin",  ev = "ybin", em = NULL,  ot = "binary",     off = NULL),
  survival    = list(y = "ysurv", ev = "ev",   em = NULL,  ot = NULL,         off = NULL),
  survival_ot = list(y = "ysurv", ev = "ev",   em = NULL,  ot = "survival",   off = NULL))
if (which_part %in% c("all", "a")) {
  for (nm in names(cases)) for (eng in c("FS", "GRF")) {
    cc <- cases[[nm]]
    o <- ext(eng, cc$y, cc$ev, cc$em, cc$ot, cc$off)
    res[[paste("fn", eng, nm, sep = ":")]] <- list(table = as.data.frame(o),
                                                   sg_hat = attr(o, "sg_hat"))
  }
}

# ---- (b) run_simulation_analysis(), 3 replicates per outcome type ----------
drop_timing <- function(x) { x$fs_minutes <- NULL; x }
if (which_part %in% c("all", "b")) {
  mk <- function(d, ot, em, off = NULL) generate_glm_dgm(
    data = d, factor_vars = c("biomarker_hi", "sex"), continuous_vars = "age",
    outcome_var = "y", treatment_var = "treat", outcome_type = ot,
    effect_measure = em, offset_var = off,
    subgroup_vars = c("biomarker_hi", "sex"),
    subgroup_cuts = list(biomarker_hi = 1L, sex = 1L),
    model = "alt", k_inter = 2.0, n_super = 2000L, seed = 8316951L,
    verbose = FALSE)
  dgms <- list(
    binary_OR      = list(dgm = mk(.make_binary_data(N = 300L), "binary", "OR"),
                          fsp = list(outcome_type = "binary", effect_measure = "OR")),
    continuous_MD  = list(dgm = mk(.make_continuous_data(N = 300L), "continuous", "MD"),
                          fsp = list(outcome_type = "continuous", effect_measure = "MD")),
    count_IRR      = list(dgm = mk(.make_count_data(N = 300L), "count", "IRR", "ftime"),
                          fsp = list(outcome_type = "count", effect_measure = "IRR",
                                     offset.name = "ftime")),
    survival       = list(dgm = setup_gbsg_dgm(model = "alt", n_super = 2000L,
                                               verbose = FALSE),
                          fsp = list()))
  for (nm in names(dgms)) {
    d <- dgms[[nm]]
    conf <- if (nm == "survival") c("v1", "v2", "v3", "v4", "v5", "v6", "v7")
            else c("biomarker_hi", "sex", "age")
    for (i in 1:3) {
      t1 <- proc.time()[["elapsed"]]
      r <- suppressWarnings(suppressMessages(run_simulation_analysis(
        sim_id = i, dgm = d$dgm, n_sample = 600, confounders_base = conf,
        run_fs = TRUE, run_fs_grf = FALSE, run_grf = TRUE,
        fs_params = c(d$fsp, list(fs.splits = 100)), grf_params = list(dmin.grf = 0), verbose = FALSE)))
      res[[paste("rsa", nm, i, sep = ":")]] <- drop_timing(as.data.frame(r))
      cat(sprintf("%-14s sim %d  %.2f s  rows %d  any.H %s\n", nm, i,
                  proc.time()[["elapsed"]] - t1, nrow(r),
                  paste(r$any.H, collapse = ",")))
    }
  }
}
saveRDS(res, out_rds)
cat("saved", length(res), "objects; elapsed", round(proc.time()[["elapsed"]] - t0, 2), "s\n")
