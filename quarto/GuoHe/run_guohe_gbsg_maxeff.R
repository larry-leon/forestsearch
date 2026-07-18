# run_guohe_gbsg_maxeff.R
#
# Driver: fit forestsearch() on GBSG with sg_focus = "maxeff" (the argmax
# configuration Guo & He's theory covers) and apply the comparator bridge.
#
# Run it from anywhere -- it locates its own directory and sources its siblings
# from there, so the working directory does not matter:
#
#   Rscript quarto/GoHe/run_guohe_gbsg_maxeff.R
#
# or open it in RStudio and Source. Requires, in the SAME directory:
#   guohe_algorithm3.R, guohe_adaptive_r.R, guohe_from_forestsearch.R
#
# ACCEPTANCE CRITERION. The self-check line must read PASS: the reconstructed
# candidate family must equal the Tier-2 gate's `n_family`. If it does not, the
# base-R enumeration in the bridge has diverged from the package's own helpers
# (generate_combination_indices / get_covs_in / get_subgroup_membership, as used
# in R/forestsearch_main.R around line 2630). Diagnose the cause -- do NOT adjust
# the enumeration until the number happens to agree.

# ---- locate this script's directory, however it was invoked ---------------
.this_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1])))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  normalizePath(getwd())
})
cat("script directory:", .this_dir, "\n")

need <- c("guohe_algorithm3.R", "guohe_adaptive_r.R", "guohe_from_forestsearch.R")
miss <- need[!file.exists(file.path(.this_dir, need))]
if (length(miss)) {
  stop("Missing alongside this script: ", paste(miss, collapse = ", "))
}
for (f in need) source(file.path(.this_dir, f))

suppressMessages({
  library(forestsearch)
  library(survival)
  library(future)
})

# ---- data, mirroring analysis_gbsg_cox_multimethod_psi_v2_2new.qmd --------
df.analysis <- gbsg
df.analysis <- within(df.analysis, {
  id          <- seq_len(nrow(df.analysis))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
  treat       <- hormon
})
confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
outcome.name <- "time_months"; event.name <- "status"
id.name      <- "id";          treat.name <- "hormon"

n_cores <- max(1L, floor(0.80 * parallel::detectCores(logical = FALSE)))
cat("cores:", n_cores, "\n\n")

# ---- fit with sg_focus = "maxeff" -----------------------------------------
# NOTE: this deliberately differs from the headline analysis, which uses
# sg_focus = "effMaxSG" with selection_rule = "neighborhood". Guo & He's
# correction is built on the argmax functional, so "maxeff" is the only
# configuration it covers. Expect a DIFFERENT selected subgroup than the
# headline fit -- that is the configuration change, not an error.
cat("=== fitting forestsearch(sg_focus = 'eff') ===\n")
t0 <- proc.time()
fs <- forestsearch(
  df.analysis,
  outcome.name = outcome.name, event.name = event.name,
  treat.name = treat.name, id.name = id.name,
  confounders.name = confounders.name,
  is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  conf_force = c("er <= 0", "pgr <= 0"),
  conf.cont_jcuts = list(er = 10, pgr = 10),
  collapse_cuts = TRUE, cont.cutoff = 4,
  subgroup_method = "consistency",
  sg_focus = "eff",                          # <-- the argmax configuration
                                             # ("eff"/"hr" is the argmax-of-effect
                                             # primitive; "maxeff" is an internal
                                             # gate reselection token, not a sg_focus)
  selection_rule = "neighborhood", effect_neighborhood = 0.10,
  debias_gate = TRUE,                        # needed for the self-check
  debias_gate_args = list(draws = 5000L),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  max_subgroups_search = 50,
  hr.threshold = 1.0, hr.consistency = 1.0,
  pconsistency.threshold = 0.90,
  consistency_method = "resample",
  fs.splits = 1000, stop_threshold = 0.95,
  use_twostage = TRUE, max.minutes = 3,
  details = FALSE, quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "multisession", workers = n_cores,
                       show_message = FALSE)
)
plan("sequential")
cat(sprintf("fit completed in %.1f s\n", (proc.time() - t0)["elapsed"]))
cat("selected subgroup:", paste(fs$sg.harm, collapse = " & "), "\n")
cat("gate n_family:", tryCatch(fs$debias_gate$n_family, error = function(e) NA), "\n\n")

# ---- apply the Guo & He comparator: sweep r + one adaptive arm ------------
# All arms reuse the SAME fitted `fs` (no forestsearch refit) and run with
# parallel = TRUE, B = 2000. Parallelism is byte-identical to serial: resample
# indices are drawn serially under `seed`, only the model fitting is distributed.
r_sweep <- c(0.03, 0.15, 0.30, 0.45)
B_comp  <- 2000L
cat("=== Guo & He comparator on the reconstructed family ===\n")
cat(sprintf("comparator parallel: TRUE | workers: %d | B: %d\n", n_cores, B_comp))
cat(sprintf("fixed-r sweep: %s | adaptive arm: r_grid = {%s}, v = 5\n\n",
            paste(sprintf("%.2f", r_sweep), collapse = ", "),
            paste(sprintf("%.2f", r_sweep), collapse = ", ")))
future::plan(future::multisession, workers = n_cores)

# pull the reporting quantities out of a guohe_a3 fit into one row
.arm_row <- function(arm, r_used, a3) {
  data.frame(
    arm          = arm,
    r            = r_used,
    debiased_HR  = a3$debiased,
    bound_side   = a3$bound_side,
    bound_1sided = a3$bound_one_sided,
    ci_lo        = a3$ci_two_sided[1],
    ci_hi        = a3$ci_two_sided[2],
    stringsAsFactors = FALSE
  )
}

arms <- list()
rows <- list()

# ---- fixed-r arms ---------------------------------------------------------
for (rr in r_sweep) {
  cat(sprintf("--- fixed r = %.2f ---\n", rr))
  a <- guohe_from_forestsearch(
    fs, r = rr, B = B_comp, level = 0.05, seed = 20260718L,
    orient = +1,          # forest-search harm convention: most HARMFUL effect
    parallel = TRUE,
    verbose = (rr == r_sweep[1])   # print the family reconstruction only once
  )
  key <- sprintf("fixed_r_%.2f", rr)
  arms[[key]] <- a
  rows[[key]] <- .arm_row(sprintf("fixed r=%.2f", rr), a$r, a)
}

# ---- adaptive-r arm (Algorithm 2, v = 5 folds over the same grid) ---------
cat("\n--- adaptive r (guohe_adaptive_r, v = 5) ---\n")
ar <- guohe_from_forestsearch(
  fs, adaptive = TRUE, r_grid = r_sweep, v = 5L,
  B = B_comp, level = 0.05, seed = 20260718L,
  orient = +1, parallel = TRUE, verbose = FALSE
)
arms[["adaptive"]] <- ar
cat(sprintf("adaptive selected r = %.2f\n", ar$r_hat))
rows[["adaptive"]] <- .arm_row(sprintf("adaptive (r*=%.2f)", ar$r_hat),
                               ar$r_hat, ar$fit)

future::plan(future::sequential)

# ---- report table ---------------------------------------------------------
tab <- do.call(rbind, rows)
rownames(tab) <- NULL
cat("\n=== SWEEP SUMMARY (de-biased HR scale) ===\n")
tab_show <- within(tab, {
  debiased_HR  <- sprintf("%.4f", debiased_HR)
  bound_1sided <- sprintf("%.4f", bound_1sided)
  two_sided_CI <- sprintf("(%.4f, %.4f)", ci_lo, ci_hi)
  r            <- sprintf("%.2f", r)
})
print(tab_show[, c("arm", "r", "debiased_HR", "bound_side",
                   "bound_1sided", "two_sided_CI")], row.names = FALSE)

# ---- acceptance check (same reconstructed family for every arm) -----------
self_check <- arms[[1]]$bridge$self_check
cat("\n=== ACCEPTANCE CHECK ===\n")
cat("self-check:", self_check, "\n")
if (!grepl("^PASS", self_check)) {
  cat("\n*** Family reconstruction did not match the gate. Diagnose against the\n",
      "*** package's own enumeration helpers before using this result.\n", sep = "")
} else {
  cat("Family reconstruction is faithful to the search's own enumeration.\n")
}

res <- list(fixed = arms[grep("^fixed", names(arms))],
            adaptive = arms[["adaptive"]], summary = tab)
saveRDS(res, file.path(.this_dir, "guohe_gbsg_maxeff_result.rds"))
cat("\nsaved:", file.path(.this_dir, "guohe_gbsg_maxeff_result.rds"), "\n")
