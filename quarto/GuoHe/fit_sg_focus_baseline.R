# fit_sg_focus_baseline.R  (analysis driver -- touches no package file)
#
# Fit forestsearch() on GBSG three times at a fixed seed, varying ONLY sg_focus
# over {"hr", "effMaxSG", "maxSG"}, mirroring the configuration of
# run_guohe_gbsg_maxeff.R (same directory).  Record, per fit, the selected
# subgroup (fs$sg.harm) and its naive HR.  Saves sg_focus_baseline.rds beside
# this script.
#
# CANDIDATE-POOL DEPTH.  run_guohe_gbsg_maxeff.R uses max_subgroups_search = 50.
# At that depth "maxSG" selects NO subgroup: the consistency pool is the top-K
# of the sg_focus-ORDERED preview list (subgroup_consistency_main.R:552-568), and
# "maxSG" orders by largest-N first, so the harmful subgroup here ({er<=0} &
# {size<=35}, n = 61, at the n.min floor) ranks in the small-N tail and is
# truncated away before consistency.  ("hr" orders by -HR and puts it at the
# top; "effMaxSG" evaluates the full pool.)  We therefore use a DEEP pool below
# so every focus screens the same, un-truncated candidate set -- an
# apples-to-apples baseline.  See MAX_SG.
#
# Run from anywhere -- it locates its own directory, so the working directory
# does not matter:
#   Rscript quarto/GuoHe/fit_sg_focus_baseline.R
# or open in RStudio and Source.

# Candidate-pool depth for the consistency screen (see header).  Deep enough to
# reach the small-N harmful subgroups that "maxSG" orders last (~90 sufficed
# here; 2000 is comfortably past the full screened family).
MAX_SG <- 2000L

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

suppressMessages({
  library(forestsearch)
  library(survival)
  library(future)
})

# ---- data, identical to run_guohe_gbsg_maxeff.R ---------------------------
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

# ---- one fit, config identical to the run script except sg_focus ----------
fit_one <- function(focus) {
  cat(sprintf("=== fitting forestsearch(sg_focus = '%s') ===\n", focus))
  t0 <- proc.time()
  future::plan(future::multisession, workers = n_cores)
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
    sg_focus = focus,                          # <-- the only thing that varies
    selection_rule = "neighborhood", effect_neighborhood = 0.10,
    debias_gate = TRUE,
    debias_gate_args = list(draws = 5000L),
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    max_subgroups_search = MAX_SG,             # deep pool (see header) -- 50 truncates maxSG to NULL
    hr.threshold = 1.0, hr.consistency = 1.0,
    pconsistency.threshold = 0.90,
    consistency_method = "resample",
    fs.splits = 1000, stop_threshold = 0.95,
    use_twostage = TRUE, max.minutes = 3,
    details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "multisession", workers = n_cores,
                         show_message = FALSE)
  )
  future::plan(future::sequential)
  elapsed <- (proc.time() - t0)["elapsed"]

  # safe scalar: NA when the field is NULL / length 0
  scalar <- function(x, na) if (is.null(x) || length(x) == 0L) na else x[[1L]]

  sg   <- fs$sg.harm
  sg_label <- if (is.null(sg) || !length(sg)) NA_character_ else paste(sg, collapse = " & ")

  # naive HR of the SELECTED subgroup (fs$sg.harm), computed on the package's
  # own estimation frame (df.est) with its own membership indicator.  df.est is
  # merge(df, temp) -- it carries the USER-FACING column names, so resolve them.
  naive_hr <- NA_real_
  if (!is.null(sg) && length(sg) &&
      !is.null(fs$grp.consistency) &&
      !is.null(fs$grp.consistency$sg.harm.id) &&
      !is.null(fs$df.est) &&
      length(fs$grp.consistency$sg.harm.id) == nrow(fs$df.est)) {
    dat <- fs$df.est
    pick <- function(cands) { hit <- cands[cands %in% names(dat)]; if (length(hit)) hit[1L] else NA_character_ }
    tcol <- pick(c("Treat", treat.name, "treat", "hormon"))
    ycol <- pick(c("Y", outcome.name, "time_months"))
    ecol <- pick(c("Event", event.name, "status"))
    memb <- fs$grp.consistency$sg.harm.id == 1L
    if (!anyNA(c(tcol, ycol, ecol))) {
      dsub <- dat[memb, , drop = FALSE]
      naive_hr <- tryCatch(
        exp(unname(coxph(Surv(dsub[[ycol]], dsub[[ecol]]) ~ dsub[[tcol]])$coefficients[[1L]])),
        error = function(e) NA_real_
      )
    }
  }

  # cross-reference: the Tier-2 gate's own naive estimate + its selected label
  gate_naive <- scalar(tryCatch(fs$debias_gate$naive$est,      error = function(e) NULL), NA_real_)
  gate_label <- scalar(tryCatch(fs$debias_gate$selected_label, error = function(e) NULL), NA_character_)
  n_sub      <- if (!is.null(fs$grp.consistency$sg.harm.id))
                  sum(fs$grp.consistency$sg.harm.id == 1L) else NA_integer_

  fnum <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) "NA" else sprintf("%.4f", x)
  fchr <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) "NA" else as.character(x)
  cat(sprintf("  selected subgroup : %s\n", fchr(sg_label)))
  cat(sprintf("  subgroup n        : %s\n", fchr(n_sub)))
  cat(sprintf("  naive HR (sg.harm): %s\n", fnum(naive_hr)))
  cat(sprintf("  [gate] naive est / label: %s / %s\n", fnum(gate_naive), fchr(gate_label)))
  cat(sprintf("  fit time          : %.1f s\n\n", elapsed))

  list(sg_focus = focus, sg.harm = sg, sg_label = sg_label,
       subgroup_n = n_sub, naive_hr = naive_hr,
       gate_naive_est = gate_naive, gate_selected_label = gate_label,
       seedit = 8316951, elapsed_sec = unname(elapsed))
}

foci <- c("hr", "effMaxSG", "maxSG")
res  <- lapply(foci, fit_one)
names(res) <- foci

# ---- summary table --------------------------------------------------------
tab <- do.call(rbind, lapply(res, function(x) data.frame(
  sg_focus       = x$sg_focus,
  selected       = x$sg_label,
  subgroup_n     = x$subgroup_n,
  naive_HR       = x$naive_hr,
  gate_naive_est = x$gate_naive_est,
  gate_label     = x$gate_selected_label,
  stringsAsFactors = FALSE
)))
rownames(tab) <- NULL

out <- list(results = res, summary = tab, seedit = 8316951,
            max_subgroups_search = MAX_SG)
saveRDS(out, file.path(.this_dir, "sg_focus_baseline.rds"))

cat("=== sg_focus baseline (seedit = 8316951) ===\n")
tab_show <- tab
tab_show$naive_HR       <- sprintf("%.4f", tab_show$naive_HR)
tab_show$gate_naive_est <- ifelse(is.na(tab$gate_naive_est), "NA",
                                  sprintf("%.4f", tab$gate_naive_est))
print(tab_show, row.names = FALSE)
cat("\nsaved:", file.path(.this_dir, "sg_focus_baseline.rds"), "\n")
