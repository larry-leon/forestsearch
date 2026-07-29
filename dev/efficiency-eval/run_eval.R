# =============================================================================
# forestsearch efficiency evaluation -- driver
# =============================================================================
#   source("run_eval.R")                       # everything
#   source("run_eval.R"); run_eval(c("G","M")) # gates + micro only
#
# Stages: G gates, M micro, P parallel, E end-to-end.
# Each writes results/<stage>.rds so long stages can be run separately.
# =============================================================================

suppressPackageStartupMessages({
  library(survival); library(data.table); library(digest)
  library(future);   library(future.apply); library(parallelly)
})

source("config.R")
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

`%||%` <- function(a, b) if (is.null(a)) b else a


#' Load the evaluation dataset (GBSG, as used throughout the package)
load_eval_data <- function() {
  df <- as.data.frame(survival::gbsg)
  df$grade3 <- as.integer(df$grade == "3")
  df$id <- seq_len(nrow(df))
  df
}


run_eval <- function(stages = c("G", "M", "P", "E"), cfg = EVAL) {

  dir.create(cfg$results_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(cfg$seed)

  if (!requireNamespace("forestsearch", quietly = TRUE)) {
    stop("forestsearch is not installed. Install it (devtools::install(), NOT ",
         "load_all() -- future workers only see installed packages), then rerun.")
  }

  # ---- read-only contract, opened -----------------------------------------
  message("== read-only guard: hashing package sources ==")
  h_before <- fs_hash_sources()
  message("   ", length(h_before), " files hashed")

  env <- fs_capture_env()
  print(env)
  saveRDS(env, file.path(cfg$results_dir, "environment.rds"))

  df <- load_eval_data()
  message("== data: GBSG, ", nrow(df), " rows ==\n")

  out <- list(environment = env)

  # ---- G: gates -----------------------------------------------------------
  if ("G" %in% stages) {
    message("== STAGE G: equivalence gates ==")
    g <- run_gates(cfg)
    print(g[, c("gate", "what", "pass", "detail")], row.names = FALSE)
    saveRDS(g, file.path(cfg$results_dir, "G.rds"))
    out$gates <- g
    if (!all(g$pass, na.rm = TRUE)) {
      message("\n!! One or more gates FAILED. Timings below are reported for ",
              "completeness but MUST NOT be used to justify incorporating a ",
              "candidate whose gate did not pass.\n")
    }
    message("")
  }

  # ---- M: micro -----------------------------------------------------------
  if ("M" %in% stages) {
    message("== STAGE M: unit-cost microbenchmarks ==")
    m <- run_micro(cfg)
    print(m[, c("id", "what", "n", "current_ms", "candidate_ms", "speedup")],
          row.names = FALSE)
    saveRDS(m, file.path(cfg$results_dir, "M.rds"))
    out$micro <- m
    message("")
  }

  # ---- P: parallel --------------------------------------------------------
  if ("P" %in% stages) {
    message("== STAGE P: parallel behaviour (the headline) ==")
    p <- tryCatch(run_parallel(cfg, df), error = function(e) {
      message("   stage P failed: ", conditionMessage(e)); NULL })
    if (!is.null(p)) {
      print(p, row.names = FALSE)
      saveRDS(p, file.path(cfg$results_dir, "P.rds"))
      out$parallel <- p
    }
    message("")
  }

  # ---- E: end-to-end ------------------------------------------------------
  if ("E" %in% stages) {
    message("== STAGE E: end-to-end profile ==")
    e <- tryCatch(run_endtoend(cfg, df, out$micro), error = function(e) {
      message("   stage E failed: ", conditionMessage(e)); NULL })
    if (!is.null(e)) {
      print(e, row.names = FALSE)
      saveRDS(e, file.path(cfg$results_dir, "E.rds"))
      out$endtoend <- e
    }
    message("")
  }

  # ---- read-only contract, closed ----------------------------------------
  message("== read-only guard: re-hashing ==")
  guard <- fs_guard_verify(h_before, fs_hash_sources())
  message("   ", guard$note)
  out$guard <- guard
  saveRDS(out, file.path(cfg$results_dir, "all.rds"))

  message("\nWrote ", normalizePath(cfg$results_dir))
  message("Render report.qmd for the decision table.")
  invisible(out)
}

if (sys.nframe() == 0L && !interactive()) run_eval()
