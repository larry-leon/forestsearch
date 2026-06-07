# =============================================================================
# fs_profile_harness.R
#
# A small, non-invasive profiling harness for forestsearch().
#
# Purpose
#   Answer the only question that decides where the "biggest efficiency gain"
#   lives: within one forestsearch() call, how is time split between
#     (a) candidate GENERATION  -- GRF (causal_forest) + DINA + (LASSO), and
#     (b) the CONSISTENCY search -- subgroup.consistency() and its splits?
#
#   forestsearch() runs these as sequential phases, and the 500-rep bootstrap
#   re-runs the whole search ~500x per simulated trial, so the per-call phase
#   split, multiplied by (1 + nb_boots), characterises the whole study budget.
#
# How it works
#   Wraps a single forestsearch() call in utils::Rprof() (statistical stack
#   sampling -- no package edits), then attributes total time to the
#   phase-entry functions. Also prints the raw top-N by total/self time so the
#   breakdown is verifiable even if an internal name differs from the buckets.
#
# Why sequential
#   Rprof samples the *main* R process; under plan = "multisession" the GRF and
#   per-split work happens in worker processes whose time the main-process
#   profiler cannot see (generation/consistency would look falsely cheap). The
#   harness therefore forces plan = "sequential", which also matches the
#   simulation/bootstrap inner-loop regime where any optimisation actually pays
#   off. details/quiet are quieted so cat() I/O does not skew the sample.
#
# Dependencies: base + utils only (Rprof/summaryRprof). No new packages.
# =============================================================================

#' Profile one forestsearch() call and decompose its phase budget
#'
#' Runs \code{forestsearch()} under \code{\link[utils]{Rprof}} and reports how
#' wall-clock time divides across data/cut prep, GRF generation, DINA
#' generation, LASSO, and the consistency search. Repeats a few times to damp
#' sampling noise and reports per-call timing plus an extrapolation to the
#' bootstrap/study budget.
#'
#' @param df Data frame passed as \code{df.analysis} to \code{forestsearch()}.
#' @param confounders Character vector passed as \code{confounders.name}.
#' @param fs_params Named list of forestsearch() arguments (e.g. the result of
#'   \code{get_fs3b_params(sg_focus)}). Copied internally; \code{parallel_args}
#'   is forced sequential and \code{details}/\code{quiet} are quieted for a
#'   clean sample. The caller's object is not modified.
#' @param reps Integer. Number of profiled calls to average. Default 3.
#' @param interval Numeric. Rprof sampling interval (s). Default 0.01.
#' @param nb_boots Integer. Bootstrap reps to extrapolate the study budget to
#'   (does not run the bootstrap). Default 500.
#' @param top_n Integer. Rows in the raw verification tables. Default 20.
#'
#' @return Invisibly, a list with \code{phase_table} (per-phase seconds and
#'   percent), \code{call_seconds} (mean per-call wall-clock), \code{raw_total}
#'   / \code{raw_self} (top-N raw Rprof tables), and \code{diagnostics}
#'   (consistency candidate counts if discoverable in the return object).
profile_forestsearch <- function(df, confounders, fs_params,
                                 reps = 3L, interval = 0.01,
                                 nb_boots = 500L, top_n = 20L) {

  if (!exists("forestsearch", mode = "function")) {
    stop("forestsearch() not found. install/attach the package first ",
         "(devtools::install(); library(forestsearch)).")
  }

  # -- Phase buckets: ordered, disjoint regexes matched against function names.
  #    First match wins, so list specific phases before shared primitives.
  phase_patterns <- list(
    "data/cut prep"   = "get_FSdata|getFSdata|collapse_redundant_cuts|collapse_cuts|exclude_cuts|tidy_cut|is_flag_drop",
    "GRF generation"  = "causal_forest|fit_causal_forest|grf\\.subg|grf_subgroup|policy_tree|regression_forest|\\bgrf\\b",
    "DINA generation" = "dina_subgroup|dina_",
    "LASSO generation"= "cv\\.glmnet|glmnet|lasso",
    "consistency search" = "subgroup\\.consistency|subgroup_consistency|evaluate_consistency|evaluate_subgroup|run_single_consistency|early_stop_decision|sort_subgroups"
  )

  # -- Prep a clean, sequential copy of the params (do not mutate caller's).
  p <- fs_params
  p$parallel_args <- list(plan = "sequential", workers = 1L,
                          show_message = FALSE)
  if (!is.null(p$details)) p$details <- FALSE
  if (!is.null(p$quiet))   p$quiet   <- TRUE
  if (!is.null(p$show_candidate_summary)) p$show_candidate_summary <- FALSE

  # -- Recursive probe for consistency diagnostics in the return object.
  find_field <- function(x, nm, depth = 0L) {
    if (depth > 6L || is.null(x)) return(NULL)
    if (is.list(x)) {
      if (!is.null(x[[nm]])) return(x[[nm]])
      for (el in x) {
        hit <- find_field(el, nm, depth + 1L)
        if (!is.null(hit)) return(hit)
      }
    }
    NULL
  }

  call_secs <- numeric(reps)
  tot_acc <- self_acc <- list()
  # per-phase sample tallies, accumulated across reps (raw-stack attribution)
  phase_samples <- setNames(numeric(length(phase_patterns)),
                            names(phase_patterns))
  other_samples <- 0; total_samples <- 0; samp_secs <- interval
  diagnostics <- NULL
  prof_file <- tempfile(fileext = ".Rprof")

  for (r in seq_len(reps)) {
    gc(verbose = FALSE)
    t0 <- proc.time()[["elapsed"]]
    utils::Rprof(prof_file, interval = interval, line.profiling = FALSE,
                 memory.profiling = FALSE)
    fs_out <- do.call(forestsearch,
                      c(list(df.analysis = df, confounders.name = confounders),
                        p))
    utils::Rprof(NULL)
    call_secs[r] <- proc.time()[["elapsed"]] - t0

    # raw top-N tables (by.total / by.self are correct for these)
    sm <- summaryRprof(prof_file)
    bt <- sm$by.total; bs <- sm$by.self
    for (fn in rownames(bt)) tot_acc[[fn]]  <- c(tot_acc[[fn]],  bt[fn, "total.time"])
    for (fn in rownames(bs)) self_acc[[fn]] <- c(self_acc[[fn]], bs[fn, "self.time"])

    # robust phase attribution: read raw sample stacks; assign each sample to
    # the FIRST phase whose pattern appears anywhere in its stack. Phases are
    # sequential siblings under forestsearch(), so >=1 entry fn appears per
    # stack -> no double counting; intra-phase nesting maps to the same phase.
    lines <- readLines(prof_file, warn = FALSE)
    hdr <- grep("^sample.interval=", lines)
    if (length(hdr)) {
      si <- as.numeric(sub("sample.interval=", "", lines[hdr[1]]))
      if (is.finite(si) && si > 0) samp_secs <- si / 1e6   # us -> s
      lines <- lines[-seq_len(hdr[1])]
    }
    lines <- lines[nzchar(lines)]
    total_samples <- total_samples + length(lines)
    for (ln in lines) {
      hit <- NA_character_
      for (ph in names(phase_patterns)) {
        if (grepl(phase_patterns[[ph]], ln, ignore.case = TRUE)) { hit <- ph; break }
      }
      if (is.na(hit)) other_samples <- other_samples + 1
      else phase_samples[hit] <- phase_samples[hit] + 1
    }

    if (is.null(diagnostics)) {
      ne <- find_field(fs_out, "n_candidates_evaluated")
      nt <- find_field(fs_out, "n_candidates_total")
      np <- find_field(fs_out, "n_passed")
      if (!is.null(ne) || !is.null(nt))
        diagnostics <- list(n_candidates_evaluated = ne,
                            n_candidates_total = nt, n_passed = np)
    }
  }
  if (file.exists(prof_file)) unlink(prof_file)

  mean_total <- vapply(tot_acc,  function(v) mean(v), numeric(1))
  mean_self  <- vapply(self_acc, function(v) mean(v), numeric(1))
  call_mean  <- mean(call_secs)

  # convert accumulated phase samples -> mean seconds per call
  phase_secs <- (phase_samples * samp_secs) / reps
  other <- max(0, (other_samples * samp_secs) / reps)
  phase_secs <- c(phase_secs, "other / overhead" = other)
  phase_table <- data.frame(
    phase   = names(phase_secs),
    seconds = round(as.numeric(phase_secs), 3),
    percent = round(100 * as.numeric(phase_secs) / call_mean, 1),
    row.names = NULL
  )

  ord_t <- head(sort(mean_total, decreasing = TRUE), top_n)
  ord_s <- head(sort(mean_self,  decreasing = TRUE), top_n)
  raw_total <- data.frame(`function` = names(ord_t),
                          total_sec = round(as.numeric(ord_t), 3),
                          check.names = FALSE, row.names = NULL)
  raw_self  <- data.frame(`function` = names(ord_s),
                          self_sec = round(as.numeric(ord_s), 3),
                          check.names = FALSE, row.names = NULL)

  # ---- report ---------------------------------------------------------------
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("forestsearch() PHASE PROFILE  (sequential, reps = ", reps,
      ", interval = ", interval, "s)\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
  cat(sprintf("Mean wall-clock per call: %.2f s  (range %.2f-%.2f)\n\n",
              call_mean, min(call_secs), max(call_secs)))
  print(phase_table, row.names = FALSE)

  cat("\nExtrapolation (per simulated trial, point fit + ", nb_boots,
      " bootstrap reruns):\n", sep = "")
  cat(sprintf("  generation phases ~ %.0f%% | consistency ~ %.0f%% of each call\n",
              sum(phase_table$percent[grepl("generation|prep", phase_table$phase)]),
              phase_table$percent[phase_table$phase == "consistency search"]))
  cat(sprintf("  est. search budget / trial: %.1f min  (= %.2f s x %d calls)\n",
              call_mean * (1 + nb_boots) / 60, call_mean, 1L + nb_boots))

  if (!is.null(diagnostics)) {
    cat("\nConsistency search candidates: evaluated = ",
        diagnostics$n_candidates_evaluated %||% NA,
        " of total = ", diagnostics$n_candidates_total %||% NA,
        " (passed = ", diagnostics$n_passed %||% NA, ")\n", sep = "")
    cat("  -> the effect-desc skip would avoid the out-of-band share of these.\n")
  }

  cat("\nRaw top-", top_n, " by TOTAL time (verify the buckets above):\n", sep = "")
  print(raw_total, row.names = FALSE)
  cat("\nRaw top-", top_n, " by SELF time:\n", sep = "")
  print(raw_self, row.names = FALSE)
  cat(strrep("=", 70), "\n", sep = "")

  invisible(list(phase_table = phase_table, call_seconds = call_mean,
                 raw_total = raw_total, raw_self = raw_self,
                 diagnostics = diagnostics))
}

# small null-coalesce used in the report
`%||%` <- function(a, b) if (is.null(a)) b else a
