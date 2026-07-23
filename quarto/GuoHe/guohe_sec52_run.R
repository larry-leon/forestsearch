# guohe_sec52_run.R
#
# Full-study driver for the Guo & He (2021) Section 5.2 / Table 7 reproduction.
# Intended for the 127-core host, NOT the chat sandbox.
#
# PREREQUISITE -- production step 0. The truth caches must exist BEFORE any
# coverage run; this driver REFUSES to proceed without them:
#
#   Rscript quarto/GuoHe/guohe_sec52_truth.R
#
# (defaults n_big = 2e6, c_step = 0.25, seeds 20260721 + 100*beta2; exits
# non-zero on any gate failure -- a failed or suspicious truth curve must
# never reach a coverage run).
#
# PILOT BEFORE PRODUCTION. Run the pilot first and REPORT THE PROJECTION TO
# THE OWNER before launching the full study:
#
#   Rscript quarto/GuoHe/guohe_sec52_run.R --pilot --B=2000 --cores=120
#
# The pilot runs beta2 = 0 (the maximal-bias case) at --reps=20 with the
# intended B and prints the projected wall-clock for 6 x 2000 at the chosen
# core count. Levers, in order: reduce B (paper silent; 2000 is the Sec 5.1
# inference; 500 mainly costs bootstrap resolution); reduce --reps (500 gives
# coverage MCSE ~ 0.010 -- report the MCSE used). NEVER thin the cutpoint
# family -- family size is a design parameter, not a tuning knob.
#
# Parallelism is over REPLICATIONS ONLY -- the outermost loop -- so no nested
# futures are created and `guohe_algorithm3(parallel = )` stays at its default
# FALSE. Each replication pins its own RNG stream via `set.seed()` inside
# `gh52_one_rep()`, so a parallel run reproduces a serial run bit-for-bit.
# Each scenario writes its own .rds on completion, so the job is restartable:
# a scenario whose output already exists is skipped unless --force is given.
#
# Usage:
#   Rscript quarto/GuoHe/guohe_sec52_run.R --pilot --B=2000 --cores=120
#   Rscript quarto/GuoHe/guohe_sec52_run.R --B=2000 --cores=120
#   Rscript quarto/GuoHe/guohe_sec52_run.R --beta2=0,0.1 --reps=2000
#
# Flags:
#   --pilot         beta2 = 0 only, reps = 20, output *_pilot.rds, projection
#   --beta2=LIST    scenario subset (default 0,0.1,0.2,0.3,0.4,0.5)
#   --cores=N       forked workers (default: 80% of physical cores)
#   --reps=N        replications per scenario (default 2000, the published value)
#   --B=N           bootstrap resamples (default 2000; the paper states no value
#                   for Section 5 -- 2000 is the Section 5.1 inference and is
#                   recorded as such in the output)
#   --out=DIR       output directory (default: alongside this script)
#   --truth-dir=DIR where the truth caches live (default: --out)
#   --force         recompute scenarios whose .rds already exists

suppressMessages(library(survival))

.gh52_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
})
source(file.path(.gh52_dir, "guohe_sec52_sim.R"))

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ---- arguments -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(args == paste0("--", nm))
opt <- function(nm, default) {
  hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
}

pilot <- flag("pilot")
b2_grid <- as.numeric(strsplit(opt("beta2", "0,0.1,0.2,0.3,0.4,0.5"), ",")[[1]])
n_rep <- as.integer(opt("reps", "2000"))
b_boot <- as.integer(opt("B", "2000"))
force <- flag("force")
out_dir <- opt("out", .gh52_dir)
truth_dir <- opt("truth-dir", out_dir)
n_cores <- as.integer(opt(
  "cores",
  as.character(max(1L, floor(0.80 * parallel::detectCores(logical = FALSE))))
))
if (pilot) {
  b2_grid <- 0
  n_rep <- as.integer(opt("reps", "20"))
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Guo & He Section 5.2 / Table 7 reproduction\n")
cat(sprintf("  mode      : %s\n", if (pilot) "PILOT (projection only)" else "production"))
cat(sprintf("  beta2     : %s\n", paste(sprintf("%.1f", b2_grid), collapse = ", ")))
cat(sprintf("  reps      : %d   B = %d\n", n_rep, b_boot))
cat(sprintf("  cores     : %d\n", n_cores))
cat(sprintf("  out       : %s\n", out_dir))
cat(sprintf("  truth dir : %s\n\n", truth_dir))

# ---- truth caches: refuse to proceed without them --------------------------
truth_file <- function(b2) {
  file.path(truth_dir, sprintf("guohe_sec52_truth_beta2_%02d.rds", round(b2 * 10)))
}
missing_tf <- b2_grid[!file.exists(vapply(b2_grid, truth_file, character(1)))]
if (length(missing_tf)) {
  stop("Truth cache(s) missing for beta2 = ",
       paste(sprintf("%.1f", missing_tf), collapse = ", "),
       ".\nRun production step 0 first:  Rscript quarto/GuoHe/guohe_sec52_truth.R\n",
       "A coverage run must never start without gate-cleared truth curves.",
       call. = FALSE)
}
truths <- lapply(b2_grid, function(b2) {
  tr <- readRDS(truth_file(b2))
  if (!inherits(tr, "gh52_truth") || !all(tr$gates$pass)) {
    stop("Truth cache for beta2 = ", b2, " is not a gate-cleared gh52_truth ",
         "object; recompute with guohe_sec52_truth.R.", call. = FALSE)
  }
  cat(sprintf("  [truth] beta2 = %.1f : n_big = %g, c_step = %.2f, seed = %d, all gates PASS\n",
              tr$beta2, tr$n_big, tr$c_step, tr$seed))
  tr
})

# ---- exact-null scoring basis (beta_exact): validate, never mutate ---------
# The truth MODULE owns the scoring basis. At beta2 = 0 the DGM makes
# beta(c) = 0 exact (b(.) == 0 everywhere, so the treatment-only Cox model is
# correctly specified within every S(c)); gh52_truth_curve() stores that
# closed form as `beta_exact` and gh52_truth_at()'s default "auto" scores
# against it, so no truth-curve Monte Carlo error enters the null scenario's
# coverage indicators. The driver only VALIDATES: a beta2 = 0 cache written
# before `beta_exact` existed would silently fall back to the simulated
# curve under "auto" -- injecting a FIXED per-scenario offset shared by all
# replicates -- so such a cache is refused here.
for (tr in truths) {
  if (tr$beta2 == 0 && is.null(tr$beta_exact)) {
    stop("The beta2 = 0 truth cache predates the exact-null scoring basis ",
         "(no `beta_exact`); regenerate the caches:\n",
         "  Rscript quarto/GuoHe/guohe_sec52_truth.R --force", call. = FALSE)
  }
  cat(sprintf("  [truth] beta2 = %.1f : scoring basis %s\n", tr$beta2,
              if (is.null(tr$beta_exact)) "simulated (isotonic)" else
                "EXACT beta(c) = 0"))
}
cat("\n")

# ---- run -------------------------------------------------------------------
run_scenario <- function(b2, truth) {
  id <- sprintf("t7_beta2_%02d", round(b2 * 10))
  f_out <- file.path(out_dir, paste0(
    "guohe_repro_", id, if (pilot) "_pilot" else "", ".rds"))
  if (file.exists(f_out) && !force) {
    cat(sprintf("[skip] %s (exists)\n", id))
    return(invisible(NULL))
  }
  # Scenario-specific seed base, so scenarios are independent but each
  # replication is individually reproducible from (id, m). NOT the Section 5.1
  # construction: its 997 multiplier puts adjacent scenario bases only 997
  # apart while the replicate-seed ranges span `reps`, so the six beta2
  # scenarios would share most of their replicate seeds. The wide prime
  # 100003 separates the ranges.
  base <- 1000000L + as.integer(sum(utf8ToInt(id)) * 100003L)
  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("[run ] %s  beta2=%.1f n=400  %d reps, B=%d ...\n",
              id, b2, n_rep, b_boot))
  utils::flush.console()

  rows <- parallel::mclapply(
    seq_len(n_rep),
    function(m) {
      t1 <- proc.time()[["elapsed"]]
      out <- try(gh52_one_rep(
        beta2 = b2, n = 400L, truth = truth, r_grid = GH52_R_GRID,
        B = b_boot, level = 0.05, min_events = 5L, seed = base + m
      ), silent = TRUE)
      if (!inherits(out, "try-error")) {
        out$rep_elapsed_sec <- proc.time()[["elapsed"]] - t1
      }
      out
    },
    mc.cores = n_cores, mc.preschedule = FALSE
  )
  bad <- vapply(rows, function(z) inherits(z, "try-error"), logical(1))
  if (any(bad)) {
    warning(sum(bad), " replication(s) errored in ", id, "; dropped. First: ",
            as.character(rows[[which(bad)[1]]]))
    rows <- rows[!bad]
  }
  res <- do.call(rbind, rows)
  el <- proc.time()[["elapsed"]] - t0

  saveRDS(
    list(
      id = id, target = "Table 7", beta2 = b2, n = 400L,
      n_rep_requested = n_rep, n_rep_used = nrow(res),
      B = b_boot, r_grid = GH52_R_GRID, orient = +1, pilot = pilot,
      seed_base = base,
      truth_identity = list(n_big = truth$n_big, c_step = truth$c_step,
                            seed = truth$seed,
                            beta_exact = !is.null(truth$beta_exact),
                            scoring_basis = if (is.null(truth$beta_exact))
                              "smooth" else "exact"),
      elapsed_sec = el,
      sessionInfo = utils::capture.output(utils::sessionInfo()),
      results = res, summary = gh52_summarise(res)
    ),
    f_out
  )
  cat(sprintf("[done] %s  %d reps in %.1f min  -> %s\n",
              id, nrow(res), el / 60, basename(f_out)))
  utils::flush.console()
  invisible(res)
}

pilot_res <- NULL
for (i in seq_along(b2_grid)) {
  r <- run_scenario(b2_grid[i], truths[[i]])
  if (pilot) pilot_res <- r
}

if (pilot && !is.null(pilot_res)) {
  per_rep <- mean(pilot_res$rep_elapsed_sec, na.rm = TRUE)
  total_core_h <- per_rep * 6 * 2000 / 3600
  cat("\n==== PILOT PROJECTION ====\n")
  cat(sprintf("  mean per-replicate (serial)          : %8.1f s  at B = %d\n",
              per_rep, b_boot))
  cat(sprintf("  full study 6 x 2000 reps, single core: %8.1f h\n", total_core_h))
  cat(sprintf("  projected wall-clock at %3d cores    : %8.1f h\n",
              n_cores, total_core_h / n_cores))
  cat("\nREPORT THIS PROJECTION TO THE OWNER BEFORE LAUNCHING PRODUCTION.\n")
  cat("Levers, in order: --B (2000 -> 500), --reps (2000 -> 500, MCSE ~0.010).\n")
  cat("Never thin the cutpoint family.\n")
} else if (!pilot) {
  cat("\nAll scenarios complete.\n")
  cat("Assemble the published-vs-reproduced report with guohe_sec52.qmd.\n")
}
