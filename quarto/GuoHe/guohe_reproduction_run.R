# guohe_reproduction_run.R
#
# Full-study driver for the Guo & He (2021) Section 5 reproduction. Intended for
# the 127-core host, NOT the chat sandbox.
#
# Parallelism is over REPLICATIONS ONLY -- the outermost loop -- so no nested
# futures are created and `guohe_algorithm3(parallel = )` stays at its default
# FALSE. Each replication pins its own RNG stream via `set.seed()` inside
# `gh_one_rep()`, so the parallel run reproduces a serial run exactly, and
# re-running a single scenario reproduces its earlier result bit-for-bit.
#
# Each scenario writes its own .rds on completion, so the job is restartable:
# a scenario whose output already exists is skipped unless --force is given.
#
# Usage:
#   Rscript quarto/GuoHe/guohe_reproduction_run.R --tables=35 --cores=120
#   Rscript quarto/GuoHe/guohe_reproduction_run.R --tables=6  --cores=120
#   Rscript quarto/GuoHe/guohe_reproduction_run.R --tables=35,6 --adaptive --B=2000
#
# Flags:
#   --tables=35,6   which published tables to target (35 = Tables 3/4/5)
#   --cores=N       forked workers (default: 80% of physical cores)
#   --reps=N        replications per scenario (default 2000, the published value)
#   --B=N           bootstrap resamples (default 2000; the paper does not state
#                   a value for Section 5 -- 2000 is the authors' MONET1 script
#                   value and is an INFERENCE, recorded as such in the output)
#   --adaptive      additionally run Algorithm 2 (dominant cost)
#   --adaptive-B=N  bootstrap resamples inside Algorithm 2 (default: --B)
#   --out=DIR       output directory (default: alongside this script)
#   --force         recompute scenarios whose .rds already exists

suppressMessages(library(survival))

.gh_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
})
source(file.path(.gh_dir, "guohe_reproduction_sim.R"))

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ---- arguments -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(args == paste0("--", nm))
opt <- function(nm, default) {
  hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
}

tables <- strsplit(opt("tables", "35,6"), ",")[[1]]
n_rep <- as.integer(opt("reps", "2000"))
b_boot <- as.integer(opt("B", "2000"))
b_adapt <- as.integer(opt("adaptive-B", as.character(b_boot)))
adaptive <- flag("adaptive")
force <- flag("force")
out_dir <- opt("out", .gh_dir)
n_cores <- as.integer(opt(
  "cores",
  as.character(max(1L, floor(0.80 * parallel::detectCores(logical = FALSE))))
))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Guo & He Section 5 reproduction\n")
cat(sprintf("  tables    : %s\n", paste(tables, collapse = ", ")))
cat(sprintf("  reps      : %d   B = %d%s\n", n_rep, b_boot,
            if (adaptive) sprintf("   adaptive B = %d", b_adapt) else ""))
cat(sprintf("  adaptive  : %s\n", if (adaptive) "yes" else "no"))
cat(sprintf("  cores     : %d\n", n_cores))
cat(sprintf("  out       : %s\n\n", out_dir))

# ---- scenario list ---------------------------------------------------------
scen <- list()
if ("35" %in% tables) {
  for (b2 in c(0, 1, 2, 3, 4, 5) / 10) {
    scen[[length(scen) + 1L]] <- list(
      id = sprintf("t35_beta2_%02d", round(b2 * 10)),
      beta = c(0, b2), n = 400L, target = "Tables 3-5"
    )
  }
}
if ("6" %in% tables) {
  for (k in c(2L, 6L, 10L, 12L)) {
    scen[[length(scen) + 1L]] <- list(
      id = sprintf("t6_k%02d", k),
      beta = rep(0, k), n = 200L * k, target = "Table 6"
    )
  }
}

# ---- run -------------------------------------------------------------------
run_scenario <- function(s) {
  f_out <- file.path(out_dir, paste0("guohe_repro_", s$id, ".rds"))
  if (file.exists(f_out) && !force) {
    cat(sprintf("[skip] %s (exists)\n", s$id))
    return(invisible(NULL))
  }
  # A scenario-specific seed base, so scenarios are independent but each
  # replication is individually reproducible from (id, m).
  base <- 1000000L + as.integer(sum(utf8ToInt(s$id)) * 997L)
  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("[run ] %s  k=%d n=%d  ...\n", s$id, length(s$beta), s$n))
  utils::flush.console()

  rows <- parallel::mclapply(
    seq_len(n_rep),
    function(m) {
      try(gh_one_rep(
        beta = s$beta, n = s$n, r_grid = GH_R_GRID, B = b_boot,
        v = 5L, adaptive = adaptive, seed = base + m
      ), silent = TRUE)
    },
    mc.cores = n_cores, mc.preschedule = FALSE
  )
  bad <- vapply(rows, function(z) inherits(z, "try-error"), logical(1))
  if (any(bad)) {
    warning(sum(bad), " replication(s) errored in ", s$id, "; dropped.")
    rows <- rows[!bad]
  }
  res <- do.call(rbind, rows)
  el <- proc.time()[["elapsed"]] - t0

  saveRDS(
    list(
      id = s$id, target = s$target, beta = s$beta, n = s$n,
      n_rep_requested = n_rep, n_rep_used = nrow(res),
      B = b_boot, adaptive_B = if (adaptive) b_adapt else NA_integer_,
      r_grid = GH_R_GRID, adaptive = adaptive, v = 5L,
      seed_base = base, elapsed_sec = el,
      sessionInfo = utils::capture.output(utils::sessionInfo()),
      results = res, summary = gh_summarise(res)
    ),
    f_out
  )
  cat(sprintf("[done] %s  %d reps in %.1f min  -> %s\n",
              s$id, nrow(res), el / 60, basename(f_out)))
  utils::flush.console()
  invisible(NULL)
}

for (s in scen) run_scenario(s)

cat("\nAll scenarios complete.\n")
cat("Assemble published-vs-reproduced tables with guohe_reproduction.qmd.\n")
