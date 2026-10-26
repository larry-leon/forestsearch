# s3_bootstrap_runs.R -- the application-shape measurements for candidates C
# (multicore vs multisession startup) and D (tail imbalance of the outer
# chunking), from the bootstrap's own recorded per-replicate timings.
#
# Runs (fixture V fit once as fs.est):
#   S   sequential,    nb = 5            -> sequential per-replicate seconds now
#   M   multisession,  nb = 4*W          -> startup + makespan decomposition
#   C   multicore,     nb = 4*W          -> same, forked workers
#   Mi  multisession,  nb = 20, 20 wrk   -> identity pair across plans
#   Ci  multicore,     nb = 20, 20 wrk
# W = min(48, parallelly::availableCores()).  Every run uses the package's
# public entry point with user-level arguments; nothing is modified.
#
# Decomposition per run, all from recorded fields:
#   outer_s      = system.time() around the call
#   inner_s      = attr(results, "timing")$total_minutes * 60  (bootstrap_results()
#                  only: worker export/package load + the replicate loop)
#   spawn_s      = outer_s - inner_s                    (plan setup, main-process fits, Ystar)
#   sum_rep_s    = sum(tmins_iteration) * 60            (total replicate CPU-seconds)
#   contig_s     = max over contiguous chunks (parallel::splitIndices(nb, W)) of the
#                  chunk's replicate seconds  -- the static-chunking makespan
#   balanced_s   = sum_rep_s / W
#   overhead_s   = inner_s - contig_s  (worker export + package load + globals discovery + collection)

source("stage0_common.R")
fx <- stage0_fixture()
fs <- run_fixture_call(fx)
stopifnot(!is.null(fs$sg.harm))
kv(boot_fixture_sg = paste(fs$sg.harm, collapse = " & "))

W  <- as.integer(min(48L, parallelly::availableCores()))
NB <- 4L * W
kv(boot_workers_W = W, boot_nb_main = NB)

volatile_cols <- c("tmins_search", "tmins_iteration")
strip_volatile <- function(b) {
  r <- as.data.frame(b$results)
  r <- r[, setdiff(names(r), volatile_cols), drop = FALSE]
  list(results = r, SG_CIs = b$SG_CIs, FSsg_tab = b$FSsg_tab, Ystar_mat = b$Ystar_mat,
       H_estimates = b$H_estimates, Hc_estimates = b$Hc_estimates,
       nb_boots = b$nb_boots, original_sg = b$original_sg)
}

run_boot <- function(plan, workers, nb, tag) {
  pa <- list(plan = plan, workers = workers, show_message = FALSE)
  t0 <- proc.time()[["elapsed"]]
  b <- suppressWarnings(suppressMessages(
    forestsearch::forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = nb,
                                                  seed = 8316951L, parallel_args = pa)))
  outer_s <- proc.time()[["elapsed"]] - t0
  future::plan("sequential")
  tm <- attr(b$results, "timing")
  inner_s <- tm$total_minutes * 60
  rep_s <- b$results$tmins_iteration * 60
  rep_s <- rep_s[order(b$results$boot_id)]
  idx <- parallel::splitIndices(nb, min(workers, nb))
  contig_s <- max(vapply(idx, function(ii) sum(rep_s[ii], na.rm = TRUE), numeric(1)))
  sum_rep_s <- sum(rep_s, na.rm = TRUE)
  stats <- list(tag = tag, plan = plan, workers = workers, nb = nb,
                outer_s = outer_s, inner_s = inner_s, spawn_s = outer_s - inner_s,
                sum_rep_s = sum_rep_s, mean_rep_s = mean(rep_s, na.rm = TRUE),
                max_rep_s = max(rep_s, na.rm = TRUE),
                contig_s = contig_s, balanced_s = sum_rep_s / min(workers, nb),
                overhead_s = inner_s - contig_s,
                n_found = sum(b$results$any_found == 1L, na.rm = TRUE),
                n_rows = nrow(b$results))
  do.call(kv, setNames(stats[-1L], paste0("boot_", tag, "_", names(stats)[-1L])))
  list(stats = stats, obj = b, rep_s = rep_s)
}

runs <- list()
runs$S  <- run_boot("sequential",   1L,  5L, "S")
runs$M  <- run_boot("multisession", W,   NB, "M")
runs$C  <- run_boot("multicore",    W,   NB, "C")
W2 <- as.integer(min(20L, W))
runs$Mi <- run_boot("multisession", W2, 20L, "Mi")
runs$Ci <- run_boot("multicore",    W2, 20L, "Ci")

# identity across plans (volatile timing columns/attributes excluded)
id_pair <- identical(strip_volatile(runs$Mi$obj), strip_volatile(runs$Ci$obj))
id_results_only <- identical(strip_volatile(runs$Mi$obj)$results, strip_volatile(runs$Ci$obj)$results)
kv(identity_multicore_vs_multisession_nb20 = id_pair,
   identity_results_table_only = id_results_only)
if (!id_pair) {
  s1 <- strip_volatile(runs$Mi$obj); s2 <- strip_volatile(runs$Ci$obj)
  cat("components differing:", paste(names(s1)[!mapply(identical, s1, s2)], collapse = ", "), "\n")
}

# per-replicate slowdown under load vs sequential
kv(rep_slowdown_M = runs$M$stats$mean_rep_s / runs$S$stats$mean_rep_s,
   rep_slowdown_C = runs$C$stats$mean_rep_s / runs$S$stats$mean_rep_s,
   startup_removed_by_multicore_s = (runs$M$stats$spawn_s + runs$M$stats$overhead_s) -
                                    (runs$C$stats$spawn_s + runs$C$stats$overhead_s),
   imbalance_M_s = runs$M$stats$contig_s - runs$M$stats$balanced_s,
   imbalance_C_s = runs$C$stats$contig_s - runs$C$stats$balanced_s)

stats_tab <- do.call(rbind, lapply(runs, function(r) as.data.frame(r$stats, stringsAsFactors = FALSE)))
utils::write.csv(stats_tab, out_path("bootstrap_runs.csv"), row.names = FALSE)
print(stats_tab, row.names = FALSE)
save_result(list(stats = stats_tab, rep_s = lapply(runs, `[[`, "rep_s"), W = W, NB = NB,
                 id_pair = id_pair), "s3_bootstrap_runs")
cat("\ns3 complete\n")
