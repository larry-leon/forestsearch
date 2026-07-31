# =============================================================================
# Precondition (brief S5): is `mr_in_replicates` RNG-neutral?
# =============================================================================
# With the outer seed fixed, run a short bootstrap twice -- once with
# mr_in_replicates = TRUE, once FALSE -- and compare Ystar_mat, the
# per-replicate selected subgroups, and the bias-corrected estimates.
#
# Read-only: nothing here writes package code, .qmd, or tests.

suppressPackageStartupMessages({
  library(forestsearch); library(survival); library(data.table)
  library(grf); library(policytree); library(doFuture)
})

OUT <- "dev/replication-check/out"
QMD <- "quarto/applications/gbsg/analysis_gbsg_cox_multimethod_psi_v2_2new.qmd"

# --- pull a chunk verbatim out of the rendered notebook ----------------------
chunk_src <- function(label, path = QMD) {
  L <- readLines(path, warn = FALSE)
  i <- grep(sprintf("^```\\{r %s\\}\\s*$", label), L)
  if (length(i) != 1L) stop("chunk '", label, "' not found uniquely")
  j <- i + which(grepl("^```\\s*$", L[(i + 1):length(L)]))[1]
  L[(i + 1):(j - 1)]
}
eval_chunk <- function(label, env = globalenv())
  eval(parse(text = chunk_src(label)), envir = env)

# --- notebook setup, verbatim ------------------------------------------------
eval_chunk("setup"); eval_chunk("parallel-setup"); eval_chunk("data-setup")

# Confirm the identification settings the brief pins (S3)
eval_chunk("setup")  # timings/opts
src_setup <- chunk_src("parallel-setup")
eval(parse(text = src_setup))

cat("\n--- config check ---\n")
cfg <- c(fs_sg_focus = "effMaxSG", fs_hr_threshold = "1", fs_hr_consistency = "1",
         fs_max_subgroups_search = "50", fs_consistency_method = "resample",
         fs_subgroup_method = "consistency", run_cv = "FALSE", run_loo = "FALSE")
for (nm in names(cfg))
  cat(sprintf("  %-24s %s\n", nm, paste(get(nm), collapse = ", ")))
cat(sprintf("  %-24s %s\n", "fs_conf_force", paste(fs_conf_force, collapse = " | ")))
cat(sprintf("  %-24s %s\n", "NB", NB))
cat(sprintf("  %-24s %s\n", "mr_draws", mr_draws))
cat(sprintf("  %-24s %s\n", "n_cores", n_cores))

# --- the main fit, verbatim from the notebook --------------------------------
# forestsearch-main chunk, minus the plotting knobs we cannot silence.
cat("\n--- fitting main forestsearch (mr_inference = TRUE) ---\n")
t0 <- proc.time()
eval_chunk("forestsearch-main")
t_fit <- (proc.time() - t0)["elapsed"]
cat("main fit elapsed:", round(t_fit, 1), "s\n")
cat("selected subgroup:", paste(fs$sg.harm, collapse = " & "), "\n")

saveRDS(list(sg.harm = fs$sg.harm, elapsed = unname(t_fit)),
        file.path(OUT, "precond_mainfit.rds"))

# --- the two short bootstraps ------------------------------------------------
NB_SHORT <- 5L
run_boot <- function(mr_in_rep) {
  cat("\n--- bootstrap nb_boots =", NB_SHORT,
      ", mr_in_replicates =", mr_in_rep, "---\n")
  t0 <- proc.time()
  r <- forestsearch_bootstrap_dofuture(
    fs.est = fs, nb_boots = NB_SHORT, show_three = FALSE, details = FALSE,
    mr_in_replicates = mr_in_rep,
    parallel_args = list(plan = "multisession", workers = NB_SHORT,
                         show_message = FALSE)
  )
  future::plan("sequential")
  attr(r, "elapsed") <- unname((proc.time() - t0)["elapsed"])
  r
}

bF <- run_boot(FALSE)
bT <- run_boot(TRUE)

# --- compare -----------------------------------------------------------------
cmp <- list()
cmp$ystar_identical <- identical(bF$Ystar_mat, bT$Ystar_mat)
cmp$ystar_allequal  <- isTRUE(all.equal(bF$Ystar_mat, bT$Ystar_mat))

# per-replicate selected subgroups: the `results` table carries the
# per-replicate subgroup definition columns.
resF <- as.data.frame(bF$results); resT <- as.data.frame(bT$results)
cmp$results_cols <- names(resF)
sg_cols <- intersect(names(resF),
                     c("sg1_name", "sg2_name", "sg_name", "subgroup", "K"))
cmp$sg_cols <- sg_cols
cmp$sg_identical <- if (length(sg_cols))
  identical(resF[, sg_cols, drop = FALSE], resT[, sg_cols, drop = FALSE]) else NA

# full per-replicate table, numeric tolerance-free
common <- intersect(names(resF), names(resT))
cmp$results_identical <- identical(resF[, common, drop = FALSE],
                                   resT[, common, drop = FALSE])
cmp$results_allequal <- isTRUE(all.equal(resF[, common, drop = FALSE],
                                         resT[, common, drop = FALSE]))
# which columns differ, if any
cmp$diff_cols <- common[!vapply(common, function(c)
  isTRUE(all.equal(resF[[c]], resT[[c]])), logical(1))]

# bias-corrected estimates
cmp$H_identical  <- isTRUE(all.equal(bF$H_estimates,  bT$H_estimates))
cmp$Hc_identical <- isTRUE(all.equal(bF$Hc_estimates, bT$Hc_estimates))
cmp$SGCIs_identical <- isTRUE(all.equal(bF$SG_CIs, bT$SG_CIs))
cmp$FSsg_identical  <- isTRUE(all.equal(bF$FSsg_tab, bT$FSsg_tab))

cmp$mr_replicates_F <- is.null(bF$mr_replicates)
cmp$mr_replicates_T_n <- if (is.null(bT$mr_replicates)) 0L else
  sum(!vapply(bT$mr_replicates, is.null, logical(1)))

cmp$elapsed_F <- attr(bF, "elapsed"); cmp$elapsed_T <- attr(bT, "elapsed")

cat("\n=================== PRECONDITION RESULT ===================\n")
str(cmp)

verdict <- isTRUE(cmp$ystar_identical) && isTRUE(cmp$results_allequal) &&
  isTRUE(cmp$H_identical) && isTRUE(cmp$Hc_identical) &&
  isTRUE(cmp$SGCIs_identical) && isTRUE(cmp$FSsg_identical)
cat("\nRNG-NEUTRAL:", verdict, "\n")

saveRDS(list(cmp = cmp, verdict = verdict,
             bF = bF[c("results", "SG_CIs", "FSsg_tab", "H_estimates", "Hc_estimates")],
             bT = bT[c("results", "SG_CIs", "FSsg_tab", "H_estimates", "Hc_estimates")]),
        file.path(OUT, "precond_rng.rds"))
cat("\nsaved:", file.path(OUT, "precond_rng.rds"), "\n")
