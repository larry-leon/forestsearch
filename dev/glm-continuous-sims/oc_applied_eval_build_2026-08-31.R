# =============================================================================
# oc_applied_eval_build_2026-08-31.R
# Applied OC evaluation, build phase (task cc_task_oc_applied_evaluation_
# 2026-08-31.md §2-§3): rebuild the anchored context from stage 0's stored
# object, gate it, build the eleven rung DGMs and enumerate their families,
# gate every rung against the first, and stash each family for the rung
# processes.  Read-only against R/.
# Outputs: oc_applied_eval_build_2026-08-31.{log,rds} here, plus one family
# file per rung under $FS_OC_SCRATCH (temporary, not committed).
# =============================================================================
suppressPackageStartupMessages({
  library(forestsearch)
  library(speff2trial)
})
cat("forestsearch version:", as.character(packageVersion("forestsearch")), "\n")
cat("R version:", R.version.string, "\n")

out_dir <- "dev/glm-continuous-sims"
stopifnot(dir.exists(out_dir))
scratch <- Sys.getenv("FS_OC_SCRATCH")
stopifnot(nzchar(scratch), dir.exists(scratch))

res <- list()

# ---------------------------------------------------------------------------
# §2 Rebuild the anchored context from stage 0's stored object, not memory
# ---------------------------------------------------------------------------
res0 <- readRDS(file.path(out_dir, "stage0_oc_applied_2026-08-31.rds"))
beta_treat <- res0$beta_treat            # -26.978725 (baseline GLM)
T_obs      <- res0$anchor$T_obs          # 87.916667 (fitted MD on Hhat)
anchor_ids <- res0$anchor$ids
n_H        <- res0$anchor$n_H
cat(sprintf("[§2] stage 0: beta_treat = %.9f;  T_obs = %.9f;  n_H = %d\n",
            beta_treat, T_obs, n_H))
stopifnot(n_H == 66L)

# Data prep, replicated exactly from stage 0 / the compare-all document
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
confounders.name <- c(cont_vars, bin_vars)
analysis_seed <- 8316951L

actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_df$id <- seq_len(nrow(actg_df))
actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$y_decline <- actg_df$cd40 - actg_df$cd420
for (v in bin_vars) actg_df[[v]] <- as.factor(actg_df[[v]])

N_analysis <- nrow(actg_df)
cat("[§2] Analysis N =", N_analysis, "\n")
if (N_analysis != 1083L || N_analysis != res0$N_analysis) {
  stop("GATE (§2): N != 1083 / stage 0's N.  Stopping.")
}

# The clause mapping (stage 0 §3.2): {age <= 37} -> age = 37;
# !{cd40 <= 507} -> cd40 = list(type = "greater", value = 507)
subgroup_vars <- c("age", "cd40")
subgroup_cuts <- list(age = 37, cd40 = list(type = "greater", value = 507))

# Re-verify once, cheaply: reconstructed flag vs the anchor membership
flag_recon <- (actg_df$age <= 37) & (actg_df$cd40 > 507)
ids_recon  <- actg_df$id[flag_recon]
match_exact <- setequal(ids_recon, anchor_ids) && sum(flag_recon) == n_H
cat(sprintf("[§2] reconstructed flag: n = %d;  matches anchor 66/66: %s\n",
            sum(flag_recon), match_exact))
if (!isTRUE(match_exact)) {
  stop("GATE (§2): reconstructed flag does not match the anchor membership ",
       "exactly.  Stopping.")
}
res$context <- list(N = N_analysis, beta_treat = beta_treat, T_obs = T_obs,
                    n_H = n_H, recon_n = sum(flag_recon), match = match_exact)

# ---------------------------------------------------------------------------
# §3 The rungs: build and gate
# ---------------------------------------------------------------------------
q_rungs <- c(0.01, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 60, T_obs)
res$q_rungs <- q_rungs
cat("[§3] rungs:", paste(format(q_rungs, digits = 9), collapse = ", "), "\n")

build_dgm <- function(q) {
  generate_glm_dgm(
    data            = actg_df,
    factor_vars     = bin_vars,
    continuous_vars = cont_vars,
    outcome_var     = "y_decline",
    treatment_var   = "treat",
    outcome_type    = "continuous",
    effect_measure  = "MD",
    subgroup_vars   = subgroup_vars,
    subgroup_cuts   = subgroup_cuts,
    model           = "alt",
    k_treat         = 1,
    k_inter         = q - beta_treat,
    adverse_outcome = TRUE,
    n_super         = 5000L,
    seed            = analysis_seed,
    verbose         = FALSE
  )
}

fs_args <- list(
  confounders.name = confounders.name,
  conf.cont_jcuts  = list(age = 10, preanti = 10, wtkg = 10,
                          karnof = 10, cd40 = 10, cd80 = 10),
  cut_type         = "default",
  cont.cutoff      = 4L,
  maxk             = 2L,
  n.min            = 60L
)

fam1 <- NULL
gate_rows <- list()
for (i in seq_along(q_rungs)) {
  q <- q_rungs[i]
  t0 <- proc.time()[["elapsed"]]
  dgm <- build_dgm(q)
  sc  <- fs_dgm_scale(dgm)
  reg <- sc$regions
  mQ  <- reg$m_tau[reg$region == "Q"]
  mQc <- reg$m_tau[reg$region == "Qc"]
  d_mQ  <- abs(abs(mQ) - q)
  d_mQc <- abs(mQc - beta_treat)
  if (d_mQ >= 1e-9) {
    stop(sprintf("GATE (§3, rung %d, q = %s): |m_tau[Q]| - q = %.3e >= 1e-9.",
                 i, format(q), d_mQ))
  }
  if (d_mQc >= 1e-9) {
    stop(sprintf("GATE (§3, rung %d, q = %s): m_tau[Qc] moved by %.3e.",
                 i, format(q), d_mQc))
  }

  fam <- fs_oc_family_enumerate(dgm, fs_args, n = N_analysis,
                                max_M = 10000L, verbose = FALSE)
  if (is.null(fam$orientation) || fam$orientation$s != 1) {
    stop(sprintf("GATE (§3, rung %d, q = %s): orientation$s != +1.",
                 i, format(q)))
  }
  if (fam$M != 4508L) {
    stop(sprintf("GATE (§3, rung %d, q = %s): M = %d != 4508.",
                 i, format(q), fam$M))
  }

  if (i == 1L) {
    fam1 <- fam
    lin_dev <- 0
  } else {
    fixed_fields <- c("lab", "Pg", "PQg", "se_g", "sens_g", "spec_g",
                      "ovl", "memb")
    bad <- fixed_fields[!vapply(fixed_fields, function(f)
      identical(fam[[f]], fam1[[f]]), logical(1))]
    if (length(bad)) {
      stop(sprintf(paste0(
        "GATE (§3, rung %d, q = %s): field(s) not identical() to the first ",
        "rung's: %s.  Something depends on k_inter that source says does ",
        "not.  Stopping."), i, format(q), paste(bad, collapse = ", ")))
    }
    lin_dev <- max(abs((fam$beta_g - fam1$beta_g) -
                         (q - q_rungs[1L]) * fam$PQg))
    if (lin_dev >= 1e-9) {
      stop(sprintf(paste0(
        "GATE (§3, rung %d, q = %s): max |beta_g(q) - beta_g(0.01) - ",
        "(q - 0.01) * PQg| = %.3e >= 1e-9.  Stopping."),
        i, format(q), lin_dev))
    }
  }

  fam_file <- file.path(scratch, sprintf("fam_rung_%02d.rds", i))
  saveRDS(fam, fam_file, compress = FALSE)
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf(
    "[§3] rung %2d  q = %-12s  m_tau[Q] ok (%.1e)  m_tau[Qc] ok (%.1e)  s = +1  M = %d  lin_dev = %.2e  %.1f s\n",
    i, format(q, digits = 9), d_mQ, d_mQc, fam$M, lin_dev, el))
  gate_rows[[i]] <- data.frame(
    rung = i, q = q, d_mQ = d_mQ, d_mQc = d_mQc, s = fam$orientation$s,
    M = fam$M, lin_dev = lin_dev, secs = el)
}
res$gates <- do.call(rbind, gate_rows)
res$fs_args <- fs_args
res$subgroup <- list(vars = subgroup_vars, cuts = subgroup_cuts)

saveRDS(res, file.path(out_dir, "oc_applied_eval_build_2026-08-31.rds"))
cat("\n[done] all §2-§3 gates PASS; families stashed under", scratch, "\n")
