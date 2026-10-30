# s1_reproduce.R -- Guo & He article Stage 0, section 3 (R): reproduction of
# the committed S4 records with the EXPORTED fs_mr_inference() and
# guohe_algorithm3() of the installed 0.4.1 build.
# Task: dev/tasks/cc_task_guohe_article_stage0_2026-10-29_v2.md
#
# Run from the repo root:
#   Rscript dev/reports/guohe_article_stage0_2026-10-29/s1_reproduce.R
# Single core. Reads quarto/GuoHe/ only; writes only under this directory.

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
suppressMessages({ library(forestsearch); library(survival) })
future::plan("sequential")
# chdir = TRUE: the sourced files resolve their siblings from commandArgs()
# --file= (this script's dir, where they are absent) and fall back to getwd();
# chdir makes that fallback quarto/GuoHe. Nothing under quarto/ is written.
source("quarto/GuoHe/mr_vs_guohe_sim.R", chdir = TRUE)

OUT <- "dev/reports/guohe_article_stage0_2026-10-29/results"
KV  <- file.path(OUT, "keyvalues.txt")
GH  <- "quarto/GuoHe"
CAP_SEC <- 15 * 60

kv <- function(name, value) {
  line <- sprintf("%s = %s", name, paste(format(value, digits = 17), collapse = ", "))
  write(line, KV, append = TRUE); cat(line, "\n")
}

main <- function() {
  setTimeLimit(elapsed = CAP_SEC, transient = TRUE)
  t_start <- proc.time()[["elapsed"]]
  stopifnot(isTRUE(all.equal(GH52_R_GRID[2], 1 / 12)))
  kv("s1_GH52_R_GRID", deparse(GH52_R_GRID))

  rows <- list()
  add <- function(id, m, q, stored, recomputed, kind) {
    stored <- unname(stored); recomputed <- unname(recomputed)
    num <- is.numeric(stored) && is.numeric(recomputed)
    rows[[length(rows) + 1L]] <<- data.frame(
      id = id, m = m, quantity = q, kind = kind,
      stored = if (num) sprintf("%.17g", stored) else as.character(stored),
      recomputed = if (num) sprintf("%.17g", recomputed) else as.character(recomputed),
      abs_diff = if (num) abs(stored - recomputed) else NA_real_,
      identical = identical(stored, recomputed),
      all_equal_1e8 = isTRUE(all.equal(stored, recomputed, tolerance = 1e-8)),
      all_equal_1e6 = isTRUE(all.equal(stored, recomputed, tolerance = 1e-6)),
      stringsAsFactors = FALSE)
  }

  last_mr <- last_fit <- NULL
  for (id in c("t7_beta2_00", "t7_beta2_03")) {
    base <- mv_gh52_base(id)
    b2 <- as.integer(sub("^t7_beta2_", "", id)) / 10
    rep_b <- readRDS(sprintf("%s/guohe_repro_%s.rds", GH, id))
    mf_b  <- readRDS(sprintf("%s/mr_field_vs_guohe_%s.rds", GH, id))
    mc_b  <- readRDS(sprintf("%s/mr_field_complement_vs_guohe_%s.rds", GH, id))
    tru   <- readRDS(sprintf("%s/guohe_sec52_truth_beta2_%02d.rds", GH, b2 * 10))
    sb_ok <- rep_b$seed_base == base
    kv(sprintf("s1_%s_seed_base_eq", id), sprintf("%s (stored %s, base %s)", sb_ok,
                                                  rep_b$seed_base, base))
    if (!isTRUE(sb_ok)) stop("GATE: seed_base != base for ", id)

    for (m in 1:2) {
      t_rep <- proc.time()[["elapsed"]]
      rep <- rep_b$results[m, ]
      mf  <- mf_b$results[m, ]
      mc  <- mc_b$results[m, ]
      if (!is.null(mf$m)) stopifnot(mf$m == m)
      if (!is.null(mc$m)) stopifnot(mc$m == m)

      # 1. data and family
      cand <- mv_gh52_regen(id, m)
      k <- length(cand$cuts)
      fits <- gh52_subgroup_fits(cand$df, cand)
      nv <- gh52_naive(fits, cand, tru)
      add(id, m, "naive_c_hat", rep$c_hat_naive, nv$c_hat, "key")
      add(id, m, "n_cand", rep$n_cand, k, "key")
      add(id, m, "naive_point", rep$naive_point, nv$point, "float")
      add(id, m, "naive_lower", rep$naive_lower, nv$lower, "float")
      add(id, m, "gamma_s_naive", rep$gamma_s_naive, nv$gamma_s, "float")

      # 2. Guo & He at r = 1/12, B = 2000 (exported)
      t0 <- proc.time()[["elapsed"]]
      fit <- guohe_algorithm3(data = cand$df, outcome = "survival",
                              treatment = "treat", candidates = cand$names,
                              time = "time", event = "event", orient = +1,
                              B = 2000L, r = 1 / 12, level = 0.05,
                              seed = base + m + 500000L, min_events = 5L,
                              diagnostics = FALSE)
      t_gh <- proc.time()[["elapsed"]] - t0
      add(id, m, "gh_c_hat", rep$c_hat_gh,
          cand$cuts[match(fit$selected, cand$names)], "key")
      add(id, m, "gh_r2_lower", rep$gamma_s - rep$r2_dist,
          log(fit$bound_one_sided), "float")
      add(id, m, "gh_r2_debiased", rep$r2_bias + rep$gamma_s,
          log(fit$debiased), "float")
      tag <- sprintf("s1_%s_m%d", id, m)
      kv(paste0(tag, "_gh_lean_cox"), fit$lean_cox)
      kv(paste0(tag, "_gh_fast"), fit$fast)
      kv(paste0(tag, "_gh_B"), fit$B)
      kv(paste0(tag, "_gh_n_estimable"), fit$n_estimable)

      # 3. MR, field and complement (exported)
      sel_lab <- cand$names[which.max(replace(fits$est, !is.finite(fits$est), -Inf))]
      cidx <- mv_cand_idx_52(cand)
      t0 <- proc.time()[["elapsed"]]
      mr <- fs_mr_inference(df = cand$df, candidates = cidx, spec = mv_spec52,
                            selected_members = cidx[[sel_lab]],
                            admission = list(effect_floor = NULL, consistency = NULL),
                            reselection = "maxeff", draws = 5000L,
                            multiplier = "poisson", ci_method = "field",
                            seed = base + m + 700000L, return_reselection = TRUE,
                            include_complement = TRUE, field_R_out = 1000L,
                            field_R_in = 500L, field_uniform = FALSE,
                            field_complement = TRUE,
                            field_scale_complement = "selected")
      t_mr <- proc.time()[["elapsed"]] - t0
      add(id, m, "mr_selected_index", mf$sel, mr$selected_index, "key")
      add(id, m, "mr_est", mf$mr_est, log(mr$debiased$est), "float")
      add(id, m, "mr_se_ij", mf$mr_se_ij, mr$debiased$se_ij, "float")
      add(id, m, "mr_lower_1s", mf$mr_lower_1s, log(mr$debiased$lower_1s), "float")
      add(id, m, "mr_bias_sel", mf$mr_bias_sel, mr$selection_bias, "float")
      add(id, m, "mr_bias_fix", mf$mr_bias_fix, mr$fixed_bias, "float")
      add(id, m, "fld_est2", mf$fld_est2, log(mr$field$est2), "float")
      add(id, m, "fld_lower_1s", mf$fld_lower_1s, log(mr$field$lower_1s), "float")
      add(id, m, "fld_lambda_mean", mf$fld_lambda_mean, mr$field$lambda_mean, "float")
      add(id, m, "p_hat_H", mf$p_hat_H,
          unname(mr$reselection$p_hat[mr$selected_index]), "float")
      # complement / joint, column names as emitted by
      # mr_field_complement_vs_guohe_run.R (.mfc_mr_cols)
      add(id, m, "cgate_est", mc$cgate_est, log(mr$complement$debiased$est), "float")
      add(id, m, "c_upper_95_s", mc$c_upper_95_s,
          log(mr$field$complement$upper_1s_s), "float")
      add(id, m, "h_lower_975 (joint_s bonf_lower_H)", mc$h_lower_975,
          log(mr$field$joint_s$bonf_lower_H), "float")
      add(id, m, "c_upper_975_s", mc$c_upper_975_s,
          log(mr$field$joint_s$bonf_upper_Hc), "float")
      js <- mr$field$joint_s
      kv(paste0(tag, "_joint_s_bonf_joint_prob"),
         sprintf("%.17g (stored joint_bonf_prob_s %.17g)", js$bonf_joint_prob,
                 mc$joint_bonf_prob_s))
      kv(paste0(tag, "_joint_s_gamma"),
         sprintf("%.17g (stored joint_gamma_s %.17g)", js$gamma, mc$joint_gamma_s))
      kv(paste0(tag, "_joint_s_corr"),
         sprintf("%.17g (stored joint_corr_s %.17g)", js$corr, mc$joint_corr_s))
      kv(paste0(tag, "_sel_lab"), sel_lab)
      kv(paste0(tag, "_sec_gh"), round(t_gh, 2))
      kv(paste0(tag, "_sec_mr"), round(t_mr, 2))
      kv(paste0(tag, "_sec_rep"), round(proc.time()[["elapsed"]] - t_rep, 2))
      last_mr <- mr; last_fit <- fit
      if (proc.time()[["elapsed"]] - t_start > CAP_SEC)
        stop("CAP: section 3 exceeded ", CAP_SEC, " s")
    }
  }

  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(OUT, "reproduction.csv"), row.names = FALSE)

  keys <- res[res$kind == "key", ]
  flts <- res[res$kind == "float", ]
  kv("s1_keys_identical", sprintf("%d / %d", sum(keys$identical), nrow(keys)))
  kv("s1_floats_identical", sprintf("%d / %d", sum(flts$identical), nrow(flts)))
  kv("s1_floats_all_equal_1e8", sprintf("%d / %d", sum(flts$all_equal_1e8), nrow(flts)))
  kv("s1_floats_all_equal_1e6", sprintf("%d / %d", sum(flts$all_equal_1e6), nrow(flts)))
  imax <- which.max(replace(flts$abs_diff, is.na(flts$abs_diff), -Inf))
  kv("s1_max_abs_diff", sprintf("%.6g (%s, %s m=%d)", flts$abs_diff[imax],
                                flts$quantity[imax], flts$id[imax], flts$m[imax]))
  mid <- flts[!flts$all_equal_1e8 & flts$all_equal_1e6, ]
  if (nrow(mid)) {
    j <- which.max(mid$abs_diff)
    kv("s1_finding_fail_1e8_pass_1e6", sprintf("%d floats; largest %.6g (%s, %s m=%d)",
       nrow(mid), mid$abs_diff[j], mid$quantity[j], mid$id[j], mid$m[j]))
  } else kv("s1_finding_fail_1e8_pass_1e6", "0 floats")
  kv("s1_elapsed_sec", round(proc.time()[["elapsed"]] - t_start, 1))

  # ---- package facts on the last mr and fit ----
  dp <- function(x) paste(deparse(x), collapse = " ")
  fm <- formals(fs_mr_inference)[c("draws", "multiplier", "ci_method",
         "include_complement", "field_R_out", "field_R_in", "field_complement",
         "field_scale_complement", "reselection")]
  fg <- formals(guohe_algorithm3)[c("orient", "B", "r", "level", "min_events",
                                    "fast", "parallel")]
  nm_line <- function(lab, x) sprintf("- `%s`: %s", lab,
                                      paste0("`", names(x), "`", collapse = ", "))
  facts <- c(
    "# Package facts (installed forestsearch 0.4.1; last replicate t7_beta2_03 m = 2)",
    "", "## Element names", "",
    nm_line("names(mr)", last_mr),
    nm_line("names(mr$debiased)", last_mr$debiased),
    nm_line("names(mr$field)", last_mr$field),
    nm_line("names(mr$field$complement)", last_mr$field$complement),
    nm_line("names(mr$field$joint_s)", last_mr$field$joint_s),
    nm_line("names(mr$complement)", last_mr$complement),
    nm_line("names(mr$reselection)", last_mr$reselection),
    nm_line("names(fit)", last_fit),
    "", "## Defaults: fs_mr_inference()", "",
    sprintf("- `%s = %s`", names(fm), vapply(fm, dp, "")),
    "", "## Defaults: guohe_algorithm3()", "",
    sprintf("- `%s = %s`", names(fg), vapply(fg, dp, "")))
  writeLines(facts, file.path(OUT, "package_facts.md"))

  gate_ok <- all(keys$identical) && all(flts$all_equal_1e6)
  kv("s1_GATE", if (gate_ok) "PASS" else "FAIL")
  gate_ok
}

ok <- main()
if (!isTRUE(ok)) quit(status = 2)
