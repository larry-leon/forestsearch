# =============================================================================
# residual_2026-08-30.R -- the two candidate explanations for the remaining
# between-rule residual (task: dev/tasks/cc_task_residual_2026-08-30.md)
# -----------------------------------------------------------------------------
# Read-only with respect to the package: no R/ change, no edit to any driver or
# document.  Reads sgdef_selection_2026-08-29.rds (the realized rules, already
# evaluated on df_super) and oc_wrapper_grid_corrected_2026-08-30.rds (the
# corrected 13-variable families and their stored fs_oc_predict() results).
#
# Parts (run from the repository root):
#   Rscript dev/glm-continuous-sims/residual_2026-08-30.R A
#       Hypothesis A -- the wrapper's se_g = seQ1000*sqrt(1000/n)*sqrt(piQ/Pg)
#       against fs_dgm_scale(dgm, regions = <the rule's own membership>).
#       Writes residual_hypA_2026-08-30.rds.
#   OC_PART_DIR=<dir> Rscript dev/glm-continuous-sims/residual_2026-08-30.R B alt500
#   OC_PART_DIR=<dir> Rscript dev/glm-continuous-sims/residual_2026-08-30.R B alt700
#   OC_PART_DIR=<dir> Rscript dev/glm-continuous-sims/residual_2026-08-30.R B null500
#       Hypothesis B -- the dedup analogue on the corrected family and the
#       fs_oc_predict() re-runs (both gates, seed 20260825, draws 2e5).
#   OC_PART_DIR=<dir> Rscript dev/glm-continuous-sims/residual_2026-08-30.R Bsum
#       Merge the three parts, the between-rule gap and the comparison against
#       measurement.  Writes residual_hypB_2026-08-30.rds.
#
# The dedup key, from source (R/subgroup_consistency_helpers.R L1102-1131):
#   remove_near_duplicate_subgroups() rounds columns 2:min(10, ncol) of
#   found.hrs to `tolerance = 0.001`, pastes them into a string key and keeps
#   the first row of each key (`!duplicated`).  Columns 2:10 are, by
#   format_search_results() (R/subgroup_search.R L878-879),
#     K, n, E, d1, m1, m0, HR, L(HR), U(HR).
#   On the GLM / continuous path (R/subgroup_search.R L640-647, L815-820):
#     E  = n0 + n1 = n,  d1 = n1 (the treated count),  m1 = m0 = NA,
#     HR = the estimate, L/U = estimate -/+ 1.96 * se.
#   found.hrs is sorted (-HR, K) before the dedup (format_search_results L883),
#   so "first" = highest estimate, then fewest factors, then enumeration order.
#
# Two population analogues are built, and labelled:
#   v1 "literal"   -- the population counterpart of every key column at the same
#                     tolerance: K, n*Pg, n*Pg (E), n*Pg/2 (d1), NA, NA, beta_g,
#                     beta_g -/+ 1.96 se_g.  Since Pg = |g|/5000 and
#                     beta_g = tauQc + bint * |g & Q| / |g|, matching to 0.001
#                     is matching (K, |g|, |g & Q|) exactly (null: (K, |g|)).
#                     Representative: the first after ordering (-beta_g, K,
#                     enumeration index), as the package orders.
#   v2 "near-twin" -- APPROXIMATE.  In a trial the sample key (n, n1, estimate
#                     to 0.001) can only coincide for two distinct population
#                     memberships when their SAMPLE memberships coincide, which
#                     for a symmetric difference of d super-population subjects
#                     happens with probability (1 - d/5000)^n.  v2 collapses,
#                     within K, pairs with d <= d_max where that probability is
#                     >= 1/2 (d_max = floor(5000 * log(2) / n)), keeping the
#                     higher beta_g.  It is a deterministic stand-in for a
#                     per-replicate random event and is reported as such.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
part <- args[1]; stopifnot(part %in% c("A", "B", "Bsum"))
PART_DIR <- Sys.getenv("OC_PART_DIR", unset = tempdir())
part_path <- function(cell) file.path(PART_DIR, sprintf("residual_hypB_part_%s.rds", cell))
D <- "dev/glm-continuous-sims"
S <- readRDS(file.path(D, "sgdef_selection_2026-08-29.rds"))
C <- readRDS(file.path(D, "oc_wrapper_grid_corrected_2026-08-30.rds"))
SEED <- C$seed; DRAWS <- C$draws
stopifnot(SEED == 20260825L, DRAWS == 2e5)
fs_args <- C$fs_args
stopifnot(length(fs_args$confounders.name) == 13L, fs_args$confounders.name[13] == "str2")
TOL <- 0.001                                  # remove_near_duplicate_subgroups() default

# ---- the DGMs, rebuilt exactly as the 08-29 / 08-30 scripts build them --------
actg_frame <- function() {
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- 1L - ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$z1  <- as.factor(ifelse(actg_df$age > 34, 1L, 0L))
  actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
  actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
  actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
  actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
  actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
  actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
  actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
  actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
  actg_df
}
build_dgm <- function(model, k_inter) generate_glm_dgm(
  data = actg_frame(), factor_vars = paste0("z", 1:12), outcome_var = "cd4_change",
  treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = model, k_treat = 1, k_inter = k_inter, n_super = 5000L,
  seed = 8316951L, verbose = FALSE)
read_pay <- function(p) { z <- readRDS(p); nm <- names(z); nm[is.na(nm)] <- ""; setNames(z, nm) }

cell_info <- list(
  alt500  = list(n = 500L, null = FALSE, fam = C$alt$families[["500"]], runs = C$alt$runs[c("n500_resample", "n500_split")], meas = C$alt$measured[["500"]], payload = C$alt$payloads[["500"]]),
  alt700  = list(n = 700L, null = FALSE, fam = C$alt$families[["700"]], runs = C$alt$runs[c("n700_resample", "n700_split")], meas = C$alt$measured[["700"]], payload = C$alt$payloads[["700"]]),
  null500 = list(n = 500L, null = TRUE,  fam = C$null$family,           runs = C$null$runs,                                  meas = C$null$measured,           payload = C$null$payload))
for (cl in names(cell_info)) names(cell_info[[cl]]$runs) <- c("resample", "split")

dgm_for <- function(cell) {
  pl <- read_pay(cell_info[[cell]]$payload); truth <- pl[["truth"]]
  if (cell_info[[cell]]$null) {
    dgm <- build_dgm("null", 0)
    stopifnot(sum(dgm$df_super$flag_harm) == 0L, is.na(dgm$hazard_ratios$harm_subgroup),
              abs(dgm$hazard_ratios$no_harm_subgroup - truth$effect_Qc) < 1e-9)
  } else {
    dgm <- build_dgm("alt", truth$beta_inter)
    stopifnot(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q) < 1e-9,
              abs(dgm$subgroup_info$proportion - truth$prevalence_Q) < 1e-9)
  }
  dgm
}

# ---- signature parser, verbatim from sgdef_selection_2026-08-29.R ------------
clause_sig <- function(cl) {
  cl <- trimws(cl); neg <- grepl("^!", cl)
  cl <- gsub("^!\\s*", "", cl); cl <- gsub("^\\{|\\}$", "", cl); cl <- trimws(cl)
  if (grepl("<=", cl, fixed = TRUE)) {
    v <- trimws(strsplit(cl, "<=", fixed = TRUE)[[1]][1]); return(paste0(v, if (neg) ">" else "<="))
  }
  if (grepl(">", cl, fixed = TRUE) && !grepl(">=", cl, fixed = TRUE)) {
    v <- trimws(strsplit(cl, ">", fixed = TRUE)[[1]][1]); return(paste0(v, if (neg) "<=" else ">"))
  }
  if (grepl("!=", cl, fixed = TRUE)) {
    v <- trimws(strsplit(cl, "!=", fixed = TRUE)[[1]][1]); return(paste0(v, if (neg) "=1" else "=0"))
  }
  if (grepl("==", cl, fixed = TRUE)) {
    v <- trimws(strsplit(cl, "==", fixed = TRUE)[[1]][1]); return(paste0(v, if (neg) "=0" else "=1"))
  }
  paste0(cl, if (neg) "=0" else "=1")
}
rule_sig <- function(rule) {
  cls <- strsplit(rule, " & ", fixed = TRUE)[[1]]
  paste(sort(vapply(cls, clause_sig, character(1))), collapse = " & ")
}

# =============================================================================
# PART A -- se_g prevalence scaling vs the candidate's own V_eff
# =============================================================================
run_A <- function() {
  out <- list()
  dgm_alt <- dgm_for("alt500"); dgm_null <- dgm_for("null500")
  fr_alt  <- fs_build_eval_frame(dgm_alt,  outcome_type = "continuous")
  fr_null <- fs_build_eval_frame(dgm_null, outcome_type = "continuous")
  for (cell in names(cell_info)) {
    ci <- cell_info[[cell]]; n <- ci$n
    dgm <- if (ci$null) dgm_null else dgm_alt
    frame <- if (ci$null) fr_null else fr_alt
    df <- dgm$df_super
    stopifnot(nrow(frame) == nrow(df),
              identical(frame$age, df$age), identical(frame$preanti, df$preanti),
              identical(frame$flag_harm, df$flag_harm))
    inQ <- df$flag_harm == 1
    # the wrapper's reference row and constants, exactly as fs_oc_family_enumerate()
    if (!ci$null) {
      sc <- fs_dgm_scale(dgm); reg <- sc$regions
      piQ <- reg$P_g[reg$region == "Q"]; V_ref <- reg$V_eff[reg$region == "Q"]; ref_lab <- "Q"
    } else {
      sc <- fs_dgm_scale(dgm, regions = list(S = rep(TRUE, nrow(df)))); reg <- sc$regions
      piQ <- 1; V_ref <- reg$V_eff[reg$region == "S"]; ref_lab <- "S"
    }
    seQ1000 <- sqrt(V_ref / (1000 * piQ))
    stopifnot(isTRUE(all.equal(seQ1000 * sqrt(1000 / n) * sqrt(piQ / ci$fam$Pg), ci$fam$se_g)))
    cat(sprintf("\n############ %s: n = %d, reference row %s: V_eff %.4f, P %.4f, seQ1000 %.5f\n",
                cell, n, ref_lab, V_ref, piQ, seQ1000))

    # ---- the realized rules, from the sgdef .rds, memberships re-resolved -----
    rules <- S$cells[[cell]]$rules
    stopifnot(all(rules$status == "ok"))
    memb <- lapply(rules$sg_def, function(r) {
      m <- forestsearch:::.fs_resolve_membership(frame, r); stopifnot(m$status == "ok"); m$in_region })
    names(memb) <- sprintf("r%04d", seq_along(memb))
    Pg_chk <- vapply(memb, mean, 1)
    stopifnot(isTRUE(all.equal(unname(Pg_chk), rules$Pg_pop, tolerance = 1e-12)))   # the stored evaluation is reproduced
    scr <- fs_dgm_scale(dgm, regions = memb)$regions
    A <- data.frame(
      sg_def = rules$sg_def, sig = rules$sig,
      K = lengths(strsplit(rules$sg_def, " & ", fixed = TRUE)),
      Pg = rules$Pg_pop, nPg = n * rules$Pg_pop, PQg = rules$PQg_pop,
      jacc_Q = if (ci$null) NA_real_ else vapply(memb, function(g) sum(g & inQ) / sum(g | inQ), 1),
      V_eff = scr$V_eff, V_mu0 = scr$V_mu0, V_mu1 = scr$V_mu1,
      se_scaled = seQ1000 * sqrt(1000 / n) * sqrt(piQ / rules$Pg_pop),
      se_direct = sqrt(scr$V_eff / (n * scr$P_g)),
      stringsAsFactors = FALSE)
    A$ratio <- A$se_scaled / A$se_direct
    stopifnot(isTRUE(all.equal(A$ratio, sqrt(V_ref / A$V_eff))))       # ratio = sqrt(V_eff[ref] / V_eff[g]); n-free
    # side column (NOT the task's ratio): the Jensen-inflated unconditional sd at the candidate's own V_eff
    jf <- vapply(seq_len(nrow(scr)), function(i) forestsearch:::.fs_jensen_factor(
      scr$P_g[i], n = n, p_treat = 0.5, V0 = scr$V_arm0[i], V1 = scr$V_arm1[i]), 1)
    A$se_direct_jensen <- A$se_direct * sqrt(jf)
    A$ratio_jensen <- A$se_scaled / A$se_direct_jensen
    # weights: how often each rule was selected
    reps <- S$cells[[cell]]$reps
    A$n_sel <- as.integer(table(factor(reps$sg_def, levels = rules$sg_def)))

    # ---- the same on the corrected analytic family (all M candidates) ----------
    fam <- ci$fam
    # membership of every family candidate: re-enumerate (10 s) and guard against the stored record
    fe <- fs_oc_family_enumerate(dgm, fs_args, n = n, max_M = 5000L, verbose = FALSE)
    stopifnot(identical(fe$lab, fam$lab), isTRUE(all.equal(fe$Pg, fam$Pg)),
              isTRUE(all.equal(fe$se_g, fam$se_g)), isTRUE(all.equal(fe$beta_g, fam$beta_g)))
    fam_memb <- lapply(seq_len(fe$M), function(j) as.logical(fe$memb[, j]))
    names(fam_memb) <- sprintf("f%04d", seq_len(fe$M))
    scf <- fs_dgm_scale(dgm, regions = fam_memb)$regions
    F <- data.frame(lab = fe$lab, K = lengths(strsplit(fe$lab, " & ", fixed = TRUE)),
                    Pg = fe$Pg, nPg = n * fe$Pg, PQg = fe$PQg,
                    jacc_Q = if (ci$null) NA_real_ else vapply(fam_memb, function(g) sum(g & inQ) / sum(g | inQ), 1),
                    V_eff = scf$V_eff, se_scaled = fe$se_g, se_direct = sqrt(scf$V_eff / (n * scf$P_g)),
                    sel_c_resample = ci$runs$resample$pred$sel_c, sel_c_split = ci$runs$split$pred$sel_c,
                    stringsAsFactors = FALSE)
    F$ratio <- F$se_scaled / F$se_direct

    # ---- summaries -------------------------------------------------------------
    qs <- function(x) round(stats::quantile(x, c(0, .01, .05, .25, .5, .75, .95, .99, 1)), 5)
    cat("ratio se_scaled / se_direct over distinct realized rules:\n"); print(qs(A$ratio))
    cat(sprintf("  selection-weighted (by replicate count) mean %.5f\n", weighted.mean(A$ratio, A$n_sel)))
    cat("ratio over the corrected analytic family (all M):\n"); print(qs(F$ratio))
    cat(sprintf("  sel_c-weighted mean (resample) %.5f\n", weighted.mean(F$ratio, F$sel_c_resample)))
    trend <- function(d, lab) {
      cat(sprintf("trend [%s]: cor(ratio, Pg) %+.3f; cor(ratio, PQg) %s; cor(ratio, Jaccard) %s\n", lab,
                  cor(d$ratio, d$Pg),
                  if (ci$null) "NA (null)" else sprintf("%+.3f", cor(d$ratio, d$PQg)),
                  if (ci$null) "NA (null)" else sprintf("%+.3f", cor(d$ratio, d$jacc_Q))))
      fit <- if (ci$null) lm(ratio ~ Pg, d) else lm(ratio ~ Pg + PQg, d)
      print(round(coef(summary(fit)), 5))
      bP <- cut(d$Pg, unique(stats::quantile(d$Pg, seq(0, 1, .25))), include.lowest = TRUE)
      cat("  by Pg quartile: mean ratio\n"); print(round(tapply(d$ratio, bP, mean), 5))
      if (!ci$null) {
        bQ <- cut(d$PQg, c(-Inf, 0.2, 0.4, 0.6, 0.8, 0.95, Inf))
        cat("  by purity PQg band: mean ratio (n)\n"); print(round(tapply(d$ratio, bQ, mean), 5)); print(table(bQ))
        big_impure <- d$Pg >= stats::quantile(d$Pg, .75) & d$PQg <= 0.4
        cat(sprintf("  large & low-purity (Pg >= Q3, PQg <= 0.4): n = %d, mean ratio %.5f, min %.5f\n",
                    sum(big_impure), mean(d$ratio[big_impure]), min(d$ratio[big_impure])))
      }
    }
    trend(A, "realized rules"); trend(F, "analytic family")
    # controls
    if (!ci$null) {
      scQ <- fs_dgm_scale(dgm, regions = list(Q = inQ))$regions
      cat(sprintf("control: Q itself -- ratio %.6f (by construction 1)\n", sqrt(V_ref / scQ$V_eff)))
      hp <- F[F$PQg >= 0.95, ]
      cat(sprintf("control: family candidates with PQg >= 0.95: n = %d, ratio mean %.5f range [%.5f, %.5f]\n",
                  nrow(hp), mean(hp$ratio), min(hp$ratio), max(hp$ratio)))
      hpA <- A[A$PQg >= 0.95, ]
      cat(sprintf("control: realized rules with PQg >= 0.95: n = %d, ratio mean %.5f range [%.5f, %.5f]\n",
                  nrow(hpA), mean(hpA$ratio), min(hpA$ratio), max(hpA$ratio)))
    }
    cat("side column, Jensen-inflated direct sd (not the task's ratio): se_scaled / se_direct_jensen over realized rules\n")
    print(qs(A$ratio_jensen)); cat(sprintf("  cor with Pg %+.3f\n", cor(A$ratio_jensen, A$Pg)))
    out[[cell]] <- list(n = n, null = ci$null, ref = ref_lab, V_ref = V_ref, piQ = piQ, seQ1000 = seQ1000,
                        realized = A, family = F, scale_regions = reg)
  }
  out$built_at <- Sys.time(); out$pkg_version <- as.character(utils::packageVersion("forestsearch"))
  out$note <- "ratio = se_scaled / se_direct = sqrt(V_eff[ref] / V_eff[g]); n-free. ratio_jensen is a side column outside the task's definition."
  saveRDS(out, file.path(D, "residual_hypA_2026-08-30.rds"))
  cat("\nwritten:", file.path(D, "residual_hypA_2026-08-30.rds"), "\n")
}

# =============================================================================
# PART B -- the dedup analogue and the fs_oc_predict() re-runs, one cell
# =============================================================================
subset_family <- function(fe, keep) {
  f <- fe
  for (nm in c("lab", "Pg", "PQg", "beta_g", "se_g", "sens_g", "spec_g")) f[[nm]] <- fe[[nm]][keep]
  f$ovl <- fe$ovl[keep, keep, drop = FALSE]; f$memb <- fe$memb[, keep, drop = FALSE]
  f$M <- sum(keep); f$counts["M"] <- f$M
  class(f) <- class(fe); f
}
mc_se_EnH <- function(pred, n, Pg) {
  w <- pred$sel_c; sqrt(sum(w * (n * Pg - pred$EnH)^2) / (DRAWS * pred$det_rate))
}

run_B <- function(cell) {
  ci <- cell_info[[cell]]; n <- ci$n; fam <- ci$fam
  dgm <- dgm_for(cell)
  fe <- fs_oc_family_enumerate(dgm, fs_args, n = n, max_M = 5000L, verbose = TRUE)
  stopifnot(identical(fe$lab, fam$lab), isTRUE(all.equal(fe$Pg, fam$Pg)), isTRUE(all.equal(fe$se_g, fam$se_g)),
            isTRUE(all.equal(fe$beta_g, fam$beta_g)), isTRUE(all.equal(fe$PQg, fam$PQg)), fe$M == fam$M)
  cat(sprintf("[%s] corrected family re-enumerated and identical to the stored record: M = %d\n", cell, fe$M))
  M <- fe$M
  inQ <- dgm$df_super$flag_harm == 1
  cnt   <- colSums(fe$memb)                       # |g|
  cntQ  <- colSums(fe$memb & inQ)                 # |g & Q|
  K     <- lengths(strsplit(fe$lab, " & ", fixed = TRUE))
  idx   <- seq_len(M)

  # ---- v1: the literal population key, the package's own rounding ----------
  # columns 2:10 of found.hrs -> (K, n, E, d1, m1, m0, HR, L, U); population counterparts:
  key_df <- data.frame(K = K, n = n * fe$Pg, E = n * fe$Pg, d1 = n * fe$Pg / 2,
                       m1 = NA_real_, m0 = NA_real_, HR = fe$beta_g,
                       L = fe$beta_g - 1.96 * fe$se_g, U = fe$beta_g + 1.96 * fe$se_g)
  kr <- key_df
  for (j in seq_len(ncol(kr))) if (is.numeric(kr[[j]])) kr[[j]] <- round(kr[[j]] / TOL) * TOL
  key1 <- apply(kr, 1, function(x) paste(x, collapse = "_"))
  # the package orders (-HR, K) before !duplicated(); enumeration order breaks exact ties
  ord <- order(-fe$beta_g, K, idx)
  keep1 <- logical(M); keep1[ord] <- !duplicated(key1[ord])
  # equivalence check: the key is (K, |g|, |g&Q|) exactly (null: (K, |g|))
  key1_eq <- if (ci$null) paste(K, cnt) else paste(K, cnt, cntQ)
  stopifnot(all(tapply(key1_eq, key1, function(v) length(unique(v))) == 1L),
            all(tapply(key1, key1_eq, function(v) length(unique(v))) == 1L))
  g1 <- split(idx, key1); g1 <- g1[lengths(g1) > 1]
  # how alike are the members of a v1 group?  Jaccard between each dropped member and its representative
  jac <- function(i, j) fe$ovl[i, j] * 5000 / (cnt[i] + cnt[j] - fe$ovl[i, j] * 5000)
  v1_pairs <- do.call(rbind, lapply(g1, function(g) {
    r <- g[keep1[g]]; stopifnot(length(r) == 1L)
    data.frame(rep = r, drop = g[!keep1[g]], jaccard = vapply(g[!keep1[g]], function(j) jac(r, j), 1))
  }))
  cat(sprintf("[%s] v1 literal key: %d groups of size > 1, %d candidates removed -> M = %d\n",
              cell, length(g1), sum(!keep1), sum(keep1)))
  if (!is.null(v1_pairs)) {
    cat(sprintf("[%s] v1: Jaccard(dropped, representative): median %.3f, mean %.3f, min %.3f, max %.3f; share with Jaccard >= 0.99: %.3f\n",
                cell, median(v1_pairs$jaccard), mean(v1_pairs$jaccard), min(v1_pairs$jaccard), max(v1_pairs$jaccard),
                mean(v1_pairs$jaccard >= 0.99)))
    ex <- head(v1_pairs[order(v1_pairs$jaccard), ], 5)
    cat(sprintf("[%s] v1 example collapsed pairs (least alike): \n", cell))
    for (i in seq_len(nrow(ex))) cat(sprintf("     keep '%s' (n*Pg %.1f, beta %.3f)  drop '%s' (Jaccard %.3f)\n",
                                             fe$lab[ex$rep[i]], n * fe$Pg[ex$rep[i]], fe$beta_g[ex$rep[i]], fe$lab[ex$drop[i]], ex$jaccard[i]))
  }

  # ---- v2: near-twins within K (APPROXIMATE) --------------------------------
  d_max <- floor(5000 * log(2) / n)
  sd_mat <- outer(cnt, cnt, "+") - 2 * fe$ovl * 5000    # symmetric difference in super-population subjects
  sameK <- outer(K, K, "==")
  near <- (sd_mat <= d_max) & sameK; diag(near) <- FALSE
  # greedy in the package's order: keep a candidate unless a kept, earlier-ordered candidate is a near-twin
  keep2 <- logical(M); kept <- integer(0)
  for (i in ord) {
    if (length(kept) && any(near[i, kept])) next
    keep2[i] <- TRUE; kept <- c(kept, i)
  }
  n_pairs2 <- sum(near[upper.tri(near)])
  cat(sprintf("[%s] v2 near-twin (d_max = %d subjects of 5000, P(coincide in a sample of %d) >= 0.5): %d pairs, %d removed -> M = %d\n",
              cell, d_max, n, n_pairs2, sum(!keep2), sum(keep2)))
  if (sum(!keep2)) {
    dr2 <- which(!keep2)
    cat(sprintf("[%s] v2 examples: \n", cell))
    for (j in head(dr2, 5)) { i <- kept[which(near[j, kept])[1]]
      cat(sprintf("     drop '%s' (n*Pg %.1f) ~ keep '%s' (n*Pg %.1f), symdiff %g\n", fe$lab[j], n * fe$Pg[j], fe$lab[i], n * fe$Pg[i], sd_mat[i, j])) }
  }
  # pairs that the sample key would coincide for with probability >= 0.05 (informational)
  d_05 <- floor(5000 * (1 - 0.05^(1 / n)))
  cat(sprintf("[%s] informational: pairs (within K) with symdiff <= %d (P(coincide) >= 0.05): %d; minimum symdiff between any two distinct candidates: %g\n",
              cell, d_05, sum(((sd_mat <= d_05) & sameK)[upper.tri(sd_mat)]), min(sd_mat[upper.tri(sd_mat)])))

  # ---- the re-runs, both gates, same seed and draws ------------------------
  versions <- list(v1_literal = keep1, v2_neartwin = keep2)
  runs <- list()
  for (v in names(versions)) {
    keep <- versions[[v]]
    if (all(keep)) { cat(sprintf("[%s] %s removes nothing: the re-run is the stored run by construction (same seed, same M)\n", cell, v)); runs[[v]] <- NULL; next }
    fsub <- subset_family(fe, keep)
    for (g in c("resample", "split")) {
      gc(); t0 <- proc.time()[["elapsed"]]
      p <- fs_oc_predict(family = fsub, n = n, c1 = fs_args$effect.threshold, c2 = fs_args$consistency.threshold,
                         consistency_method = g, pconsistency = fs_args$pconsistency.threshold,
                         draws = DRAWS, seed = SEED)
      secs <- proc.time()[["elapsed"]] - t0
      p$family <- NULL
      runs[[v]][[g]] <- list(pred = p, secs = secs, keep = keep, M = fsub$M)
      st <- ci$runs[[g]]$pred
      cat(sprintf("[%s] %-12s %-8s M %d -> %d | det %.5f -> %.5f | EnH %.3f -> %.3f | Esens %.4f -> %.4f | Eppv %.4f -> %.4f | naive %.3f -> %.3f | %.0f s\n",
                  cell, v, g, M, fsub$M, st$det_rate, p$det_rate, st$EnH, p$EnH, st$Esens, p$Esens, st$Eppv, p$Eppv, st$Enaive_bias, p$Enaive_bias, secs))
    }
  }
  out <- list(cell = cell, n = n, null = ci$null, M = M, lab = fe$lab, Pg = fe$Pg, PQg = fe$PQg, beta_g = fe$beta_g, K = K,
              key1 = key1, keep1 = keep1, v1_groups = g1, v1_pairs = v1_pairs,
              keep2 = keep2, d_max = d_max, n_pairs2 = n_pairs2, min_symdiff = min(sd_mat[upper.tri(sd_mat)]),
              runs = runs, stored = ci$runs, seed = SEED, draws = DRAWS, tol = TOL)
  saveRDS(out, part_path(cell)); cat("[", cell, "] part written: ", part_path(cell), "\n", sep = "")
}

# =============================================================================
# Bsum -- merge, between-rule gap, comparison against measurement
# =============================================================================
run_Bsum <- function() {
  parts <- lapply(names(cell_info), function(cl) readRDS(part_path(cl))); names(parts) <- names(cell_info)
  summ <- list(); tabs <- list()
  for (cl in names(parts)) {
    P <- parts[[cl]]; ci <- cell_info[[cl]]; n <- ci$n; d <- S$cells[[cl]]$reps
    pop_realized <- mean(d$nPg_pop); meas_EnH <- mean(d$n_harm)
    size_sig <- tapply(d$nPg_pop, d$sig, mean); s <- names(size_sig)
    rw <- function(sel_c, lab) { A <- tapply(sel_c, vapply(lab, rule_sig, character(1), USE.NAMES = FALSE), sum)
                                 w <- A[s]; w[is.na(w)] <- 0; sum(w * size_sig) / sum(w) }
    cat(sprintf("\n################ %s (n = %d): M = %d; v1 removes %d, v2 removes %d\n", cl, n, P$M, sum(!P$keep1), sum(!P$keep2)))
    rows <- list()
    for (g in c("resample", "split")) {
      st <- P$stored[[g]]$pred
      rows[[length(rows) + 1]] <- data.frame(cell = cl, version = "corrected (stored)", gate = g, M = P$M,
        det_rate = st$det_rate, EnH = st$EnH, EnH_mcse = mc_se_EnH(st, n, P$Pg), Esens = st$Esens, Espec = st$Espec, Eppv = st$Eppv,
        Enaive_bias = st$Enaive_bias, between = pop_realized - st$EnH, reweight = rw(st$sel_c, P$lab), stringsAsFactors = FALSE)
      for (v in names(P$runs)) {
        r <- P$runs[[v]][[g]]; p <- r$pred; keep <- r$keep
        rows[[length(rows) + 1]] <- data.frame(cell = cl, version = v, gate = g, M = r$M,
          det_rate = p$det_rate, EnH = p$EnH, EnH_mcse = mc_se_EnH(p, n, P$Pg[keep]), Esens = p$Esens, Espec = p$Espec, Eppv = p$Eppv,
          Enaive_bias = p$Enaive_bias, between = pop_realized - p$EnH, reweight = rw(p$sel_c, P$lab[keep]), stringsAsFactors = FALSE)
      }
    }
    tb <- do.call(rbind, rows)
    tb$measured_EnH <- meas_EnH; tb$pop_realized <- pop_realized
    tb$gap_EnH_meas <- meas_EnH - tb$EnH
    tb$gap_sens_meas <- if (ci$null) NA else ci$meas$Esens - tb$Esens
    tb$gap_ppv_meas  <- if (ci$null) NA else ci$meas$Eppv - tb$Eppv
    tb$gap_naive_meas <- ci$meas$Enaive_bias - tb$Enaive_bias
    print(tb[, c("version", "gate", "M", "det_rate", "EnH", "EnH_mcse", "between", "reweight", "Esens", "Eppv", "Enaive_bias",
                 "gap_EnH_meas", "gap_sens_meas", "gap_ppv_meas", "gap_naive_meas")], digits = 5, row.names = FALSE)
    cat(sprintf("measured: EnH %.2f (SE %.2f), sens %s, ppv %s, naive bias %.2f (SE %.2f); pop. size of realized rules %.2f\n",
                ci$meas$EnH, ci$meas$EnH_se, if (ci$null) "NA" else sprintf("%.4f", ci$meas$Esens),
                if (ci$null) "NA" else sprintf("%.4f", ci$meas$Eppv), ci$meas$Enaive_bias, ci$meas$Enaive_bias_se, pop_realized))
    # where did the v1-removed candidates' mass sit in the stored run?
    for (g in c("resample", "split")) {
      st <- P$stored[[g]]$pred
      cat(sprintf("  stored %-8s sel_c mass on the v1-removed candidates: %.4f (on v2-removed: %.4f); sel_c mass on v1 representatives of groups: %.4f\n",
                  g, sum(st$sel_c[!P$keep1]), sum(st$sel_c[!P$keep2]),
                  sum(st$sel_c[unique(P$v1_pairs$rep)])))
    }
    tabs[[cl]] <- tb
  }
  all_tab <- do.call(rbind, tabs); rownames(all_tab) <- NULL
  cat("\n==== between-rule size gap, corrected (stored) vs deduped families (resample gate):\n")
  print(all_tab[all_tab$gate == "resample", c("cell", "version", "M", "EnH", "between", "reweight", "gap_EnH_meas", "gap_naive_meas")], digits = 4, row.names = FALSE)
  out <- list(table = all_tab, parts = lapply(parts, function(P) P[setdiff(names(P), c("stored"))]),
              built_at = Sys.time(), pkg_version = as.character(utils::packageVersion("forestsearch")),
              key_source = "remove_near_duplicate_subgroups(): columns 2:10 of found.hrs (K, n, E, d1, m1, m0, HR, L(HR), U(HR)) rounded to 0.001, first kept; found.hrs ordered (-HR, K)")
  saveRDS(out, file.path(D, "residual_hypB_2026-08-30.rds"))
  cat("\nwritten:", file.path(D, "residual_hypB_2026-08-30.rds"), "\n")
}

switch(part, A = run_A(), B = run_B(args[2]), Bsum = run_Bsum())
