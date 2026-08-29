# =============================================================================
# residual_quantiles_2026-08-30.R -- population-quantile cuts against
# replicate-quantile cuts (task: dev/tasks/cc_task_residual_quantiles_2026-08-30.md)
# -----------------------------------------------------------------------------
# Read-only with respect to the package.  Reads the corrected 13-variable
# families and their stored fs_oc_predict() results
# (oc_wrapper_grid_corrected_2026-08-30.rds) and the realized rules
# (sgdef_selection_2026-08-29.rds).
#
# Parts (repository root):
#   Rscript dev/glm-continuous-sims/residual_quantiles_2026-08-30.R stage1 <cell>   # S2, one cell (alt500 | alt700 | null500)
#   Rscript dev/glm-continuous-sims/residual_quantiles_2026-08-30.R gate            # S2 tables + the gate, all cells
#   Rscript dev/glm-continuous-sims/residual_quantiles_2026-08-30.R stage2 <cell>   # S3, only if the gate says continue
#   Rscript dev/glm-continuous-sims/residual_quantiles_2026-08-30.R summary         # S3 comparison
#
# Stage 1.  R replicates are drawn from the cell's DGM at its n with
# simulate_from_glm_dgm(dgm, n, seed = 20260825 + r), r = 1..R -- the corrected
# run's seed family.  For each replicate the search's cut matrix is built with
# get_FSdata() under the drivers' arguments exactly as fs_oc_family_enumerate()
# builds the population one (same call, the replicate frame in place of
# df_super), and each corrected-family candidate is located in the replicate
# family by its CLAUSE SPECIFICATION: (variable, operator, rank of the
# threshold within that variable+operator's cut list) for continuous clauses,
# (variable, level) for binaries -- never by threshold value.  A candidate
# whose clause has no counterpart (the replicate's collapse produced fewer cuts
# for that variable) is recorded as "disappears" for that replicate; replicate
# columns with no population counterpart are recorded as "appears".
# Per replicate and candidate: the replicate sample prevalence minus the
# population Pg.  Running sums of prevalence, P(g & Q), and the M x M overlap
# are kept for Stage 2 so the draws are not repeated.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE); part <- args[1]
stopifnot(part %in% c("stage1", "gate", "stage2", "summary"))
D <- "dev/glm-continuous-sims"
C <- readRDS(file.path(D, "oc_wrapper_grid_corrected_2026-08-30.rds"))
S <- readRDS(file.path(D, "sgdef_selection_2026-08-29.rds"))
fs_args <- C$fs_args; stopifnot(length(fs_args$confounders.name) == 13L)
SEED_FAMILY <- 20260825L; R_REPS <- 200L
DRAWS <- C$draws; SEED_OC <- C$seed; stopifnot(SEED_OC == 20260825L, DRAWS == 2e5)
stage1_path <- function(cell) file.path(D, sprintf("residual_quantiles_stage1_%s_2026-08-30.rds", cell))
stage2_path <- function(cell, variant) file.path(D, sprintf("residual_quantiles_stage2_%s_%s_2026-08-30.rds", cell, variant))

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
    stopifnot(sum(dgm$df_super$flag_harm) == 0L, abs(dgm$hazard_ratios$no_harm_subgroup - truth$effect_Qc) < 1e-9)
  } else {
    dgm <- build_dgm("alt", truth$beta_inter)
    stopifnot(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q) < 1e-9, abs(dgm$subgroup_info$proportion - truth$prevalence_Q) < 1e-9)
  }
  dgm
}

# ---- the search's cut matrix on a frame, exactly as fs_oc_family_enumerate() builds it ----
fs_formals <- formals(forestsearch)
.arg <- function(nm) if (nm %in% names(fs_args)) fs_args[[nm]] else eval(fs_formals[[nm]])
cut_matrix <- function(frame) {
  df_cut <- frame; df_cut[[".fs_oc_y"]] <- 0; df_cut[[".fs_oc_event"]] <- 1
  FSdata <- get_FSdata(
    df.analysis = df_cut, use_lasso = FALSE, use_grf = FALSE, grf_cuts = NULL,
    confounders.name = .arg("confounders.name"), cont.cutoff = .arg("cont.cutoff"),
    conf_force = .arg("conf_force"), conf.cont_medians = .arg("conf.cont_medians"),
    conf.cont_medians_force = .arg("conf.cont_medians_force"), conf.cont_jcuts = .arg("conf.cont_jcuts"),
    dina_cuts = NULL, collapse_cuts = .arg("collapse_cuts"), collapse_cuts_args = .arg("collapse_cuts_args"),
    defaultcut_names = .arg("defaultcut_names"), cut_type = .arg("cut_type"), exclude_cuts = .arg("exclude_cuts"),
    outcome.name = ".fs_oc_y", event.name = ".fs_oc_event", details = FALSE, outcome_type = "continuous")
  Zdf <- dummy(FSdata$df[, FSdata$confs_names, drop = FALSE])
  Z <- as.matrix(Zdf); storage.mode(Z) <- "numeric"; colnames(Z) <- names(Zdf)
  lab <- forestsearch:::.fs_oc_column_labels(colnames(Z), FSdata$confs_names, FSdata$confs)
  list(Z = Z, lab = lab, cuts = FSdata$confs)
}
# clause label -> specification key: continuous "var op value" -> var|op|rank (rank of value within the
# var+op cut list, ascending); "var > v" is the negation of cut "var <= v", so its rank is that of v in the "<=" list
PAT <- "^\\s*([A-Za-z0-9_.]+)\\s*(<=|>=|<|>|!=|==)\\s*(-?[0-9.]+([eE][-+]?[0-9]+)?)\\s*$"
parse_clause <- function(lab) {
  m <- regmatches(lab, regexec(PAT, lab))[[1]]
  if (length(m)) return(list(var = m[2], op = m[3], val = as.numeric(m[4]), binary = FALSE))
  if (grepl("!= 1$", lab)) return(list(var = trimws(sub("!= 1$", "", lab)), op = "=0", val = NA_real_, binary = TRUE))
  list(var = trimws(lab), op = "=1", val = NA_real_, binary = TRUE)
}
NEG <- c("<=" = ">", ">" = "<=", "<" = ">=", ">=" = "<")
spec_keys <- function(labs) {
  cl <- lapply(labs, parse_clause)
  # threshold lists per (var, base operator): base op is the cut's own operator ("<=" for a ">" negation)
  base_op <- vapply(cl, function(z) if (z$binary) z$op else if (z$op %in% c(">", ">=")) NEG[[z$op]] else z$op, character(1))
  var <- vapply(cl, `[[`, character(1), "var"); val <- vapply(cl, `[[`, numeric(1), "val")
  key <- character(length(labs))
  for (i in seq_along(labs)) {
    if (cl[[i]]$binary) { key[i] <- paste(var[i], cl[[i]]$op, sep = "|"); next }
    g <- var == var[i] & base_op == base_op[i]
    vals <- sort(unique(val[g]))
    key[i] <- paste(var[i], cl[[i]]$op, match(val[i], vals), length(vals), sep = "|")   # var|op|rank|count
  }
  key
}
cand_clauses <- function(lab) strsplit(lab, " & ", fixed = TRUE)

run_stage1 <- function(cell) {
  ci <- cell_info[[cell]]; n <- ci$n; fam <- ci$fam; M <- fam$M
  dgm <- dgm_for(cell)
  # population: re-enumerate to get memberships and assert identity to the stored record
  fe <- fs_oc_family_enumerate(dgm, fs_args, n = n, max_M = 5000L, verbose = FALSE)
  stopifnot(identical(fe$lab, fam$lab), isTRUE(all.equal(fe$Pg, fam$Pg)), fe$M == M)
  pop_cm <- cut_matrix(dgm$df_super)
  stopifnot(identical(pop_cm$cuts, fam$cuts))
  pop_key <- spec_keys(pop_cm$lab); names(pop_key) <- pop_cm$lab
  cand_cl <- cand_clauses(fe$lab)
  cand_keys <- lapply(cand_cl, function(cls) unname(pop_key[cls]))
  stopifnot(!anyNA(unlist(cand_keys)))
  K <- lengths(cand_cl)
  inQ_fun <- function(frame) frame$flag_harm == 1
  cat(sprintf("[%s] n = %d, M = %d, population cut columns %d (%d cuts); R = %d replicates, seeds %d + 1..%d\n",
              cell, n, M, ncol(pop_cm$Z), length(pop_cm$cuts), R_REPS, SEED_FAMILY, R_REPS))
  # per-variable population cut counts
  pk <- do.call(rbind, lapply(strsplit(unique(pop_key), "|", fixed = TRUE), function(z) data.frame(var = z[1], op = z[2], count = if (length(z) == 4) as.integer(z[4]) else NA_integer_)))
  pop_counts <- unique(pk[!is.na(pk$count) & pk$op %in% c("<=", "<"), c("var", "count")])
  cat("population continuous cut counts per variable: "); cat(sprintf("%s=%d", pop_counts$var, pop_counts$count), "\n")

  dP <- matrix(NA_real_, R_REPS, M)                     # replicate prevalence - Pg
  sumP <- numeric(M); sumPQ <- numeric(M); sumPQc <- numeric(M); nvalid <- integer(M)
  sum_ovl <- matrix(0, M, M); n_ovl <- 0L; cnt_ovl <- matrix(0L, M, M)
  PQmat <- matrix(NA_real_, R_REPS, M); PQcmat <- matrix(NA_real_, R_REPS, M)   # per replicate P(g & Q), P(g & Qc)
  appear <- list(); disappear <- integer(M); rep_counts <- list()
  for (r in seq_len(R_REPS)) {
    df <- simulate_from_glm_dgm(dgm, n = n, seed = SEED_FAMILY + r)
    cm <- cut_matrix(df)
    rk <- spec_keys(cm$lab)
    # replicate cut counts per variable (continuous "<=" lists)
    rc <- do.call(rbind, lapply(strsplit(unique(rk), "|", fixed = TRUE), function(z) if (length(z) == 4 && z[2] == "<=") data.frame(var = z[1], count = as.integer(z[4])) else NULL))
    rep_counts[[r]] <- rc
    col_of <- match(unlist(lapply(cand_keys, function(k) k)), rk)   # flattened
    # candidate membership: AND over its clauses' replicate columns
    memb <- matrix(FALSE, n, M); ok <- logical(M)
    pos <- 0L
    for (m in seq_len(M)) {
      kk <- cand_keys[[m]]; idx <- match(kk, rk)
      if (anyNA(idx)) { ok[m] <- FALSE; disappear[m] <- disappear[m] + 1L; next }
      v <- cm$Z[, idx[1]] == 1
      if (length(idx) > 1) for (j in idx[-1]) v <- v & (cm$Z[, j] == 1)
      memb[, m] <- v; ok[m] <- TRUE
    }
    appear[[r]] <- setdiff(rk, unlist(cand_keys))
    prev <- colMeans(memb); prev[!ok] <- NA
    dP[r, ] <- prev - fam$Pg
    inQ <- inQ_fun(df)
    sumP[ok] <- sumP[ok] + prev[ok]; nvalid[ok] <- nvalid[ok] + 1L
    PQmat[r, ok] <- colMeans(memb[, ok, drop = FALSE] & inQ); PQcmat[r, ok] <- colMeans(memb[, ok, drop = FALSE] & !inQ)
    sumPQ[ok] <- sumPQ[ok] + PQmat[r, ok]; sumPQc[ok] <- sumPQc[ok] + PQcmat[r, ok]
    # overlap: pairwise-valid running sum (both candidates located) and its count; unlocated columns are all FALSE
    okm <- outer(ok, ok)
    sum_ovl <- sum_ovl + crossprod(memb * 1) / n; cnt_ovl <- cnt_ovl + okm
    if (all(ok)) n_ovl <- n_ovl + 1L
    if (r %% 25 == 0) cat(sprintf("[%s] replicate %d: cut columns %d, candidates located %d / %d, new replicate keys %d\n", cell, r, ncol(cm$Z), sum(ok), M, length(appear[[r]])))
  }
  # ---- summaries ----
  selc <- ci$runs$resample$pred$sel_c; selc_s <- ci$runs$split$pred$sel_c
  EdP <- colMeans(dP, na.rm = TRUE)                       # per-candidate mean shift
  located_all <- colSums(!is.na(dP)) == R_REPS
  cat(sprintf("[%s] label correspondence: candidates located in every replicate %d / %d; located in none %d; mean located per replicate %.1f\n",
              cell, sum(located_all), M, sum(colSums(!is.na(dP)) == 0), mean(rowSums(!is.na(dP)))))
  disc <- which(!located_all)
  if (length(disc)) { cat(sprintf("[%s] candidates that disappear in some replicates (n = %d); replicates missing, top 10:\n", cell, length(disc)))
    o <- order(-disappear[disc])[1:min(10, length(disc))]; for (i in disc[o]) cat(sprintf("     %-40s missing in %d / %d\n", fe$lab[i], disappear[i], R_REPS)) }
  ap <- table(unlist(appear)); if (length(ap)) { cat(sprintf("[%s] replicate keys with no population counterpart (key: var|op|rank|count, replicates):\n", cell)); print(head(sort(ap, decreasing = TRUE), 12)) }
  rcc <- do.call(rbind, lapply(seq_along(rep_counts), function(r) { z <- rep_counts[[r]]; if (!is.null(z)) { z$r <- r; z } }))
  cat(sprintf("[%s] replicate continuous cut counts per variable (mode, range) vs population:\n", cell))
  for (v in pop_counts$var) { x <- rcc$count[rcc$var == v]; cat(sprintf("     %-8s pop %d | rep mode %d, range %d..%d, share equal to pop %.3f\n", v, pop_counts$count[pop_counts$var == v], as.integer(names(which.max(table(x)))), min(x), max(x), mean(x == pop_counts$count[pop_counts$var == v]))) }
  cat(sprintf("[%s] delta prevalence (replicate - population), per candidate mean over replicates: mean %+.5f, median %+.5f, sd %.5f, range [%+.4f, %+.4f]\n",
              cell, mean(EdP, na.rm = TRUE), median(EdP, na.rm = TRUE), sd(EdP, na.rm = TRUE), min(EdP, na.rm = TRUE), max(EdP, na.rm = TRUE)))
  cat(sprintf("[%s] in subjects (x n): unweighted mean %+.3f\n", cell, n * mean(EdP, na.rm = TRUE)))
  byK <- tapply(EdP, K, mean, na.rm = TRUE); cat(sprintf("[%s] by K: ", cell)); cat(sprintf("K=%s %+.5f (n=%d)", names(byK), byK, as.integer(table(K))), sep = " | "); cat("\n")
  qP <- cut(fam$Pg, unique(quantile(fam$Pg, seq(0, 1, .25))), include.lowest = TRUE)
  byP <- tapply(EdP, qP, mean, na.rm = TRUE); cat(sprintf("[%s] by Pg quartile: ", cell)); cat(sprintf("%s %+.5f", names(byP), byP), sep = " | "); cat("\n")
  cat(sprintf("[%s] cor(E[dP], Pg) = %+.3f; lm slope of E[dP] on Pg = %+.4f (SE %.4f)\n", cell, cor(EdP, fam$Pg, use = "complete"), coef(summary(lm(EdP ~ fam$Pg)))[2, 1], coef(summary(lm(EdP ~ fam$Pg)))[2, 2]))
  w <- function(sc) { ok <- !is.na(EdP); sum(sc[ok] * EdP[ok]) / sum(sc[ok]) }
  shift_r <- w(selc); shift_s <- w(selc_s)
  gaps <- c(alt500 = 2.11, alt700 = 0.61, null500 = 1.65)
  cat(sprintf("[%s] sel_c-weighted E[dP]: resample %+.5f -> %+.3f subjects; split %+.5f -> %+.3f subjects | between-rule gap %+.2f | share %.1f%% / %.1f%%\n",
              cell, shift_r, n * shift_r, shift_s, n * shift_s, gaps[[cell]], 100 * n * shift_r / gaps[[cell]], 100 * n * shift_s / gaps[[cell]]))
  # the selection-weighted shift restricted to the candidates the search actually picked (sgdef signatures), informational
  out <- list(cell = cell, n = n, M = M, R = R_REPS, seed_family = SEED_FAMILY, lab = fe$lab, K = K, Pg = fam$Pg, PQg = fam$PQg,
              sel_c_resample = selc, sel_c_split = selc_s, dP = dP, EdP = EdP, located_all = located_all, disappear = disappear,
              appear = ap, rep_counts = rcc, pop_counts = pop_counts,
              sums = list(P = sumP, PQ = sumPQ, PQc = sumPQc, nvalid = nvalid, ovl = sum_ovl, n_ovl = n_ovl, cnt_ovl = cnt_ovl),
              PQmat = PQmat, PQcmat = PQcmat,
              shift = c(resample = n * shift_r, split = n * shift_s), gap = gaps[[cell]],
              byK = byK, byP = byP, when = Sys.time())
  saveRDS(out, stage1_path(cell)); cat("written:", stage1_path(cell), "\n")
}

run_gate <- function() {
  tab <- do.call(rbind, lapply(names(cell_info), function(cl) {
    s <- readRDS(stage1_path(cl))
    data.frame(cell = cl, n = s$n, M = s$M, located_all = sum(s$located_all),
               shift_resample = s$shift[["resample"]], shift_split = s$shift[["split"]], gap = s$gap,
               share_resample = s$shift[["resample"]] / s$gap, share_split = s$shift[["split"]] / s$gap,
               cor_dP_Pg = cor(s$EdP, s$Pg, use = "complete"),
               byK1 = s$byK[["1"]], byK2 = s$byK[["2"]], byP_Q1 = s$byP[[1]], byP_Q4 = s$byP[[length(s$byP)]])
  }))
  print(tab, digits = 4, row.names = FALSE)
  small <- abs(tab$share_resample) < 0.20 & abs(tab$share_split) < 0.20
  cat(sprintf("\nGATE: sel_c-weighted shift under 20%% of the gap in every cell: %s\n", if (all(small)) "YES -> skip S3" else "NO -> run S3"))
  saveRDS(tab, file.path(D, "residual_quantiles_gate_2026-08-30.rds"))
}

run_stage2 <- function(cell, variant) {
  # variant "located": every field is the mean over the replicates in which that candidate (pair) was located --
  #                    Pg, P(g&Q), P(g&Qc) over nvalid[m] replicates, ovl[i, j] over cnt_ovl[i, j] replicates.
  # variant "full":    every field is the mean over the fully-corresponding replicates only (all M located) --
  #                    internally consistent, fewer replicates.
  stopifnot(variant %in% c("located", "full"))
  ci <- cell_info[[cell]]; n <- ci$n; fam <- ci$fam; M <- fam$M
  s <- readRDS(stage1_path(cell))
  dgm <- dgm_for(cell)
  fe <- fs_oc_family_enumerate(dgm, fs_args, n = n, max_M = 5000L, verbose = FALSE)
  stopifnot(identical(fe$lab, fam$lab))
  Pmat <- sweep(s$dP, 2, fam$Pg, "+")                    # replicate prevalence (NA where not located)
  if (variant == "located") {
    Pg <- colMeans(Pmat, na.rm = TRUE); PgQ <- colMeans(s$PQmat, na.rm = TRUE); PgQc <- colMeans(s$PQcmat, na.rm = TRUE)
    ovl <- s$sums$ovl / pmax(s$sums$cnt_ovl, 1L); ovl[s$sums$cnt_ovl == 0L] <- NA
    stopifnot(!anyNA(ovl), !anyNA(Pg))
    R_used <- range(s$sums$nvalid)
  } else {
    full <- rowSums(is.na(s$dP)) == 0L; stopifnot(sum(full) >= 10L)
    Pg <- colMeans(Pmat[full, , drop = FALSE]); PgQ <- colMeans(s$PQmat[full, , drop = FALSE]); PgQc <- colMeans(s$PQcmat[full, , drop = FALSE])
    # the overlap sum over the full replicates only is not stored separately; recover it by re-drawing those replicates
    sum_ovl <- matrix(0, M, M)
    pop_key <- spec_keys(cut_matrix(dgm$df_super)$lab); names(pop_key) <- cut_matrix(dgm$df_super)$lab
    cand_keys <- lapply(cand_clauses(fe$lab), function(cls) unname(pop_key[cls]))
    for (r in which(full)) {
      df <- simulate_from_glm_dgm(dgm, n = n, seed = SEED_FAMILY + r); cm <- cut_matrix(df); rk <- spec_keys(cm$lab)
      memb <- matrix(FALSE, n, M)
      for (m in seq_len(M)) { idx <- match(cand_keys[[m]], rk); v <- cm$Z[, idx[1]] == 1; if (length(idx) > 1) for (j in idx[-1]) v <- v & (cm$Z[, j] == 1); memb[, m] <- v }
      stopifnot(isTRUE(all.equal(colMeans(memb), Pmat[r, ])))
      sum_ovl <- sum_ovl + crossprod(memb * 1) / n
    }
    ovl <- sum_ovl / sum(full); R_used <- c(sum(full), sum(full))
  }
  cat(sprintf("[%s] variant %s: replicates per field %d..%d; Pg mean shift %+.5f; min eigenvalue of Rho %.4f\n", cell, variant, R_used[1], R_used[2],
              mean(Pg - fam$Pg), min(eigen((ovl / sqrt(outer(Pg, Pg)) + t(ovl / sqrt(outer(Pg, Pg)))) / 2, symmetric = TRUE, only.values = TRUE)$values)))
  reg <- fe$scale$regions
  if (!ci$null) {
    piQ <- reg$P_g[reg$region == "Q"]; V_ref <- reg$V_eff[reg$region == "Q"]
    tauQc <- abs(reg$m_tau[reg$region == "Qc"]); bint <- abs(reg$m_tau[reg$region == "Q"]) - tauQc
    PQ <- fe$PQ; PQg <- PgQ / Pg; beta_g <- tauQc + bint * PQg
    sens_g <- PgQ / PQ; spec_g <- 1 - PgQc / (1 - PQ)
  } else {
    piQ <- 1; V_ref <- reg$V_eff[reg$region == "S"]; tauQc <- abs(reg$m_tau[reg$region == "S"])
    PQ <- 0; PQg <- rep(0, M); beta_g <- rep(tauQc, M); sens_g <- rep(NA_real_, M); spec_g <- 1 - Pg
  }
  seQ1000 <- sqrt(V_ref / (1000 * piQ))
  se_g <- seQ1000 * sqrt(1000 / n) * sqrt(piQ / Pg)      # the corrected family's formula at the replicate-mean Pg
  favg <- fe; favg$Pg <- unname(Pg); favg$PQg <- unname(PQg); favg$beta_g <- unname(beta_g); favg$se_g <- unname(se_g)
  favg$sens_g <- unname(sens_g); favg$spec_g <- unname(spec_g); favg$ovl <- ovl; favg$memb <- NULL
  forestsearch:::.fs_oc_check_family(favg)
  runs <- list()
  for (g in c("resample", "split")) {
    gc(); t0 <- proc.time()[["elapsed"]]
    p <- fs_oc_predict(family = favg, n = n, c1 = fs_args$effect.threshold, c2 = fs_args$consistency.threshold,
                       consistency_method = g, pconsistency = fs_args$pconsistency.threshold, draws = DRAWS, seed = SEED_OC)
    p$family <- NULL; runs[[g]] <- list(pred = p, secs = proc.time()[["elapsed"]] - t0)
    st <- ci$runs[[g]]$pred
    cat(sprintf("[%s] %-8s det %.5f -> %.5f | EnH %.3f -> %.3f | Esens %s -> %s | Eppv %.4f -> %.4f | naive %.3f -> %.3f | %.0f s\n",
                cell, g, st$det_rate, p$det_rate, st$EnH, p$EnH, format(st$Esens, digits = 4), format(p$Esens, digits = 4), st$Eppv, p$Eppv, st$Enaive_bias, p$Enaive_bias, runs[[g]]$secs))
  }
  # inverted c1 at the corrected targets, resample gate only (the visible "anything else moved" check)
  targets <- if (ci$null) C$null$invert$targets else C$alt$invert$targets
  iv <- fs_oc_invert(family = favg, n = n, target = targets, solve_for = "c1", c2 = fs_args$consistency.threshold,
                     consistency_method = "resample", pconsistency = fs_args$pconsistency.threshold, draws = DRAWS, seed = SEED_OC)
  inv <- attr(iv, "table")
  saveRDS(list(cell = cell, n = n, variant = variant, R_used = R_used, family = favg[c("lab", "Pg", "PQg", "beta_g", "se_g", "sens_g", "spec_g", "M", "PQ")],
               runs = runs, invert = inv, targets = targets, when = Sys.time()), stage2_path(cell, variant))
  cat("written:", stage2_path(cell, variant), "\n")
}

run_summary <- function() {
  for (cl in names(cell_info)) for (variant in c("located", "full")) {
    if (!file.exists(stage2_path(cl, variant))) next
    ci <- cell_info[[cl]]; s2 <- readRDS(stage2_path(cl, variant)); d <- S$cells[[cl]]$reps
    cat(sprintf("\n======== variant %s (replicates per field %d..%d)", variant, s2$R_used[1], s2$R_used[2]))
    pop_realized <- mean(d$nPg_pop)
    cat(sprintf("\n######## %s (n = %d): pop. size of realized rules %.2f; measured EnH %.2f (SE %.2f)\n", cl, ci$n, pop_realized, ci$meas$EnH, ci$meas$EnH_se))
    for (g in c("resample", "split")) {
      st <- ci$runs[[g]]$pred; p <- s2$runs[[g]]$pred
      cat(sprintf("  %-8s between-rule %+.3f -> %+.3f | EnH gap to measured %+.3f -> %+.3f | sens %s -> %s (meas %s) | ppv %s -> %s (meas %s) | naive %.2f -> %.2f (meas %.2f) | det %.5f -> %.5f\n",
                  g, pop_realized - st$EnH, pop_realized - p$EnH, ci$meas$EnH - st$EnH, ci$meas$EnH - p$EnH,
                  format(st$Esens, digits = 4), format(p$Esens, digits = 4), format(ci$meas$Esens, digits = 4),
                  format(st$Eppv, digits = 4), format(p$Eppv, digits = 4), format(ci$meas$Eppv, digits = 4),
                  st$Enaive_bias, p$Enaive_bias, ci$meas$Enaive_bias, st$det_rate, p$det_rate))
    }
    old_inv <- if (ci$null) C$null$invert$table else C$alt$invert$table[C$alt$invert$table$n == ci$n, ]
    old_inv <- old_inv[old_inv$consistency_method == "resample", ]
    cat("  inverted c1 (resample): corrected -> replicate-averaged\n")
    print(data.frame(target = s2$targets, c1_corrected = old_inv$value, c1_avg = s2$invert$value), digits = 5, row.names = FALSE)
  }
}

switch(part, stage1 = run_stage1(args[2]), gate = run_gate(), stage2 = run_stage2(args[2], args[3]), summary = run_summary())
