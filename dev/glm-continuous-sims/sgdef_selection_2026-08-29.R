# =============================================================================
# sgdef_selection_2026-08-29.R -- the realized selection distribution
# -----------------------------------------------------------------------------
# Read-only diagnostic.  For the three tracked payloads (alt n = 500, alt
# n = 700, null n = 500):
#   §2 tabulate results$sg_def over detected replicates -- verbatim strings and
#      covariate signatures (thresholds stripped), K, per-covariate usage;
#   §3 evaluate every distinct verbatim rule on the DGM's df_super with the
#      package's own resolver (.fs_resolve_membership, the eval-frame path of
#      betaHhat_truth.R) and compute its population quantities;
#   §4 within-signature realized-vs-population comparison (the discriminating
#      test) and the naive bias by signature;
#   §5 measured vs analytic selection distribution over signatures, analytic
#      p_sel read from the stored oc_wrapper_grid / oc_wrapper_null .rds.
# No fs_oc_* re-runs.  Writes dev/glm-continuous-sims/sgdef_selection_2026-08-29.rds
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "quarto/simulations/actg175/continuous/mr_md_harm"
PAY <- c(
  alt500  = file.path(D, "fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000",   "fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds"),
  alt700  = file.path(D, "fs_maxeffCons_mr_md40_knoise0_n700_s1000_d5000",   "fs_maxeffCons_mr_md40_knoise0_n700_res_1_1000.rds"),
  null500 = file.path(D, "fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000", "fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds"))
OUT <- "dev/glm-continuous-sims/sgdef_selection_2026-08-29.rds"
G  <- readRDS("dev/glm-continuous-sims/oc_wrapper_grid_2026-08-29.rds")
N0 <- readRDS("dev/glm-continuous-sims/oc_wrapper_null_2026-08-29.rds")

read_pay <- function(p) { z <- readRDS(p); nm <- names(z); nm[is.na(nm)] <- ""; setNames(z, nm) }
pl <- lapply(PAY, read_pay)

# ---- the DGMs, rebuilt and gated ----------------------------------------------
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
build <- function(model, k_inter) generate_glm_dgm(
  data = actg_df, factor_vars = paste0("z", 1:12), outcome_var = "cd4_change",
  treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = model, k_treat = 1, k_inter = k_inter, n_super = 5000L,
  seed = 8316951L, verbose = FALSE)
tr_alt <- pl$alt500[["truth"]]; tr_null <- pl$null500[["truth"]]
dgm_alt  <- build("alt", tr_alt$beta_inter)
dgm_null <- build("null", 0)
stopifnot(abs(dgm_alt$hazard_ratios$harm_subgroup - tr_alt$effect_Q) < 1e-9,
          abs(dgm_alt$subgroup_info$proportion - tr_alt$prevalence_Q) < 1e-9,
          abs(dgm_alt$hazard_ratios$harm_subgroup - pl$alt700[["truth"]]$effect_Q) < 1e-9,
          sum(dgm_null$df_super$flag_harm) == 0L, is.na(dgm_null$hazard_ratios$harm_subgroup),
          abs(dgm_null$hazard_ratios$no_harm_subgroup - tr_null$effect_Qc) < 1e-9)
frame_alt  <- fs_build_eval_frame(dgm_alt,  outcome_type = "continuous")
frame_null <- fs_build_eval_frame(dgm_null, outcome_type = "continuous")
sc_alt <- fs_dgm_scale(dgm_alt); ra <- sc_alt$regions
tauQc <- abs(ra$m_tau[ra$region == "Qc"]); bint <- abs(ra$m_tau[ra$region == "Q"]) - tauQc
sc_null <- fs_dgm_scale(dgm_null, regions = list(S = rep(TRUE, nrow(frame_null))))
beta_null <- abs(sc_null$regions$m_tau)
cat(sprintf("DGMs gated. alt: tauQc %.4f bint %.4f; null common effect %.4f\n", tauQc, bint, beta_null))

# ---- signature parsing (thresholds stripped) -----------------------------------
# realized clause forms: "{var <= v}", "!{var <= v}", "{var}", "!{var}"
# analytic label forms:  "var <= v", "var > v", "var", "var != 1"
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
  paste0(cl, if (neg) "=0" else "=1")            # bare binary
}
rule_sig <- function(rule) {
  cls <- strsplit(rule, " & ", fixed = TRUE)[[1]]
  paste(sort(vapply(cls, clause_sig, character(1))), collapse = " & ")
}
rule_vars <- function(rule) unique(gsub("[<=>]+.*$|=[01]$", "", strsplit(rule_sig(rule), " & ", fixed = TRUE)[[1]]))

# ---- per cell ------------------------------------------------------------------
cells <- list()
for (cell in names(pl)) {
  p <- pl[[cell]]; res <- p[["results"]]; oc <- p[["oc"]]; meta <- p[["meta"]]
  n_trial <- meta$n_sample; is_null <- grepl("null", cell)
  frame <- if (is_null) frame_null else frame_alt
  inQ <- frame$flag_harm == 1; PQ <- mean(inQ)
  d <- res[res$detected %in% TRUE, ]
  orient <- oc$targets$orient
  d$beta_or <- orient * d$betaHhat_H                # oriented as the OC summary orients
  d$naive_or <- d$nv_H_est                          # already oriented (avg matches oc$estimation naive)
  d$naive_bias <- d$naive_or - d$beta_or
  d$sig <- vapply(d$sg_def, rule_sig, character(1), USE.NAMES = FALSE)
  d$K   <- lengths(strsplit(d$sg_def, " & ", fixed = TRUE))
  cat(sprintf("\n############ %s: n = %d, detected %d of %d\n", cell, n_trial, nrow(d), nrow(res)))
  # §2(a) verbatim
  tv <- sort(table(d$sg_def), decreasing = TRUE)
  cat(sprintf("distinct verbatim sg_def: %d; top 15:\n", length(tv))); print(head(tv, 15))
  # §2(b) signatures
  ts <- sort(table(d$sig), decreasing = TRUE)
  cat(sprintf("distinct signatures: %d; top 20:\n", length(ts))); print(head(ts, 20))
  cat("K distribution:\n"); print(table(d$K))
  vars_all <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80", "hemo", "homo", "drugs", "race", "gender", "symptom", "str2")
  use <- sapply(vars_all, function(v) mean(vapply(d$sg_def, function(r) v %in% rule_vars(r), logical(1))))
  cat("fraction of winners using each covariate:\n"); print(round(use, 3))
  cat(sprintf("winners using age or preanti (the document's two axes): %.3f; using neither: %.3f\n",
              mean(sapply(d$sg_def, function(r) any(c("age", "preanti") %in% rule_vars(r)))),
              mean(sapply(d$sg_def, function(r) !any(c("age", "preanti") %in% rule_vars(r))))))
  # §3 evaluate distinct verbatim rules on the population
  u <- unique(d$sg_def)
  ev <- do.call(rbind, lapply(u, function(r) {
    m <- forestsearch:::.fs_resolve_membership(frame, r)
    if (!identical(m$status, "ok") && !identical(m$status, "empty")) {
      return(data.frame(sg_def = r, status = m$status, missing = paste(m$missing, collapse = ","),
                        Pg_pop = NA_real_, nPg_pop = NA_real_, PQg_pop = NA_real_, sens_pop = NA_real_,
                        spec_pop = NA_real_, beta_pop = NA_real_, stringsAsFactors = FALSE))
    }
    g <- m$in_region; Pg <- mean(g); PgQ <- mean(g & inQ); PgQc <- mean(g & !inQ)
    PQg <- if (Pg > 0) PgQ / Pg else NA_real_
    data.frame(sg_def = r, status = m$status, missing = "",
               Pg_pop = Pg, nPg_pop = n_trial * Pg,
               PQg_pop = if (is_null) 0 else PQg,
               sens_pop = if (is_null) NA_real_ else PgQ / PQ,
               spec_pop = 1 - PgQc / (1 - PQ),
               beta_pop = if (is_null) beta_null else tauQc + bint * PQg,
               stringsAsFactors = FALSE)
  }))
  ev$sig <- vapply(ev$sg_def, rule_sig, character(1), USE.NAMES = FALSE)
  cat(sprintf("§3: %d distinct rules; resolved ok %d, empty %d, unresolved %d\n",
              nrow(ev), sum(ev$status == "ok"), sum(ev$status == "empty"), sum(ev$status == "unresolved")))
  if (any(ev$status == "unresolved")) print(ev[ev$status == "unresolved", c("sg_def", "missing")])
  # join population values onto replicates
  evj <- ev[, c("sg_def", "status", "Pg_pop", "nPg_pop", "PQg_pop", "sens_pop", "spec_pop", "beta_pop")]
  names(evj)[2] <- "pop_status"
  d <- merge(d, evj, by = "sg_def", all.x = TRUE, sort = FALSE)
  # §4 within-signature comparison
  pair <- function(x, y) { z <- x - y; z <- z[is.finite(z)]; c(diff = mean(z), se = stats::sd(z) / sqrt(length(z)), n = length(z)) }
  sigs20 <- names(ts)[ts >= 20]
  w <- lapply(sigs20, function(s) {
    r <- d[d$sig == s, ]
    data.frame(signature = s, n = nrow(r),
               size_real = mean(r$n_harm), size_pop = mean(r$nPg_pop),
               size_diff = pair(r$n_harm, r$nPg_pop)[["diff"]], size_se = pair(r$n_harm, r$nPg_pop)[["se"]],
               ppv_real = mean(r$ppv), ppv_pop = mean(r$PQg_pop),
               ppv_diff = pair(r$ppv, r$PQg_pop)[["diff"]], ppv_se = pair(r$ppv, r$PQg_pop)[["se"]],
               sens_real = mean(r$sens), sens_pop = mean(r$sens_pop),
               sens_diff = pair(r$sens, r$sens_pop)[["diff"]], sens_se = pair(r$sens, r$sens_pop)[["se"]],
               spec_real = mean(r$spec), spec_pop = mean(r$spec_pop),
               spec_diff = pair(r$spec, r$spec_pop)[["diff"]], spec_se = pair(r$spec, r$spec_pop)[["se"]],
               beta_real = mean(r$beta_or), beta_pop = mean(r$beta_pop),
               beta_diff = pair(r$beta_or, r$beta_pop)[["diff"]], beta_se = pair(r$beta_or, r$beta_pop)[["se"]],
               naive_bias_real = mean(r$naive_bias), naive_bias_se = stats::sd(r$naive_bias) / sqrt(nrow(r)),
               stringsAsFactors = FALSE)
  })
  w <- do.call(rbind, w)
  # pooled over ALL replicates (paired), and over the >=20 signatures weighted by count
  pooled_all <- data.frame(
    quantity = c("size", "ppv", "sens", "spec", "beta"),
    real = c(mean(d$n_harm), mean(d$ppv), mean(d$sens), mean(d$spec), mean(d$beta_or)),
    pop  = c(mean(d$nPg_pop), mean(d$PQg_pop), mean(d$sens_pop), mean(d$spec_pop), mean(d$beta_pop)),
    diff = c(pair(d$n_harm, d$nPg_pop)[["diff"]], pair(d$ppv, d$PQg_pop)[["diff"]], pair(d$sens, d$sens_pop)[["diff"]],
             pair(d$spec, d$spec_pop)[["diff"]], pair(d$beta_or, d$beta_pop)[["diff"]]),
    se   = c(pair(d$n_harm, d$nPg_pop)[["se"]], pair(d$ppv, d$PQg_pop)[["se"]], pair(d$sens, d$sens_pop)[["se"]],
             pair(d$spec, d$spec_pop)[["se"]], pair(d$beta_or, d$beta_pop)[["se"]]))
  wsum <- function(col) if (nrow(w)) sum(w[[col]] * w$n) / sum(w$n) else NA_real_
  pooled_sig <- data.frame(quantity = c("size", "ppv", "sens", "spec", "beta"),
                           diff_weighted = c(wsum("size_diff"), wsum("ppv_diff"), wsum("sens_diff"), wsum("spec_diff"), wsum("beta_diff")),
                           n_reps = sum(w$n), n_sigs = nrow(w))
  cat(sprintf("\n§4: signatures with >= 20 replicates: %d (covering %d of %d detected replicates)\n", nrow(w), sum(w$n), nrow(d)))
  print(w[order(-w$n), c("signature", "n", "size_real", "size_pop", "size_diff", "size_se", "ppv_real", "ppv_pop", "ppv_diff", "ppv_se")], digits = 4, row.names = FALSE)
  print(w[order(-w$n), c("signature", "n", "sens_diff", "sens_se", "spec_diff", "spec_se", "beta_real", "beta_pop", "beta_diff", "beta_se", "naive_bias_real", "naive_bias_se")], digits = 4, row.names = FALSE)
  cat("pooled over ALL detected replicates (paired realized - population):\n"); print(pooled_all, digits = 4, row.names = FALSE)
  cat("pooled over the >= 20 signatures, weighted by count:\n"); print(pooled_sig, digits = 4, row.names = FALSE)
  # naive bias: aggregate realized vs analytic, and within-signature sign
  an <- if (is_null) N0$runs$resample$pred else G$runs[[sprintf("n%d_resample", n_trial)]]$pred
  an_s <- if (is_null) N0$runs$split$pred else G$runs[[sprintf("n%d_split", n_trial)]]$pred
  cat(sprintf("naive bias: realized aggregate %.3f (SE %.3f); analytic resample %.3f, split %.3f; within-signature realized minus analytic(resample): min %.2f, max %.2f, share of signatures above analytic %.2f\n",
              mean(d$naive_bias), stats::sd(d$naive_bias) / sqrt(nrow(d)), an$Enaive_bias, an_s$Enaive_bias,
              min(w$naive_bias_real - an$Enaive_bias), max(w$naive_bias_real - an$Enaive_bias),
              mean(w$naive_bias_real > an$Enaive_bias)))
  # by K (a mix axis the aggregate cannot see)
  byK <- do.call(rbind, lapply(sort(unique(d$K)), function(k) { r <- d[d$K == k, ]
    data.frame(K = k, n = nrow(r), size_real = mean(r$n_harm), size_pop = mean(r$nPg_pop), naive_bias = mean(r$naive_bias),
               ppv_real = mean(r$ppv), ppv_pop = mean(r$PQg_pop)) }))
  cat("by K (factors in the winner):\n"); print(byK, digits = 4, row.names = FALSE)
  # §5 analytic selection distribution over signatures
  an_lab <- an$lab; an_sig <- vapply(an_lab, rule_sig, character(1), USE.NAMES = FALSE)
  agg <- function(p) tapply(p, an_sig, sum)
  a_res <- agg(an$sel_c); a_spl <- agg(an_s$sel_c)
  m_freq <- ts / sum(ts)
  all_sigs <- union(names(a_res), names(m_freq))
  A <- setNames(rep(0, length(all_sigs)), all_sigs); A[names(a_res)] <- a_res
  Asp <- setNames(rep(0, length(all_sigs)), all_sigs); Asp[names(a_spl)] <- a_spl
  Mf <- setNames(rep(0, length(all_sigs)), all_sigs); Mf[names(m_freq)] <- m_freq
  top_an <- names(sort(A, decreasing = TRUE))[1:15]
  cmp <- data.frame(signature = top_an, analytic_resample = round(A[top_an], 4), analytic_split = round(Asp[top_an], 4),
                    measured = round(Mf[top_an], 4), measured_n = as.integer(round(Mf[top_an] * sum(ts))))
  cat("\n§5 top 15 analytic signatures (sel_c) beside measured frequency:\n"); print(cmp, row.names = FALSE)
  never_sel <- sum(A[Mf == 0]); absent_an <- sum(Mf[A == 0])
  absent_sigs <- names(Mf)[A == 0 & Mf > 0]
  uses_str2 <- grepl("str2", absent_sigs)
  cat(sprintf("analytic mass on signatures the search never selected: resample %.4f, split %.4f\n", never_sel, sum(Asp[Mf == 0])))
  cat(sprintf("measured frequency on signatures absent from the analytic family: %.4f (%d signatures; %d of them use str2, carrying %.4f; %d others carry %.4f)\n",
              absent_an, length(absent_sigs), sum(uses_str2), sum(Mf[absent_sigs][uses_str2]), sum(!uses_str2), sum(Mf[absent_sigs][!uses_str2])))
  if (any(!uses_str2)) print(sort(Mf[absent_sigs][!uses_str2], decreasing = TRUE)[1:min(10, sum(!uses_str2))])
  tvd <- 0.5 * sum(abs(A - Mf)); tvd_s <- 0.5 * sum(abs(Asp - Mf))
  keep <- !grepl("str2", all_sigs)
  Mf2 <- Mf[keep] / sum(Mf[keep]); A2 <- A[keep] / sum(A[keep])
  tvd2 <- 0.5 * sum(abs(A2 - Mf2))
  cat(sprintf("total variation distance over signatures: resample %.4f, split %.4f; excluding str2 signatures and renormalising: %.4f\n", tvd, tvd_s, tvd2))
  # mass on the top-15 measured signatures
  top_m <- names(sort(Mf, decreasing = TRUE))[1:15]
  cat("top 15 measured signatures beside analytic:\n")
  print(data.frame(signature = top_m, measured = round(Mf[top_m], 4), measured_n = as.integer(round(Mf[top_m] * sum(ts))),
                   analytic_resample = round(A[top_m], 4)), row.names = FALSE)
  cells[[cell]] <- list(n = n_trial, null = is_null, detected = nrow(d), n_reps = nrow(res),
                        verbatim = tv, signatures = ts, K = table(d$K), covariate_use = use,
                        rules = ev, reps = d[, c("sim_id", "sg_def", "sig", "K", "n_harm", "ppv", "sens", "spec", "beta_or", "naive_or", "naive_bias",
                                                  "pop_status", "Pg_pop", "nPg_pop", "PQg_pop", "sens_pop", "spec_pop", "beta_pop")],
                        within = w, pooled_all = pooled_all, pooled_sig = pooled_sig, byK = byK,
                        analytic = list(lab = an_lab, sig = an_sig, sel_c_resample = an$sel_c, sel_c_split = an_s$sel_c,
                                        Enaive_bias_resample = an$Enaive_bias, Enaive_bias_split = an_s$Enaive_bias,
                                        EnH_resample = an$EnH, Esens_resample = an$Esens, Eppv_resample = an$Eppv),
                        dist = list(A_resample = A, A_split = Asp, M = Mf, tvd_resample = tvd, tvd_split = tvd_s, tvd_excl_str2 = tvd2,
                                    never_selected_mass = never_sel, absent_from_analytic = absent_an, absent_sigs = absent_sigs))
}
saveRDS(list(cells = cells, tauQc = tauQc, bint = bint, beta_null = beta_null, built_at = Sys.time(),
             pkg_version = as.character(utils::packageVersion("forestsearch"))), OUT)
cat("\nwritten:", OUT, "\n")
