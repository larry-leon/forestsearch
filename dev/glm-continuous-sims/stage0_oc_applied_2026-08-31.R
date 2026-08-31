# =============================================================================
# stage0_oc_applied_2026-08-31.R
# Applied OC evaluation, stage 0: feasibility gates for anchoring the wrapper
# on the ACTG175 continuous analysis (cc_task_oc_applied_stage0_2026-08-31.md)
# Read-only against R/.  Outputs: this log + stage0_oc_applied_2026-08-31.rds
# =============================================================================
suppressPackageStartupMessages({
  library(forestsearch)
  library(speff2trial)
})
cat("forestsearch version:", as.character(packageVersion("forestsearch")), "\n")
cat("R version:", R.version.string, "\n")
res <- list()  # accumulated results, saved at the end

out_dir <- "dev/glm-continuous-sims"
stopifnot(dir.exists(out_dir))

# ---------------------------------------------------------------------------
# §2.2 Data prep, replicated exactly from the compare-all document
# ---------------------------------------------------------------------------
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
cat("\n[§2.2] Analysis N =", N_analysis, "\n")
res$N_analysis <- N_analysis

# ---------------------------------------------------------------------------
# §3.1 ITT MD on y_decline
# ---------------------------------------------------------------------------
itt <- lm(y_decline ~ treat, data = actg_df)
itt_row <- summary(itt)$coefficients["treat", , drop = FALSE]
cat("\n[§3.1] ITT MD on y_decline (adverse scale):\n")
print(itt_row, digits = 8)
res$itt <- itt_row

# ---------------------------------------------------------------------------
# §2.3 The anchor.  The committed payloads do not hold the fixed-family
# maxeffCons result (selected_subgroups_continuous.rds row 8 is the
# comparison run with use_lasso/use_grf ON; comparison_continuous.rds holds
# the same eight runs and no fs_anchor), so the single authorized anchor
# call runs here: the MR-anchor chunk minus mr_inference, sequential plan.
# ---------------------------------------------------------------------------
cat("\n[§2.3] Running the fixed-family maxeffCons anchor (sequential)...\n")
t_anchor <- proc.time()
fs_anchor <- forestsearch(
  df.analysis      = actg_df,
  confounders.name = confounders.name,
  outcome.name     = "y_decline",
  treat.name       = "treat",
  id.name          = "id",
  outcome_type     = "continuous",
  effect_measure   = "MD",
  adverse_outcome  = TRUE,
  seedit           = analysis_seed,

  sg_focus         = "maxeffCons",
  selection_rule   = "neighborhood",
  consistency_method = "resample",

  effect.threshold       = 10,
  consistency.threshold  = 5,
  pconsistency.threshold = 0.90,
  use_twostage     = TRUE,
  conf.cont_jcuts  = list(age = 10, preanti = 10, wtkg = 10,
                          karnof = 10, cd40 = 10, cd80 = 10),
  cut_type         = "default",
  maxk             = 2L,
  n.min            = 60L,
  d0.min           = 10L,
  d1.min           = 10L,
  fs.splits        = 500L,

  use_lasso        = FALSE,
  use_grf          = FALSE,
  use_dina         = FALSE,

  is.RCT           = TRUE,
  parallel_args    = list(plan = "sequential"),
  details          = FALSE,
  quiet            = FALSE
)
t_anchor <- (proc.time() - t_anchor)["elapsed"]
cat(sprintf("\n[§2.3] anchor wall-clock (sequential): %.1f s\n", t_anchor))
res$t_anchor <- t_anchor

if (is.null(fs_anchor$sg.harm) || length(fs_anchor$sg.harm) == 0L) {
  stop("GATE (§2.3): the anchor declares nothing (no subgroup); ",
       "the design's Q needs a found subgroup.  Stopping.")
}
H_def <- paste(fs_anchor$sg.harm, collapse = " & ")
cat("[§2.3] anchor Ĥ definition:", H_def, "\n")

stopifnot(!is.null(fs_anchor$df.est),
          "treat.recommend" %in% names(fs_anchor$df.est))
# treat.recommend == 0 marks membership in the harm subgroup Ĥ
memb_anchor <- fs_anchor$df.est[
  fs_anchor$df.est$treat.recommend == 0L, "id", drop = TRUE]
n_H <- length(memb_anchor)

win <- fs_anchor$grp.consistency$out_sg$result[1L, ]
T_obs  <- as.numeric(win$hr)      # MD scale: hr column holds the MD
p_cons <- as.numeric(win$Pcons)
cat(sprintf("[§2.3] n(Ĥ) = %d;  T̂_obs (fitted MD) = %.6f;  p.consistency = %.4f\n",
            n_H, T_obs, p_cons))
res$anchor <- list(def = H_def, sg.harm = fs_anchor$sg.harm, n_H = n_H,
                   T_obs = T_obs, p_cons = p_cons, ids = memb_anchor)

# ---------------------------------------------------------------------------
# §3.2 Ĥ's clauses -> subgroup_vars / subgroup_cuts, verified subject-for-
# subject on the analysis frame, then the anchored DGM at q = 20
# ---------------------------------------------------------------------------
# Parse each clause "{var <= c}" / "!{var <= c}" (the fixed-family cut forms).
parse_clause <- function(cl) {
  m <- regmatches(cl, regexec(
    "^(!?)\\{\\s*([A-Za-z0-9_.]+)\\s*(<=|>=|==|!=|<|>)\\s*([-0-9.eE]+)\\s*\\}$",
    cl))[[1]]
  if (length(m) == 0L) stop("unparseable clause: ", cl)
  list(neg = identical(m[2], "!"), var = m[3], op = m[4],
       val = as.numeric(m[5]))
}
clauses <- lapply(fs_anchor$sg.harm, parse_clause)

subgroup_vars <- vapply(clauses, `[[`, character(1), "var")
subgroup_cuts <- lapply(clauses, function(cl) {
  if (!cl$neg && cl$op == "<=") cl$val                       # var <= c
  else if (cl$neg && cl$op == "<=")                          # !(var <= c): var > c
    list(type = "greater", value = cl$val)
  else stop("no subgroup_cuts mapping for clause form: ",
            if (cl$neg) "!" else "", "{", cl$var, " ", cl$op, " ", cl$val, "}")
})
names(subgroup_cuts) <- subgroup_vars
cat("\n[§3.2] clause -> subgroup_cuts mapping:\n")
for (i in seq_along(clauses)) {
  cl <- clauses[[i]]
  cat(sprintf("  %s%s %s %s  ->  %s = %s\n",
              if (cl$neg) "NOT " else "", cl$var, cl$op, format(cl$val),
              cl$var,
              if (is.list(subgroup_cuts[[i]]))
                sprintf("list(type='greater', value=%s)", format(cl$val))
              else format(cl$val)))
}

# Reconstructed flag on the analysis frame, from the same semantics
# process_cutpoint() applies: numeric -> (var <= c); greater -> (var > c)
flag_recon <- rep(TRUE, nrow(actg_df))
for (i in seq_along(clauses)) {
  cl <- clauses[[i]]
  v  <- actg_df[[cl$var]]
  ind <- if (is.list(subgroup_cuts[[i]])) v > subgroup_cuts[[i]]$value
         else v <= subgroup_cuts[[i]]
  flag_recon <- flag_recon & ind
}
ids_recon <- actg_df$id[flag_recon]
match_exact <- setequal(ids_recon, memb_anchor) &&
  sum(flag_recon) == n_H
cat(sprintf("[§3.2] reconstructed flag: n = %d;  matches anchor membership subject-for-subject: %s\n",
            sum(flag_recon), match_exact))
if (!match_exact) {
  stop("GATE (§3.2): reconstructed flag does not match the anchor's ",
       "membership exactly.  Stopping.")
}
res$recon <- list(n = sum(flag_recon), match = match_exact)

# The linear relation for k_inter: on the identity link,
# m_tau[Q] = k_treat * beta_treat + k_inter, so k_inter(q) = q - beta_treat.
fml_base <- as.formula(paste("y_decline ~ treat +",
                             paste(c(bin_vars, cont_vars), collapse = " + ")))
fit_base <- glm(fml_base, data = actg_df, family = gaussian())
beta_treat <- unname(coef(fit_base)["treat"])
cat(sprintf("\n[§3.2] baseline GLM beta_treat = %.6f;  k_inter(q) = q - beta_treat\n",
            beta_treat))
res$beta_treat <- beta_treat

build_dgm <- function(q, model = "alt") {
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
    model           = model,
    k_treat         = 1,
    k_inter         = if (model == "alt") q - beta_treat else 0,
    adverse_outcome = TRUE,
    n_super         = 5000L,
    seed            = analysis_seed,
    verbose         = FALSE
  )
}

dgm20 <- build_dgm(20)
sc20  <- fs_dgm_scale(dgm20)
reg20 <- sc20$regions
mQ20  <- reg20$m_tau[reg20$region == "Q"]
cat(sprintf("[§3.2] DGM at q = 20: fs_dgm_scale readback m_tau[Q] = %.10f (target 20)\n",
            mQ20))
stopifnot(abs(mQ20 - 20) < 1e-8)
res$dgm20_readback <- mQ20

# ---------------------------------------------------------------------------
# §3.3 The signed orientation table, one DGM build per rung
# ---------------------------------------------------------------------------
q_rungs <- c(2.5, 5, 7.5, 10, 15, 20, 30, 40, T_obs)
tab <- do.call(rbind, lapply(q_rungs, function(q) {
  d <- build_dgm(q)
  r <- fs_dgm_scale(d)$regions
  mQ  <- r$m_tau[r$region == "Q"]
  mQc <- r$m_tau[r$region == "Qc"]
  data.frame(q = q, m_tau_Q = mQ, m_tau_Qc = mQc,
             same_sign = sign(mQ) == sign(mQc))
}))
rownames(tab) <- NULL
cat("\n[§3.3] signed orientation table (last rung: q = T̂_obs):\n")
print(tab, digits = 8)
res$sign_table <- tab

# ---------------------------------------------------------------------------
# §3.4 One enumeration attempt at a rung where the signs differ
# ---------------------------------------------------------------------------
fs_args <- list(
  confounders.name = confounders.name,
  conf.cont_jcuts  = list(age = 10, preanti = 10, wtkg = 10,
                          karnof = 10, cd40 = 10, cd80 = 10),
  cut_type         = "default",
  cont.cutoff      = 4L,
  maxk             = 2L,
  n.min            = 60L
)
diff_rung <- tab$q[!tab$same_sign]
if (length(diff_rung) > 0L) {
  qd <- diff_rung[which(abs(diff_rung - 20) == min(abs(diff_rung - 20)))][1]
  cat(sprintf("\n[§3.4] enumeration attempt at differing-sign rung q = %s:\n", qd))
  msg <- tryCatch({
    fs_oc_family_enumerate(build_dgm(qd), fs_args, n = N_analysis,
                           max_M = 10000L, verbose = TRUE)
    "(no stop: enumeration succeeded)"
  }, error = function(e) conditionMessage(e))
  cat("[§3.4] captured message: ", msg, "\n", sep = "")
  res$stop_message <- msg
  res$stop_rung <- qd
} else {
  cat("\n[§3.4] no rung with differing signs; no enumeration attempt.\n")
  res$stop_message <- NA_character_
}

# §3.5 verdict
binds <- any(!tab$same_sign)
res$verdict_binds <- binds
res$binding_rungs <- tab$q[!tab$same_sign]
cat(sprintf("\n[§3.5] VERDICT: the sign gate %s%s\n",
            if (binds) "BINDS" else "does NOT bind",
            if (binds) paste0(" (rungs: ",
                              paste(format(res$binding_rungs), collapse = ", "),
                              ")") else ""))

# ---------------------------------------------------------------------------
# §4 The family via the null DGM, at the analyst's spec and n = N_analysis
# ---------------------------------------------------------------------------
cat("\n[§4] building the null DGM and enumerating the family...\n")
dgm_null <- build_dgm(0, model = "null")
stopifnot(sum(dgm_null$df_super$flag_harm) == 0L)

t_enum <- proc.time()
fam <- fs_oc_family_enumerate(dgm_null, fs_args, n = N_analysis,
                              max_M = 10000L, verbose = TRUE)
t_enum <- (proc.time() - t_enum)["elapsed"]
cat(sprintf("[§4] enumeration wall-clock: %.2f s\n", t_enum))
print(fam)
res$fam_counts <- fam$counts
res$t_enum <- t_enum
res$M <- fam$M

# Ĥ on df_super via the reconstructed flag (memberships depend only on the
# covariates, identical across null/alt builds at the same seed)
ds <- dgm_null$df_super
H_super <- rep(TRUE, nrow(ds))
for (i in seq_along(clauses)) {
  cl <- clauses[[i]]
  v  <- ds[[cl$var]]
  ind <- if (is.list(subgroup_cuts[[i]])) v > subgroup_cuts[[i]]$value
         else v <= subgroup_cuts[[i]]
  H_super <- H_super & ind
}
P_H <- mean(H_super)
cat(sprintf("\n[§4] Ĥ on df_super: %d of %d (P = %.4f)\n",
            sum(H_super), nrow(ds), P_H))

purity  <- colMeans(fam$memb & H_super) / fam$Pg          # P(g & Ĥ)/P(g)
inter   <- colMeans(fam$memb & H_super)
uni     <- fam$Pg + P_H - inter
jaccard <- inter / uni
i_pur <- which.max(purity); i_jac <- which.max(jaccard)
cat(sprintf("[§4] max purity  P(g∩Ĥ)/P(g) = %.4f  at member %d: %s (Pg = %.4f)\n",
            purity[i_pur], i_pur, fam$lab[i_pur], fam$Pg[i_pur]))
cat(sprintf("[§4] max Jaccard              = %.4f  at member %d: %s (Pg = %.4f)\n",
            jaccard[i_jac], i_jac, fam$lab[i_jac], fam$Pg[i_jac]))
identical_member <- which(inter == fam$Pg & inter == P_H)  # g == Ĥ exactly
cat(sprintf("[§4] member identical to Ĥ: %s\n",
            if (length(identical_member)) paste0("YES — member ",
              identical_member[1], ": ", fam$lab[identical_member[1]])
            else "none"))
res$H_super <- list(P_H = P_H, n_H = sum(H_super))
res$purity <- list(max_purity = purity[i_pur], max_jaccard = jaccard[i_jac],
                   nearest_lab = fam$lab[i_jac], nearest_Pg = fam$Pg[i_jac],
                   nearest_purity = purity[i_jac],
                   identical = length(identical_member) > 0L)

# ---------------------------------------------------------------------------
# §5 Cost: one fs_oc_predict at draws = 2e4, extrapolated to 2e5
# ---------------------------------------------------------------------------
cat("\n[§5] timing fs_oc_predict on the null family, draws = 2e4 ...\n")
invisible(gc(reset = TRUE, verbose = FALSE))
t_pred <- proc.time()
pred <- fs_oc_predict(family = fam, n = N_analysis, c1 = 10, c2 = 5,
                      consistency_method = "resample", pconsistency = 0.90,
                      draws = 2e4, seed = analysis_seed)
t_pred <- (proc.time() - t_pred)["elapsed"]
g <- gc(verbose = FALSE)
peak_mb <- sum(g[, "max used"] * c(56, 8)) / 2^20   # Ncells*56B + Vcells*8B
cat(sprintf("[§5] fs_oc_predict wall-clock: %.2f s;  peak memory (R gc): %.0f MB\n",
            t_pred, peak_mb))
print(pred)
res$t_pred_2e4 <- t_pred
res$peak_mb <- peak_mb

t_2e5 <- t_pred * 10          # linear in draws
n_sets_ladder <- 11           # 10 rungs + the null
n_sets_sweeps <- 2 * 10       # two sweeps at ~10 grid points each (assumed)
n_sets_invert <- 2 * 10       # inversions, ~10 evaluations each (assumed)
total_h <- (n_sets_ladder + n_sets_sweeps + n_sets_invert) * t_2e5 / 3600
cat(sprintf("\n[§5] extrapolation: 2e5 draws ≈ %.1f s per set (linear).\n", t_2e5))
cat(sprintf("[§5] stage-2 projection: (%d ladder + %d sweep + %d inversion) sets ≈ %.2f h\n",
            n_sets_ladder, n_sets_sweeps, n_sets_invert, total_h))
cat(sprintf("[§5] M×M memory: M = %d -> %.1f MB per M×M matrix (8B doubles)\n",
            fam$M, 8 * fam$M^2 / 2^20))
res$projection <- list(t_2e5 = t_2e5, total_h = total_h,
                       mxm_mb = 8 * fam$M^2 / 2^20)

# ---------------------------------------------------------------------------
# §6 The scale table from the q = 20 anchored DGM
# ---------------------------------------------------------------------------
cat("\n[§6] scale table at q = 20:\n")
cat(sprintf("  sigma = %.6f (sigma^2 = %.4f)\n", sc20$sigma, sc20$sigma^2))
print(reg20[, c("region", "n_g", "P_g", "m_tau", "V_eff")], digits = 8)
ratio <- reg20$V_eff[reg20$region == "Q"] / reg20$V_eff[reg20$region == "S"]
cat(sprintf("  V_eff[Q] / V_eff[S] = %.6f\n", ratio))
res$scale <- list(sigma = sc20$sigma,
                  regions = reg20[, c("region", "n_g", "P_g", "m_tau", "V_eff")],
                  veff_ratio = ratio)

saveRDS(res, file.path(out_dir, "stage0_oc_applied_2026-08-31.rds"))
cat("\n[done] results saved to", file.path(out_dir, "stage0_oc_applied_2026-08-31.rds"), "\n")
