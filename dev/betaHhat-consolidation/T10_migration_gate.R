# T10_migration_gate.R
# ---------------------------------------------------------------------------
# The migration gate for consolidating betaHhat into package functions.
#
# Compares each simulation module's betaHhat_one*() against the package
# fs_betaHhat_one() on a FIXED fixture per outcome family.  The spec changes no
# estimand, so any movement is a defect in the consolidation until shown
# otherwise.
#
# GATE SEMANTICS -- read this before editing the summary.
#
# An NA difference is a FAILURE, never a skip.  The first version of this
# script computed max(abs(diff), na.rm = TRUE), which silently dropped the one
# row where the old module returned NA and the new function returned a finite
# value -- and then printed "PASS -- every value bitwise identical".  It was
# not.  na.rm = TRUE in a gate summary is the same vacuous-check class as a
# pgrep that matches its own command line: the check cannot fail, so it proves
# nothing.  Three categories are now reported separately and only the first is
# silent:
#
#   identical      both finite and diff == 0
#   moved          both finite and diff != 0            -> FAIL
#   na_mismatch    exactly one side NA                  -> FAIL unless
#                                                          explicitly sanctioned
#
# SANCTIONED movements are listed in .SANCTIONED below, each with the reason.
# A row not on that list that is not `identical` fails the gate.
#
# Result of record: dev/betaHhat-consolidation/T10_GATE_RESULT.md
# ---------------------------------------------------------------------------

setwd("/home/larryleon/Documents/GitHub/forestsearch")
suppressMessages(devtools::load_all(quiet = TRUE))
suppressPackageStartupMessages(library(speff2trial))
`%||%` <- function(a, b) if (is.null(a)) b else a

# Rows allowed to move, and why.  Anything else that moves fails.
.SANCTIONED <- data.frame(
  family = "survival", rule = "disjunction",
  reason = paste("D1 fix: betaHhat_truth.R splits on \" & \" before testing",
                 "for \"|\", shredding the rule.  The old module cannot",
                 "express a disjunction at all; binary/continuous are exactly",
                 "0 here because e6f6024 already fixed D1 for them."),
  stringsAsFactors = FALSE)

cmp <- function(family, rule_id, old_H, old_Hc, old_nH, old_nHc, new) {
  na_mis <- xor(is.na(old_H), is.na(new$betaHhat_H)) ||
            xor(is.na(old_Hc), is.na(new$betaHhat_Hc)) ||
            xor(is.na(old_nH), is.na(new$nH_eval))
  dH  <- if (is.na(old_H)  && is.na(new$betaHhat_H))  0 else
           new$betaHhat_H  - old_H
  dHc <- if (is.na(old_Hc) && is.na(new$betaHhat_Hc)) 0 else
           new$betaHhat_Hc - old_Hc
  data.frame(family = family, rule = rule_id,
             old_H = old_H, new_H = new$betaHhat_H, diff_H = dH,
             old_Hc = old_Hc, new_Hc = new$betaHhat_Hc, diff_Hc = dHc,
             old_nH = old_nH, new_nH = new$nH_eval,
             old_nHc = old_nHc, new_nHc = new$nHc_eval,
             na_mismatch = na_mis,
             stringsAsFactors = FALSE)
}
rows <- list()

# --- CONTINUOUS: calibrated ACTG175 DGM from acceptance_betaHhat_md.qmd -----
cat("=== building the continuous fixture ===\n")
abv <- c("hemo","homo","drugs","race","gender","symptom")
a <- subset(speff2trial::ACTG175, arms %in% c(1L,3L)); a$id <- seq_len(nrow(a))
a$treat_orig <- ifelse(a$arms==1L,1L,0L); a$treat <- 1L-a$treat_orig
a$cd4_change <- a$cd420-a$cd40; a <- a[!is.na(a$cd420),]
a$z1<-as.factor(ifelse(a$age>34,1L,0L)); a$z2<-as.factor(ifelse(a$preanti<=744.5,1L,0L))
a$z3<-as.factor(ifelse(a$wtkg<=75,1L,0L)); a$z4<-as.factor(ifelse(a$karnof<=median(a$karnof),1L,0L))
a$z5<-as.factor(ifelse(a$cd40<=median(a$cd40),1L,0L)); a$z6<-as.factor(ifelse(a$cd80<=median(a$cd80),1L,0L))
a$z7<-as.factor(a$hemo);a$z8<-as.factor(a$homo);a$z9<-as.factor(a$drugs)
a$z10<-as.factor(a$race);a$z11<-as.factor(a$gender);a$z12<-as.factor(a$symptom)
for (v in abv) a[[v]] <- as.factor(a[[v]])
dgm_c <- calibrate_glm_interaction(data=a, factor_vars=paste0("z",1:12),
  outcome_var="cd4_change", treatment_var="treat", target_effect=40,
  outcome_type="continuous", effect_measure="MD", subgroup_vars=c("z1","z2"),
  subgroup_cuts=list(z1=1L,z2=1L), k_inter_range=c(0,120), grid_step=2,
  n_super=5000L, seed=8316951L, verbose=FALSE)
S_c <- dgm_c$df_super
cat(sprintf("  df_super n=%d  MD(Q)=%.6f\n", nrow(S_c), dgm_c$hazard_ratios$harm_subgroup))
env_md <- new.env()
source("quarto/simulations/actg175/continuous/betaHhat_truth_md.R", local = env_md)

RULES_C <- list(
  exact_Q     = c("{age > 34}", "{preanti <= 744.5}"),
  narrow      = c("{age > 34}", "{preanti <= 744.5}", "{wtkg <= 75}"),
  one_cut     = c("{age > 34}"),
  disjunction = "(age > 34 & preanti <= 744.5) | (wtkg <= 60)",
  negation    = c("!{cd40 <= 400}"),
  disjoint_Q  = c("{age <= 34}"))
for (rn in names(RULES_C)) {
  o <- env_md$betaHhat_one_md(RULES_C[[rn]], S_c, focus = "harm")
  n <- fs_betaHhat_one(RULES_C[[rn]], S_c, focus = "harm",
                       outcome_type = "continuous", effect_measure = "MD")
  rows[[length(rows)+1]] <- cmp("continuous", rn, unname(o[["betaHhat_H"]]),
    unname(o[["betaHhat_Hc"]]), as.integer(o[["nH_eval"]]),
    as.integer(o[["nHc_eval"]]), n)
}

# --- BINARY: the DGM behind the actg175_or075_* mr_sweep bundles ------------
cat("=== building the binary fixture ===\n")
cont_vars <- c("age","preanti","wtkg","karnof","cd40","cd80")
bin_vars  <- c("hemo","homo","drugs","race","gender","symptom")
b <- subset(speff2trial::ACTG175, arms %in% c(1L,3L))
b$id <- seq_len(nrow(b)); b$treat <- ifelse(b$arms==1L,1L,0L)
b <- b[!is.na(b$cd420),]
b$y_neg <- 1L - as.integer(b$cd420 > b$cd40)
dgm_b <- calibrate_glm_interaction(data=b, factor_vars=bin_vars,
  continuous_vars=cont_vars, outcome_var="y_neg", treatment_var="treat",
  target_effect=0.75, outcome_type="binary", effect_measure="OR",
  subgroup_vars=c("wtkg","cd40"),
  subgroup_cuts=list(wtkg=list(type="greater", quantile=0.70),
                     cd40=list(type="greater", quantile=0.70)),
  k_treat=1, adverse_outcome=FALSE, k_inter_range=c(0.3,1.5), grid_step=0.025,
  n_super=100000L, seed=8316951L)
env_or <- new.env()
source("quarto/simulations/actg175/binary/betaHhat_truth_glm.R", local = env_or)
E_b <- env_or$build_eval_frame_glm(dgm_b, eval_seed = 20260628L)
cat(sprintf("  eval frame n=%d  OR(H)=%.6f\n", nrow(E_b), dgm_b$hazard_ratios$AHR_harm))

qw <- unname(quantile(b$wtkg, 0.70)); qc <- unname(quantile(b$cd40, 0.70))
RULES_B <- list(
  exact_Q     = c(sprintf("{wtkg > %s}", qw), sprintf("{cd40 > %s}", qc)),
  one_cut     = c(sprintf("{wtkg > %s}", qw)),
  disjunction = sprintf("(wtkg > %s & cd40 > %s) | (age <= 30)", qw, qc),
  negation    = c(sprintf("!{cd80 <= %s}", unname(quantile(b$cd80, 0.5)))),
  disjoint_Q  = c(sprintf("{wtkg <= %s}", qw)))
for (rn in names(RULES_B)) {
  o <- env_or$betaHhat_one_or(RULES_B[[rn]], E_b)
  n <- fs_betaHhat_one(RULES_B[[rn]], E_b, focus = "harm",
                       outcome_type = "binary", outcome.name = "y_sim",
                       treat.name = "treat_sim")
  rows[[length(rows)+1]] <- cmp("binary", rn, unname(o[["betaHhat_H"]]),
    unname(o[["betaHhat_Hc"]]), as.integer(o[["nH_eval"]]),
    as.integer(o[["nHc_eval"]]), n)
}

# --- SURVIVAL: setup_gbsg_dgm() + build_eval_frame(eval_seed = 20260628) ----
cat("=== building the survival fixture ===\n")
env_s <- new.env()
source("quarto/simulations/gbsg_redux/betaHhat_truth.R", local = env_s)
k_inter <- calibrate_k_inter(target_hr_harm = 1.0, model = "alt", use_ahr = FALSE)
dgm_s <- setup_gbsg_dgm(model = "alt", k_inter = k_inter, n_super = 100000L,
                        seed = 8316951L)
E_s <- env_s$build_eval_frame(dgm_s, analysis_time = 84, cens_adjust = log(1.5),
                              eval_seed = 20260628L)
cat(sprintf("  eval frame n=%d\n", nrow(E_s)))

RULES_S <- list(
  conjunction = c("{er > 125}", "{size > 20}"),
  one_cut     = c("{nodes > 2}"),
  disjunction = "(er > 125 & size > 20) | (nodes > 5)",
  negation    = c("!{age <= 50}"),
  small       = c("{er > 400}", "{nodes > 8}"))
for (rn in names(RULES_S)) {
  # The old module's warnings on the disjunction ARE the D1 evidence; show them.
  o <- withCallingHandlers(env_s$betaHhat_one(RULES_S[[rn]], E_s),
        warning = function(w) {
          cat(sprintf("  [old module, %s] WARNING: %s\n", rn,
                      conditionMessage(w)))
          invokeRestart("muffleWarning")
        })
  n <- fs_betaHhat_one(RULES_S[[rn]], E_s, focus = "harm",
                       outcome_type = "survival", outcome.name = "y_sim",
                       event.name = "event_sim", treat.name = "treat_sim")
  rows[[length(rows)+1]] <- cmp("survival", rn, unname(o[["betaHhat_H"]]),
    unname(o[["betaHhat_Hc"]]), as.integer(o[["nH_eval"]]),
    as.integer(o[["nHc_eval"]]), n)
}

T10 <- do.call(rbind, rows)

# --- classification.  NO na.rm ANYWHERE IN THE VERDICT. ---------------------
T10$sanctioned <- paste(T10$family, T10$rule) %in%
                  paste(.SANCTIONED$family, .SANCTIONED$rule)
T10$class <- ifelse(T10$na_mismatch, "na_mismatch",
              ifelse(T10$diff_H == 0 & T10$diff_Hc == 0, "identical", "moved"))
T10$verdict <- ifelse(T10$class == "identical", "ok",
                ifelse(T10$sanctioned, "sanctioned", "FAIL"))

cat("\n\n================ T10 MIGRATION GATE ================\n")
for (fam in c("continuous","binary","survival")) {
  d <- T10[T10$family == fam, ]
  cat(sprintf("\n--- %s ---\n", toupper(fam)))
  print(d[, c("rule","old_H","new_H","diff_H","old_Hc","new_Hc","diff_Hc",
              "class","verdict")], row.names = FALSE, digits = 12)
  fin <- d$class == "identical"
  cat(sprintf("  identical: %d   moved: %d   na_mismatch: %d\n",
      sum(d$class=="identical"), sum(d$class=="moved"),
      sum(d$class=="na_mismatch")))
  if (any(fin))
    cat(sprintf("  max |diff| over identical+moved rows: %.3e\n",
        max(abs(c(d$diff_H[!d$na_mismatch], d$diff_Hc[!d$na_mismatch])))))
  cat(sprintf("  counts identical on comparable rows: %s\n",
      all(d$old_nH[fin] == d$new_nH[fin]) && all(d$old_nHc[fin] == d$new_nHc[fin])))
}

cat("\n================ VERDICT ================\n")
cat(sprintf("rows compared : %d\n", nrow(T10)))
cat(sprintf("  identical   : %d\n", sum(T10$class == "identical")))
cat(sprintf("  moved       : %d\n", sum(T10$class == "moved")))
cat(sprintf("  na_mismatch : %d\n", sum(T10$class == "na_mismatch")))
cat(sprintf("  sanctioned  : %d\n", sum(T10$verdict == "sanctioned")))
cat(sprintf("  FAIL        : %d\n", sum(T10$verdict == "FAIL")))
comparable <- !T10$na_mismatch
cat(sprintf("max |diff| over comparable rows: %.6e\n",
    max(abs(c(T10$diff_H[comparable], T10$diff_Hc[comparable])))))
if (any(T10$verdict == "sanctioned")) {
  cat("\nSANCTIONED MOVEMENTS:\n")
  for (i in which(T10$verdict == "sanctioned")) {
    s <- .SANCTIONED[.SANCTIONED$family == T10$family[i] &
                     .SANCTIONED$rule   == T10$rule[i], ]
    cat(sprintf("  %s / %s : old=%s new=%.12f  partition %d + %d = %d\n",
        T10$family[i], T10$rule[i], format(T10$old_H[i]), T10$new_H[i],
        T10$new_nH[i], T10$new_nHc[i], T10$new_nH[i] + T10$new_nHc[i]))
    cat(sprintf("    reason: %s\n", s$reason[1]))
  }
}
cat(sprintf("\nGATE: %s\n", if (any(T10$verdict == "FAIL")) "FAIL -- STOP"
            else "PASS -- no unsanctioned movement"))
saveRDS(T10, "dev/betaHhat-consolidation/T10_result.rds")
