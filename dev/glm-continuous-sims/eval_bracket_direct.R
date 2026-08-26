# =============================================================================
# eval_bracket_direct.R -- direct evaluation of the Level 4 bracket on df_super
# =============================================================================
#
# Queue item 2 of the continuous/MD workstream.
#
# PURPOSE
#   The scale anchor of
#     quarto/simulations/actg175/continuous/analytic_verification_and_prediction_md_harm.qmd
#   is currently obtained backwards: one measured sigma[betahat_1000(Q)] = 13.6786
#   is inverted to V_eff = 64,476, hence a Level 4 bracket of 16,119.  This
#   script computes the bracket FORWARDS, by enumeration, from quantities the
#   DGM object already stores.
#
# WHAT IS BEING COMPUTED
#   For a region g, with mu_w the arm-mean surfaces and tau = mu_1 - mu_0,
#
#     B(g) = sigma^2 + V_g[mu_0] + C_g[mu_0, tau] + (1/2) V_g[tau]
#     V_eff(g) = 4 * B(g)
#     sigma[betahat_n(g)] = sqrt( V_eff(g) / (n * P(g)) )
#
#   V_g[.] and C_g[.,.] are FINITE-POPULATION functionals over the rows of
#   df_super in g (divisor N, not N-1), matching the document's indexed region
#   moments m_g[a], V_g[a], C_g[a,b].  The N vs N-1 difference is under 0.06%
#   at |Q| = 1723 and is reported for completeness.
#
# INPUTS, ALL ALREADY STORED (verified against R/generate_glm_dgm.R)
#   dgm$model_params$sigma  -- L489-493, 534.  NOT an auxiliary estimate: this
#                              is the literal SD passed to rnorm() when outcomes
#                              are drawn (R/simulate_from_glm_dgm.R:118-120), so
#                              it is a DGM constant with no sampling error.
#   dgm$df_super$mu0, $mu1  -- L464-465.
#   dgm$df_super$flag_harm  -- true region membership.
#
# ORIENTATION
#   mu0, mu1 and hence tau are on the ORIGINAL coding (tau = -40 on Q).  The
#   bracket is orientation-invariant provided mu_0 and tau are taken from the
#   same coding, which they are here.  Do not flip one without the other:
#   V[.] terms are invariant but C[mu_0, tau] would change sign.
#
# REGIONS EVALUATED
#   Q, Q^c and S (everyone).  These require no candidate enumeration.  The
#   per-candidate bracket vector is deliberately OUT OF SCOPE -- it needs the
#   package's candidate-generation path, and it is only meaningful once the
#   formula itself is confirmed.
#
# VALIDATED DIAGNOSTIC
#   Six exact algebraic identities are checked BEFORE any bracket is reported.
#   If any fails, the bracket numbers are not printed.  This follows the
#   session lesson that a diagnostic which has not been validated against a
#   known positive is not evidence.
#
# COMPUTE
#   Zero.  Enumeration and arithmetic on a stored object.  No simulation, no
#   model fitting, no RNG draw of any kind, so no .Random.seed handling is
#   required and the caller's stream is untouched.
#
# USAGE
#   Rscript dev/glm-continuous-sims/eval_bracket_direct.R <path-to-dgm.rds>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("usage: Rscript eval_bracket_direct.R <path-to-dgm.rds>", call. = FALSE)
}
dgm_path <- args[[1L]]
if (!file.exists(dgm_path)) {
  stop("no such file: ", dgm_path, call. = FALSE)
}

dgm <- readRDS(dgm_path)

cat("=====================================================================\n")
cat("Direct bracket evaluation -- queue item 2\n")
cat("=====================================================================\n")
cat("file  : ", normalizePath(dgm_path), "\n", sep = "")
cat("mtime : ", format(file.info(dgm_path)$mtime, "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")


# -- 1. Object identity -------------------------------------------------------

cat("--- 1. object identity ---\n")
cat(sprintf("  class              : %s\n", paste(class(dgm), collapse = ", ")))
cat(sprintf("  outcome_type       : %s\n", dgm$outcome_type))
cat(sprintf("  effect_measure     : %s\n", dgm$effect_measure))
cat(sprintf("  model_type         : %s\n", dgm$model_type))
cat(sprintf("  n_super            : %d\n", nrow(dgm$df_super)))
cat(sprintf("  prevalence P(Q)    : %.6f\n", dgm$subgroup_info$proportion))
cat(sprintf("  |Q|                : %d\n", dgm$subgroup_info$size))
cat(sprintf("  effect_Q           : %.9f\n", dgm$hazard_ratios$harm_subgroup))
cat(sprintf("  effect_Qc          : %.9f\n", dgm$hazard_ratios$no_harm_subgroup))
cat(sprintf("  effect_ITT         : %.9f\n", dgm$hazard_ratios$overall))
cat(sprintf("  k_treat / k_inter  : %.6f / %.6f\n",
            dgm$model_params$k_treat, dgm$model_params$k_inter))
cat(sprintf("  beta_inter         : %.9f\n", dgm$model_params$beta_inter))
cat(sprintf("  adverse_outcome    : %s\n", dgm$adverse_outcome))
cat(sprintf("  noise_scheme       : %s\n", dgm$noise_scheme))
cat(sprintf("  sigma (residual SD): %.9f\n", dgm$model_params$sigma))
cat("\n")

stopifnot(
  "object is not a glm_dgm"       = inherits(dgm, "glm_dgm"),
  "outcome_type is not continuous" = identical(dgm$outcome_type, "continuous"),
  "effect_measure is not MD"       = identical(dgm$effect_measure, "MD"),
  "model_params$sigma is missing"  = !is.null(dgm$model_params$sigma),
  "df_super lacks mu0/mu1"         =
    all(c("mu0", "mu1", "flag_harm") %in% names(dgm$df_super))
)


# -- 2. Extract the three inputs ---------------------------------------------

sigma <- dgm$model_params$sigma
mu0   <- dgm$df_super$mu0
mu1   <- dgm$df_super$mu1
tau   <- mu1 - mu0

in_Q  <- dgm$df_super$flag_harm == 1
N     <- length(mu0)
piQ   <- mean(in_Q)

regions <- list(
  "Q"   = in_Q,
  "Qc"  = !in_Q,
  "S"   = rep(TRUE, N)
)


# -- 3. Finite-population region moments (the document's m_g, V_g, C_g) ------

m_g <- function(a, g) mean(a[g])
v_g <- function(a, g) { z <- a[g]; mean((z - mean(z))^2) }
c_g <- function(a, b, g) {
  x <- a[g]; y <- b[g]
  mean((x - mean(x)) * (y - mean(y)))
}


# -- 4. Validated diagnostic: six exact identities ---------------------------
#
# tau takes exactly two values by construction (R/generate_glm_dgm.R:459-462):
#   tau = k_treat * treat_effect + beta_inter * 1{x in Q}
# so it is CONSTANT on Q and on Q^c, and two-valued on S.  Each identity below
# is exact algebra, not an approximation, and must hold to numerical tolerance.

cat("--- 2. validated diagnostic (must pass before any bracket is reported) ---\n")

tau_vals <- sort(unique(round(tau, 9)))
b_int    <- dgm$model_params$beta_inter
tol      <- 1e-6

checks <- list()

checks[["tau takes exactly two values"]] <-
  c(length(tau_vals), 2)

checks[["tau on Q     == effect_Q"]] <-
  c(m_g(tau, in_Q), dgm$hazard_ratios$harm_subgroup)

checks[["tau on Qc    == effect_Qc"]] <-
  c(m_g(tau, !in_Q), dgm$hazard_ratios$no_harm_subgroup)

checks[["V_Q[tau]     == 0"]] <-
  c(v_g(tau, in_Q), 0)

checks[["V_Qc[tau]    == 0"]] <-
  c(v_g(tau, !in_Q), 0)

checks[["C_Q[mu0,tau] == 0"]] <-
  c(c_g(mu0, tau, in_Q), 0)

checks[["V_S[tau]     == b_int^2 pi (1-pi)"]] <-
  c(v_g(tau, regions$S), b_int^2 * piQ * (1 - piQ))

checks[["C_S[mu0,tau] == b_int pi (m_Q[mu0] - m_S[mu0])"]] <-
  c(c_g(mu0, tau, regions$S),
    b_int * piQ * (m_g(mu0, in_Q) - m_g(mu0, regions$S)))

ok <- TRUE
for (nm in names(checks)) {
  got <- checks[[nm]][1]; want <- checks[[nm]][2]
  pass <- abs(got - want) < max(tol, tol * abs(want))
  ok <- ok && pass
  cat(sprintf("  [%s] %-42s got %18.9f   expected %18.9f\n",
              if (pass) "PASS" else "FAIL", nm, got, want))
}
cat("\n")

if (!ok) {
  cat("DIAGNOSTIC FAILED -- bracket values are NOT reported.\n")
  cat("Either the object is not the maxeffCons MD40 fixture, or the\n")
  cat("construction of tau differs from R/generate_glm_dgm.R:459-462.\n")
  quit(status = 1L)
}


# -- 5. The bracket, per region ----------------------------------------------

cat("--- 3. Level 4 bracket, by region ---\n")
cat(sprintf("  sigma   = %.9f      sigma^2 = %.4f\n", sigma, sigma^2))
cat("\n")

out <- data.frame()
for (nm in names(regions)) {
  g   <- regions[[nm]]
  Pg  <- mean(g)
  vm  <- v_g(mu0, g)
  cm  <- c_g(mu0, tau, g)
  vt  <- v_g(tau, g)
  brk <- sigma^2 + vm + cm + 0.5 * vt

  cat(sprintf("  region %-3s  |g| = %5d   P(g) = %.6f\n", nm, sum(g), Pg))
  cat(sprintf("     m_g[mu_0]        = %16.6f\n", m_g(mu0, g)))
  cat(sprintf("     m_g[tau]         = %16.6f\n", m_g(tau, g)))
  cat(sprintf("     sigma^2          = %16.6f\n", sigma^2))
  cat(sprintf("     V_g[mu_0]        = %16.6f\n", vm))
  cat(sprintf("     C_g[mu_0, tau]   = %16.6f\n", cm))
  cat(sprintf("     (1/2) V_g[tau]   = %16.6f\n", 0.5 * vt))
  cat(sprintf("     ------------------------------------\n"))
  cat(sprintf("     bracket B(g)     = %16.6f\n", brk))
  cat(sprintf("     V_eff(g) = 4B    = %16.6f\n", 4 * brk))
  cat(sprintf("     B(g) / sigma^2   = %16.6f\n", brk / sigma^2))
  cat("\n")

  out <- rbind(out, data.frame(region = nm, Pg = Pg, brk = brk))
}


# -- 6. Comparison against the committed anchor ------------------------------

BRK_COMMITTED <- 16119          # from sigma[betahat_1000(Q)] = 13.6786
VEF_COMMITTED <- 64476
SD_COMMITTED  <- 13.6786

brkQ <- out$brk[out$region == "Q"]

cat("--- 4. comparison against the committed anchor ---\n")
cat(sprintf("  bracket(Q) direct        : %12.2f\n", brkQ))
cat(sprintf("  bracket    committed     : %12.2f\n", BRK_COMMITTED))
cat(sprintf("  ratio                    : %12.6f   (%+.2f%% in variance, %+.2f%% in SD)\n",
            brkQ / BRK_COMMITTED,
            100 * (brkQ / BRK_COMMITTED - 1),
            100 * (sqrt(brkQ / BRK_COMMITTED) - 1)))
cat("\n")
cat(sprintf("  V_eff(Q) direct          : %12.2f   (committed %d)\n",
            4 * brkQ, VEF_COMMITTED))
cat(sprintf("  sigma[betahat_1000(Q)]   : %12.6f   (committed %.4f)\n",
            sqrt(4 * brkQ / (1000 * piQ)), SD_COMMITTED))
cat(sprintf("  per-subject effective SD : %12.4f   (committed %.1f)\n",
            sqrt(brkQ), sqrt(BRK_COMMITTED)))
cat("\n")

cat("  Predicted sigma[betahat_n(Q)] at the two measured designs\n")
cat("  (idealized scale; multiply by sqrt(Jensen) for the unconditional value):\n")
for (n in c(500, 700)) {
  cat(sprintf("    n = %3d : direct %8.4f   committed-anchor %8.4f\n",
              n, sqrt(4 * brkQ / (n * piQ)), sqrt(VEF_COMMITTED / (n * piQ))))
}
cat("\n")

cat("  BRANCH READOUT\n")
cat("    ~16,777  -> the anchor was mismeasured; replace it and recompute\n")
cat("                the downstream rates.  Variance formula vindicated.\n")
cat("    ~16,119  -> the anchor is right and the difference-in-means variance\n")
cat("                formula is missing a term.  Larger finding; stop and\n")
cat("                report before any further work.\n")
cat("    neither  -> report verbatim; do not interpret.\n\n")


# -- 7. Divisor sensitivity (N vs N-1), for completeness ---------------------

cat("--- 5. divisor sensitivity (finite-population N vs sample N-1) ---\n")
for (nm in names(regions)) {
  g  <- regions[[nm]]
  ng <- sum(g)
  vm_pop  <- v_g(mu0, g)
  vm_samp <- vm_pop * ng / (ng - 1)
  cat(sprintf("  %-3s : V_g[mu_0] = %14.6f (N)  vs %14.6f (N-1)   diff %+.4f%%\n",
              nm, vm_pop, vm_samp, 100 * (vm_samp / vm_pop - 1)))
}
cat("\n")
cat("=====================================================================\n")
cat("end\n")
cat("=====================================================================\n")
