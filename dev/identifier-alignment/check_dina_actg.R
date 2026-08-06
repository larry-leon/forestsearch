# check_dina_actg.R  (revised 2026-08-05, after reading the vignette)
#
# WHICH RUN THIS APPLIES TO.  There are three distinct DINA analyses in
# analysis_actg175_binary_multimethod, and they must not be conflated:
#
#   (A) sec 7.1  fs_dina_screen -- ForestSearch with DINA supplying screening
#                cuts.  FS is the identifier.  8 non-dominated Pareto
#                candidates; selects !{wtkg >= 83.4624} & !{cd40 >= 413}, OR 6.34.
#   (B) sec 7.2  standalone dina_subgroup() -- selects cd40 >= 500, N = 100.
#   (C) sec 7.3  subgroup_method = "dina" -- DINA is the identifier.  At
#                or_threshold = 1.0 the harm floor is m_diff = 0: 10284 searched,
#                *147* qualifying, selects {preanti >= 849.4} & {cd40 >= 338},
#                n = 98, mean tau-hat = 0.0098, naive OR 1.936.
#
# CC's diagnostic was (C) run at hr.threshold = 1.25, NOT at 1.0.  That is a
# different configuration and it is the only one in which the family collapses
# to a single candidate.  The check below applies to CC's configuration alone.
#
# The observed family-size sequence, across the three floors reported:
#
#     floor = log(1.00) = 0.0000  ->  147 qualifying   (vignette sec 7.3)
#     floor = log(1.05) = 0.0488  ->   89 qualifying   (CC)
#     floor = log(1.25) = 0.2231  ->    1 qualifying   (CC)
#
# consistent with a tau-bar distribution concentrated well below 0.223 -- the
# selected candidate at floor 0 carries tau-bar = 0.0098.  "Family of one" is
# therefore a property of the threshold, not of DINA on ACTG175.

cat("=== M = 1 active-floor closed form vs CC's reported numbers ===\n")
cat("    (configuration: subgroup_method = \"dina\", hr.threshold = 1.25)\n\n")

beta_hat <- 0        # CC reported 1.7e-15, i.e. numerically zero
t_g      <- log(1.25)
p        <- 0.338    # CC's selection_rate
reported <- 0.538    # CC's selection_bias, divide-by-winners as coded

c_abs <- t_g - beta_hat
z     <- qnorm(1 - p)          # c / sigma implied by the selection rate
sigma <- c_abs / z

cf_byB <- sigma * dnorm(z)     # Eq. 9 denominator B
cf_byW <- cf_byB / p           # denominator = winners, as coded

cat(sprintf("  implied sigma_D            %.4f\n", sigma))
cat(sprintf("  implied c/sigma            %.4f\n", z))
cat(sprintf("  closed form, divide-by-B   %.4f\n", cf_byB))
cat(sprintf("  closed form, by winners    %.4f\n", cf_byW))
cat(sprintf("  CC reported                %.4f   ratio %.3f\n\n",
            reported, reported / cf_byW))

# Plausibility cross-check on the implied sigma_D, from a quantity that was
# measured rather than inferred.  The payload's DINA naive 95% CI for the
# n = 98 subgroup is (0.8374, 4.4748) on the OR scale.  That is a DIFFERENT
# subgroup from CC's, so it bounds the order of magnitude only.
ci_lo <- 0.8374208; ci_hi <- 4.4748487
sigma_payload <- (log(ci_hi) - log(ci_lo)) / (2 * qnorm(0.975))
cat(sprintf("  sigma_D from payload CI (n = 98 subgroup)  %.4f\n", sigma_payload))
cat(sprintf("  inferred sigma_D for CC's subgroup         %.4f\n", sigma))
cat("  -> same order; the inference is not obviously off-scale.\n\n")

cat(sprintf("  M = 1 IJ consequence: se_ij -> 2*sigma_D = %.4f vs robust %.4f\n",
            2 * sigma, sigma))

cat("\nCAVEATS, in order of size:\n")
cat("  1. beta_hat = 0 is CC's reported value, not one recomputed here.\n")
cat("  2. sigma_D is INFERRED from the selection rate under a Gaussian law;\n")
cat("     the package default multiplier is Poisson, so D is not exactly normal.\n")
cat("  3. The 7% gap is consistent with (2) plus Monte Carlo error at B = 5000,\n")
cat("     but this is corroboration, not verification.  Only a fixture with\n")
cat("     multiplier = \"gaussian\" and a known sigma_D settles it.\n")
