get_dgm_stratified <- function(df, knot = 5, kappa = 10, 
                                log.hrs = log(c(0.75, 0.75, 0.75)),
                                details = FALSE, strata_tte = NULL, 
                                masked = TRUE) {
  
  # === 1. Data preparation (single pass) ===
  dfa2 <- df
  if (masked) {
    dfa2 <- dfa2[dfa2$mflag == 1, ]
    dfa2$event <- dfa2$event - dfa2$en
    dfa2$tte <- dfa2$tte - dfa2$yn
  }
  
  # Create interaction terms (vectorized)
  dfa2$z.treat <- dfa2$z * dfa2$treat
  dfa2$z.k <- pmax(dfa2$z - knot, 0)  # Faster than ifelse
  dfa2$z.k.treat <- dfa2$z.k * dfa2$treat
  
  # === 2. Build formula (simplified) ===
  base_terms <- "treat + z + z.treat + z.k + z.k.treat"
  
  if (!is.null(strata_tte)) {
    weib.formula <- as.formula(
      paste0("Surv(tte, event) ~ ", base_terms, " + strata(", strata_tte, ")")
    )
  } else {
    weib.formula <- as.formula(paste0("Surv(tte, event) ~ ", base_terms))
  }
  
  # === 3. Fit models ===
  fit.weibk <- survreg(weib.formula, dist = 'weibull', data = dfa2)
  fitC.weib <- survreg(Surv(tte, 1 - event) ~ 1, dist = 'weibull', data = dfa2)
  
  # Extract parameters
  mu <- coef(fit.weibk)[1]
  gamma <- coef(fit.weibk)[-1]
  tau <- fit.weibk$scale
  muC <- coef(fitC.weib)[1]
  tauC <- fitC.weib$scale
  
  # === 4. Handle stratification (vectorized) ===
  if (!is.null(strata_tte)) {
    strataO <- dfa2[[strata_tte]]
    
    # Vectorized stratum matching
    strata_names <- sub("^strata\\(.*\\)=", "", names(tau))
    tau.strataO <- tau[match(strataO, strata_names)]
    
    if (any(is.na(tau.strataO))) {
      stop("strata_tte not uniquely identified via matching")
    }
    
    tau.approx <- median(tau.strataO)
  } else {
    strataO <- rep("All", nrow(dfa2))
    tau.strataO <- rep(tau, nrow(dfa2))
    tau.approx <- tau
  }
  
  dfa2$tau.strataO <- tau.strataO
  
  # === 5. Compute hazard ratio parameters ===
  # Unpack log hazard ratios
  loghr.0 <- log.hrs[1]
  loghr.knot <- log.hrs[2]
  loghr.kappa <- log.hrs[3]
  
  # Compute piecewise linear parameters
  b0 <- numeric(length(gamma))
  b0[1] <- loghr.0
  b0[3] <- (loghr.knot - loghr.0) / knot
  b0[5] <- (loghr.kappa - loghr.0 - kappa * b0[3]) / (kappa - knot)
  
  # Convert to AFT parameterization
  gamma.true <- -b0 * tau.approx
  
  # === 6. Sort and return ===
  # Use data.table for efficient sorting
  if (!inherits(dfa2, "data.table")) {
    dfa2 <- data.table::as.data.table(dfa2)
  }
  data.table::setorder(dfa2, z)
  
  # === 7. Optional diagnostic plot ===
  if (details) {
    z_seq <- seq(0, max(dfa2$z), by = 1)
    z_k_seq <- pmax(z_seq - knot, 0)
    loghr_seq <- b0[1] + b0[3] * z_seq + b0[5] * z_k_seq
    
    plot(z_seq, loghr_seq, type = "s", lty = 1,
         xlab = "z", ylab = "psi(z)")
    rug(dfa2$z)
    abline(h = log.hrs, col = "gray60")
    abline(v = c(0, knot, kappa), lwd = 0.5, col = "blue", lty = 2)
    abline(h = 0, lwd = 0.25, col = "red", lty = 1)
  }
  
  # === 8. Return results ===
  list(
    df_super = dfa2,
    gamma.true = gamma.true,
    mu = mu,
    tau = tau,
    muC = muC,
    tauC = tauC,
    strata_tte = strata_tte,
    tau.approx = tau.approx
  )
}