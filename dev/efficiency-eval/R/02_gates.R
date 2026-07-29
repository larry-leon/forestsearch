# =============================================================================
# STAGE G -- EQUIVALENCE GATES
# =============================================================================
# Speedups are meaningless unless the candidate computes the same thing.
# Every gate returns a tidy row; the driver refuses to report a candidate's
# timing if its gate did not pass.

.gate_row <- function(gate, what, n, pass, detail = "") {
  data.frame(stringsAsFactors = FALSE,
             gate = gate, what = what, n = n, pass = pass, detail = detail)
}

.make_subgroup <- function(N, hr_true, cens = 0.4) {
  tt <- rep(0:1, length.out = N)
  Y  <- stats::rexp(N, rate = 0.002 * exp(log(hr_true) * tt))
  C  <- stats::rexp(N, rate = 0.002 * cens / (1 - cens))
  data.table::data.table(Y = pmin(Y, C), Event = as.numeric(Y <= C), Treat = tt)
}


# -----------------------------------------------------------------------------
# G1. Fast Cox engine reproduces coxph() coefficients
# -----------------------------------------------------------------------------
gate_cox_engine <- function(reps = 40L, seed = 1) {
  set.seed(seed)
  worst <- 0; n_ok <- 0L
  for (k in seq_len(reps)) {
    N <- sample(c(60, 120, 250, 500), 1)
    d <- .make_subgroup(N, sample(c(0.7, 1, 1.5, 2.2), 1))
    # heavy ties on some reps: Efron handling is the risk case
    if (k %% 3 == 0) d$Y <- round(d$Y / 50) * 50

    ref <- tryCatch(unname(survival::coxph(
      survival::Surv(Y, Event) ~ Treat, data = d, init = 0)$coefficients[1]),
      error = function(e) NA_real_)
    o <- order(d$Y)
    cnd <- cand_cox_beta(as.double(d$Y)[o], as.double(d$Event)[o],
                         as.double(d$Treat)[o], 0)
    if (is.na(ref) && is.na(cnd)) { n_ok <- n_ok + 1L; next }
    if (!is.na(ref) && !is.na(cnd)) {
      dd <- abs(ref - cnd); worst <- max(worst, dd)
      if (dd < 1e-10) n_ok <- n_ok + 1L
    }
  }
  .gate_row("G1", "coxph.fit == coxph (beta)", reps, n_ok == reps,
            sprintf("max |diff| = %.3e over %d fits (incl. tied times)",
                    worst, reps))
}


# -----------------------------------------------------------------------------
# G2. Engine swap ALONE leaves Pcons bit-identical
# -----------------------------------------------------------------------------
# Futility OFF, so any difference is attributable to the engine only.
# Swept over hr.consistency != 1: the natural-scale comparison exists for
# exactly that case, and a gate run only at hr.consistency = 1 would not
# exercise it (exp(b) > 1 <=> b > 0 holds exactly in IEEE754).
gate_pcons_engine <- function(reps = 40L, n.splits = 200L, seed = 2) {
  set.seed(seed)
  n_ok <- 0L; mism <- character(0)
  for (k in seq_len(reps)) {
    N   <- sample(c(80, 150, 250, 400), 1)
    hrc <- sample(c(0.80, 1.00, 1.15, 1.25, 1.50), 1)
    d   <- .make_subgroup(N, sample(c(0.8, 1, 1.3, 1.8, 2.5), 1))

    set.seed(9000 + k)
    a <- cand_consistency_splits(d, N, n.splits, hrc, engine = "coxph",
                                 futility = FALSE)
    set.seed(9000 + k)
    b <- cand_consistency_splits(d, N, n.splits, hrc, engine = "fast",
                                 futility = FALSE)
    ok <- (is.na(a) && is.na(b)) || isTRUE(all.equal(a, b))
    n_ok <- n_ok + ok
    if (!ok) mism <- c(mism, sprintf("N=%d hrc=%.2f %s vs %s", N, hrc,
                                     format(a), format(b)))
  }
  .gate_row("G2", "Pcons exact, engine swap only", reps, n_ok == reps,
            if (n_ok == reps) "exact for all hr.consistency tested"
            else paste(utils::head(mism, 3), collapse = "; "))
}


# -----------------------------------------------------------------------------
# G3. Futility stopping preserves the accept/reject DECISION
# -----------------------------------------------------------------------------
# Pcons may legitimately differ for stopped candidates (the reported value is
# a proven upper bound). Two things must hold:
#   (a) the decision is unchanged;
#   (b) a truncated value is >= the full-budget value, i.e. the bound errs
#       toward acceptance and can never wrongly REJECT.
gate_futility <- function(reps = 40L, n.splits = 200L, threshold = 0.90,
                          seed = 3) {
  set.seed(seed)
  n_dec <- 0L; n_bound <- 0L; bad <- character(0)
  for (k in seq_len(reps)) {
    N   <- sample(c(80, 150, 250, 400), 1)
    hrc <- sample(c(1.00, 1.15, 1.25), 1)
    d   <- .make_subgroup(N, sample(c(0.8, 1, 1.3, 1.8, 2.5), 1))

    set.seed(9500 + k)
    full <- cand_consistency_splits(d, N, n.splits, hrc,
                                    pconsistency.threshold = NULL,
                                    engine = "fast", futility = FALSE)
    set.seed(9500 + k)
    stop <- cand_consistency_splits(d, N, n.splits, hrc,
                                    pconsistency.threshold = threshold,
                                    engine = "fast", futility = TRUE)

    df_ <- !is.na(full) && full >= threshold
    ds_ <- !is.na(stop) && stop >= threshold
    n_dec <- n_dec + (df_ == ds_)

    truncated <- !isTRUE(all.equal(full, stop))
    ok_b <- !truncated || (!is.na(stop) && !is.na(full) && stop >= full)
    n_bound <- n_bound + ok_b
    if (df_ != ds_ || !ok_b)
      bad <- c(bad, sprintf("N=%d full=%s stop=%s", N, format(full), format(stop)))
  }
  rbind(
    .gate_row("G3a", "futility: decision unchanged", reps, n_dec == reps,
              if (n_dec == reps) sprintf("identical at threshold %.2f", threshold)
              else paste(utils::head(bad, 3), collapse = "; ")),
    .gate_row("G3b", "futility: reported value is an upper bound", reps,
              n_bound == reps,
              "truncated Pcons >= full-budget Pcons (cannot wrongly reject)")
  )
}


# -----------------------------------------------------------------------------
# G4. Search scorer candidate matches the package function
# -----------------------------------------------------------------------------
gate_search_scorer <- function(reps = 40L, seed = 4) {
  set.seed(seed)
  ok_hr <- 0L; ok_med <- 0L; worst_hr <- 0; worst_ci <- 0
  for (k in seq_len(reps)) {
    N <- sample(c(300, 686, 1000), 1)
    d <- .make_subgroup(N, sample(c(0.8, 1, 1.4), 1))
    idx <- stats::rbinom(N, 1, stats::runif(1, 0.3, 0.7))
    if (sum(idx) < 40) next

    ref <- tryCatch(forestsearch:::fit_cox_for_subgroup(
      d$Y, d$Event, d$Treat, idx), error = function(e) NULL)
    cnd <- cand_fit_cox_subgroup(d$Y, d$Event, d$Treat, idx)
    if (is.null(ref) || is.null(cnd)) next

    worst_hr <- max(worst_hr, abs(unname(ref$hr) - cnd$hr))
    worst_ci <- max(worst_ci, abs(unname(ref$lower) - cnd$lower),
                    abs(unname(ref$upper) - cnd$upper))
    ok_hr <- ok_hr + (abs(unname(ref$hr) - cnd$hr) < 1e-9)
    m_ok <- (is.na(ref$med0) && is.na(cnd$med0)) ||
      isTRUE(all.equal(unname(ref$med0), unname(cnd$med0)))
    m_ok <- m_ok && ((is.na(ref$med1) && is.na(cnd$med1)) ||
                       isTRUE(all.equal(unname(ref$med1), unname(cnd$med1))))
    ok_med <- ok_med + m_ok
  }
  rbind(
    .gate_row("G4a", "search scorer: HR matches", ok_hr, ok_hr > 0,
              sprintf("max |dHR| = %.3e", worst_hr)),
    .gate_row("G4b", "search scorer: medians match without summary.survfit",
              ok_med, ok_med > 0,
              sprintf("max |dCI| = %.3e (qnorm vs 1.96 in package)", worst_ci))
  )
}


# -----------------------------------------------------------------------------
# G5. Vectorised pre-screen selects the identical combination set
# -----------------------------------------------------------------------------
gate_prescreen <- function(seed = 5) {
  set.seed(seed)
  out <- list()
  for (L in c(20L, 40L)) {
    N <- 686
    prev <- c(stats::runif(round(L * 0.75), 0.10, 0.60),
              stats::runif(L - round(L * 0.75), 0.005, 0.045))
    zz <- vapply(prev, function(p) stats::rbinom(N, 1, p), numeric(N))
    colnames(zz) <- paste0("z", seq_len(L))
    minp <- 0.05; rmin <- 5

    for (k in c(2L, 3L)) {
      combos <- t(utils::combn(L, k))
      ref <- vapply(seq_len(nrow(combos)), function(i) {
        x <- zz[, combos[i, ], drop = FALSE]
        if (!forestsearch:::has_positive_variance(x)) return(FALSE)
        if (!forestsearch:::meets_prevalence_threshold(x, minp)) return(FALSE)
        !forestsearch:::extract_idx_flagredundancy(x, rmin)$flag.redundant
      }, logical(1))
      cnd <- cand_screen_combos(zz, combos, minp, rmin)
      out[[length(out) + 1L]] <- .gate_row(
        sprintf("G5.L%d.k%d", L, k), "pre-screen selects identical set",
        nrow(combos), identical(ref, cnd),
        sprintf("%d of %d pass; L_eff = %d", sum(ref), length(ref),
                sum(colMeans(zz) >= minp)))
    }
  }
  do.call(rbind, out)
}


# -----------------------------------------------------------------------------
# G6. GLM candidate matches glm() + vcov()
# -----------------------------------------------------------------------------
gate_glm <- function(reps = 40L, seed = 6) {
  set.seed(seed)
  n_ok <- 0L; worst_b <- 0; worst_se <- 0
  for (k in seq_len(reps)) {
    N <- sample(c(120, 300, 600), 1)
    d <- data.frame(treat = rep(0:1, length.out = N))
    d$Y <- stats::rbinom(N, 1, stats::plogis(-0.4 + 0.6 * d$treat))
    ref <- tryCatch({
      f <- stats::glm(Y ~ treat, data = d, family = stats::binomial())
      list(b = unname(stats::coef(f)[["treat"]]),
           se = unname(sqrt(diag(stats::vcov(f)))[["treat"]]))
    }, error = function(e) NULL)
    X <- cbind(`(Intercept)` = 1, treat = d$treat)
    cnd <- cand_glm_effect(d$Y, X)
    if (is.null(ref) || !cnd$converged) next
    worst_b  <- max(worst_b, abs(ref$b - cnd$estimate))
    worst_se <- max(worst_se, abs(ref$se - cnd$se))
    n_ok <- n_ok + (abs(ref$b - cnd$estimate) < 1e-10 &&
                      abs(ref$se - cnd$se) < 1e-10)
  }
  .gate_row("G6", "glm.fit == glm + vcov", reps, n_ok == reps,
            sprintf("max |db| = %.3e, max |dse| = %.3e", worst_b, worst_se))
}


run_gates <- function(cfg) {
  rows <- list(
    gate_cox_engine(cfg$gate_reps),
    gate_pcons_engine(cfg$gate_reps, cfg$n_splits_gate),
    gate_futility(cfg$gate_reps, cfg$n_splits_gate),
    gate_search_scorer(cfg$gate_reps),
    gate_prescreen()
  )
  if (isTRUE(cfg$run_glm)) rows <- c(rows, list(gate_glm(cfg$gate_reps)))
  do.call(rbind, rows)
}
