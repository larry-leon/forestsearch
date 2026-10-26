# b2_lean.R -- the candidates sized and their bit-identity tested, per fixture,
# with NO change to the package.  One traced forestsearch() call per fixture
# records (in this session only, untraced before any timing):
#   * the consistency-stage evaluator's arguments
#     (evaluate_consistency_twostage / evaluate_subgroup_consistency): the
#     shared objects once, the candidate index m per call;
#   * the consistency "pieces" inputs per evaluated candidate
#     (.consistency_glm_pieces / .consistency_cox_pieces): the needed columns
#     of the subgroup frame, and cox_init for survival;
#   * the search-stage fit inputs per candidate (binary:
#     fit_glm_for_subgroup; continuous: fit_glm_for_subgroup_fast).
# Then, per candidate, the package's own function (unmodified) is run on the
# recorded arguments and timed against a lean equivalent, and the outputs are
# compared with identical().
#
# Candidates measured here (letters as in the survey / Stage 0):
#   B  survival consistency pieces: coxph.fit + residuals.coxph on a constructed
#      object instead of two formula coxph() fits + a model-frame re-evaluation
#   E  continuous consistency pieces: lm.fit + QR hat values instead of
#      lm() + lm.influence()
#   F  binary search-stage fit: glm.fit on vectors instead of glm() on a
#      data-frame slice, SE from the fit's QR instead of vcov()/summary.glm()
#   G  binary consistency pieces: glm.fit + QR hat values instead of
#      glm() + lm.influence()
#   H  consistency-stage member extraction: a logical index over the cut
#      columns instead of as.data.table(df) + sequential row filters, per
#      candidate; measured end to end as the whole evaluator (package) vs a
#      lean evaluator (H + the lean pieces of B/E/G).
# A (survival search fit) is Stage 0's and is quoted, not re-measured.

source("stage0b_common.R")

# ------------------------------------------------------------------------------
# lean pieces (vectors in; beta_hat, sigma_D out; NULL where the package's
# function returns NULL).  The dfbeta arithmetic is .dfbeta_glm()'s, line for
# line, with lm.influence()'s hat values recomputed from the fit's own QR.
# ------------------------------------------------------------------------------
hat_from_qr <- function(qr, n, k) {
  h <- rowSums(base::qr.qy(qr, diag(1, n, k))^2)
  h[h > 1 - 10 * .Machine$double.eps] <- 1
  h
}
lean_glm_pieces <- function(y, z, effect_measure, adverse_outcome, outcome_type, treat.name) {
  if (!adverse_outcome) { if (outcome_type == "continuous") y <- -y else if (outcome_type == "binary") y <- 1L - y }
  n0 <- sum(z == 0L, na.rm = TRUE); n1 <- sum(z == 1L, na.rm = TRUE)
  if (n0 < 3L || n1 < 3L) return(NULL)
  X <- cbind("(Intercept)" = 1, as.double(z)); colnames(X)[2L] <- treat.name
  if (effect_measure == "MD") {
    fit <- stats::lm.fit(X, as.double(y)); w <- rep(1, nrow(X))
  } else if (effect_measure == "OR") {
    fit <- tryCatch(suppressWarnings(stats::glm.fit(X, as.double(y), family = stats::binomial("logit"))),
                    error = function(e) NULL)
    if (is.null(fit) || isFALSE(fit$converged)) return(NULL)
    w <- fit$weights
  } else return(NULL)
  cf <- fit$coefficients
  if (any(is.na(cf))) return(NULL)
  beta_hat <- unname(cf[[treat.name]])
  h <- hat_from_qr(fit$qr, nrow(X), fit$rank)
  A <- solve(crossprod(X, w * X))
  r <- sqrt(w) * fit$residuals
  dfb <- t(A %*% t(X * (sqrt(w) * r / (1 - h))))
  dfb <- as.numeric(dfb[, 2L])
  sigma_D <- sqrt(sum(dfb^2))
  if (!is.finite(beta_hat) || !is.finite(sigma_D) || sigma_D <= 0) return(NULL)
  list(beta_hat = beta_hat, sigma_D = sigma_D)
}

# survival: the conventions coxph() passes to coxph.fit() are read at run time
# (as in Stage 0); `init` is taken from the captured second fit (the warm start).
.coxfit_conv <- NULL
.cox_terms <- stats::terms(survival::Surv(Y, E) ~ Treat, specials = c("strata", "cluster", "tt"))
lean_cox_fit <- function(Y, E, Treat, init) {
  a <- .coxfit_conv
  n <- length(Y)
  ys <- survival::Surv(Y, E)
  if (isTRUE(a$control$timefix)) ys <- survival::aeqSurv(ys)
  args <- c(list(x = matrix(as.double(Treat), ncol = 1L, dimnames = list(NULL, "Treat")), y = ys,
                 offset = if (is.null(a$offset)) NULL else rep(0, n),
                 rownames = if (is.null(a$rownames)) NULL else as.character(seq_len(n)),
                 init = init),
            a[setdiff(names(a), c("x", "y", "offset", "rownames", "init"))])
  do.call(survival::coxph.fit, args)
}
lean_cox_pieces <- function(Y, E, Treat) {
  n <- length(Y)
  if (n < 2L || sum(E) < 2L) return(NULL)
  fit0 <- tryCatch(suppressWarnings(lean_cox_fit(Y, E, Treat, NULL)), error = function(e) NULL)
  cox_init <- if (is.null(fit0)) 0 else fit0$coefficients[1]
  fit <- tryCatch(suppressWarnings(lean_cox_fit(Y, E, Treat, cox_init)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- fit$coefficients
  if (length(cf) == 0L || any(is.na(cf))) return(NULL)
  beta_hat <- unname(cf[[1L]])
  obj <- fit
  obj$x <- matrix(as.double(Treat), ncol = 1L, dimnames = list(NULL, "Treat"))
  ys <- survival::Surv(Y, E); if (isTRUE(.coxfit_conv$control$timefix)) ys <- survival::aeqSurv(ys)
  obj$y <- ys; obj$strata <- NULL; obj$method <- .coxfit_conv$method; obj$weights <- NULL
  obj$terms <- .cox_terms      # residuals.coxph() requires a terms component
  class(obj) <- "coxph"
  dfbeta <- tryCatch(stats::residuals(obj, type = "dfbeta"), error = function(e) NULL)
  if (is.null(dfbeta)) return(NULL)
  dfbeta <- as.numeric(dfbeta)
  sigma_D <- sqrt(sum(dfbeta^2))
  if (!is.finite(beta_hat) || !is.finite(sigma_D) || sigma_D <= 0) return(NULL)
  list(beta_hat = beta_hat, sigma_D = sigma_D, cox_init = cox_init)
}

# ------------------------------------------------------------------------------
# per fixture
# ------------------------------------------------------------------------------
summ <- list()
fixture_sel <- strsplit(Sys.getenv("STAGE0B_FIXTURES", "survival,continuous,binary"), ",")[[1L]]
for (tag in intersect(names(FIXTURES), fixture_sel)) {
  cat("\n==========", tag, "==========\n")
  fx <- FIXTURES[[tag]]()
  is_surv <- tag == "survival"
  out_nm <- fx$args$outcome.name; trt_nm <- fx$args$treat.name

  # ---- 1. traced call: capture -------------------------------------------------
  .stage <- new.env(); .stage$args <- NULL; .stage$m <- integer(0); .stage$fun <- NA_character_
  assign(".stage", .stage, envir = globalenv())
  stage_tracer <- function(fname) bquote({
    .fr <- environment()
    if (is.null(.GlobalEnv$.stage$args)) {
      .nm <- setdiff(names(formals(get(.(fname), envir = asNamespace("forestsearch")))), c("m", "..."))
      .GlobalEnv$.stage$args <- lapply(setNames(.nm, .nm), function(.a) get(.a, envir = .fr, inherits = FALSE))
      .GlobalEnv$.stage$fun <- .(fname)
    }
    .GlobalEnv$.stage$m <- c(.GlobalEnv$.stage$m, m)
  })
  invisible(capture.output({
    trace("evaluate_consistency_twostage", where = asNamespace("forestsearch"), print = FALSE,
          tracer = stage_tracer("evaluate_consistency_twostage"))
    trace("evaluate_subgroup_consistency", where = asNamespace("forestsearch"), print = FALSE,
          tracer = stage_tracer("evaluate_subgroup_consistency"))
  }))
  .pieces <- new_cap(); assign(".pieces", .pieces, envir = globalenv())
  if (is_surv) {
    invisible(capture.output(trace(".consistency_cox_pieces", where = asNamespace("forestsearch"), print = FALSE,
      tracer = quote({ .d <- as.data.frame(df)
        .GlobalEnv$.pieces$push(list(Y = .d[[tte.name]], E = .d[[event.name]], Treat = .d[[treat.name]],
                                     cox_init = cox_init, adj = adjust_covariates, ids = .d[["id"]])) }))))
  } else {
    invisible(capture.output(trace(".consistency_glm_pieces", where = asNamespace("forestsearch"), print = FALSE,
      tracer = quote({ .d <- as.data.frame(df)
        .GlobalEnv$.pieces$push(list(y = .d[[outcome.name]], z = .d[[treat.name]], ids = .d[["id"]],
                                     offset = if (!is.null(offset.name)) .d[[offset.name]] else NULL))
        if (is.null(.GlobalEnv$.pieces$conv))
          .GlobalEnv$.pieces$conv <- list(outcome_type = outcome_type, effect_measure = effect_measure,
            treat.name = treat.name, outcome.name = outcome.name, offset.name = offset.name,
            adjust_covariates = adjust_covariates, adverse_outcome = adverse_outcome) }))))
  }
  .sfit <- new_cap(); assign(".sfit", .sfit, envir = globalenv())
  if (tag == "binary") {
    invisible(capture.output(trace("fit_glm_for_subgroup", where = asNamespace("forestsearch"), print = FALSE,
      tracer = quote({ .GlobalEnv$.sfit$push(list(id.x = id.x))
        if (is.null(.GlobalEnv$.sfit$df_clean)) { .GlobalEnv$.sfit$df_clean <- df_clean; .GlobalEnv$.sfit$estimator_fn <- estimator_fn } }))))
  } else if (tag == "continuous") {
    invisible(capture.output(trace("fit_glm_for_subgroup_fast", where = asNamespace("forestsearch"), print = FALSE,
      tracer = quote({ .GlobalEnv$.sfit$push(list(id.x = id.x))
        if (is.null(.GlobalEnv$.sfit$yy)) { .GlobalEnv$.sfit$yy <- yy; .GlobalEnv$.sfit$tt <- tt; .GlobalEnv$.sfit$estimator_fn <- estimator_fn } }))))
  }
  t_traced <- system.time(fs <- run_fixture_call(fx))[["elapsed"]]
  invisible(capture.output({
    untrace("evaluate_consistency_twostage", where = asNamespace("forestsearch"))
    untrace("evaluate_subgroup_consistency", where = asNamespace("forestsearch"))
    if (is_surv) untrace(".consistency_cox_pieces", where = asNamespace("forestsearch"))
    else untrace(".consistency_glm_pieces", where = asNamespace("forestsearch"))
    if (tag == "binary") untrace("fit_glm_for_subgroup", where = asNamespace("forestsearch"))
    if (tag == "continuous") untrace("fit_glm_for_subgroup_fast", where = asNamespace("forestsearch"))
  }))
  cc <- candidate_counts(fs)
  P <- cap_items(.pieces); S <- cap_items(.sfit); M <- .stage$m; SA <- .stage$args
  kv_tag <- function(...) { v <- list(...); do.call(kv, setNames(v, paste0(tag, "_", names(v)))) }
  kv_tag(capture_call_s = t_traced, capture_stage_fun = .stage$fun,
         capture_n_stage_calls = length(M), capture_n_pieces_calls = length(P),
         capture_n_searchfit_calls = length(S),
         capture_stage_calls_match_counts = identical(length(M), as.integer(cc$n_consistency_evaluated)))
  if (length(M) == 0L || is.null(SA)) {
    cat("GATE FAILED: no consistency-stage call captured for", tag, "\n"); quit(status = 4)
  }
  if (is_surv) {
    # the conventions coxph() passes to coxph.fit(): captured during ONE package pieces call
    .capfit <- new.env(); .capfit$calls <- list(); assign(".capfit", .capfit, envir = globalenv())
    invisible(capture.output(trace("coxph.fit", where = asNamespace("survival"), print = FALSE, tracer = quote({
      .fr <- environment(); .nm <- setdiff(names(formals(survival::coxph.fit)), "...")
      .GlobalEnv$.capfit$calls[[length(.GlobalEnv$.capfit$calls) + 1L]] <-
        lapply(setNames(.nm, .nm), function(.a) tryCatch(get(.a, envir = .fr, inherits = FALSE), error = function(e) NULL)) }))))
    p1 <- P[[1L]]
    d1 <- data.frame(Y = p1$Y, Event = p1$E, Treat = p1$Treat)
    invisible(forestsearch:::.consistency_cox_pieces(d1, "Y", "Event", "Treat", cox_init = p1$cox_init, adjust_covariates = NULL))
    invisible(capture.output(untrace("coxph.fit", where = asNamespace("survival"))))
    cf <- .capfit$calls
    kv_tag(coxfit_calls_per_pieces = length(cf),
           coxfit_init_as_passed = if (length(cf)) paste(deparse(cf[[length(cf)]]$init), collapse = "") else "<none>",
           coxfit_init_equals_cox_init = if (length(cf)) isTRUE(all.equal(unname(cf[[length(cf)]]$init), unname(p1$cox_init))) else NA,
           coxfit_nocenter = if (length(cf)) paste(deparse(cf[[1]]$nocenter), collapse = "") else "<none>",
           coxfit_timefix = isTRUE(cf[[1]]$control$timefix))
    .coxfit_conv <<- cf[[length(cf)]]
  }

  # ---- 2. search-stage fit (binary: F; continuous: bucket only) --------------
  if (tag == "binary") {
    dfc <- .sfit$df_clean; efn <- .sfit$estimator_fn
    yy <- dfc[[out_nm]]; tt <- dfc[[trt_nm]]
    fn_env <- environment(efn)
    adv <- isTRUE(fn_env$adverse_outcome)
    pkg_fit <- forestsearch:::fit_glm_for_subgroup
    lean_fit <- function(id.x) {
      idx <- id.x == 1
      if (sum(idx) < 6L) return(NULL)
      y <- yy[idx]; z <- tt[idx]
      if (!adv) y <- 1L - y
      n0 <- sum(z == 0L, na.rm = TRUE); n1 <- sum(z == 1L, na.rm = TRUE)
      why <- forestsearch:::.fs_existence_reason("OR", forestsearch:::.fs_binary_cells(y, z))
      if (!is.null(why)) return(list(nonestimable = TRUE, reason = why))
      X <- cbind("(Intercept)" = 1, as.double(z)); colnames(X)[2L] <- trt_nm
      fit <- tryCatch(suppressWarnings(stats::glm.fit(X, as.double(y), family = stats::binomial("logit"))),
                      error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      est <- fit$coefficients[[trt_nm]]
      if (isFALSE(fit$converged) && !is.na(est)) return(list(nonestimable = TRUE, reason = "non-convergent fit"))
      if (is.na(est)) return(NULL)
      p1 <- seq_len(fit$rank)
      se <- unname(sqrt(diag(1 * chol2inv(fit$qr$qr[p1, p1, drop = FALSE])))[2L])
      z_crit <- 1.96
      list(hr = est, lower = est - z_crit * se, upper = est + z_crit * se, med0 = NA_real_, med1 = NA_real_)
    }
    res_f <- lapply(S, function(s) pkg_fit(dfc, s$id.x, efn))
    res_l <- lapply(S, function(s) lean_fit(s$id.x))
    same <- mapply(function(u, v) identical(u, v), res_f, res_l)
    kind <- function(u) if (is.null(u)) "NULL" else if (isTRUE(u$nonestimable)) "nonestimable" else "fit"
    kf <- vapply(res_f, kind, ""); kl <- vapply(res_l, kind, "")
    tp <- time_pair(function() for (s in S) pkg_fit(dfc, s$id.x, efn), function() for (s in S) lean_fit(s$id.x))
    kv_tag(F_n = length(S), F_identical = sum(same), F_kinds_agree = identical(kf, kl),
           F_n_fit = sum(kf == "fit"), F_n_nonestimable = sum(kf == "nonestimable"), F_n_null = sum(kf == "NULL"),
           F_formula_s_reps = tp$a, F_lean_s_reps = tp$b, F_formula_s = tp$a_med, F_lean_s = tp$b_med,
           F_machinery_s = tp$a_med - tp$b_med, F_machinery_share = (tp$a_med - tp$b_med) / tp$a_med)
    if (any(!same)) { k <- which(!same)[1]; cat("first F mismatch at", k, "\n"); str(res_f[[k]]); str(res_l[[k]]) }
  } else if (tag == "continuous") {
    yy <- .sfit$yy; tt <- .sfit$tt; efn <- .sfit$estimator_fn
    t_fast <- vapply(1:3, function(i) system.time(for (s in S) forestsearch:::fit_glm_for_subgroup_fast(yy, tt, s$id.x, efn))[["elapsed"]], numeric(1))
    kv_tag(searchfit_fast_n = length(S), searchfit_fast_s = median(t_fast))
  }

  # ---- 3. pieces alone (B / E / G) -------------------------------------------
  if (is_surv) {
    pkg_pieces <- function(p) forestsearch:::.consistency_cox_pieces(
      data.frame(Y = p$Y, Event = p$E, Treat = p$Treat), "Y", "Event", "Treat",
      cox_init = p$cox_init, adjust_covariates = p$adj)
    lean_pieces <- function(p) lean_cox_pieces(p$Y, p$E, p$Treat)
    c_cmp <- log(SA$hr.consistency)
  } else {
    cv <- .pieces$conv
    pkg_pieces <- function(p) {
      d <- data.frame(y = p$y, z = p$z); names(d) <- c(cv$outcome.name, cv$treat.name)
      if (!is.null(p$offset)) d[[cv$offset.name]] <- p$offset
      forestsearch:::.consistency_glm_pieces(d, cv$outcome_type, cv$effect_measure, cv$treat.name,
                                             cv$outcome.name, cv$offset.name, cv$adjust_covariates, cv$adverse_outcome)
    }
    lean_pieces <- function(p) lean_glm_pieces(p$y, p$z, cv$effect_measure, cv$adverse_outcome, cv$outcome_type, cv$treat.name)
    c_cmp <- SA$glm_resample_spec$comparison_threshold
  }
  res_f <- lapply(P, pkg_pieces); res_l <- lapply(P, lean_pieces)
  nf <- vapply(res_f, is.null, logical(1)); nl <- vapply(res_l, is.null, logical(1))
  id_beta  <- mapply(function(u, v) !is.null(u) && !is.null(v) && identical(u$beta_hat, v$beta_hat), res_f, res_l)
  id_sigma <- mapply(function(u, v) !is.null(u) && !is.null(v) && identical(u$sigma_D, v$sigma_D), res_f, res_l)
  rate_of <- function(r) if (is.null(r)) NA_real_ else max(0, 2 * stats::pnorm((r$beta_hat - c_cmp) / r$sigma_D) - 1)
  id_rate <- mapply(function(u, v) identical(rate_of(u), rate_of(v)), res_f, res_l)
  tp <- time_pair(function() for (p in P) pkg_pieces(p), function() for (p in P) lean_pieces(p))
  L <- if (is_surv) "B" else if (tag == "continuous") "E" else "G"
  kv_tag(pieces_letter = L, pieces_n = length(P), pieces_null_formula = sum(nf), pieces_null_lean = sum(nl),
         pieces_null_sets_agree = identical(nf, nl), pieces_beta_identical = sum(id_beta),
         pieces_sigma_identical = sum(id_sigma), pieces_rate_identical = sum(id_rate),
         pieces_formula_s_reps = tp$a, pieces_lean_s_reps = tp$b,
         pieces_formula_s = tp$a_med, pieces_lean_s = tp$b_med,
         pieces_machinery_s = tp$a_med - tp$b_med, pieces_machinery_share = (tp$a_med - tp$b_med) / tp$a_med,
         pieces_lean_identical_all = all(id_rate))
  if (any(!id_rate)) { k <- which(!id_rate)[1]; cat("first pieces mismatch at", k, ":\n"); str(res_f[[k]]); str(res_l[[k]]) }
  pieces_identical_all <- all(id_rate)

  # ---- 4. the whole consistency evaluator: package vs lean (H + pieces) -------
  stage_fun <- get(.stage$fun, envir = asNamespace("forestsearch"))
  pkg_stage <- function(m) do.call(stage_fun, c(list(m = m), SA))
  df_all <- as.data.frame(SA$df)
  index_mat <- as.matrix(SA$index.Z); storage.mode(index_mat) <- "double"
  names.Z <- SA$names.Z
  col_cache <- lapply(setNames(names.Z, names.Z), function(v) df_all[[v]] == 1)
  pdig <- SA$pconsistency.digits; pthr <- SA$pconsistency.threshold
  if (is_surv) { Y_all <- df_all$Y; E_all <- df_all$Event; T_all <- df_all$Treat } else {
    y_all <- df_all[[cv$outcome.name]]; z_all <- df_all[[cv$treat.name]]
    off_all <- if (!is.null(cv$offset.name)) df_all[[cv$offset.name]] else NULL }
  lean_stage <- function(m, use_lean_pieces = TRUE) {
    this.m <- names.Z[index_mat[m, ] == 1]
    keep <- col_cache[[this.m[1L]]]
    for (v in this.m[-1L]) keep <- keep & col_cache[[v]]
    N.x <- sum(keep)
    if (N.x < 10) return(NULL)
    rr <- if (is_surv) {
      if (use_lean_pieces) lean_cox_pieces(Y_all[keep], E_all[keep], T_all[keep])
      else pkg_pieces(list(Y = Y_all[keep], E = E_all[keep], Treat = T_all[keep],
                           cox_init = lean_cox_fit(Y_all[keep], E_all[keep], T_all[keep], NULL)$coefficients[1], adj = NULL))
    } else {
      if (use_lean_pieces) lean_glm_pieces(y_all[keep], z_all[keep], cv$effect_measure, cv$adverse_outcome, cv$outcome_type, cv$treat.name)
      else pkg_pieces(list(y = y_all[keep], z = z_all[keep], offset = if (is.null(off_all)) NULL else off_all[keep]))
    }
    if (is.null(rr)) return(structure(list(), class = "stage0b_fallback"))  # package would fall through to splitting
    p <- round(rate_of(rr), pdig)
    if (p < pthr) return(NULL)
    p
  }
  res_pkg <- lapply(M, pkg_stage)
  res_lean <- lapply(M, lean_stage)
  pcons_pkg <- vapply(res_pkg, function(r) if (is.null(r)) NA_real_ else as.numeric(r[["Pcons"]]), numeric(1))
  pcons_lean <- vapply(res_lean, function(r) if (is.null(r) || inherits(r, "stage0b_fallback")) NA_real_ else r, numeric(1))
  n_fallback <- sum(vapply(res_lean, inherits, logical(1), "stage0b_fallback"))
  stage_ident <- identical(pcons_pkg, pcons_lean)
  tp <- time_pair(function() for (m in M) pkg_stage(m), function() for (m in M) lean_stage(m))
  # variant: lean extraction with the package's own pieces (what H alone gives)
  tH <- vapply(1:3, function(i) system.time(for (m in M) lean_stage(m, use_lean_pieces = FALSE))[["elapsed"]], numeric(1))
  res_H <- lapply(M, function(m) lean_stage(m, use_lean_pieces = FALSE))
  pcons_H <- vapply(res_H, function(r) if (is.null(r) || inherits(r, "stage0b_fallback")) NA_real_ else r, numeric(1))
  kv_tag(stage_n = length(M), stage_n_passed_pkg = sum(!is.na(pcons_pkg)), stage_n_passed_lean = sum(!is.na(pcons_lean)),
         stage_pcons_identical_all = stage_ident, stage_n_lean_fallback = n_fallback,
         stage_pcons_identical_H_only = identical(pcons_pkg, pcons_H),
         stage_formula_s_reps = tp$a, stage_lean_s_reps = tp$b, stage_formula_s = tp$a_med, stage_lean_s = tp$b_med,
         stage_lean_Honly_s = median(tH),
         stage_machinery_s = tp$a_med - tp$b_med, stage_machinery_share = (tp$a_med - tp$b_med) / tp$a_med,
         stage_H_only_machinery_s = tp$a_med - median(tH))
  # extraction overhead = whole package evaluator minus the package pieces alone (both on the same candidates)
  kv_tag(stage_extraction_overhead_s = tp$a_med - median(vapply(1:3, function(i) system.time(for (p in P) pkg_pieces(p))[["elapsed"]], numeric(1))))
  if (!stage_ident) { k <- which(!(is.na(pcons_pkg) == is.na(pcons_lean)) | (pcons_pkg != pcons_lean))[1]; cat("first stage mismatch at m =", M[k], ": pkg", pcons_pkg[k], "lean", pcons_lean[k], "\n") }

  summ[[tag]] <- list(cc = cc, P_n = length(P), S_n = length(S), M_n = length(M),
                      pieces_identical_all = pieces_identical_all, stage_identical = stage_ident)
  save_result(list(cc = cc, M = M, pcons_pkg = pcons_pkg, pcons_lean = pcons_lean,
                   stage_args_names = names(SA), coxfit_conv = if (is_surv) .coxfit_conv else NULL), paste0("b2_", tag))
}
save_result(summ, "b2_summary")
cat("\nb2 complete\n")
