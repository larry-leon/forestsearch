# =============================================================================
# fs_oc_family_enumerate() -- the population enumeration of the candidate
#                             family under the forestsearch() cut specification
# -----------------------------------------------------------------------------
# Built for the operating-characteristics wrapper fs_oc_predict().  The family
# is enumerated on the DGM's population frame (dgm$df_super), not on a sample:
# cuts land at POPULATION quantiles, every membership proportion is a
# population proportion, and the enumeration is deterministic in (dgm, args, n).
#
# The package's own machinery is reused end to end:
#   - cuts:          get_FSdata() (cut_var_jq / get_conf_force / collapse)
#   - directions:    dummy() -> both indicator columns of every 2-level factor
#   - combinations:  generate_combination_indices() / get_covs_in() /
#                    get_subgroup_membership(), exactly as forestsearch_main.R
#                    composes the MR family
#   - floors:        meets_prevalence_threshold() (minp, per constituent
#                    factor), the rmin redundancy rule, and the n.min size
#                    floor expressed as Pg >= n.min / n
# The GRF pre-screen, the effect screen, the consistency screen and the
# statistics-keyed near-duplicate removal are deliberately NOT applied: those
# are the gate, which fs_oc_predict() evaluates on the family.  Candidates with
# identical population membership are collapsed to one.
# =============================================================================


#' Enumerate the candidate family on a DGM's population frame
#'
#' Builds the family of candidate subgroups that \code{\link{forestsearch}}
#' would search over, but on the DGM's super-population frame
#' \code{dgm$df_super} rather than on a sampled trial, so that every cut lands
#' at a population quantile and every prevalence, purity and overlap is a
#' population proportion.  The result carries the nine quantities the
#' operating-characteristics prediction of \code{\link{fs_oc_predict}} consumes.
#'
#' @details
#' \strong{Cuts.} The cut-related entries of \code{forestsearch_args}
#' (\code{confounders.name}, \code{conf.cont_jcuts}, \code{cut_type},
#' \code{cont.cutoff}, \code{conf.cont_medians}, \code{conf.cont_medians_force},
#' \code{conf_force}, \code{defaultcut_names}, \code{exclude_cuts},
#' \code{collapse_cuts}, \code{collapse_cuts_args}) are handed to the package's
#' own \code{get_FSdata()} on \code{df_super}.  Both directions of every cut are
#' generated, as the search does.  LASSO, GRF and DINA cut sources are off: they
#' are sample-fitted screens, not part of the cut specification.
#'
#' \strong{Combinations.} All combinations of up to \code{maxk} indicator
#' columns are enumerated with \code{generate_combination_indices()}, in the
#' order \code{forestsearch()} uses when it composes the MR family.
#'
#' \strong{Structural floors} (all on the population frame, in the order of
#' \code{evaluate_combination_with_status()}):
#' \enumerate{
#'   \item combinations with a constant column or an empty pairwise
#'     intersection are skipped (the search's status 0);
#'   \item \code{minp}: every constituent factor must have population
#'     prevalence \code{>= minp};
#'   \item \code{rmin}: each added factor must shrink the membership by more
#'     than \code{rmin} subjects \emph{of a trial of size n}, i.e. by more than
#'     \code{rmin / n} in population proportion;
#'   \item size: \code{Pg >= n.min / n}, where \code{n.min} is resolved as
#'     \code{forestsearch()} resolves it -- the supplied value (default 60), or
#'     \code{max(60, ceiling(n.min.frac * n))} when \code{n.min = NULL}.
#' }
#' The GRF variable-importance pre-screen, the effect screen, the consistency
#' screen and the near-duplicate removal are \emph{not} applied.  Candidates
#' whose population membership vectors are identical are collapsed to the
#' first one enumerated.
#'
#' \strong{Null DGMs.} A DGM whose \code{df_super$flag_harm} has no member
#' (Q empty; \code{generate_glm_dgm(model = "null")}) is detected
#' structurally, cross-checked against \code{dgm$model} when present.  Under
#' the null every candidate has the same true effect, so the fields become:
#' \code{beta_g} = the common effect (the DGM's \code{effect_Qc} =
#' \code{effect_ITT}, oriented as the alternative path orients);
#' \code{se_g} from the whole-population effective variance
#' (\code{fs_dgm_scale(dgm, regions = list(S = ...))}, the S row) at
#' \code{(n, Pg)} with the same prevalence scaling; \code{PQg = 0};
#' \code{sens_g = NA} (0/0, undefined -- not zero); \code{spec_g = 1 - Pg};
#' \code{PQ = 0} (from which NPV = 1 follows downstream).  Enumeration, floors,
#' \code{ovl} and the covariance are unchanged.  The element \code{null}
#' records the branch taken.
#'
#' \strong{Scale.} Every mean and standard error is derived from
#' \code{\link{fs_dgm_scale}(dgm)}.  The orientation is the harm direction
#' \code{s = sign(m_tau[Q])}, so the planted effect is oriented positive:
#' with \code{Q} the true harm region, \code{tauQc = s * m_tau[Qc]},
#' \code{bint = s * (m_tau[Q] - m_tau[Qc])},
#' \code{seQ1000 = sqrt(V_eff[Q] / (1000 * P(Q)))}, and for each candidate
#' \code{beta_g = tauQc + bint * PQg} -- the signed mixture
#' \code{s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)} -- and
#' \code{se_g = seQ1000 * sqrt(1000 / n) * sqrt(P(Q) / Pg)} (sign-free, from
#' \code{V_eff[Q]}).  Opposite-sign families (\code{sign(m_tau[Qc]) != s})
#' are supported: benefit-direction candidates carry oriented-negative
#' \code{beta_g}, and \code{tauQc} may be negative for such families.  When
#' both region effects share a sign the values coincide exactly with the
#' former oriented-absolute reading.  A DGM with \code{m_tau[Q]} exactly zero
#' is rejected (no harm direction to orient by): plant a nonzero Q effect or
#' use the null path.
#'
#' @param dgm An object of class \code{"glm_dgm"} from
#'   \code{\link{generate_glm_dgm}}, with the true-region indicator in
#'   \code{df_super$flag_harm}.
#' @param forestsearch_args Named list of \code{\link{forestsearch}} arguments.
#'   \code{confounders.name} is required; the cut, \code{maxk}, \code{minp},
#'   \code{rmin}, \code{n.min} and \code{n.min.frac} entries are read, and any
#'   omitted entry takes \code{forestsearch()}'s default.
#' @param n Integer.  Trial size at which the family is evaluated: sets the
#'   size floor \code{n.min / n} and the \code{sqrt(1000 / n)} rescaling of the
#'   anchored standard errors.
#' @param max_M Integer.  Size guard: if more than \code{max_M} candidates
#'   survive the floors, the function stops before allocating the M-by-M
#'   overlap matrix (the downstream symmetric root is O(M^3)).
#' @param verbose Logical.  Print the enumeration counts at each stage.
#'
#' @return An object of class \code{c("fs_oc_family", "list")} with elements
#'   \describe{
#'     \item{\code{lab}}{Character, length M: the rule of each candidate.}
#'     \item{\code{Pg}}{Population prevalence \code{P(g)}.}
#'     \item{\code{PQg}}{Purity \code{P(g & Q) / P(g)}.}
#'     \item{\code{beta_g}}{Oriented mixture mean
#'       \code{s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)}, with
#'       \code{s = sign(m_tau[Q])} the harm direction.}
#'     \item{\code{se_g}}{Anchored standard error at \code{n}.}
#'     \item{\code{sens_g}}{\code{P(g & Q) / P(Q)}.}
#'     \item{\code{spec_g}}{\code{1 - P(g & Qc) / P(Qc)}.}
#'     \item{\code{ovl}}{M-by-M matrix of \code{P(g_i & g_j)}.}
#'     \item{\code{M}}{Number of candidates.}
#'     \item{\code{PQ}}{\code{P(Q)}, the prevalence of the true region.}
#'     \item{\code{memb}}{Logical matrix, \code{nrow(df_super)} by M, the
#'       population membership of each candidate.}
#'     \item{\code{null}}{Logical: \code{TRUE} when the null branch was
#'       taken (Q empty).}
#'     \item{\code{orientation}}{Alternative branch only: list with the
#'       harm-direction sign \code{s}, the signed region effects
#'       \code{m_tau_Q} and \code{m_tau_Qc}, and the oriented mixture
#'       coefficients \code{tauQc = s * m_tau_Qc} (may be negative for
#'       opposite-sign families) and \code{bint = s * (m_tau_Q - m_tau_Qc)}.
#'       Absent on the null branch.}
#'     \item{\code{scale}}{The \code{fs_dgm_scale} object used.}
#'     \item{\code{n}}{The trial size.}
#'     \item{\code{args_used}}{The \code{forestsearch_args} entries consumed,
#'       with the resolved \code{n.min}.}
#'     \item{\code{cuts}}{The cut expressions \code{get_FSdata()} produced.}
#'     \item{\code{counts}}{Named integer vector: candidates at each stage.}
#'   }
#'
#' @seealso \code{\link{fs_oc_predict}}, \code{\link{fs_dgm_scale}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 2000
#' age <- round(rnorm(N, 35, 9)); pre <- round(rexp(N, 1 / 500))
#' V   <- factor(rbinom(N, 1, 0.42), levels = 0:1)
#' inQ <- as.integer(age > 34 & pre <= 745)
#' mu0 <- 40 + 0.2 * age
#' dgm <- structure(list(
#'   df_super = data.frame(age = age, preanti = pre, V = V,
#'                         mu0 = mu0, mu1 = mu0 - 26 - 14 * inQ,
#'                         flag_harm = inQ),
#'   outcome_type = "continuous", effect_measure = "MD",
#'   model_params = list(sigma = 127.5)), class = c("glm_dgm", "list"))
#' fam <- fs_oc_family_enumerate(
#'   dgm, list(confounders.name = c("age", "preanti", "V"),
#'             conf.cont_jcuts = list(age = 4, preanti = 4), n.min = 60),
#'   n = 500)
#' fam
#' }
#'
#' @export
fs_oc_family_enumerate <- function(dgm, forestsearch_args, n,
                                   max_M = 2000L, verbose = FALSE) {

  # ---- validation -----------------------------------------------------------
  if (!inherits(dgm, "glm_dgm")) {
    stop("'dgm' must be an object of class 'glm_dgm'.", call. = FALSE)
  }
  if (!is.data.frame(dgm$df_super)) {
    stop("'dgm$df_super' must be a data frame.", call. = FALSE)
  }
  if (!"flag_harm" %in% names(dgm$df_super)) {
    stop("'dgm$df_super' has no 'flag_harm' column.", call. = FALSE)
  }
  if (!is.list(forestsearch_args)) {
    stop("'forestsearch_args' must be a named list.", call. = FALSE)
  }
  if (is.null(forestsearch_args$confounders.name)) {
    stop("'forestsearch_args$confounders.name' is required.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0) {
    stop("'n' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(max_M) || length(max_M) != 1L || is.na(max_M) || max_M < 1) {
    stop("'max_M' must be a single positive number.", call. = FALSE)
  }

  df_super <- dgm$df_super
  n_super  <- nrow(df_super)
  inQ      <- df_super$flag_harm == 1
  # ---- null detection: structural, from the harm flag itself ---------------
  # A null DGM has Q = {} -- generate_glm_dgm() zeroes flag_harm under
  # model = "null" (R/generate_glm_dgm.R L316-319), and the driver asserts
  # sum(inQ) == 0L, is.na(effect_Q), beta_inter == 0 for the null cell.  An
  # alternative DGM always has a non-empty calibrated Q, so sum(flag_harm) == 0
  # cannot misfire on one.  When the object declares dgm$model (exact match --
  # `$` would partial-match model_params), it must agree.
  is_null <- !any(inQ)
  if (!is.null(dgm[["model"]]) && identical(dgm[["model"]], "null") && !is_null) {
    stop("dgm$model is \"null\" but df_super$flag_harm has ", sum(inQ),
         " members in Q; the object is inconsistent.", call. = FALSE)
  }
  if (!is.null(dgm[["model"]]) && !identical(dgm[["model"]], "null") && is_null) {
    stop("df_super$flag_harm has no members in Q but dgm$model is \"",
         dgm[["model"]], "\"; the object is inconsistent.", call. = FALSE)
  }

  # ---- the forestsearch() arguments consumed, with its defaults -------------
  fs_formals <- formals(forestsearch)
  .arg <- function(nm) {
    if (nm %in% names(forestsearch_args)) forestsearch_args[[nm]]
    else eval(fs_formals[[nm]])
  }
  cut_arg_names <- c("confounders.name", "conf.cont_jcuts", "cut_type",
                     "cont.cutoff", "conf.cont_medians",
                     "conf.cont_medians_force", "conf_force",
                     "defaultcut_names", "exclude_cuts", "collapse_cuts",
                     "collapse_cuts_args")
  cut_args <- lapply(cut_arg_names, .arg)
  names(cut_args) <- cut_arg_names

  maxk <- as.integer(.arg("maxk"))
  minp <- .arg("minp")
  # rmin is not a forestsearch() formal: subgroup.search()'s default (5) is
  # what the search receives unless the caller overrides it.
  rmin <- if ("rmin" %in% names(forestsearch_args)) forestsearch_args[["rmin"]]
          else eval(formals(subgroup.search)[["rmin"]])
  # sg_focus = "maxeff" relaxes both floors to their most permissive setting
  # (forestsearch_main.R: minp -> 0, search_overrides$rmin <- 0).
  if (identical(.arg("sg_focus"), "maxeff")) { minp <- 0; rmin <- 0 }
  n.min.frac <- .arg("n.min.frac")
  # n.min: forestsearch_main.R SECTION 1A2 -- the supplied value (default 60)
  # or, when n.min = NULL is passed, max(60, ceiling(n.min.frac * N)).
  if ("n.min" %in% names(forestsearch_args) &&
      is.null(forestsearch_args[["n.min"]])) {
    if (!is.numeric(n.min.frac) || length(n.min.frac) != 1L ||
        is.na(n.min.frac) || n.min.frac <= 0 || n.min.frac >= 1) {
      stop("n.min.frac must be a single number in (0, 1).", call. = FALSE)
    }
    n.min <- max(60L, as.integer(ceiling(n.min.frac * n)))
  } else {
    n.min <- .arg("n.min")
  }
  if (!is.numeric(maxk) || maxk < 1L || maxk > 3L) {
    stop("'maxk' must be 1, 2 or 3.", call. = FALSE)
  }

  # ---- 1. scale: every mean and se from fs_dgm_scale(dgm) -------------------
  if (!is_null) {
    scale <- fs_dgm_scale(dgm)
    reg   <- scale$regions
    iQ  <- match("Q",  reg$region)
    iQc <- match("Qc", reg$region)
    if (is.na(iQ) || is.na(iQc)) {
      stop("fs_dgm_scale(dgm) did not return both 'Q' and 'Qc' regions.",
           call. = FALSE)
    }
    # Signed orientation: s is the harm direction sign(m_tau[Q]), so the
    # planted effect is oriented positive and the general signed mixture
    #   beta_g = s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)
    # applies.  When both region effects share the sign this is algebraically
    # (and in IEEE arithmetic bit-) identical to the previous |.| reading
    # (s*m_Qc = |m_Qc|, s*(m_Q - m_Qc) = |m_Q| - |m_Qc|); when they differ,
    # benefit-direction candidates carry oriented-negative means, so tauQc
    # may be negative.
    m_Q_signed  <- reg$m_tau[iQ]
    m_Qc_signed <- reg$m_tau[iQc]
    s <- sign(m_Q_signed)
    if (s == 0) {
      stop("m_tau[Q] is exactly zero: there is no harm direction to orient ",
           "by.  Plant a nonzero Q effect, or use the null path ",
           "(generate_glm_dgm(model = \"null\")) for a homogeneous-effect ",
           "family.", call. = FALSE)
    }
    piQ     <- reg$P_g[iQ]
    tauQc   <- s * m_Qc_signed
    bint    <- s * (m_Q_signed - m_Qc_signed)
    seQ1000 <- sqrt(reg$V_eff[iQ] / (1000 * piQ))
    PQ      <- mean(inQ)
  } else {
    # ---- null branch: no Q/Qc partition; the whole-population row only ----
    # fs_dgm_scale()'s default regions need a non-empty Q (it stops with
    # "region 'Q' is empty"); its public `regions` argument computes the same
    # columns on any region, so the S row is obtained without touching that
    # file.  Every candidate has the same true effect (the DGM's effect_Qc =
    # effect_ITT), oriented as the alternative path orients (abs of a common
    # sign), and its se is the whole-population effective variance at (n, Pg):
    #   se_g = sqrt(V_eff[S] / (1000 * 1)) * sqrt(1000 / n) * sqrt(1 / Pg),
    # the alternative's prevalence scaling with P(Q) -> P(S) = 1.
    scale <- fs_dgm_scale(dgm, regions = list(S = rep(TRUE, n_super)))
    reg   <- scale$regions
    iS    <- match("S", reg$region)
    piQ     <- 1                      # the reference region is S, P(S) = 1
    tauQc   <- abs(reg$m_tau[iS])     # the common effect, oriented
    bint    <- 0
    seQ1000 <- sqrt(reg$V_eff[iS] / (1000 * piQ))
    PQ      <- 0
  }

  # ---- 2. cuts on the population frame via get_FSdata() ---------------------
  # get_FSdata() validates an outcome and an event column; neither is used
  # with LASSO off, so two constant numeric columns satisfy the contract.
  df_cut <- df_super
  df_cut[[".fs_oc_y"]]     <- 0
  df_cut[[".fs_oc_event"]] <- 1
  FSdata <- get_FSdata(
    df.analysis             = df_cut,
    use_lasso               = FALSE,
    use_grf                 = FALSE,
    grf_cuts                = NULL,
    confounders.name        = cut_args$confounders.name,
    cont.cutoff             = cut_args$cont.cutoff,
    conf_force              = cut_args$conf_force,
    conf.cont_medians       = cut_args$conf.cont_medians,
    conf.cont_medians_force = cut_args$conf.cont_medians_force,
    conf.cont_jcuts         = cut_args$conf.cont_jcuts,
    dina_cuts               = NULL,
    collapse_cuts           = cut_args$collapse_cuts,
    collapse_cuts_args      = cut_args$collapse_cuts_args,
    defaultcut_names        = cut_args$defaultcut_names,
    cut_type                = cut_args$cut_type,
    exclude_cuts            = cut_args$exclude_cuts,
    outcome.name            = ".fs_oc_y",
    event.name              = ".fs_oc_event",
    details                 = FALSE,
    outcome_type            = dgm$outcome_type %||% "continuous"
  )
  cuts <- FSdata$confs
  # Both directions: dummy() expands each 2-level factor into two indicator
  # columns (forestsearch_main.R, "dummy() expands each 2-level factor ...").
  Zdf <- dummy(FSdata$df[, FSdata$confs_names, drop = FALSE])
  Z   <- as.matrix(Zdf)
  storage.mode(Z) <- "numeric"
  colnames(Z) <- names(Zdf)
  L <- ncol(Z)
  col_lab <- .fs_oc_column_labels(colnames(Z), FSdata$confs_names, cuts)

  # ---- 3. enumerate all <= maxk combinations (forestsearch_main.R MR family)
  combo <- generate_combination_indices(L, maxk)
  tot   <- calculate_max_combinations(L, maxk)

  counts <- c(cut_columns = L, enumerated = 0L, empty = 0L,
              minp = 0L, rmin = 0L, size = 0L, kept = 0L,
              duplicate = 0L, M = 0L)
  memb_list <- list()
  lab_list  <- character(0)
  rmin_prop <- rmin / n

  for (kk in seq_len(tot)) {
    covs.in <- get_covs_in(kk, maxk, L,
                           combo$counts_1, combo$indices_1,
                           combo$counts_2, combo$indices_2,
                           combo$counts_3, combo$indices_3)
    k_sel <- sum(covs.in)
    if (k_sel < 1L || k_sel > maxk) next
    counts["enumerated"] <- counts["enumerated"] + 1L
    sel <- which(covs.in == 1)
    x   <- Z[, sel, drop = FALSE]

    # status 0: a constant column or an empty pairwise intersection
    # (has_positive_variance: every entry of x'x positive)
    if (!has_positive_variance(x)) {
      counts["empty"] <- counts["empty"] + 1L; next
    }
    # status 1: per-factor population prevalence floor
    if (!meets_prevalence_threshold(x, minp)) {
      counts["minp"] <- counts["minp"] + 1L; next
    }
    # status 2: redundancy -- each added factor must shrink the membership by
    # more than rmin subjects of a trial of size n, i.e. rmin / n in proportion
    if (.fs_oc_redundant(x, rmin_prop)) {
      counts["rmin"] <- counts["rmin"] + 1L; next
    }
    # status 4: subgroup size floor as a population prevalence
    memb <- get_subgroup_membership(Z, covs.in)
    Pg_k <- mean(memb)
    if (Pg_k < n.min / n) {
      counts["size"] <- counts["size"] + 1L; next
    }
    counts["kept"] <- counts["kept"] + 1L
    memb_list[[length(memb_list) + 1L]] <- memb
    lab_list <- c(lab_list, paste(col_lab[sel], collapse = " & "))
  }

  if (length(memb_list) == 0L) {
    stop("no candidate survives the structural floors (n.min = ", n.min,
         ", n = ", n, ").", call. = FALSE)
  }

  # ---- 4. collapse identical population memberships -------------------------
  memb_mat <- do.call(cbind, memb_list)
  key <- apply(memb_mat, 2L, function(v) paste(which(v), collapse = ","))
  first <- !duplicated(key)
  counts["duplicate"] <- sum(!first)
  memb_mat <- memb_mat[, first, drop = FALSE]
  lab <- lab_list[first]
  M <- ncol(memb_mat)
  counts["M"] <- M
  colnames(memb_mat) <- NULL

  if (verbose) {
    cat("fs_oc_family_enumerate: cut columns", L, "| enumerated",
        counts["enumerated"], "| dropped: empty", counts["empty"],
        "minp", counts["minp"], "rmin", counts["rmin"], "size",
        counts["size"], "| kept", counts["kept"], "| duplicates",
        counts["duplicate"], "| M =", M, "\n")
  }

  # ---- size guard, before any M x M allocation ------------------------------
  if (M > max_M) {
    stop(sprintf(paste0(
      "fs_oc_family_enumerate: %d candidates survive the floors, above ",
      "max_M = %d.  The M x M overlap/covariance matrices would take %.1f MB ",
      "each and the symmetric root is O(M^3) (~%.2e flops).  Raise max_M ",
      "deliberately, or tighten the cut specification."),
      M, as.integer(max_M), 8 * M^2 / 2^20, M^3), call. = FALSE)
  }

  # ---- 5. the interface fields, all population proportions ------------------
  Pg   <- colMeans(memb_mat)
  PgQ  <- colMeans(memb_mat & inQ)          # P(g & Q)
  PgQc <- colMeans(memb_mat & !inQ)         # P(g & Qc)
  if (!is_null) {
    PQg    <- PgQ / Pg
    beta_g <- tauQc + bint * PQg
    se_g   <- seQ1000 * sqrt(1000 / n) * sqrt(piQ / Pg)
    sens_g <- PgQ / PQ
    spec_g <- 1 - PgQc / (1 - PQ)
  } else {
    # null: P(g & Q) = 0 for every g, P(Q) = 0
    PQg    <- rep(0, M)                 # purity 0
    beta_g <- rep(tauQc, M)             # the common effect
    se_g   <- seQ1000 * sqrt(1000 / n) * sqrt(piQ / Pg)   # piQ = 1
    sens_g <- rep(NA_real_, M)          # 0 / 0: undefined, not zero
    spec_g <- 1 - Pg                    # Qc is the whole population
    # npv = (1 - Pg - (PQ - PQg*Pg)) / (1 - Pg) = 1 follows from PQ = 0 in
    # fs_oc_predict(); nothing to store.
  }
  ovl    <- crossprod(memb_mat * 1) / n_super
  dimnames(ovl) <- NULL

  out <- list(
    lab = lab, Pg = unname(Pg), PQg = unname(PQg), beta_g = unname(beta_g),
    se_g = unname(se_g), sens_g = unname(sens_g), spec_g = unname(spec_g),
    ovl = ovl, M = M, PQ = PQ,
    memb = memb_mat,
    null = is_null,
    scale = scale, n = n,
    args_used = c(cut_args, list(maxk = maxk, minp = minp, rmin = rmin,
                                 n.min = n.min, n.min.frac = n.min.frac)),
    cuts = cuts,
    counts = counts
  )
  if (!is_null) {
    # Orientation provenance (alternative branch only; the null branch is
    # unchanged in structure): the harm-direction sign, the signed region
    # effects, and the oriented mixture coefficients actually used.
    out$orientation <- list(s = s, m_tau_Q = m_Q_signed,
                            m_tau_Qc = m_Qc_signed,
                            tauQc = tauQc, bint = bint)
  }
  class(out) <- c("fs_oc_family", "list")
  out
}


#' @export
print.fs_oc_family <- function(x, ...) {
  cat("Population-enumerated candidate family (fs_oc_family)\n")
  cat(sprintf("  M            : %d candidates (from %d cut columns, maxk = %d)\n",
              x$M, x$counts[["cut_columns"]], x$args_used$maxk))
  cat(sprintf("  n            : %s\n", format(x$n)))
  cat(sprintf("  floors       : n.min = %s (Pg >= %.4f), minp = %s, rmin = %s\n",
              format(x$args_used$n.min), x$args_used$n.min / x$n,
              format(x$args_used$minp), format(x$args_used$rmin)))
  cat(sprintf("  prevalence   : %.4f to %.4f\n", min(x$Pg), max(x$Pg)))
  if (isTRUE(x$null)) {
    cat(sprintf("  NULL DGM     : Q empty; common effect %.4f; se from the S-row V_eff\n",
                x$beta_g[1]))
  } else {
    cat(sprintf("  P(Q)         : %.4f;  purity range %.3f to %.3f\n",
                x$PQ, min(x$PQg), max(x$PQg)))
  }
  cat(sprintf("  dropped      : empty %d, minp %d, rmin %d, size %d, duplicate %d\n",
              x$counts[["empty"]], x$counts[["minp"]], x$counts[["rmin"]],
              x$counts[["size"]], x$counts[["duplicate"]]))
  invisible(x)
}


# -----------------------------------------------------------------------------
# internal helpers
# -----------------------------------------------------------------------------

# Population analogue of extract_idx_flagredundancy(): TRUE when some added
# factor shrinks the membership PROPORTION by at most rmin_prop.
.fs_oc_redundant <- function(x, rmin_prop) {
  id.x <- rep(1, nrow(x))
  p.last <- 1
  for (m in seq_len(ncol(x))) {
    id.x <- id.x * x[, m]
    p.this <- mean(id.x)
    if (p.last - p.this <= rmin_prop) return(TRUE)
    p.last <- p.this
  }
  FALSE
}

# Readable rule for each dummy() column "qK.<level>": level 1 is the cut
# expression itself; level 0 is its complement.
.fs_oc_column_labels <- function(zcols, q_names, cuts) {
  vapply(zcols, function(z) {
    dot <- regexpr("\\.[^.]*$", z)
    qn  <- substr(z, 1L, dot - 1L)
    lev <- substr(z, dot + 1L, nchar(z))
    k   <- match(qn, q_names)
    if (is.na(k)) return(z)
    expr <- cuts[k]
    if (identical(lev, "1")) return(expr)
    .fs_oc_negate(expr)
  }, character(1), USE.NAMES = FALSE)
}

.fs_oc_negate <- function(expr) {
  ops <- c("<=", ">=", "!=", "==", "<", ">")
  neg <- c(">", "<", "==", "!=", ">=", "<=")
  for (i in seq_along(ops)) {
    if (grepl(ops[i], expr, fixed = TRUE)) {
      parts <- strsplit(expr, ops[i], fixed = TRUE)[[1L]]
      if (length(parts) == 2L) {
        return(paste0(trimws(parts[1L]), " ", neg[i], " ", trimws(parts[2L])))
      }
    }
  }
  # bare variable: evaluate_comparison() reads it as `var == 1`
  paste0(trimws(expr), " != 1")
}
