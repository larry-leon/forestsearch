# verify_loo.R
# One-off check that forestsearch_Kfold(fs, Kfolds = n) implements the
# LOO procedure as expected:
#   For each subject k = 1, ..., n:
#     - Refit forestsearch() on the remaining n - 1 subjects.
#     - Predict subgroup status for the deleted subject k using the
#       training-fold subgroup definition.
#
# This file is intentionally lightweight and project-local; it is not a
# package addition.  Source it from your analysis qmd / R session,
# run once, then delete.
#
# Usage (after a forestsearch run `fs_result` and a corresponding
# `forestsearch_Kfold(fs_result, Kfolds = n)` result `fs_loo_result`):
#
#   source("verify_loo.R")
#   chk <- verify_loo(fs_result, fs_loo_result)
#   chk$summary
#   chk$disagreements   # rows where manual LOO != forestsearch_Kfold output
#
# Optional: skip the package's K-fold result and let the helper run
# forestsearch_Kfold() internally:
#
#   chk <- verify_loo(fs_result, kfold_subjects = 1:25)
#
# `kfold_subjects` lets you spot-check a subset rather than all n —
# useful when n is large and you just want a sanity sample.

verify_loo <- function(fs,
                       fs_kfold       = NULL,
                       kfold_subjects = NULL,
                       seedit         = 8316951L,
                       verbose        = TRUE) {

  # --- 0. Inputs and setup -------------------------------------------------
  stopifnot(!is.null(fs$df.est),
            !is.null(fs$args_call_all))

  df_est  <- fs$df.est
  n       <- nrow(df_est)
  id_name <- fs$args_call_all$id.name

  subjects <- if (is.null(kfold_subjects)) seq_len(n) else kfold_subjects
  stopifnot(all(subjects %in% seq_len(n)))

  # --- 1. Obtain the forestsearch_Kfold LOO result ------------------------
  if (is.null(fs_kfold)) {
    if (verbose) cat("verify_loo(): running forestsearch_Kfold(Kfolds = ",
                     n, ") ...\n", sep = "")
    fs_kfold <- forestsearch_Kfold(fs, Kfolds = n, seedit = seedit,
                                   details = FALSE)
  }

  # The per-subject K-fold predictions live on fs_kfold$resCV, keyed by
  # id.name.  For LOO (Kfolds == n) the random shuffle at
  # forestsearch_cross_validation.R:189 is disabled, so each fold
  # corresponds to a single subject; we key on id.name to be safe
  # regardless of fold size.
  if (is.null(fs_kfold$resCV)) {
    stop("verify_loo(): fs_kfold$resCV is missing.  ",
         "Is `fs_kfold` a forestsearch_Kfold() result?")
  }
  if (!id_name %in% names(fs_kfold$resCV)) {
    stop("verify_loo(): id column '", id_name,
         "' not found in fs_kfold$resCV.")
  }
  kfold_lookup <- setNames(fs_kfold$resCV$treat.recommend,
                           fs_kfold$resCV[[id_name]])

  # --- 2. Build the cv_args template (matches lines 208-217 of
  #         forestsearch_cross_validation.R) ---------------------------------
  cv_args_template <- fs$args_call_all
  cv_args_template$parallel_args <- list(plan = "sequential", workers = 1L,
                                         show_message = FALSE)
  cv_args_template$details <- FALSE
  cv_args_template$quiet   <- TRUE
  cv_args_template$plot.sg <- FALSE
  cv_args_template$ps_hat  <- NULL

  # --- 3. Manual LOO loop --------------------------------------------------
  result_rows <- vector("list", length(subjects))
  if (verbose) {
    cat("verify_loo(): manually replicating LOO for ", length(subjects),
        " subject(s).\n", sep = "")
    pb <- txtProgressBar(min = 0, max = length(subjects), style = 3)
  }

  for (i in seq_along(subjects)) {
    k <- subjects[i]

    cv_args <- cv_args_template
    cv_args$df.analysis <- df_est[-k, , drop = FALSE]

    fs_train <- suppressWarnings(try(do.call(forestsearch, cv_args),
                                     silent = TRUE))

    if (inherits(fs_train, "try-error")) {
      manual_recommend <- NA_integer_
      manual_sg_factors <- "<train error>"
      manual_sg_found   <- FALSE
    } else if (is.null(fs_train$sg.harm)) {
      manual_recommend  <- 1L                                  # ITT fallback
      manual_sg_factors <- ""
      manual_sg_found   <- FALSE
    } else {
      pred_k <- get_dfpred(df.predict = df_est[k, , drop = FALSE],
                           sg.harm    = fs_train$sg.harm,
                           version    = 2)
      manual_recommend  <- as.integer(pred_k$treat.recommend)
      manual_sg_factors <- paste(fs_train$sg.harm, collapse = " & ")
      manual_sg_found   <- TRUE
    }

    subject_id <- df_est[[id_name]][k]
    kfold_recommend <- kfold_lookup[[as.character(subject_id)]]

    result_rows[[i]] <- data.frame(
      subject_idx     = k,
      subject_id      = subject_id,
      manual_recommend = manual_recommend,
      kfold_recommend  = if (is.null(kfold_recommend)) NA_integer_
                         else as.integer(kfold_recommend),
      agree            = !is.na(manual_recommend) &&
                         !is.null(kfold_recommend) &&
                         manual_recommend == as.integer(kfold_recommend),
      manual_sg_found  = manual_sg_found,
      manual_sg        = manual_sg_factors,
      stringsAsFactors = FALSE
    )

    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) { close(pb); cat("\n") }

  per_subject <- do.call(rbind, result_rows)
  rownames(per_subject) <- NULL

  # --- 4. Summary ----------------------------------------------------------
  n_chk        <- nrow(per_subject)
  n_agree      <- sum(per_subject$agree, na.rm = TRUE)
  n_disagree   <- sum(!per_subject$agree, na.rm = TRUE)
  n_sg_found   <- sum(per_subject$manual_sg_found)
  sg_find_rate <- 100 * n_sg_found / n_chk

  summary_tbl <- data.frame(
    metric = c("Subjects checked",
               "Manual LOO matches forestsearch_Kfold",
               "Disagreements",
               "Manual LOO found a subgroup",
               "Subgroup-find rate (%)"),
    value  = c(n_chk, n_agree, n_disagree, n_sg_found,
               sprintf("%.1f", sg_find_rate)),
    stringsAsFactors = FALSE
  )

  disagreements <- per_subject[!per_subject$agree, , drop = FALSE]

  if (verbose) {
    cat("\n=== verify_loo() summary ===\n")
    print(summary_tbl, row.names = FALSE)
    if (nrow(disagreements) > 0L) {
      cat("\nDisagreements (manual vs forestsearch_Kfold):\n")
      print(disagreements[, c("subject_idx", "subject_id",
                              "manual_recommend", "kfold_recommend",
                              "manual_sg")],
            row.names = FALSE)
    } else {
      cat("\nAll subjects agree -- forestsearch_Kfold(Kfolds = n) ",
          "implements LOO as expected.\n", sep = "")
    }
  }

  invisible(list(
    summary       = summary_tbl,
    per_subject   = per_subject,
    disagreements = disagreements
  ))
}
