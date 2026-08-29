# =============================================================================
# oc_wrapper_confs_compare_2026-08-30.R -- does the confounder correction close
# the residuals?  (§4 of dev/tasks/cc_task_oc_wrapper_confs_2026-08-30.md)
# -----------------------------------------------------------------------------
# Reads: oc_wrapper_grid_corrected_2026-08-30.rds (13-variable family),
#        oc_wrapper_grid_2026-08-29.rds / oc_wrapper_null_2026-08-29.rds
#        (12-variable family, superseded), sgdef_selection_2026-08-29.rds (the
#        realized selection distribution and per-rule population values --
#        NOT recomputed).
# Read-only diagnostic; no fs_oc_* re-runs.  Writes
#   dev/glm-continuous-sims/oc_wrapper_confs_compare_2026-08-30.rds
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
D <- "dev/glm-continuous-sims"
C  <- readRDS(file.path(D, "oc_wrapper_grid_corrected_2026-08-30.rds"))
G  <- readRDS(file.path(D, "oc_wrapper_grid_2026-08-29.rds"))
N0 <- readRDS(file.path(D, "oc_wrapper_null_2026-08-29.rds"))
S  <- readRDS(file.path(D, "sgdef_selection_2026-08-29.rds"))
OUT <- file.path(D, "oc_wrapper_confs_compare_2026-08-30.rds")
set.seed(20260830L)   # the multinomial noise-floor simulation only

# ---- signature parser: sgdef_selection_2026-08-29.R, verbatim ----------------
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
  paste0(cl, if (neg) "=0" else "=1")
}
rule_sig <- function(rule) {
  cls <- strsplit(rule, " & ", fixed = TRUE)[[1]]
  paste(sort(vapply(cls, clause_sig, character(1))), collapse = " & ")
}

# ---- the three cells: old / corrected / measured / sgdef ---------------------
cell_spec <- list(
  alt500  = list(old = list(resample = G$runs$n500_resample$pred, split = G$runs$n500_split$pred),
                 new = list(resample = C$alt$runs$n500_resample$pred, split = C$alt$runs$n500_split$pred),
                 meas = C$alt$measured[["500"]], n = 500, null = FALSE,
                 M_old = G$families[["500"]]$M, M_new = C$alt$families[["500"]]$M,
                 fam_old = G$families[["500"]], fam_new = C$alt$families[["500"]]),
  alt700  = list(old = list(resample = G$runs$n700_resample$pred, split = G$runs$n700_split$pred),
                 new = list(resample = C$alt$runs$n700_resample$pred, split = C$alt$runs$n700_split$pred),
                 meas = C$alt$measured[["700"]], n = 700, null = FALSE,
                 M_old = G$families[["700"]]$M, M_new = C$alt$families[["700"]]$M,
                 fam_old = G$families[["700"]], fam_new = C$alt$families[["700"]]),
  null500 = list(old = list(resample = N0$runs$resample$pred, split = N0$runs$split$pred),
                 new = list(resample = C$null$runs$resample$pred, split = C$null$runs$split$pred),
                 meas = C$null$measured, n = 500, null = TRUE,
                 M_old = N0$family$M, M_new = C$null$family$M,
                 fam_old = N0$family, fam_new = C$null$family))

qs <- c("det_rate", "EnH", "Esens", "Espec", "Eppv", "Enpv", "EbetaH", "Enaive_bias", "mass_below")
tvd <- function(a, b) 0.5 * sum(abs(a - b))
align <- function(x, keys) { z <- setNames(rep(0, length(keys)), keys); z[names(x)] <- x; z }

res <- list()
for (cl in names(cell_spec)) {
  sp <- cell_spec[[cl]]; sg <- S$cells[[cl]]; d <- sg$reps
  cat(sprintf("\n################ %s (n = %d)  M: %d -> %d\n", cl, sp$n, sp$M_old, sp$M_new))

  # -- 1. the functionals, old / corrected / delta / measured, both gates --------
  fun <- do.call(rbind, lapply(c("resample", "split"), function(g) {
    o <- sp$old[[g]]; nw <- sp$new[[g]]
    data.frame(gate = g, quantity = qs,
               old = vapply(qs, function(q) as.numeric(o[[q]]), 1),
               old_se = vapply(qs, function(q) { v <- o[[paste0(q, "_se")]]; if (is.null(v)) NA_real_ else as.numeric(v) }, 1),
               corrected = vapply(qs, function(q) as.numeric(nw[[q]]), 1),
               corrected_se = vapply(qs, function(q) { v <- nw[[paste0(q, "_se")]]; if (is.null(v)) NA_real_ else as.numeric(v) }, 1),
               measured = vapply(qs, function(q) { v <- sp$meas[[q]]; if (is.null(v)) NA_real_ else as.numeric(v) }, 1),
               measured_se = vapply(qs, function(q) { v <- sp$meas[[paste0(q, "_se")]]; if (is.null(v)) NA_real_ else as.numeric(v) }, 1),
               stringsAsFactors = FALSE)
  }))
  fun$delta <- fun$corrected - fun$old
  fun$gap_old <- fun$measured - fun$old
  fun$gap_corrected <- fun$measured - fun$corrected
  cat("functionals (old = 12-variable family, corrected = 13-variable, measured = payload):\n")
  print(fun, digits = 5, row.names = FALSE)

  # -- 2. signatures: corrected analytic sel_c aggregated by signature ------------
  lab_new <- sp$new$resample$lab; sig_new <- vapply(lab_new, rule_sig, character(1), USE.NAMES = FALSE)
  lab_old <- sg$analytic$lab;     sig_old <- sg$analytic$sig
  stopifnot(length(lab_old) == sp$M_old, length(lab_new) == sp$M_new)
  A_new_r <- tapply(sp$new$resample$sel_c, sig_new, sum); A_new_s <- tapply(sp$new$split$sel_c, sig_new, sum)
  A_old_r <- sg$dist$A_resample; Mf <- sg$dist$M
  keys <- union(union(names(A_old_r), names(Mf)), names(A_new_r))
  A_old <- align(A_old_r, keys); A_new <- align(A_new_r, keys); A_new_sp <- align(A_new_s, keys); Mf <- align(Mf, keys)
  n_det <- sg$detected

  # str2 mass
  str2_new_r <- sum(sp$new$resample$sel_c[grepl("str2", lab_new)])
  str2_new_s <- sum(sp$new$split$sel_c[grepl("str2", lab_new)])
  str2_meas  <- sum(Mf[grepl("str2", keys)])
  str2_meas_n <- sum(grepl("str2", d$sg_def))
  str2_meas_se <- sqrt(str2_meas * (1 - str2_meas) / n_det)
  cat(sprintf("\nstr2 mass: corrected analytic resample %.4f, split %.4f | measured %.4f (%d of %d; binomial SE %.4f) | old analytic 0 by construction\n",
              str2_new_r, str2_new_s, str2_meas, str2_meas_n, n_det, str2_meas_se))
  # per-str2-signature
  s2 <- keys[grepl("str2", keys) & (Mf > 0 | A_new > 0)]
  str2_tab <- data.frame(signature = s2, analytic_corrected = round(A_new[s2], 4), measured = round(Mf[s2], 4),
                         measured_n = as.integer(round(Mf[s2] * n_det)))
  str2_tab <- str2_tab[order(-str2_tab$analytic_corrected), ]
  cat("str2 signatures (corrected analytic resample beside measured):\n"); print(head(str2_tab, 12), row.names = FALSE)
  # the population size and purity of the str2 candidates the family now contains
  i2 <- grepl("str2", lab_new)
  cat(sprintf("str2 candidates in the corrected family: %d of %d; sel_c-weighted mean n*Pg %.2f; family mean n*Pg %.2f; PQg (str2 cands, sel_c-weighted) %.3f\n",
              sum(i2), sp$M_new, sum(sp$new$resample$sel_c[i2] * sp$fam_new$Pg[i2]) / sum(sp$new$resample$sel_c[i2]) * sp$n,
              sum(sp$new$resample$sel_c * sp$fam_new$Pg) * sp$n,
              if (sp$null) 0 else sum(sp$new$resample$sel_c[i2] * sp$fam_new$PQg[i2]) / sum(sp$new$resample$sel_c[i2])))

  # -- 3. the between-rule size gap, sgdef §4 definition -------------------------
  size_sig <- tapply(d$nPg_pop, d$sig, mean)           # population size of the measured rules, per signature
  s <- names(size_sig)
  rw <- function(A) { w <- A[s]; w[is.na(w)] <- 0; sum(w * size_sig) / sum(w) }
  pop_realized <- mean(d$nPg_pop)                       # measured-frequency-weighted population size
  between <- data.frame(
    gate = c("resample", "split"),
    analytic_EnH_old = c(sp$old$resample$EnH, sp$old$split$EnH),
    analytic_EnH_corrected = c(sp$new$resample$EnH, sp$new$split$EnH),
    pop_of_realized = pop_realized,
    between_old = pop_realized - c(sp$old$resample$EnH, sp$old$split$EnH),
    between_corrected = pop_realized - c(sp$new$resample$EnH, sp$new$split$EnH),
    reweight_old = c(rw(A_old), rw(align(sg$dist$A_split, keys))),
    reweight_corrected = c(rw(A_new), rw(A_new_sp)),
    within = mean(d$n_harm) - pop_realized,
    measured_EnH = mean(d$n_harm))
  cat("\nbetween-rule size gap (pop. size of the realized rules minus analytic sel_c-weighted expectation):\n")
  print(between, digits = 5, row.names = FALSE)

  # -- 4. distribution comparison on the corrected family ------------------------
  top_new <- names(sort(A_new, decreasing = TRUE))[1:15]
  top_tab <- data.frame(signature = top_new, analytic_corrected = round(A_new[top_new], 4),
                        analytic_old = round(A_old[top_new], 4), measured = round(Mf[top_new], 4),
                        measured_n = as.integer(round(Mf[top_new] * n_det)))
  cat("\ntop 15 corrected analytic signatures (resample sel_c) beside old analytic and measured:\n"); print(top_tab, row.names = FALSE)
  never_new <- sum(A_new[Mf == 0]); never_old <- sum(A_old[Mf == 0])
  absent_new <- sum(Mf[A_new == 0]); absent_old <- sum(Mf[A_old == 0])
  absent_sigs_new <- keys[A_new == 0 & Mf > 0]
  cat(sprintf("analytic mass on never-selected signatures: old %.4f -> corrected %.4f (split %.4f)\n", never_old, never_new, sum(A_new_sp[Mf == 0])))
  cat(sprintf("measured mass absent from the analytic family: old %.4f (%d sigs) -> corrected %.4f (%d sigs)\n",
              absent_old, sum(A_old == 0 & Mf > 0), absent_new, length(absent_sigs_new)))
  if (length(absent_sigs_new)) { cat("still absent:\n"); print(sort(Mf[absent_sigs_new], decreasing = TRUE)) }
  # why still absent: evaluate against the family's floors.  The realized rules for these
  # signatures have population sizes in sgdef $rules; compare to n.min/n.
  if (length(absent_sigs_new)) {
    ru <- sg$rules[sg$rules$sig %in% absent_sigs_new, c("sg_def", "sig", "Pg_pop", "nPg_pop")]
    ru$floor_nPg <- sp$fam_new$n.min; ru$below_floor <- ru$nPg_pop < sp$fam_new$n.min
    cat("realized rules on the still-absent signatures, population size vs the family floor n.min:\n"); print(ru, digits = 4, row.names = FALSE)
  }
  tvd_old <- tvd(A_old, Mf); tvd_new <- tvd(A_new, Mf); tvd_new_s <- tvd(A_new_sp, Mf)
  keep <- !grepl("str2", keys)
  tvd_old_x <- tvd(A_old[keep] / sum(A_old[keep]), Mf[keep] / sum(Mf[keep]))
  tvd_new_x <- tvd(A_new[keep] / sum(A_new[keep]), Mf[keep] / sum(Mf[keep]))
  # multinomial noise floor: 1000-replicate draws from the analytic distribution itself
  floor_sim <- function(A, B = 2000L) {
    p <- A[A > 0]; p <- p / sum(p)
    v <- replicate(B, { x <- stats::rmultinom(1, n_det, p)[, 1] / n_det; tvd(x, p) })
    c(mean = mean(v), q025 = unname(stats::quantile(v, 0.025)), q975 = unname(stats::quantile(v, 0.975)))
  }
  fl_old <- floor_sim(A_old); fl_new <- floor_sim(A_new)
  dist_tab <- data.frame(
    family = c("old (12)", "corrected (13)"),
    tvd_resample = c(tvd_old, tvd_new), tvd_split = c(sg$dist$tvd_split, tvd_new_s),
    tvd_excl_str2 = c(tvd_old_x, tvd_new_x),
    floor_mean = c(fl_old[["mean"]], fl_new[["mean"]]), floor_q025 = c(fl_old[["q025"]], fl_new[["q025"]]),
    floor_q975 = c(fl_old[["q975"]], fl_new[["q975"]]),
    excess = c(tvd_old - fl_old[["mean"]], tvd_new - fl_new[["mean"]]),
    support_analytic = c(sum(A_old > 0), sum(A_new > 0)), support_measured = sum(Mf > 0))
  cat("\ntotal variation distance, analytic vs measured over signatures, and the multinomial floor:\n")
  print(dist_tab, digits = 4, row.names = FALSE)
  # where does the corrected family move mass?  signatures with the largest |A_new - A_old|
  mv <- sort(A_new - A_old); mv_tab <- data.frame(signature = c(names(head(mv, 8)), names(tail(mv, 8))),
                                                  delta = round(c(head(mv, 8), tail(mv, 8)), 4),
                                                  analytic_old = round(A_old[c(names(head(mv, 8)), names(tail(mv, 8)))], 4),
                                                  analytic_corrected = round(A_new[c(names(head(mv, 8)), names(tail(mv, 8)))], 4),
                                                  measured = round(Mf[c(names(head(mv, 8)), names(tail(mv, 8)))], 4))
  cat("largest shifts of analytic mass, old -> corrected (8 down, 8 up):\n"); print(mv_tab, row.names = FALSE)
  # the tilt: mass on age-band and age x {cd40, cd80, wtkg} vs preanti<= & preanti> and cd40/cd80 x wtkg/preanti
  grp <- function(A, pat) sum(A[grepl(pat, keys)])
  tilt <- data.frame(group = c("age<= & age>", "age> & (cd40|cd80|wtkg)", "preanti<= & preanti>", "(cd40|cd80).*(wtkg|preanti)"),
                     old = c(grp(A_old, "^age<= & age>$"), grp(A_old, "^age> & (cd40|cd80|wtkg)"), grp(A_old, "^preanti<= & preanti>$"), grp(A_old, "^(cd40|cd80)[<>=]+ & (preanti|wtkg)")),
                     corrected = c(grp(A_new, "^age<= & age>$"), grp(A_new, "^age> & (cd40|cd80|wtkg)"), grp(A_new, "^preanti<= & preanti>$"), grp(A_new, "^(cd40|cd80)[<>=]+ & (preanti|wtkg)")),
                     measured = c(grp(Mf, "^age<= & age>$"), grp(Mf, "^age> & (cd40|cd80|wtkg)"), grp(Mf, "^preanti<= & preanti>$"), grp(Mf, "^(cd40|cd80)[<>=]+ & (preanti|wtkg)")))
  cat("the sgdef report's tilt groups:\n"); print(tilt, digits = 4, row.names = FALSE)

  # -- 5. sweep and inversion side by side ---------------------------------------
  if (sp$null) { sw_old <- N0$sweep$table; sw_new <- C$null$sweep$table; iv_old <- N0$invert$table; iv_new <- C$null$invert$table
  } else { sw_old <- G$sweep$table[G$sweep$table$n == sp$n, ]; sw_new <- C$alt$sweep$table[C$alt$sweep$table$n == sp$n, ]
           iv_old <- G$invert$table[G$invert$table$n == sp$n, ]; iv_new <- C$alt$invert$table[C$alt$invert$table$n == sp$n, ] }
  sw <- merge(sw_old[, c("consistency_method", "c1", "det_rate", "det_rate_se", "EnH", "Enaive_bias")],
              sw_new[, c("consistency_method", "c1", "det_rate", "det_rate_se", "EnH", "Enaive_bias")],
              by = c("consistency_method", "c1"), suffixes = c("_old", "_new"))
  sw$delta_det <- sw$det_rate_new - sw$det_rate_old; sw <- sw[order(sw$consistency_method, sw$c1), ]
  cat("\nc1 sweep, det_rate old vs corrected:\n"); print(sw, digits = 5, row.names = FALSE)
  old_c1 <- if (sp$null) iv_old$value else iv_old$c1
  new_c1 <- iv_new$value
  iv <- data.frame(gate = iv_new$consistency_method, target = iv_new$target, c1_old = old_c1, c1_corrected = new_c1,
                   delta = new_c1 - old_c1, achieved = iv_new$achieved, ceiling_old = iv_old$ceiling, ceiling_new = iv_new$ceiling)
  cat("inversions, c1 at target, old vs corrected:\n"); print(iv, digits = 5, row.names = FALSE)

  res[[cl]] <- list(n = sp$n, null = sp$null, M_old = sp$M_old, M_new = sp$M_new, functionals = fun,
                    str2 = list(analytic_resample = str2_new_r, analytic_split = str2_new_s, measured = str2_meas,
                                measured_n = str2_meas_n, measured_se = str2_meas_se, table = str2_tab),
                    between = between, top15 = top_tab, never_selected = c(old = never_old, corrected = never_new),
                    absent = list(old = absent_old, corrected = absent_new, sigs = absent_sigs_new),
                    dist = dist_tab, shifts = mv_tab, tilt = tilt, sweep = sw, invert = iv,
                    A_old = A_old, A_new = A_new, A_new_split = A_new_sp, Mf = Mf, sig_new = sig_new, lab_new = lab_new)
}
saveRDS(list(cells = res, built_at = Sys.time(), pkg_version = as.character(utils::packageVersion("forestsearch"))), OUT)
cat("\nwritten:", OUT, "\n")
