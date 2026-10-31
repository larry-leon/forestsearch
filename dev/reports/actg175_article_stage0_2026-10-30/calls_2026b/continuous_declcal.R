# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: fs-glms-interpretable/quarto/actg175/analysis_actg175_continuous_mr.qmd
# fs-glms-interpretable HEAD 19c6ac41c1198dbc513f5c6978d0d6c5df1fc157 (file last changed d97ca1f / 935ae87, 2026-09-23)
#   chunk 'declcal': fence lines 592-605
#   chunk 'declcal-table': fence lines 607-662
#   chunk 'declcal-observed': fence lines 681-690
#   chunk 'declcal-decision': fence lines 694-699
# ----------------------------------------------------------

## ==== chunk 'declcal' ====
c0_grid <- c(-5, -10, -15); c0_primary <- -10
c1_nat <- 10; c2_nat <- 5; fmt_c0 <- function(x) sprintf("%g", x)   # the claim, CD4 cells/mm³
stopifnot(!is.null(g$declaration_field), identical(g$declaration_field$meta$c0, c0_grid))   # the capture fired
cal    <- fs_declaration_calibration(fit, alpha = 0.10, c0 = c0_grid)
cal_c2 <- fs_declaration_calibration(fit, alpha = 0.10)                  # c0 = NULL: the c2 boundary
cal_05 <- fs_declaration_calibration(fit, alpha = 0.05, c0 = c0_grid)    # stored in the payload, not reported
stopifnot(identical(cal$c0$source, "capture"), cal$n_family_prereduction == g$n_family,
          isTRUE(cal$reduction$replay_check), !is.null(cal$screened),
          length(cal$screened) == fit$grp.consistency$n_candidates_evaluated)
z_scr <- cal$z_pstar   # the p* = 0.90 screen's cutoff as the package applies it (its print's "T >= ...")
stopifnot(length(z_scr) == 1L, is.finite(z_scr), identical(z_scr, cal_c2$z_pstar))
cal

## ==== chunk 'declcal-table' ====
# The family's identifiers are the search's factor codes: q<i>.1 is the i-th evaluated factor, q<i>.0 its negation.
fam_def <- function(id) vapply(strsplit(id, " & ", fixed = TRUE), function(p) {
  i <- as.integer(sub("^q(\\d+)\\.[01]$", "\\1", p))
  paste0(ifelse(endsWith(p, ".0"), "!", ""), "{", fit$confounders.evaluated[i], "}", collapse = " & ")
}, "")
fld  <- g$declaration_field
id_H <- fld$family_id[fld$meta$selected_index]
stopifnot(identical(fam_def(id_H), label_fs))
T_by_candidate <- data.frame(id = names(cal$T_hat), definition = fam_def(names(cal$T_hat)),
                             N = unname(fld$meta$sizes[names(cal$T_hat)]),
                             beta_hat = unname(cal$beta_hat), sigma_D = unname(cal$sigma_D),
                             T = unname(cal$T_hat), screened = names(cal$T_hat) %in% cal$screened)
i_max <- which.max(T_by_candidate$T)
max_T <- as.list(T_by_candidate[i_max, c("id", "definition", "N", "T", "screened")])
T_H   <- cal$T_hat[[id_H]]

# The applied admitted set: the calibrated set over the pre-reduction family, intersected with the
# candidates this document's own screen evaluated. Intersecting can only remove declarations.
adm     <- c(cal$c0$admitted_calibrated, list(c2 = cal_c2$admitted_calibrated))
applied <- lapply(adm, intersect, cal$screened)
tb      <- cal$c0$table
kappa   <- c(tb$kappa_hat, cal_c2$kappa_hat)
# The settable-pair columns, all of them, under the package's names and in its order, for the c0 rows. The c2 row
# comes from the c0 = NULL call, which returns none of them, so its cells are left blank.
st_nm <- c("pstar_settable", "pstar_achievable", "pcons_eff_settable", "z_eff_settable", "z_gap",
           "digits_fine", "pstar_fine", "z_gap_fine")
stopifnot(identical(names(tb)[names(tb) %in% st_nm], st_nm))
st     <- tb[st_nm]
st_fmt <- data.frame(
  pstar_settable     = formatC(st$pstar_settable, format = "f", digits = cal$digits),
  pstar_achievable   = as.character(st$pstar_achievable),
  pcons_eff_settable = formatC(st$pcons_eff_settable, format = "f", digits = cal$digits + 1L),
  z_eff_settable     = sprintf("%.3f", st$z_eff_settable),
  z_gap              = sprintf("%.3f", st$z_gap),
  digits_fine        = sprintf("%d", st$digits_fine),
  pstar_fine         = vapply(seq_len(nrow(st)), function(i) if (is.na(st$digits_fine[i])) "NA" else
                                formatC(st$pstar_fine[i], format = "f", digits = st$digits_fine[i]), ""),
  z_gap_fine         = sprintf("%.4f", st$z_gap_fine))
st_fmt <- rbind(st_fmt, as.data.frame(as.list(setNames(rep("", length(st_nm)), st_nm))))
stopifnot(identical(names(st_fmt), st_nm))
declcal_table <- data.frame(
  c0       = c(paste0(fmt_c0(tb$c0), ifelse(tb$c0 == c0_primary, " (primary)", "")),
               paste0(fmt_c0(c2_nat), ", the c2 boundary (c0 = c2)")),
  kappa    = sprintf("%.3f", kappa),
  st_fmt,
  fw       = sprintf("%.3f", c(tb$fw_size, cal_c2$fw_size)),
  admitted = sprintf("%d (%d)", lengths(applied), lengths(adm)),
  H_cal    = vapply(applied, function(a) if (length(a)) paste(describe_subgroup(fam_def(a)), collapse = "; ")
                    else "none", ""))
knitr::kable(declcal_table, row.names = FALSE, align = paste0("lr", strrep("r", length(st_nm)), "rrl"),
             col.names = c("c0 (MD, cells/mm³)", "κ̂₀.₁₀", paste0("`", st_nm, "`"), "FŴ₀.₁₀ of the p\\* = 0.90 screen at c0",
                           "candidates the calibrated rule admits", "calibrated Ĥ"),
             caption = sprintf("Declaration calibration at α = 0.10. Admitted, k (m): k in the applied set, the calibrated admitted set intersected with the %d candidates this document's screen evaluated; m in the calibrated admitted set over the pre-reduction family of %d.",
                               length(cal$screened), cal$n_family_prereduction))

## ==== chunk 'declcal-observed' ====
obs_tab <- data.frame(c0 = declcal_table$c0, kappa = declcal_table$kappa,
                      gap_max = sprintf("%.3f", kappa - max_T$T), gap_H = sprintf("%.3f", kappa - T_H))
knitr::kable(obs_tab, row.names = FALSE, align = "lrrr",
             col.names = c("c0 (MD, cells/mm³)", "κ̂₀.₁₀", "κ̂ − max T", "κ̂ − T(Ĥ)"),
             caption = sprintf("Maximum over the family: T = %.3f, %s (N %d)%s. The declared subgroup, %s (N %d): T = %.3f. The p* = 0.90 screen's cutoff is %.3f.",
                               max_T$T, check_description(max_T$definition, max_T$N, df), max_T$N,
                               if (max_T$screened) ", among the candidates the screen evaluated" else ", a candidate the screen did not evaluate",
                               describe_subgroup(label_fs), n_H, T_H, z_scr))

## ==== chunk 'declcal-decision' ====
# The readings below are written for an empty applied set at every c0; the document stops otherwise.
if (any(lengths(applied) > 0)) stop("the calibrated rule admits a candidate: the calibrated-H section must be rewritten")
stopifnot(max_T$T < min(kappa), identical(max_T$id, names(which.max(cal$T_hat[cal$screened]))))
i_p <- which(tb$c0 == c0_primary)

