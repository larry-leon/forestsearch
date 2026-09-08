# Re-selection map check (Stage 1b): the observed Hhat is reproduced by the gate's own
# selection map applied to the UNPERTURBED candidate effects.
suppressPackageStartupMessages(library(forestsearch))
x <- readRDS(file.path(Sys.getenv("SP"), "rep1_field.rds")); fs <- x$fs.est; g <- fs$mr_inference
cat("args_call_all names:", paste(names(fs$args_call_all), collapse = " "), "\n")
A <- fs$args_call_all
# Locate the cut matrix Z: a 0/1 matrix or frame, nrow = n, whose column names include the
# selected label's terms (e.g. "q25.0").  Walk the object recursively.
.terms_sel <- strsplit(g$selected_label, " & ", fixed = TRUE)[[1]]
# Z exactly as forestsearch_main.R:2933-2940 builds it: dummy() over the screened
# two-level cut factors q1..qL of the estimation frame (q<k>.0 / q<k>.1 columns).
qcols <- grep("^q[0-9]+$", names(fs$df.est), value = TRUE)
dfc <- fs$df.est[, qcols, drop = FALSE]
for (k in qcols) dfc[[k]] <- factor(dfc[[k]], levels = c(0, 1))
dfc <- forestsearch:::dummy(dfc)
Z <- as.matrix(dfc); colnames(Z) <- names(dfc); storage.mode(Z) <- "integer"
fz <- list(Z = Z, path = sprintf("dummy(df.est[, q1..q%d])", length(qcols)))
cat("Z found at", fz$path, ":", paste(dim(Z), collapse = "x"), "| n.min:", A$n.min, "| maxk:", A$maxk, "\n")
cat("df.est:", paste(dim(fs$df.est), collapse = "x"), "| has y_sim/treat_sim:", all(c("y_sim", "treat_sim") %in% names(fs$df.est)), "\n")
cat("admission: "); str(fs$admission, max.level = 3)
cat("settings: "); str(g$settings)
# --- rebuild the family exactly as forestsearch_main.R:3351-3365 does ------------------
L <- ncol(Z); maxk <- A$maxk; n.min <- A$n.min
combo <- forestsearch:::generate_combination_indices(L, maxk)
tot   <- forestsearch:::calculate_max_combinations(L, maxk)
fam <- list()
for (kk in seq_len(tot)) {
  covs.in <- forestsearch:::get_covs_in(kk, maxk, L, combo$counts_1, combo$indices_1,
                                        combo$counts_2, combo$indices_2, combo$counts_3, combo$indices_3)
  k_sel <- sum(covs.in); if (k_sel < 1L || k_sel > maxk) next
  mem <- which(forestsearch:::get_subgroup_membership(Z, covs.in))
  if (length(mem) >= n.min) fam[[paste(colnames(Z)[covs.in == 1], collapse = " & ")]] <- mem
}
cat(sprintf("family rebuilt: %d candidates (gate reported n_family = %d)\n", length(fam), g$n_family))
obs <- which(fs$grp.consistency$sg.harm.id == 1L)
cat(sprintf("observed Hhat: |Hhat| = %d; selected_label '%s' in family: %s; its members == observed (as sets): %s\n",
            length(obs), g$selected_label, g$selected_label %in% names(fam),
            setequal(fam[[g$selected_label]], obs)))
# --- assemble the effects on the gate's frame and spec ---------------------------------
df_g <- if (all(c("y_sim", "treat_sim") %in% names(fs$df.est))) fs$df.est else x$df
gspec <- list(outcome_type = "continuous", effect_measure = "MD", treat.name = "treat_sim",
              outcome.name = "y_sim", event.name = NULL, offset.name = NULL,
              adjust_covariates = NULL, adverse_outcome = FALSE)
asm <- forestsearch:::.fs_mr_assemble(df_g, fam, gspec)
cat(sprintf("assembled: %d estimable candidates; log_scale = %s\n", length(asm$names), asm$log_scale))
bh <- asm$beta_hat; sdv <- asm$sigma_D; sz <- asm$sizes
sel <- match(g$selected_label, asm$names)
cat(sprintf("selected index in assembly: %d (gate: %d) | beta_hat[sel] = %.10f vs gate naive est %.10f | identical to 1e-10: %s\n",
            sel, g$selected_index, bh[sel], g$naive$est, abs(bh[sel] - g$naive$est) < 1e-10))
# --- the gate's admission set and selection rule on the UNPERTURBED effects -------------
adm <- fs$admission
has_e <- !is.null(adm$effect_floor); has_c <- !is.null(adm$consistency)
c_cons <- if (has_c) adm$consistency$c_cons else NULL
t_g <- if (has_e && has_c) pmax(adm$effect_floor, c_cons + qnorm((1 + adm$consistency$p_star) / 2) * sdv) else if (has_e) adm$effect_floor else NULL
pass <- if (is.null(t_g)) seq_along(bh) else which(bh >= t_g)
zc <- if (has_c) (bh - c_cons) / sdv else NULL
st <- g$settings
s0 <- forestsearch:::.fs_mr_select(bh, zc, sz, pass, st$reselection, 0.10, st$selection_rule, asm$log_scale)
cat(sprintf("admission: effect_floor %s, consistency floor %s (p_star %s) -> %d passers of %d\n",
            format(adm$effect_floor), format(c_cons), format(adm$consistency$p_star), length(pass), length(bh)))
cat(sprintf("RE-SELECTION MAP on the unperturbed effects: S(beta_hat) = %d ('%s'); gate's selected_index = %d ('%s') -> %s\n",
            s0, asm$names[s0], g$selected_index, g$selected_label, if (identical(s0, sel) && identical(asm$names[s0], g$selected_label)) "REPRODUCED" else "MISMATCH"))
cat(sprintf("members of S(beta_hat) == observed Hhat (as sets): %s\n", setequal(fam[[asm$names[s0]]], obs)))
cat(sprintf("p_hat(Hhat) = %.4f | top-3: %s\n", g$reselection$p_hat[g$selected_index],
            paste(sprintf("%s=%.3f", names(sort(g$reselection$p_hat, decreasing = TRUE))[1:3], sort(g$reselection$p_hat, decreasing = TRUE)[1:3]), collapse = ", ")))
