# Does GRF's proposed family depend on the outcome?  Compares the matched prevalence
# pairs (12.4% vs 31% at n 500 and n 1500) in the Part B probe bundles.
# The candidate LISTS are not stored per replicate -- the bundles hold only n_family,
# the selected rule and p_hat_top_labels -- so this reports the strongest condition the
# bundles support (per-replicate size equality) beside the selected rule, and does NOT
# add a recorder change.  The source settles the rest: .grf_dr_candidates()
# (grf_subgroup_labels.R:255-277) admits on X and n_min alone.
setwd("/Users/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020")
g <- function(n, z) readRDS(sprintf(
  "results/grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n%d%s_nb20_grfprobe_res_1_36.rds",
  n, if (z) "_z1q60" else ""))$results
for (n in c(500L, 1500L)) {
  a <- g(n, FALSE); b <- g(n, TRUE)
  a <- a[order(a$sim_id), ]; b <- b[order(b$sim_id), ]
  cat(sprintf("\n===== matched pair at n = %d : 12.4%% vs 31%% =====\n", n))
  cat(sprintf("  sim_id vectors identical                : %s\n", identical(a$sim_id, b$sim_id)))
  cat(sprintf("  n_true (same draws?) identical          : %s\n", identical(a$n_true, b$n_true)))
  cat(sprintf("  n_family identical on all %d replicates : %s  (agree on %d of %d)\n",
      nrow(a), identical(a$n_family, b$n_family), sum(a$n_family == b$n_family), nrow(a)))
  d <- a$n_family - b$n_family
  cat(sprintf("  n_family max |difference|               : %d\n", max(abs(d))))
  cat(sprintf("  SELECTED rule identical                 : %s  (agree on %d of %d)\n",
      identical(a$sg_def, b$sg_def), sum(a$sg_def == b$sg_def, na.rm = TRUE), nrow(a)))
  cat(sprintf("  selected |Hhat| (n_sel) identical       : %s\n", identical(a$n_sel, b$n_sel)))
  cat(sprintf("  p_hat_top_labels identical              : %s  (agree on %d of %d)\n",
      identical(a$p_hat_top_labels, b$p_hat_top_labels),
      sum(a$p_hat_top_labels == b$p_hat_top_labels, na.rm = TRUE), nrow(a)))
  cat("  first four replicates, selected rule side by side:\n")
  for (i in 1:4) cat(sprintf("    sim %2d | K %4d vs %4d | 12.4%%: %-46s | 31%%: %s\n",
      a$sim_id[i], a$n_family[i], b$n_family[i], a$sg_def[i], b$sg_def[i]))
}
