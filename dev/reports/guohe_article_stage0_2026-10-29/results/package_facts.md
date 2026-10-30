# Package facts (installed forestsearch 0.4.1; last replicate t7_beta2_03 m = 2)

## Element names

- `names(mr)`: `selected_index`, `selected_label`, `measure`, `log_scale`, `ci_method`, `naive`, `debiased`, `selection_bias`, `fixed_bias`, `selection_rate`, `mean_r`, `mean_r_c`, `complement`, `settings`, `harm_flag`, `n_family`, `n_selected`, `timing_seconds`, `reselection`, `field`, `ij_residual`
- `names(mr$debiased)`: `est`, `lower`, `upper`, `lower_1s`, `se`, `se_ij`, `se_wald`, `var_ij`, `ij_source`, `ij_draws`, `se_ij_two_term`, `se_ij_winner`, `ij_source_winner`, `se_ij_winner_floor`, `ij_source_winner_floor`, `lower_w`, `upper_w`, `lower_1s_w`, `lower_wf`, `upper_wf`, `lower_1s_wf`
- `names(mr$field)`: `lambda_mean`, `lambda_sd`, `q05`, `q25`, `q50`, `q75`, `q95`, `q025`, `q975`, `n_out_used`, `n_in_used_mean`, `est2`, `lower_1s`, `lower_2s`, `upper_2s`, `se_field`, `lower_se`, `upper_se`, `R_out`, `R_in`, `seed_offset`, `timing_seconds`, `complement`, `joint`, `joint_s`
- `names(mr$field$complement)`: `lambda_mean`, `lambda_sd`, `q05`, `q25`, `q50`, `q75`, `q95`, `q025`, `q975`, `n_in_used_mean`, `est2`, `upper_1s`, `lower_1s`, `lower_2s`, `upper_2s`, `se_field`, `lower_se`, `upper_se`, `n_out_used`, `n_out_dropped_unfit`, `n_complement_fits`, `n_new_fits`, `share_draws_new_fit`, `R_out`, `R_in`, `timing_seconds`, `lambda_mean_s`, `se_field_s`, `est2_s`, `upper_1s_s`, `lower_1s_s`, `lower_2s_s`, `upper_2s_s`, `lower_se_s`, `upper_se_s`
- `names(mr$field$joint_s)`: `gamma`, `joint_prob`, `alpha`, `lower_H`, `upper_Hc`, `bonf_gamma`, `bonf_lower_H`, `bonf_upper_Hc`, `bonf_joint_prob`, `corr`, `n_joint_draws`, `grid_gamma`, `grid_joint_prob`
- `names(mr$complement)`: `naive`, `debiased`, `selection_bias`, `fixed_bias`, `n`
- `names(mr$reselection)`: `winner`, `p_hat`
- `names(fit)`: `outcome`, `scale`, `orient`, `selected`, `naive`, `debiased`, `bound_one_sided`, `bound_side`, `ci_two_sided`, `ci_wald`, `n`, `B`, `B_requested`, `r`, `level`, `n_candidates`, `n_estimable`, `fast`, `lean_cox`, `adjust_covariates`, `drop_rate`

## Defaults: fs_mr_inference()

- `draws = 2000L`
- `multiplier = c("poisson", "gaussian", "rademacher")`
- `ci_method = c("field", "ij", "wald")`
- `include_complement = FALSE`
- `field_R_out = 1000L`
- `field_R_in = 500L`
- `field_complement = TRUE`
- `field_scale_complement = c("selected", "none")`
- `reselection = c("maxcons", "maxeff", "maxSG", "minSG", "effMaxSG", "effMinSG")`

## Defaults: guohe_algorithm3()

- `orient = -1`
- `B = 2000L`
- `r = 0.03`
- `level = 0.05`
- `min_events = 5L`
- `fast = NULL`
- `parallel = FALSE`
