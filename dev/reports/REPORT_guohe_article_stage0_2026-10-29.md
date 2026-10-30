# REPORT — Guo & He (S4) article, Stage 0: reproduction gate and timing ladder

Task: `dev/tasks/cc_task_guohe_article_stage0_2026-10-29_v2.md` (supersedes v1). Run 2026-10-30, single core.

## 1. Provenance

- MASTER_SHA: `c2967362e7ee9a6361cb1b07fdf83c0098e712db`
- Branch: `feature/articles-guohe`
- v1 stop: `e751c215` (v1 task doc). v1 halted at its §1.5 because devtools 2.5.2 rejects `upgrade = "never"` (`upgrade` must be TRUE, FALSE or NA). v1's logs `logs/install.log` and `logs/s0_env.log` are kept.
- This document's task commit: `462ca74c` (v2 task doc, alone). Filename received `~/Downloads/cc_task_guohe_article_stage0_2026-10-29_v2.md`; committed as `dev/tasks/cc_task_guohe_article_stage0_2026-10-29_v2.md`.
- Install route: `Rscript -e 'devtools::install(".", upgrade = FALSE, build_vignettes = FALSE)'` → `logs/install_v2.log`, exit 0, `* DONE (forestsearch)`. The `R CMD INSTALL` fallback was not needed.
- §1.2 `git status --porcelain` (verbatim): `?? dev/reports/guohe_article_stage0_2026-10-29/` (expected; v1's folder). No other dirt.

From `results/keyvalues.txt` (verbatim):

```
task_start_utc = 2026-10-30 16:08:27 
pd_Version = 0.4.1
pd_Built = R 4.6.1; ; 2026-10-30 16:08:45 UTC; unix
pd_Packaged = 2026-10-30 16:08:42 UTC; larryleon
library_path = /home/larryleon/R/x86_64-pc-linux-gnu-library/4.6
survival_version = 3.8.12
R_version = R version 4.6.1 (2026-06-24)
sysname = Linux
machine = x86_64
availableCores = 128
gate_1_6_version_0.4.1 = TRUE
gate_1_6_built_after_start = TRUE (Built 2026-10-30 16:08:45 UTC vs start 2026-10-30 16:08:27)
gate_1_7_exports = TRUE
```

Gate outcomes:
- §1.1 branch/history: PASS (`e751c215` directly on `c2967362`).
- §1.6 installed build: PASS (`pd$Version == "0.4.1"`; Built 16:08:45 UTC is later than task_start_utc 16:08:27).
- §1.7 exports: PASS (`fs_mr_inference`, `guohe_algorithm3` both exported).
- §2 call-identity gate: PASS. `mv_mr()` calls `forestsearch:::fs_mr_inference` with `reselection = "maxeff"` and `admission = list(effect_floor = NULL, consistency = NULL)` (lines 135–136, 151). `gh52_one_rep()` calls `guohe_algorithm3()` with `orient = +1`, `diagnostics = FALSE`, `seed = boot_seed`, and `boot_seed <- as.integer(seed) + 500000L` (lines 199, 209–210). `level` and `min_events` are passed through from `gh52_one_rep()` formals that default to `0.05` and `5L` (lines 190–191). The production driver `guohe_sec52_run.R:167` passes `level = 0.05, min_events = 5L` explicitly.
- §3 reproduction: PASS (see §2 below).
- §4 truth curve at n_big = 2e5: gates PASS. At 5e4: gates PASS.
- §6 post-conditions: PASS (see closing block).

### §2 source quotes (current tree, line-numbered)

N1 standard as the complement driver states it: `MFC_DISCRETE` (lines 322–333) is "identical() REQUIRED"; floats use `all.equal()` at `MFC_TOL <- 1e-8` (line 350, applied at line 384).

`R/fs_mr_inference.R` lines 584–609:

```r
 584  fs_mr_inference <- function(df, candidates, spec, selected_members,
 585                             admission,
 586                             t_confirm = NULL, confirm_rule = c("point", "ci"),
 587                             reselection = c("maxcons", "maxeff", "maxSG",
 588                                             "minSG", "effMaxSG", "effMinSG"),
 589                             effect_neighborhood = 0.10,
 590                             selection_rule = c("neighborhood", "pareto", "both"),
 591                             draws = 2000L,
 592                             multiplier = c("poisson", "gaussian", "rademacher"),
 593                             include_complement = FALSE,
 594                             ci_method = c("field", "ij", "wald"),
 595                             seed = NULL,
 596                             return_reselection = TRUE,
 597                             field_R_out = 1000L,
 598                             field_R_in = 500L,
 599                             field_uniform = FALSE,
 600                             field_M_cap = NULL,
 601                             field_complement = TRUE,
 602                             field_decompose = FALSE,
 603                             field_scale_complement = c("selected", "none"),
 604                             ij_residual = c("two_term", "winner", "winner_floor"),
 605                             field_recovery = FALSE,
 606                             keep_declaration_field = FALSE,
 607                             keep_field_matrix = FALSE,
 608                             declaration_c0 = NULL,
 609                             pconsistency.digits = NULL) {
```

`R/fs_mr_inference.R` lines 661–679:

```r
 661    .has_effect <- !is.null(admission$effect_floor)
 662    .has_cons   <- !is.null(admission$consistency)
 663    c_cons <- if (.has_cons) admission$consistency$c_cons else NULL
 664  
 665    if (.has_effect && .has_cons) {
 666      # The search admits on round(Pcons, digits) >= p_star, i.e. on the
 667      # effective threshold .fs_pcons_eff(p_star, digits), not on p_star itself;
 668      # MR's floor uses the same rule so its re-selection map is the one that ran.
 669      # Applies on both consistency paths (under "split" the closed form was
 670      # already the approximation; only the rounding is corrected).
 671      digits <- pconsistency.digits
 672      if (is.null(digits))
 673        digits <- eval(formals(subgroup.consistency)$pconsistency.digits)
 674      z   <- stats::qnorm((1 + .fs_pcons_eff(admission$consistency$p_star,
 675                                             as.integer(digits))) / 2)
 676      t_g <- pmax(admission$effect_floor, c_cons + z * sdv)
 677      .admit <- function(bs) which(bs >= t_g)
 678    } else if (.has_effect) {
 679      t_g <- admission$effect_floor
```

`R/guohe_algorithm3.R` lines 232–248:

```r
 232  guohe_algorithm3 <- function(data,
 233                               outcome = c("survival", "binary", "continuous"),
 234                               treatment,
 235                               candidates,
 236                               time = NULL,
 237                               event = NULL,
 238                               y = NULL,
 239                               orient = -1,
 240                               B = 2000L,
 241                               r = 0.03,
 242                               level = 0.05,
 243                               seed = NULL,
 244                               min_events = 5L,
 245                               diagnostics = TRUE,
 246                               parallel = FALSE,
 247                               adjust_covariates = NULL,
 248                               fast = NULL) {
```

`quarto/GuoHe/mr_vs_guohe_sim.R` lines 60–62:

```r
  60  MV_DRAWS      <- 5000L        # D3 default
  61  MV_MULTIPLIER <- "poisson"    # D3 default (centred Poisson)
  62  MV_SEED_MR    <- 700000L      # MR seed offset from the replicate data seed
```

`quarto/GuoHe/mr_vs_guohe_sim.R` lines 79–79:

```r
  79  mv_gh52_base <- function(id) 1000000L + as.integer(sum(utf8ToInt(id))) * 100003L
```

`quarto/GuoHe/mr_vs_guohe_sim.R` lines 93–98:

```r
  93  mv_gh52_regen <- function(id, m) {
  94    b2 <- as.integer(sub("^t7_beta2_", "", id)) / 10
  95    set.seed(mv_gh52_base(id) + m)
  96    df0 <- gh52_sim_data(b2, 400L)
  97    gh52_candidates(df0)   # list(df, cuts, names)
  98  }
```

`quarto/GuoHe/mr_vs_guohe_sim.R` lines 107–110:

```r
 107  mv_cand_idx_52 <- function(cand) {
 108    stats::setNames(lapply(cand$names, function(nm) which(cand$df[[nm]] == 1L)),
 109                    cand$names)
 110  }
```

`quarto/GuoHe/mr_vs_guohe_sim.R` lines 114–117:

```r
 114  mv_spec52 <- list(outcome_type = "survival", effect_measure = "HR",
 115                    treat.name = "treat", outcome.name = "time",
 116                    event.name = "event", offset.name = NULL,
 117                    adjust_covariates = NULL, adverse_outcome = NULL)
```

`quarto/GuoHe/mr_vs_guohe_sim.R` lines 132–151:

```r
 132    args <- list(
 133      df = df, candidates = cands, spec = spec,
 134      selected_members = cands[[sel_label]],
 135      admission = list(effect_floor = NULL, consistency = NULL),
 136      reselection = "maxeff",
 137      draws = draws, multiplier = multiplier,
 138      ci_method = ci_method, seed = seed, return_reselection = TRUE,
 139      include_complement = include_complement,
 140      # winner-only IJ variants (TASK_complement_refinements_2026-09-06);
 141      # "two_term" keeps the reported interval byte-identical.
 142      ij_residual = ij_residual)
 143    if (identical(ci_method, "field"))
 144      args <- c(args, list(field_R_out = field_R_out, field_R_in = field_R_in,
 145                           # kappa(Sigma-hat) sweep (TASK_mr_field_uniform_2026-09-05);
 146                           # FALSE keeps the 2026-09-05 field output byte-identical.
 147                           field_uniform = field_uniform,
 148                           # complement field (TASK_mr_field_complement_2026-09-06);
 149                           # FALSE keeps the field output byte-identical.
 150                           field_complement = field_complement))
 151    do.call(forestsearch:::fs_mr_inference, args)
```

`quarto/GuoHe/guohe_sec52_sim.R` lines 85–85:

```r
  85  GH52_R_GRID <- c(1 / 3, 1 / 12, 1 / 21, 1 / 30)
```

`quarto/GuoHe/guohe_sec52_sim.R` lines 116–123:

```r
 116  gh52_candidates <- function(df) {
 117    cuts <- c(GH52_C_LO, sort(df$w[df$w > GH52_C_LO & df$w <= GH52_C_HI]))
 118    nm <- sprintf("cand_%03d", seq_along(cuts))
 119    memb <- vapply(cuts, function(cc) as.integer(df$w <= cc),
 120                   integer(nrow(df)))
 121    colnames(memb) <- nm
 122    list(df = cbind(df, as.data.frame(memb)), cuts = cuts, names = nm)
 123  }
```

`quarto/GuoHe/guohe_sec52_sim.R` lines 132–146:

```r
 132  gh52_subgroup_fits <- function(df, cand) {
 133    k <- length(cand$names)
 134    est <- se <- rep(NA_real_, k)
 135    for (j in seq_len(k)) {
 136      sub <- df[df[[cand$names[j]]] == 1L, , drop = FALSE]
 137      fit <- try(suppressWarnings(
 138        survival::coxph(survival::Surv(time, event) ~ treat, data = sub)
 139      ), silent = TRUE)
 140      if (!inherits(fit, "try-error")) {
 141        est[j] <- unname(fit$coefficients[[1L]])
 142        se[j] <- unname(sqrt(diag(stats::vcov(fit)))[1L])
 143      }
 144    }
 145    list(est = est, se = se)
 146  }
```

`quarto/GuoHe/guohe_sec52_sim.R` lines 149–162:

```r
 149  gh52_naive <- function(fits, cand, truth, level = 0.05) {
 150    ok <- is.finite(fits$est)
 151    sel <- which.max(replace(fits$est, !ok, -Inf))
 152    c_hat <- cand$cuts[sel]
 153    gamma_s <- gh52_truth_at(truth, c_hat)
 154    point <- fits$est[sel]
 155    lower <- point - stats::qnorm(1 - level) * fits$se[sel]
 156    list(
 157      sel = sel, c_hat = c_hat, gamma_s = gamma_s,
 158      point = point, lower = lower,
 159      cover = as.integer(lower <= gamma_s),
 160      dist = gamma_s - lower, bias = point - gamma_s
 161    )
 162  }
```

`quarto/GuoHe/guohe_sec52_sim.R` lines 189–191:

```r
 189  gh52_one_rep <- function(beta2, n = 400L, truth, r_grid = GH52_R_GRID,
 190                           B = 2000L, level = 0.05, min_events = 5L,
 191                           seed = NULL) {
```

`quarto/GuoHe/guohe_sec52_sim.R` lines 199–211:

```r
 199    boot_seed <- if (is.null(seed)) NULL else as.integer(seed) + 500000L
 200  
 201    sel_name <- NA_character_
 202    c_hat_gh <- gamma_s <- NA_real_
 203    n_sel <- NA_integer_
 204    eng <- list()
 205    for (i in seq_along(r_grid)) {
 206      fit <- try(suppressWarnings(guohe_algorithm3(
 207        data = cand$df, outcome = "survival", treatment = "treat",
 208        candidates = cand$names, time = "time", event = "event",
 209        orient = +1, B = B, r = r_grid[i], level = level, seed = boot_seed,
 210        min_events = min_events, diagnostics = FALSE
 211      )), silent = TRUE)
```

`quarto/GuoHe/guohe_sec52_truth.R` lines 72–84:

```r
  72  gh52_sim_data <- function(beta2, n) {
  73    w <- stats::runif(n, 0, GH52_W_MAX)
  74    d <- stats::rbinom(n, 1L, 0.5)
  75    b <- ifelse(w <= GH52_C_LO, beta2, 0)
  76    t_event <- stats::rexp(n, rate = exp(b * d))
  77    t_cens <- exp(stats::runif(n, -1.25, 1.00))
  78    data.frame(
  79      time = pmin(t_event, t_cens),
  80      event = as.integer(t_event <= t_cens),
  81      treat = d,
  82      w = w
  83    )
  84  }
```

`quarto/GuoHe/guohe_sec52_truth.R` lines 147–151:

```r
 147  gh52_truth_curve <- function(beta2, n_big = 2e6, c_step = 0.25, seed,
 148                               anchor_k = 4, null_k = 4,
 149                               tol_mono = 0.02, tol_total_rise = 0.01,
 150                               drop_frac = 0.15, cens_tol = 0.005,
 151                               verbose = TRUE) {
```

`quarto/GuoHe/guohe_sec52_truth.R` lines 300–300:

```r
 300  gh52_truth_at <- function(truth, c_hat, use = c("auto", "exact", "smooth", "raw")) {
```

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R` lines 108–133:

```r
 108  mfc_mr <- function(df, cands, sel_label, spec, draws = MV_DRAWS,
 109                     multiplier = MV_MULTIPLIER, seed = NULL,
 110                     ci_method = "field", field_R_out = MF_R_OUT,
 111                     field_R_in = MF_R_IN,
 112                     field_uniform = FALSE, ij_residual = "two_term") {
 113    args <- list(
 114      df = df, candidates = cands, spec = spec,
 115      selected_members = cands[[sel_label]],
 116      admission = list(effect_floor = NULL, consistency = NULL),
 117      reselection = "maxeff",
 118      draws = draws, multiplier = multiplier,
 119      ci_method = ci_method, seed = seed, return_reselection = TRUE,
 120      include_complement = TRUE,                    # A3 addition (was FALSE)
 121      # winner-only IJ variants (TASK_complement_refinements_2026-09-06);
 122      # "two_term" keeps the reported interval byte-identical.
 123      ij_residual = ij_residual)
 124    if (identical(ci_method, "field"))
 125      args <- c(args, list(field_R_out = field_R_out, field_R_in = field_R_in,
 126                           # kappa(Sigma-hat) sweep (TASK_mr_field_uniform_2026-09-05);
 127                           # FALSE keeps the 2026-09-05 field output byte-identical.
 128                           field_uniform = field_uniform,
 129                           # complement field (TASK_mr_field_complement_2026-09-06);
 130                           field_complement = TRUE,              # A3 addition
 131                           field_scale_complement = "selected")) # A3 addition
 132    do.call(forestsearch:::fs_mr_inference, args)
 133  }
```

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R` lines 147–166:

```r
 147    f  <- mr$field
 148    fc <- f$complement
 149    jt <- f$joint
 150    js <- f$joint_s
 151    gc <- mr$complement
 152    e2 <- log(f$est2); flo1 <- log(f$lower_1s)
 153    flo2 <- log(f$lower_2s); fup2 <- log(f$upper_2s)
 154  
 155    # --- the two 0.975 one-sided bounds and the both-correct indicator (A2) ---
 156    # Engine's own Bonferroni pair at gamma = alpha/2 = 0.025. The harm side is
 157    # identical in `joint` and `joint_s` (only the complement draws differ), so
 158    # bonf_lower_H is quoted once.
 159    h_lo_975  <- .lg(jt$bonf_lower_H)
 160    c_up_975  <- .lg(jt$bonf_upper_Hc)
 161    c_up_975s <- .lg(js$bonf_upper_Hc)
 162    h_lo_975s <- .lg(js$bonf_lower_H)
 163    # 0.95 one-sided: harm lower (field) and complement upper (field / field-s).
 164    h_lo_95   <- flo1
 165    c_up_95   <- .lg(fc$upper_1s)
 166    c_up_95s  <- .lg(fc$upper_1s_s)
```

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R` lines 197–201:

```r
 197      # ---- NEW: complement upper at both levels, unscaled and field-s --------
 198      # Truth for the complement is exactly 0 (see header), so coverage is
 199      # upper >= 0.
 200      c_upper_95 = c_up_95, c_upper_975 = c_up_975,
 201      c_upper_95_s = c_up_95s, c_upper_975_s = c_up_975s,
```

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R` lines 212–216:

```r
 212      joint_gamma = .nz(jt$gamma), joint_gamma_s = .nz(js$gamma),
 213      joint_prob = .nz(jt$joint_prob), joint_prob_s = .nz(js$joint_prob),
 214      joint_bonf_prob = .nz(jt$bonf_joint_prob),
 215      joint_bonf_prob_s = .nz(js$bonf_joint_prob),
 216      joint_corr = .nz(jt$corr), joint_corr_s = .nz(js$corr),
```

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R` lines 228–229:

```r
 228      cgate_naive_est = .lg(gc$naive$est),
 229      cgate_est = .lg(gc$debiased$est),
```

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R` lines 322–333:

```r
 322  # ---- N1 pairing standard (v4; replaces A3's plain identical()) -------------
 323  # Integer / selection / flag / seed columns: identical() REQUIRED, every
 324  # replicate; any failure is a STOP. Floating-point columns: all.equal() at
 325  # 1e-8, with the worst absolute and relative deviation and the bit-identical
 326  # fraction recorded per cell. The cross-machine comparison is a provenance
 327  # measurement, not a gate: stored bundles are x86_64 / R 4.6.1 / reference
 328  # BLAS, this Mac is arm64 / R 4.5.2 / Accelerate.
 329  MFC_DISCRETE <- c("sel", "c_hat", "n_sel", "seed_data", "seed_mr", "m",
 330                    "k_family", "naive_cover", "mr_cover", "fld_cover_1s",
 331                    "fld_cover_2s", "fld_n_out_used", "mr_ij_source",
 332                    "mr_ij_draws",
 333                    paste0("gh_r", 1:4, "_cov"))
```

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R` lines 350–350:

```r
 350  MFC_TOL <- 1e-8
```

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R` lines 380–386:

```r
 380      if (max(rl) > worst_rel) {
 381        worst_rel <- max(rl)
 382        worst_rel_at <- sprintf("%s / m=%s", nm, res$m[ok][which.max(rl)])
 383      }
 384      if (!isTRUE(all.equal(a[ok], b[ok], tolerance = MFC_TOL)))
 385        ae_bad <- c(ae_bad, nm)
 386    }
```

## 2. R — the reproduction

Four replicates: `t7_beta2_00` and `t7_beta2_03`, m ∈ {1, 2}. Exported `guohe_algorithm3()` and `fs_mr_inference()` of the installed 0.4.1 build were compared against `guohe_repro_<id>.rds`, `mr_field_vs_guohe_<id>.rds` and `mr_field_complement_vs_guohe_<id>.rds`. Values shown to 12 significant digits; `results/reproduction.csv` has them at %.17g.

| id | m | quantity | kind | stored | recomputed | abs_diff | identical | all.equal 1e-8 | all.equal 1e-6 |
|---|---:|---|---|---:|---:|---:|---|---|---|
| t7_beta2_00 | 1 | naive_c_hat | key | 46.3580730371 | 46.3580730371 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | n_cand | key | 154 | 154 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | naive_point | float | 0.206942588869 | 0.206942588869 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | naive_lower | float | -0.0762030721295 | -0.0762030721295 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | gamma_s_naive | float | 0 | 0 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | gh_c_hat | key | 46.3580730371 | 46.3580730371 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | gh_r2_lower | float | -0.211230214267 | -0.211230214267 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | gh_r2_debiased | float | 0.0899929195506 | 0.0899929195506 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | mr_selected_index | key | 100 | 100 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | mr_est | float | 0.150710680413 | 0.150710680413 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | mr_se_ij | float | 0.327928729572 | 0.327928729572 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | mr_lower_1s | float | -0.388684079804 | -0.388684079804 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | mr_bias_sel | float | 0.0563307391407 | 0.0563307391407 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | mr_bias_fix | float | -9.88306851643e-05 | -9.88306851643e-05 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | fld_est2 | float | 0.156312037323 | 0.156312037323 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | fld_lower_1s | float | -0.130778193742 | -0.130778193742 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | fld_lambda_mean | float | -0.00560135690942 | -0.00560135690942 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | p_hat_H | float | 0.299 | 0.299 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | cgate_est | float | 0.0739390113232 | 0.0739390113232 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | c_upper_95_s | float | 0.443090992444 | 0.443090992444 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 1 | h_lower_975 (joint_s bonf_lower_H) | float | -0.19541111837 | -0.19541111837 | 2.78e-16 | FALSE | TRUE | TRUE |
| t7_beta2_00 | 1 | c_upper_975_s | float | 0.495227025415 | 0.495227025415 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | naive_c_hat | key | 30 | 30 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | n_cand | key | 132 | 132 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | naive_point | float | -0.0925560810303 | -0.0925560810303 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | naive_lower | float | -0.449403859508 | -0.449403859508 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | gamma_s_naive | float | 0 | 0 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | gh_c_hat | key | 30 | 30 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | gh_r2_lower | float | -0.506899812697 | -0.506899812697 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | gh_r2_debiased | float | -0.196383920547 | -0.196383920547 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | mr_selected_index | key | 1 | 1 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | mr_est | float | -0.19359344757 | -0.19359344757 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | mr_se_ij | float | 0.39382212575 | 0.39382212575 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | mr_lower_1s | float | -0.841373199483 | -0.841373199483 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | mr_bias_sel | float | 0.0977488285828 | 0.0977488285828 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | mr_bias_fix | float | 0.00328853795685 | 0.00328853795685 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | fld_est2 | float | -0.20777164501 | -0.20777164501 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | fld_lower_1s | float | -0.527785486747 | -0.527785486747 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | fld_lambda_mean | float | 0.0141781974401 | 0.0141781974401 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | p_hat_H | float | 0.2602 | 0.2602 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | cgate_est | float | 0.278022433664 | 0.278022433664 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_00 | 2 | c_upper_95_s | float | 0.613025778864 | 0.613025778864 | 2.22e-16 | FALSE | TRUE | TRUE |
| t7_beta2_00 | 2 | h_lower_975 (joint_s bonf_lower_H) | float | -0.602437037244 | -0.602437037244 | 3.33e-16 | FALSE | TRUE | TRUE |
| t7_beta2_00 | 2 | c_upper_975_s | float | 0.68917770913 | 0.68917770913 | 2.22e-16 | FALSE | TRUE | TRUE |
| t7_beta2_03 | 1 | naive_c_hat | key | 35.4605453089 | 35.4605453089 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | n_cand | key | 141 | 141 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | naive_point | float | 0.454507982456 | 0.454507982456 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | naive_lower | float | 0.123973588852 | 0.123973588852 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | gamma_s_naive | float | 0.251996750552 | 0.251996750552 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | gh_c_hat | key | 35.4605453089 | 35.4605453089 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | gh_r2_lower | float | 0.0236576354158 | 0.0236576354158 | 1.04e-17 | FALSE | TRUE | TRUE |
| t7_beta2_03 | 1 | gh_r2_debiased | float | 0.346222049348 | 0.346222049348 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | mr_selected_index | key | 17 | 17 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | mr_est | float | 0.394571009733 | 0.394571009733 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | mr_se_ij | float | 0.378324686202 | 0.378324686202 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | mr_lower_1s | float | -0.227717722531 | -0.227717722531 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | mr_bias_sel | float | 0.0627417496589 | 0.0627417496589 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | mr_bias_fix | float | -0.00280477693631 | -0.00280477693631 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | fld_est2 | float | 0.387244270049 | 0.387244270049 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | fld_lower_1s | float | 0.0793365466388 | 0.0793365466388 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | fld_lambda_mean | float | 0.00732673968363 | 0.00732673968363 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | p_hat_H | float | 0.2178 | 0.2178 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | cgate_est | float | -0.0185285384729 | -0.0185285384729 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | c_upper_95_s | float | 0.270239548533 | 0.270239548533 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 1 | h_lower_975 (joint_s bonf_lower_H) | float | 0.0142758914443 | 0.0142758914443 | 2.2e-16 | FALSE | TRUE | TRUE |
| t7_beta2_03 | 1 | c_upper_975_s | float | 0.297325895249 | 0.297325895249 | 1.67e-16 | FALSE | TRUE | TRUE |
| t7_beta2_03 | 2 | naive_c_hat | key | 49.9390554056 | 49.9390554056 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | n_cand | key | 153 | 153 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | naive_point | float | 0.414934900397 | 0.414934900397 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | naive_lower | float | 0.140839188712 | 0.140839188712 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | gamma_s_naive | float | 0.176908146399 | 0.176908146399 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | gh_c_hat | key | 49.9390554056 | 49.9390554056 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | gh_r2_lower | float | 0.000732496644086 | 0.000732496644086 | 2.17e-18 | FALSE | TRUE | TRUE |
| t7_beta2_03 | 2 | gh_r2_debiased | float | 0.307543149494 | 0.307543149494 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | mr_selected_index | key | 93 | 93 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | mr_est | float | 0.329397198936 | 0.329397198936 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | mr_se_ij | float | 0.318570180303 | 0.318570180303 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | mr_lower_1s | float | -0.194604117575 | -0.194604117575 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | mr_bias_sel | float | 0.0849145470079 | 0.0849145470079 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | mr_bias_fix | float | 0.000623154453164 | 0.000623154453164 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | fld_est2 | float | 0.315815623384 | 0.315815623384 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | fld_lower_1s | float | 0.0314575559452 | 0.0314575559452 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | fld_lambda_mean | float | 0.013581575552 | 0.013581575552 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | p_hat_H | float | 0.0686 | 0.0686 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | cgate_est | float | 0.284381318774 | 0.284381318774 | 0 | TRUE | TRUE | TRUE |
| t7_beta2_03 | 2 | c_upper_95_s | float | 0.652043799124 | 0.652043799124 | 3.33e-16 | FALSE | TRUE | TRUE |
| t7_beta2_03 | 2 | h_lower_975 (joint_s bonf_lower_H) | float | -0.0410060878714 | -0.0410060878714 | 2.36e-16 | FALSE | TRUE | TRUE |
| t7_beta2_03 | 2 | c_upper_975_s | float | 0.710598039161 | 0.710598039161 | 4.44e-16 | FALSE | TRUE | TRUE |

Gate lines (verbatim from `results/keyvalues.txt`):

```
s1_t7_beta2_00_seed_base_eq = TRUE (stored 92902757, base 92902757)
s1_t7_beta2_03_seed_base_eq = TRUE (stored 93202766, base 93202766)
s1_keys_identical = 16 / 16
s1_floats_identical = 61 / 72
s1_floats_all_equal_1e8 = 72 / 72
s1_floats_all_equal_1e6 = 72 / 72
s1_max_abs_diff = 4.44089e-16 (c_upper_975_s, t7_beta2_03 m=2)
s1_finding_fail_1e8_pass_1e6 = 0 floats
s1_GATE = PASS
```

**Verdict on R.** The 0.4.1 exported functions reproduce the S4 records to 1e-8 (all.equal), not bit-for-bit. All 16 selection keys are `identical()`. 61 of 72 floats are bit-identical. The 11 that are not differ by at most 4.44e-16, all in complement/joint quantities and the Guo & He r = 1/12 lower bound at β₂ = 0.3. No float fails at 1e-8, so there is no 1e-8/1e-6 finding. The article's code path is the manuscript's code path.

Recorded per replicate (verbatim):

```
s1_t7_beta2_00_m1_gh_lean_cox = TRUE
s1_t7_beta2_00_m1_gh_fast = FALSE
s1_t7_beta2_00_m1_gh_B = 2000
s1_t7_beta2_00_m1_gh_n_estimable = 154
s1_t7_beta2_00_m1_joint_s_bonf_joint_prob = 0.95099999999999996 (stored joint_bonf_prob_s 0.95099999999999996)
s1_t7_beta2_00_m1_joint_s_gamma = 0.025000000000000001 (stored joint_gamma_s 0.025000000000000001)
s1_t7_beta2_00_m1_joint_s_corr = -0.018804201961662412 (stored joint_corr_s -0.018804201961662444)
s1_t7_beta2_00_m1_sel_lab = cand_100
s1_t7_beta2_00_m1_sec_gh = 67.079999999999998
s1_t7_beta2_00_m1_sec_mr = 2.4500000000000002
s1_t7_beta2_00_m1_sec_rep = 69.989999999999995
s1_t7_beta2_00_m2_gh_lean_cox = TRUE
s1_t7_beta2_00_m2_gh_fast = FALSE
s1_t7_beta2_00_m2_gh_B = 2000
s1_t7_beta2_00_m2_gh_n_estimable = 132
s1_t7_beta2_00_m2_joint_s_bonf_joint_prob = 0.95099999999999996 (stored joint_bonf_prob_s 0.95099999999999996)
s1_t7_beta2_00_m2_joint_s_gamma = 0.025000000000000001 (stored joint_gamma_s 0.025000000000000001)
s1_t7_beta2_00_m2_joint_s_corr = 0.027578739914417862 (stored joint_corr_s 0.027578739914417873)
s1_t7_beta2_00_m2_sel_lab = cand_001
s1_t7_beta2_00_m2_sec_gh = 57.909999999999997
s1_t7_beta2_00_m2_sec_mr = 1.8100000000000001
s1_t7_beta2_00_m2_sec_rep = 60.039999999999999
s1_t7_beta2_03_m1_gh_lean_cox = TRUE
s1_t7_beta2_03_m1_gh_fast = FALSE
s1_t7_beta2_03_m1_gh_B = 2000
s1_t7_beta2_03_m1_gh_n_estimable = 141
s1_t7_beta2_03_m1_joint_s_bonf_joint_prob = 0.95199999999999996 (stored joint_bonf_prob_s 0.95199999999999996)
s1_t7_beta2_03_m1_joint_s_gamma = 0.026000000000000002 (stored joint_gamma_s 0.026000000000000002)
s1_t7_beta2_03_m1_joint_s_corr = -0.018346431214076576 (stored joint_corr_s -0.018346431214076545)
s1_t7_beta2_03_m1_sel_lab = cand_017
s1_t7_beta2_03_m1_sec_gh = 61.659999999999997
s1_t7_beta2_03_m1_sec_mr = 1.95
s1_t7_beta2_03_m1_sec_rep = 63.920000000000002
s1_t7_beta2_03_m2_gh_lean_cox = TRUE
s1_t7_beta2_03_m2_gh_fast = FALSE
s1_t7_beta2_03_m2_gh_B = 2000
s1_t7_beta2_03_m2_gh_n_estimable = 153
s1_t7_beta2_03_m2_joint_s_bonf_joint_prob = 0.95199999999999996 (stored joint_bonf_prob_s 0.95199999999999996)
s1_t7_beta2_03_m2_joint_s_gamma = 0.026000000000000002 (stored joint_gamma_s 0.026000000000000002)
s1_t7_beta2_03_m2_joint_s_corr = -0.040053356029274373 (stored joint_corr_s -0.040053356029274401)
s1_t7_beta2_03_m2_sel_lab = cand_093
s1_t7_beta2_03_m2_sec_gh = 66.439999999999998
s1_t7_beta2_03_m2_sec_mr = 2.1400000000000001
s1_t7_beta2_03_m2_sec_rep = 68.939999999999998
```

## 3. T — the ladder

Single core, cell `t7_beta2_03`, data seeds base + 101, 102, 103 (base = 93202766).

Per-replicate measurements (`results/timing.csv`):

| rep | component | rung | setting | sec |
|---:|---|---|---|---:|
| 101 | family | r1 | k = 163 | 0.47 |
| 101 | guohe | r1 | B = 2000 | 71.20 |
| 101 | guohe | r2 | B = 500 | 17.70 |
| 101 | guohe | r3 | B = 200 | 7.22 |
| 101 | mr | r1 | draws 5000, 1000/500 | 2.44 |
| 101 | mr | r2 | draws 2000, 500/250 | 1.35 |
| 101 | mr | r3 | draws 1000, 200/100 | 1.05 |
| 102 | family | r1 | k = 156 | 0.38 |
| 102 | guohe | r1 | B = 2000 | 67.45 |
| 102 | guohe | r2 | B = 500 | 16.83 |
| 102 | guohe | r3 | B = 200 | 6.75 |
| 102 | mr | r1 | draws 5000, 1000/500 | 2.12 |
| 102 | mr | r2 | draws 2000, 500/250 | 1.28 |
| 102 | mr | r3 | draws 1000, 200/100 | 1.02 |
| 103 | family | r1 | k = 162 | 0.38 |
| 103 | guohe | r1 | B = 2000 | 70.12 |
| 103 | guohe | r2 | B = 500 | 17.44 |
| 103 | guohe | r3 | B = 200 | 7.16 |
| 103 | mr | r1 | draws 5000, 1000/500 | 2.19 |
| 103 | mr | r2 | draws 2000, 500/250 | 1.38 |
| 103 | mr | r3 | draws 1000, 200/100 | 1.09 |

Component summary (verbatim):

```
s2_sec_family_r1 = mean 0.41, max 0.47 (k = 163)
s2_sec_guohe_r1 = mean 69.59, max 71.20 (B = 2000)
s2_sec_guohe_r2 = mean 17.32, max 17.70 (B = 500)
s2_sec_guohe_r3 = mean 7.04, max 7.22 (B = 200)
s2_sec_mr_r1 = mean 2.25, max 2.44 (draws 5000, 1000/500)
s2_sec_mr_r2 = mean 1.34, max 1.38 (draws 2000, 500/250)
s2_sec_mr_r3 = mean 1.05, max 1.09 (draws 1000, 200/100)
s2_per_rep_sec = r1 72.25, r2 19.07, r3 8.51
s2_elapsed_sec = 304.2
```

| component | rung 1 (manuscript) | rung 2 | rung 3 |
|---|---|---|---|
| data + family + naive fits | mean 0.41 s, max 0.47 s | — | — |
| `guohe_algorithm3()`, r = 1/12 | B = 2000: mean 69.59 s, max 71.20 s | B = 500: mean 17.32 s, max 17.70 s | B = 200: mean 7.04 s, max 7.22 s |
| `fs_mr_inference()`, complement on | draws 5000, 1000/500: mean 2.25 s, max 2.44 s | draws 2000, 500/250: mean 1.34 s, max 1.38 s | draws 1000, 200/100: mean 1.05 s, max 1.09 s |

**Bound movement across MR rungs** (log scale; informational, no gate; `results/rung_bounds.csv`):

| rep | rung | setting | log field lower_1s (H) | log MR lower_1s (H) | log complement upper_1s_s |
|---:|---|---|---:|---:|---:|
| 101 | r1 | draws 5000, 1000/500 | 0.0496 | -0.2590 | 0.0446 |
| 101 | r2 | draws 2000, 500/250 | 0.0505 | -0.2326 | 0.0371 |
| 101 | r3 | draws 1000, 200/100 | 0.0652 | -0.2416 | 0.0341 |
| 102 | r1 | draws 5000, 1000/500 | 0.0197 | -0.2707 | 0.0362 |
| 102 | r2 | draws 2000, 500/250 | 0.0870 | -0.2771 | 0.0381 |
| 102 | r3 | draws 1000, 200/100 | 0.0092 | -0.2611 | 0.0484 |
| 103 | r1 | draws 5000, 1000/500 | 0.0395 | -0.2655 | 0.4042 |
| 103 | r2 | draws 2000, 500/250 | 0.0474 | -0.2744 | 0.3979 |
| 103 | r3 | draws 1000, 200/100 | 0.0729 | -0.2503 | 0.4292 |

Largest spread across the three rungs on any one data set: field lower_1s (H) 0.078 (rep 102: 0.0197 / 0.0870 / 0.0092). MR IJ lower_1s 0.026 (rep 101). Complement field-s upper 0.031 (rep 103).

**Truth curve** (β₂ = 0.3, c_step = 1, seed 20260751; verbatim):

```
s2_cache_n_big = 2e+06
s2_cache_c_step = 0.25
s2_truth_nbig_200000_sec = 5.91
s2_truth_nbig_200000_gates = all PASS
s2_truth_nbig_200000_gate_table = G1 censoring value 0.0019015 <= 0.005; G2 anchor beta(30) value 0.0033106 <= 0.036809; G3 max single rise value -0.0011011 <= 0.02; G4 end-to-end rise value -0.14831 <= 0.018405; G6 dilution drop value 0.14831 >= 0.045
s2_truth_nbig_200000_max_abs_dbeta_vs_cache = 0.0063013
s2_truth_nbig_200000_beta30 = 0.29669 (cache 0.29919)
s2_truth_nbig_200000_beta60 = 0.14838 (cache 0.14547)
s2_truth_nbig_50000_sec = 1.20
s2_truth_nbig_50000_gates = all PASS
s2_truth_nbig_50000_gate_table = G1 censoring value 0.00067845 <= 0.005; G2 anchor beta(30) value 0.0067885 <= 0.073438; G3 max single rise value -0.00047547 <= 0.02; G4 end-to-end rise value -0.1571 <= 0.036719; G6 dilution drop value 0.1571 >= 0.045
s2_truth_nbig_50000_max_abs_dbeta_vs_cache = 0.019017
s2_truth_nbig_50000_beta30 = 0.29321 (cache 0.29919)
s2_truth_nbig_50000_beta60 = 0.13612 (cache 0.14547)
```

Cache `str()` top level (`results/truth_cache_str.txt`):

```
List of 19
 $ beta2        : num 0.3
 $ n_big        : num 2e+06
 $ c_step       : num 0.25
 $ seed         : int 20260751
 $ c_grid       : num [1:121] 30 30.2 30.5 30.8 31 ...
 $ beta_raw     : num [1:121] 0.299 0.296 0.294 0.291 0.289 ...
 $ beta_smooth  : num [1:121] 0.299 0.296 0.294 0.291 0.289 ...
 $ beta_exact   : NULL
 $ scoring_basis: chr "smooth"
 $ offset_anchor: num -0.000813
 $ offset_mean  : num NA
 $ se           : num [1:121] 0.00291 0.0029 0.00289 0.00288 0.00287 ...
 $ n_sub        : int [1:121] 750263 756538 762875 769170 775361 781673 787909 794216 800355 806418 ...
 $ n_events     : int [1:121] 476910 480597 484295 488031 491677 495350 498995 502606 506242 509845 ...
 $ cens_emp     : num 0.392
 $ cens_analytic: num 0.393
 $ gates        :'data.frame':	5 obs. of  5 variables:
 $ elapsed_sec  : num 335
 $ r_version    : chr "R version 4.6.1 (2026-06-24)"
 - attr(*, "class")= chr "gh52_truth"
```

### Projection (`results/projection.md`, transplanted)


Per-replicate seconds (family + Guo-He + MR, means over 3 replicates): rung 1 72.2, rung 2 19.1, rung 3 8.5.
Truth curve beta2 = 0.3, c_step = 1: n_big 2e5 5.9 s; n_big 5e4 1.2 s. Added once for the single beta2 > 0 cell.

Cells {beta2 = 0, 0.3}; 2 cells x n_rep replicates. Truth at n_big = 2e5 included.

| n_rep | rung | 1 core | 2 workers (assumed linear scaling) | 1 core, truth at 5e4 |
|---:|---|---:|---:|---:|
| 50 | 1 | 2.0 h (121 min) | 1.0 h (60 min) | 2.0 h (120 min) |
| 50 | 2 | 0.5 h (32 min) | 0.3 h (16 min) | 0.5 h (32 min) |
| 50 | 3 | 0.2 h (14 min) | 0.1 h (7 min) | 0.2 h (14 min) |
| 100 | 1 | 4.0 h (241 min) | 2.0 h (121 min) | 4.0 h (241 min) |
| 100 | 2 | 1.1 h (64 min) | 0.5 h (32 min) | 1.1 h (64 min) |
| 100 | 3 | 0.5 h (28 min) | 0.2 h (14 min) | 0.5 h (28 min) |
| 200 | 1 | 8.0 h (482 min) | 4.0 h (241 min) | 8.0 h (482 min) |
| 200 | 2 | 2.1 h (127 min) | 1.1 h (64 min) | 2.1 h (127 min) |
| 200 | 3 | 0.9 h (57 min) | 0.5 h (28 min) | 0.9 h (57 min) |

2-worker column: ASSUMPTION -- replicate-level parallelism scaling linearly; the truth curve is not divided (one serial call).

Reference, not an article size: the manuscript's shape, 6 cells x 2000 replicates at rung 1, single core:
240.8 h (14450 min) (truth curves excluded; the manuscript's caches were n_big = 2e6).

## 4. Package facts (`results/package_facts.md`, transplanted)


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

## 5. Script notes

- **Install route:** `devtools::install(".", upgrade = FALSE, build_vignettes = FALSE)`, exit 0. The `R CMD INSTALL` fallback was not used. §1.6–1.7 ran in a fresh Rscript as `s0_install_gate.R` (in this record directory; log `logs/s0_install_gate.log`).
- **Deviation, sourcing (`s1_reproduce.R`, `s2_timing.R`):** `source("quarto/GuoHe/mr_vs_guohe_sim.R", chdir = TRUE)` instead of a plain `source()`.
  - Why: `mr_vs_guohe_sim.R:50–56` and `guohe_sec52_sim.R:70–77` resolve their sibling files from `commandArgs()` `--file=` (this script's directory, where the siblings are absent) and then fall back to `getwd()`. From the repo root that fallback is the repo root, where `guohe_reproduction_sim.R` does not exist, so a plain `source()` would fail.
  - `chdir = TRUE` makes the fallback `quarto/GuoHe` for the duration of the source. Nothing is written there. Side effect: two harmless `normalizePath()` warnings in each log ("No such file or directory" for the relative `--file=` path after the chdir).
- **Single core:** `future::plan("sequential")` as specified. Both scripts were also launched with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` to pin BLAS/OpenMP threads (an addition, so the measured seconds are one core). `guohe_algorithm3()` was called at its default `parallel = FALSE`, which is not passed explicitly, matching the §3.2 call as written.
- **Caps:** `setTimeLimit(elapsed = cap, transient = TRUE)` inside `main()` (one top-level call), plus a between-replicate elapsed check in s1. §3 took 262.9 s (cap 900 s); §4 took 304.2 s (cap 1500 s).
- **Row alignment:** `results[m, ]` assumed to be replicate m. `stopifnot(mf$m == m)` and `stopifnot(mc$m == m)` were asserted where the column exists.
- **seed_base:** asserted with `==` as specified (a stop if false), recorded in keyvalues rather than in `reproduction.csv`.
- **Complement/joint mapping** (column names from `mr_field_complement_vs_guohe_run.R` `.mfc_mr_cols`):
  - `log(mr$complement$debiased$est)` ↔ `cgate_est` (line 229).
  - `log(mr$field$complement$upper_1s_s)` ↔ `c_upper_95_s` (lines 166, 201).
  - `log(mr$field$joint_s$bonf_upper_Hc)` ↔ `c_upper_975_s` (lines 161, 201).
  - `log(mr$field$joint_s$bonf_lower_H)` ↔ `h_lower_975`. The driver emits no column for `joint_s$bonf_lower_H`. `h_lower_975` is taken from `joint$bonf_lower_H` (line 159), and the driver states the harm side is identical in `joint` and `joint_s` (lines 156–158). Measured agreement ≤ 3.33e-16, not bit-identical.
  - `bonf_joint_prob`, `gamma`, `corr` from `joint_s` were recorded next to the stored `joint_*_s` values (keyvalues), not gated.
- **Additions beyond the spec, informational only:**
  - Per-replicate component seconds in s1 (`*_sec_gh`, `*_sec_mr`, `*_sec_rep`) and `sel_lab` in keyvalues.
  - `stopifnot(all.equal(GH52_R_GRID[2], 1/12))`, confirming r2 = 1/12 from the sourced grid.
- **Keyvalues cosmetics:**
  - s1's `kv()` formats numbers with `format(digits = 17)` (hence `0.95199999999999996`, `61.659999999999997`).
  - The `s1_GH52_R_GRID` line shows `", , "` and trailing spaces because `deparse()` wrapped the vector into two strings. The grid is `c(1/3, 1/12, 1/21, 1/30)` (`guohe_sec52_sim.R:85`), and the r2 check above passed.
- **Truth curve:** `gh52_truth_curve()` hard-stops on its own gate failures (`guohe_sec52_truth.R:239–245`). The script wraps it in `tryCatch`, so a 5e4 failure would be recorded and the run continued, and a 2e5 failure stops the run. Neither failed. The comparison with the cache uses `gh52_truth_at(·, 30:60)` at its `"auto"` basis (`beta_smooth` for β₂ = 0.3) for both curves.
- **Projection arithmetic:**
  - Per replicate = mean(family) + mean(Guo & He at the rung) + mean(MR at the rung).
  - The article shape is 2 cells × n_rep. Truth cost (2e5) is added once, for the β₂ = 0.3 cell only; β₂ = 0 has an exact curve (`guohe_sec52_truth.R:197`). A 5e4-truth column is added alongside.
  - 2-worker column: replicate time is halved and the serial truth call is not; this is labelled an assumption.
  - The manuscript shape excludes truth curves.
- **No §3 replicate re-run; no write under `quarto/`; no `R/` change.**

## 6. OPEN ITEMS

1. **Per-replicate Guo & He cost is higher here than in the record.**
   - This machine (Linux x86_64, single core, load avg ≈ 1): `guohe_algorithm3()` at B = 2000 took 57.9–71.2 s per replicate (§3 and §4 together).
   - `REPORT_guohe_fixedr_t7_2026-09-09.md` §5 records 39.1 core-s per replicate on its 12-worker run. That machine is not identified in the lines read here.
   - The projection is a rate on this machine. The article's actual build machine may differ by roughly 1.8×.
2. **Defaults the article must override explicitly.**
   - `fs_mr_inference()` defaults: `draws = 2000L` (manuscript 5000), `include_complement = FALSE`, `reselection` defaults to `"maxcons"` (manuscript `"maxeff"`). `ci_method` defaults to `"field"`.
   - `guohe_algorithm3()` defaults: `orient = -1` and `r = 0.03` (manuscript `+1` and `1/12`). `level`, `B` and `min_events` defaults match the manuscript.
3. **Element names.**
   - `mr$selected_label` exists in 0.4.1 (v1 used it; v2 switched to `selected_index`). Both are available.
   - `joint_s` carries `gamma` *and* `bonf_gamma`. `gamma` was 0.025 at β₂ = 0 and 0.026 at β₂ = 0.3 (stored values agree), so it is not the fixed α/2. The article should name `bonf_*` for the Bonferroni pair.
4. **Complement-driver gap.** No emitted column for `joint_s$bonf_lower_H`. The reproduction used `h_lower_975` (`joint`) per the driver's equality statement.
5. **Truth-curve size.** Cache is n_big = 2e6, c_step = 0.25, and took 335 s.
   - Reduced curves, max |Δβ(c)| on c = 30…60 vs cache: 0.0063 at 2e5, 0.019 at 5e4.
   - Both pass their own gates. The article's scoring basis (cache, or a regenerated smaller curve) is open.
6. **Rung movement.** On single data sets the field lower bound moved by up to 0.078 (log scale) between rungs (rep 102). This is a Monte Carlo-noise question for the article's settings.

## 7. Verdict (numbers only)

1. Selection keys identical: 16 / 16.
2. Floats bit-identical: 61 / 72.
3. Floats all.equal 1e-8: 72 / 72. At 1e-6: 72 / 72.
4. Max |diff|: 4.44089e-16 (c_upper_975_s, t7_beta2_03 m=2).
5. Reproduction: to 1e-8, not bit-for-bit. Gate PASS.
6. Per-replicate seconds, single core: r1 72.25, r2 19.07, r3 8.51.
7. Guo & He B = 2000/500/200: mean 69.59 / 17.32 / 7.04 s. MR rungs: 2.25 / 1.34 / 1.05 s.
8. Truth curve: 2e5 5.91 s, max |Δβ| 0.0063. 5e4 1.20 s, max |Δβ| 0.019. Gates PASS both.
9. 2 cells × 100 reps, single core: rung 1 4.0 h, rung 2 1.1 h, rung 3 0.5 h.
10. Manuscript shape 6 × 2000 at rung 1: 240.8 h single core (not an article size).

## Closing block (§6 post-conditions, after commit `7adb9c16`)

All three commands print nothing (no pre-existing dirt was listed in §1.2).

```
$ git diff --name-only c2967362 HEAD | grep -v '^dev/'
$ git status --porcelain | grep -v '^?? ' | grep -v '^.. dev/'
$ git status --porcelain | grep '^?? ' | grep -v ' dev/'
```

Commits on `feature/articles-guohe` above `c2967362`: `e751c215` (v1 task doc), `462ca74c` (v2 task doc), `7adb9c16` (record + report), and this closing-block commit. No push, no fetch, no merge.
