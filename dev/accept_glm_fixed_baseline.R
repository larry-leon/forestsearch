# accept_glm_fixed_baseline.R ----------------------------------------------
# [delivery sentinel: c1r1-640c8ccc]
# Arc C' step C1 gate: baseline = "fixed" for glm_dgm.
# Resample-path byte-identity is delegated to the standing phase gates
# (accept_phase45/44/41), re-run in the same battery.
suppressPackageStartupMessages({library(forestsearch); library(survival); library(future)})
fail <- function(...) stop(sprintf(...), call. = FALSE)
gate <- function(id, ok, d = "") { cat(sprintf("[%s] %s %s\n", id, if (ok) "PASS" else "FAIL", d)); if (!ok) fail("Gate %s failed. %s", id, d) }
msg_of <- function(e) tryCatch({ e; NA_character_ }, error = function(x) conditionMessage(x))
gb <- survival::gbsg
gc2 <- data.frame(y = as.numeric(gb$age), yb = as.integer(gb$status),
                  ycnt = as.integer(pmin(gb$nodes, 20L)), t = gb$rfstime/365.25,
                  w = as.integer(gb$hormon), meno = as.integer(gb$meno),
                  grade3 = as.integer(gb$grade == 3), size = gb$size)
mk <- function(otype, ovar, off = NULL) generate_glm_dgm(
  data = gc2, factor_vars = c("meno","grade3"), outcome_var = ovar,
  treatment_var = "w", outcome_type = otype, offset_var = off,
  subgroup_vars = c("meno","grade3"), subgroup_cuts = list(meno=1L, grade3=1L),
  model = "null", n_super = 2000L, seed = 301L)
dg_c <- mk("continuous","y"); dg_b <- mk("binary","yb"); dg_k <- mk("count","ycnt","t")
cat("== G1: df_source stored, consistent with df_super ==========\n")
for (nm in c("c","b","k")) {
  d <- get(paste0("dg_", nm))
  gate(paste0("G1/", nm, "/present"), !is.null(d$df_source) &&
         nrow(d$df_source) == nrow(gc2), sprintf("nrow=%s", nrow(d$df_source)))
  po <- intersect(c("p0","p1","mu0","mu1"), names(d$df_source))
  gate(paste0("G1/", nm, "/po-cols"), length(po) == 2L, paste(po, collapse=","))
}
# internal consistency: df_super rows are draws of df_source rows; PO values
# are deterministic per row, so matching covariate rows must match POs.
d <- dg_c; key <- paste(d$df_source$y, d$df_source$w, d$df_source$meno, d$df_source$grade3, d$df_source$size)
ks <- paste(d$df_super$y, d$df_super$w, d$df_super$meno, d$df_super$grade3, d$df_super$size)
ix <- match(ks, key)
gate("G1/consistency", !anyNA(ix) &&
       isTRUE(all.equal(d$df_super$mu0, d$df_source$mu0[ix], tolerance = 0)) &&
       isTRUE(all.equal(d$df_super$mu1, d$df_source$mu1[ix], tolerance = 0)),
     "df_super PO columns == df_source PO columns at matched rows (exact)")
cat("== G2: simulator fixed semantics ============================\n")
s1 <- simulate_from_glm_dgm(dg_c, baseline = "fixed", seed = 1)
s2 <- simulate_from_glm_dgm(dg_c, baseline = "fixed", seed = 2)
cov_cols <- c("y","meno","grade3","size","flag_harm")
gate("G2/panel-frozen", identical(s1[cov_cols], s2[cov_cols]) &&
       nrow(s1) == nrow(gc2), "covariates+flag_harm identical across seeds")
gate("G2/redraws", !identical(s1$treat_sim, s2$treat_sim) &&
       !identical(s1$y_sim, s2$y_sim), "treatment and outcomes redrawn")
gate("G2/n-rule-sim", is.na(msg_of(simulate_from_glm_dgm(dg_c, n = nrow(gc2), baseline = "fixed", seed = 3))) &&
       grepl("every source patient exactly once",
             msg_of(simulate_from_glm_dgm(dg_c, n = 50, baseline = "fixed", seed = 3))),
     "simulator: NULL-or-full accepted, other n stops (survival mirror)")
old <- dg_c; old$df_source <- NULL
gate("G2/actionable", grepl("rebuild it with generate_glm_dgm()",
       msg_of(simulate_from_glm_dgm(old, baseline = "fixed")), fixed = TRUE), "")
cat("== G3: fixed parity vs hand loop through the runner =========\n")
sg <- list(list(id="flag_itt == 1", name="All", grp="ITT"),
           list(id="meno == 1", name="Post", grp="Clinical"),
           list(id="grade3 == 1", name="G3", grp="Clinical"))
plan(sequential)
res_f <- run_subgroup_sims(dg_c, sg, n_sims = 6, benchmarks = NULL,
                           baseline = "fixed", workers = 1, seed_base = 0L,
                           hr_true = dg_c$hazard_ratios$overall)
hand <- sapply(1:6, function(ss) {
  df <- simulate_from_glm_dgm(dg_c, baseline = "fixed", seed = 0L + ss)
  df$flag_itt <- 1L
  d2 <- df[df$meno == 1, ]
  fit <- stats::lm(y_sim ~ treat_sim, data = d2)
  co <- summary(fit)$coefficients["treat_sim", ]
  co[["Estimate"]]
})
gate("G3/parity", isTRUE(all.equal(unname(res_f$sim_hrs[, "Post"]), unname(hand), tolerance = 0)),
     "runner fixed estimates == hand lm loop on frozen panel (exact)")
gate("G3/ns-constant", all(res_f$sim_ns[, "All"] == nrow(gc2)),
     "per-trial N == nrow(df_source) every replicate")
gate("G3/wrapper-n-rule", grepl("must be NULL when baseline",
       msg_of(run_subgroup_sims(dg_c, sg, 2, n = nrow(gc2), benchmarks = NULL,
                                baseline = "fixed", workers = 1, seed_base = 0L))),
     "wrapper: strict n = NULL (stricter than the simulator; both families)")
gate("G3/wrapper-actionable", grepl("Rebuild it with generate_glm_dgm()",
       msg_of(run_subgroup_sims(old, sg, 2, benchmarks = NULL,
                                baseline = "fixed", workers = 1, seed_base = 0L)), fixed = TRUE), "")
cat("== G4: per-type fixed smoke + serialization ================\n")
for (nm in c("b","k")) {
  d <- get(paste0("dg_", nm))
  r <- run_subgroup_sims(d, sg, n_sims = 4, benchmarks = NULL,
                         baseline = "fixed", workers = 1, seed_base = 0L,
                         hr_true = d$hazard_ratios$overall)
  gate(paste0("G4/", nm), all(r$sim_ns[, "All"] == nrow(gc2)) &&
         any(is.finite(r$sim_hrs)), sprintf("measure=%s", r$effect$measure))
}
fx <- subgroup_glm()
gate("G4/512KB", length(serialize(fx, NULL)) < 512*1024, "")
cat("\nALL GATES PASS -- GLM fixed-baseline (Arc C' step C1) accepted.\n")
