# s2_timing.R -- Guo & He article Stage 0, section 4 (T): timing ladder,
# single core, cell t7_beta2_03, data seeds base + 101..103; truth curve at
# reduced n_big; projection table (pure arithmetic from measured seconds).
# Task: dev/tasks/cc_task_guohe_article_stage0_2026-10-29_v2.md
#
# Run from the repo root:
#   Rscript dev/reports/guohe_article_stage0_2026-10-29/s2_timing.R

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
suppressMessages({ library(forestsearch); library(survival) })
future::plan("sequential")
source("quarto/GuoHe/mr_vs_guohe_sim.R", chdir = TRUE)   # see s1_reproduce.R

OUT <- "dev/reports/guohe_article_stage0_2026-10-29/results"
KV  <- file.path(OUT, "keyvalues.txt")
CAP_SEC <- 25 * 60
ID <- "t7_beta2_03"

kv <- function(name, value) {
  line <- sprintf("%s = %s", name, paste(value, collapse = ", "))
  write(line, KV, append = TRUE); cat(line, "\n")
}
el <- function() proc.time()[["elapsed"]]

GH_RUNGS <- list(r1 = 2000L, r2 = 500L, r3 = 200L)
MR_RUNGS <- list(r1 = c(draws = 5000L, R_out = 1000L, R_in = 500L),
                 r2 = c(draws = 2000L, R_out = 500L,  R_in = 250L),
                 r3 = c(draws = 1000L, R_out = 200L,  R_in = 100L))

main <- function() {
  setTimeLimit(elapsed = CAP_SEC, transient = TRUE)
  t_start <- el()
  base <- mv_gh52_base(ID)
  tim <- list(); bnd <- list()

  for (j in 101:103) {
    sd <- base + j
    t0 <- el()
    cand <- mv_gh52_regen(ID, j)
    fits <- gh52_subgroup_fits(cand$df, cand)
    sel_lab <- cand$names[which.max(replace(fits$est, !is.finite(fits$est), -Inf))]
    cidx <- mv_cand_idx_52(cand)
    tim[[length(tim) + 1L]] <- data.frame(rep = j, component = "family",
      rung = "r1", setting = sprintf("k = %d", length(cand$cuts)), sec = el() - t0)

    for (rg in names(GH_RUNGS)) {
      t0 <- el()
      fit <- guohe_algorithm3(data = cand$df, outcome = "survival",
                              treatment = "treat", candidates = cand$names,
                              time = "time", event = "event", orient = +1,
                              B = GH_RUNGS[[rg]], r = 1 / 12, level = 0.05,
                              seed = sd + 500000L, min_events = 5L,
                              diagnostics = FALSE)
      tim[[length(tim) + 1L]] <- data.frame(rep = j, component = "guohe",
        rung = rg, setting = sprintf("B = %d", GH_RUNGS[[rg]]), sec = el() - t0)
    }

    for (rg in names(MR_RUNGS)) {
      s <- MR_RUNGS[[rg]]
      t0 <- el()
      mr <- fs_mr_inference(df = cand$df, candidates = cidx, spec = mv_spec52,
                            selected_members = cidx[[sel_lab]],
                            admission = list(effect_floor = NULL, consistency = NULL),
                            reselection = "maxeff", draws = s[["draws"]],
                            multiplier = "poisson", ci_method = "field",
                            seed = sd + 700000L, return_reselection = TRUE,
                            include_complement = TRUE,
                            field_R_out = s[["R_out"]], field_R_in = s[["R_in"]],
                            field_uniform = FALSE, field_complement = TRUE,
                            field_scale_complement = "selected")
      setting <- sprintf("draws %d, %d/%d", s[["draws"]], s[["R_out"]], s[["R_in"]])
      tim[[length(tim) + 1L]] <- data.frame(rep = j, component = "mr",
        rung = rg, setting = setting, sec = el() - t0)
      bnd[[length(bnd) + 1L]] <- data.frame(rep = j, rung = rg, setting = setting,
        fld_lower_1s = log(mr$field$lower_1s),
        mr_lower_1s = log(mr$debiased$lower_1s),
        c_upper_95_s = log(mr$field$complement$upper_1s_s))
    }
    cat(sprintf("rep %d done at %.1f s\n", j, el() - t_start))
  }
  tim <- do.call(rbind, tim); bnd <- do.call(rbind, bnd)
  utils::write.csv(tim, file.path(OUT, "timing.csv"), row.names = FALSE)
  utils::write.csv(bnd, file.path(OUT, "rung_bounds.csv"), row.names = FALSE)

  agg <- do.call(rbind, lapply(split(tim, list(tim$component, tim$rung), drop = TRUE),
    function(d) data.frame(component = d$component[1], rung = d$rung[1],
                           setting = d$setting[1], mean_sec = mean(d$sec),
                           max_sec = max(d$sec))))
  agg <- agg[order(match(agg$component, c("family", "guohe", "mr")), agg$rung), ]
  rownames(agg) <- NULL
  utils::write.csv(agg, file.path(OUT, "timing_summary.csv"), row.names = FALSE)
  for (i in seq_len(nrow(agg)))
    kv(sprintf("s2_sec_%s_%s", agg$component[i], agg$rung[i]),
       sprintf("mean %.2f, max %.2f (%s)", agg$mean_sec[i], agg$max_sec[i],
               agg$setting[i]))

  # ---- truth curves (beta2 = 0.3) ----
  cache <- readRDS("quarto/GuoHe/guohe_sec52_truth_beta2_03.rds")
  writeLines(utils::capture.output(str(cache, max.level = 1)),
             file.path(OUT, "truth_cache_str.txt"))
  kv("s2_cache_n_big", format(cache$n_big)); kv("s2_cache_c_step", cache$c_step)
  cc <- 30:60
  y_cache <- gh52_truth_at(cache, cc)
  truth_sec <- c()
  for (nb in c(2e5, 5e4)) {
    lab <- sprintf("s2_truth_nbig_%g", nb)
    t0 <- el()
    tr <- tryCatch(gh52_truth_curve(0.3, n_big = nb, c_step = 1,
                                    seed = 20260721L + 30L, verbose = FALSE),
                   error = function(e) e)
    sec <- el() - t0
    truth_sec[as.character(nb)] <- sec
    kv(paste0(lab, "_sec"), sprintf("%.2f", sec))
    if (inherits(tr, "error")) {
      kv(paste0(lab, "_gates"), paste("FAIL:", conditionMessage(tr)))
      if (nb == 2e5) stop("GATE: truth curve gates failed at n_big = 2e5")
      next
    }
    kv(paste0(lab, "_gates"), "all PASS")
    kv(paste0(lab, "_gate_table"), paste(sprintf("%s value %.5g %s %.5g",
       tr$gates$gate, tr$gates$value, tr$gates$direction, tr$gates$threshold),
       collapse = "; "))
    y <- gh52_truth_at(tr, cc)
    kv(paste0(lab, "_max_abs_dbeta_vs_cache"), sprintf("%.5g", max(abs(y - y_cache))))
    kv(paste0(lab, "_beta30"), sprintf("%.5f (cache %.5f)", y[1], y_cache[1]))
    kv(paste0(lab, "_beta60"), sprintf("%.5f (cache %.5f)", y[length(cc)],
                                       y_cache[length(cc)]))
  }

  # ---- projection (pure arithmetic) ----
  fam <- agg$mean_sec[agg$component == "family"]
  per_rep <- vapply(c("r1", "r2", "r3"), function(rg)
    fam + agg$mean_sec[agg$component == "guohe" & agg$rung == rg] +
      agg$mean_sec[agg$component == "mr" & agg$rung == rg], 0)
  tr_sec <- truth_sec[["2e+05"]]; tr_sec5 <- truth_sec[["50000"]]
  hm <- function(s) sprintf("%.1f h (%.0f min)", s / 3600, s / 60)
  lines <- c(
    "# Projection (pure arithmetic from measured single-core seconds)", "",
    sprintf("Per-replicate seconds (family + Guo-He + MR, means over 3 replicates): rung 1 %.1f, rung 2 %.1f, rung 3 %.1f.",
            per_rep[1], per_rep[2], per_rep[3]),
    sprintf("Truth curve beta2 = 0.3, c_step = 1: n_big 2e5 %.1f s; n_big 5e4 %.1f s. Added once for the single beta2 > 0 cell.",
            tr_sec, tr_sec5), "",
    "Cells {beta2 = 0, 0.3}; 2 cells x n_rep replicates. Truth at n_big = 2e5 included.", "",
    "| n_rep | rung | 1 core | 2 workers (assumed linear scaling) | 1 core, truth at 5e4 |",
    "|---:|---|---:|---:|---:|")
  for (nr in c(50, 100, 200)) for (i in 1:3) {
    s1 <- 2 * nr * per_rep[i] + tr_sec
    s2 <- 2 * nr * per_rep[i] / 2 + tr_sec
    s5 <- 2 * nr * per_rep[i] + tr_sec5
    lines <- c(lines, sprintf("| %d | %d | %s | %s | %s |", nr, i, hm(s1), hm(s2), hm(s5)))
  }
  lines <- c(lines, "",
    "2-worker column: ASSUMPTION -- replicate-level parallelism scaling linearly; the truth curve is not divided (one serial call).",
    "", "Reference, not an article size: the manuscript's shape, 6 cells x 2000 replicates at rung 1, single core:",
    sprintf("%s (truth curves excluded; the manuscript's caches were n_big = 2e6).",
            hm(6 * 2000 * per_rep[1])))
  writeLines(lines, file.path(OUT, "projection.md"))
  kv("s2_per_rep_sec", sprintf("r1 %.2f, r2 %.2f, r3 %.2f", per_rep[1], per_rep[2], per_rep[3]))
  kv("s2_elapsed_sec", sprintf("%.1f", el() - t_start))
  TRUE
}

ok <- main()
if (!isTRUE(ok)) quit(status = 2)
