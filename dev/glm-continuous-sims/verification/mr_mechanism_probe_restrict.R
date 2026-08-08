src <- readLines("dev/glm-continuous-sims/verification/mr_mechanism_A1.R")
cut <- grep("^## --- run ---", src)
eval(parse(text = paste(src[1:(cut-1)], collapse="\n")))
library(doFuture); library(foreach); library(future)

probe <- function(i, n_sample) {
  sd_i <- SEED_TABLE[i]; RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)
  df <- simulate_from_glm_dgm(dgm, n=n_sample, seed=sd_i); df$id <- seq_len(nrow(df))
  fs <- tryCatch(suppressWarnings(forestsearch(df.analysis=df,
    confounders.name=confounders_analysis, outcome.name=outcome_name,
    treat.name=treat_name, id.name=id_name, outcome_type="continuous",
    effect_measure="MD", effect.threshold=md_threshold,
    consistency.threshold=md_consistency, pconsistency.threshold=pconsistency,
    fs.splits=fs_splits, n.min=n_min, d0.min=d0_min, d1.min=d1_min, maxk=maxk,
    vi.grf.min=vi_grf_min, consistency_method="resample",
    use_twostage=use_twostage, is.RCT=is_rct, adverse_outcome=adverse_outcome,
    details=FALSE, quiet=TRUE, seedit=sd_i,
    parallel_args=list(plan="sequential"), mr_inference=FALSE)), error=function(e) NULL)
  if (is.null(fs) || is.null(fs$sg.harm)) return(NULL)
  Z <- build_Z(df, fs$confounders.evaluated, fs$confounders.candidate)
  fam <- build_family(Z, maxk, n_min)
  spec <- list(outcome_type="continuous", effect_measure="MD", treat.name=treat_name,
    outcome.name=outcome_name, event.name=NULL, offset.name=NULL,
    adjust_covariates=NULL, adverse_outcome=adverse_outcome)
  asm <- get(".fs_mr_assemble", envir=asNamespace("forestsearch"))(df, fam, spec)
  adm <- fs$admission; zq <- qnorm((1+adm$consistency$p_star)/2)
  t_g <- pmax(adm$effect_floor, adm$consistency$c_cons + zq*asm$sigma_D)
  pass <- which(asm$beta_hat >= t_g)
  if (!length(pass)) return(NULL)
  z <- (asm$beta_hat - adm$consistency$c_cons)/asm$sigma_D
  win <- pass[which.max(z[pass])]
  win_lab <- asm$names[win]

  # The identifier's QUALIFYING set: rows of find.grps$out.found$hr.subgroups,
  # relabelled into the same "qJ.d & qK.e" vocabulary as MR's family names.
  tbl <- as.data.frame(fs$find.grps$out.found$hr.subgroups)
  ind <- setdiff(names(tbl), c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)"))
  qual <- vapply(seq_len(nrow(tbl)), function(k) {
    on <- ind[which(as.numeric(tbl[k, ind]) == 1)]
    paste(on, collapse=" & ") }, character(1))
  sel_members <- which(fs$grp.consistency$sg.harm.id == 1L)
  hit <- which(vapply(fam, function(ix) setequal(ix, sel_members), logical(1)))
  sel_lab <- if (length(hit)) names(fam)[hit[1]] else NA_character_
  # argmax z RESTRICTED to the identifier's own qualifying set
  in_q <- pass[asm$names[pass] %in% qual]
  win_r <- if (length(in_q)) in_q[which.max(z[in_q])] else NA_integer_
  restricted_match <- !is.na(win_r) && !is.na(sel_lab) && identical(asm$names[win_r], sel_lab)
  data.frame(sim_id=i, n=n_sample,
    restricted_match = restricted_match,
    restricted_size = if (is.na(win_r)) NA_integer_ else asm$sizes[win_r],
    n_in_q = length(in_q),
    mr_family=length(asm$names), qualifying=nrow(tbl), admitted=length(pass),
    mr_argmax_in_qualifying = win_lab %in% qual,
    mr_argmax_size = asm$sizes[win],
    realized_size = length(sel_members),
    realized_in_mr_family = length(hit) > 0,
    stringsAsFactors=FALSE)
}
plan("multisession", workers=40)
for (cl in list(list(n=1000L,k=40L), list(n=4000L,k=20L))) {
  res <- foreach(i=seq_len(cl$k), .options.future=list(packages="forestsearch", seed=TRUE)) %dofuture% probe(i, cl$n)
  r <- do.call(rbind, Filter(Negate(is.null), res))
  cat(sprintf("\n=== n = %d, %d replicates ===\n", cl$n, nrow(r)))
  cat(sprintf("  MR family size          median %.0f\n", median(r$mr_family)))
  cat(sprintf("  identifier qualifying   median %.0f\n", median(r$qualifying)))
  cat(sprintf("  MR admitted             median %.0f\n", median(r$admitted)))
  cat(sprintf("  MR argmax-z winner IS in the identifier's qualifying set : %d/%d = %.3f\n",
      sum(r$mr_argmax_in_qualifying), nrow(r), mean(r$mr_argmax_in_qualifying)))
  cat(sprintf("  realized winner is in MR's family                        : %d/%d = %.3f\n",
      sum(r$realized_in_mr_family), nrow(r), mean(r$realized_in_mr_family)))
  cat(sprintf("  MR argmax size median %.0f  vs realized median %.0f\n",
      median(r$mr_argmax_size), median(r$realized_size)))
  cat(sprintf("  RESTRICTED to the identifier's qualifying set:\n"))
  cat(sprintf("    argmax-z == realized winner : %d/%d = %.3f\n",
      sum(r$restricted_match), nrow(r), mean(r$restricted_match)))
  cat(sprintf("    restricted argmax size median %.0f   (candidates in intersection: median %.0f)\n",
      median(r$restricted_size, na.rm=TRUE), median(r$n_in_q)))
}
plan("sequential")
