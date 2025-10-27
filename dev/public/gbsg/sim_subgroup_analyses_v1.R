get_SGanalyses <- function(dgm, subgroups_name,subgroups_id, sims, wname, bw, Ndraw = nrow(dgm$df_super),bmcut=log(1),
outfile=c("results/bm1_bw=log5_v1.Rdata")){

res <- list()

if(length(subgroups_id) != length(subgroups_name)) stop("Subgroup ids and names do not match") 
  
# Template for storing estimates Include analysis name
temp <- matrix(nrow = sims, ncol = length(subgroups_id))
colnames(temp) <- c(subgroups_name)

# Conduct 6 analyses with different approaches for baseline adjustment
# Store the point and upper bound estimates for the Cox HRs
subgroup_hrs1 <- subgroup_hrs2 <- subgroup_hrs3 <- subgroup_hrs4 <- subgroup_hrs5 <- subgroup_hrs6 <- temp
subgroup_ubs1 <- subgroup_ubs2 <- subgroup_ubs3 <- subgroup_ubs4 <- subgroup_ubs5 <- subgroup_ubs6 <- temp
subgroup_ns <- temp

rm("temp")

for (ss in 1:sims) {
  
  df_sim <- draw_sim_stratified(dgm = dgm, ss = ss, Ndraw = Ndraw, wname = wname, bw = bw,
                                strata_rand = "stratum")
  # ITT flag
  df_sim$itt <- c("Y")
  
  # ER positive
  df_sim$cut_opt <- with(df_sim, ifelse(z < bmcut, 0, 1))
  
  if (ss == 1) cat("Optimal biomarker cut=",c(bmcut), "\n")
  
  for (gg in 1:length(subgroups_id)) {
    
    df_sg <- subset(df_sim, eval(parse(text = c(subgroups_id[gg]))))
    
    subgroup_ns[ss, gg] <- c(nrow(df_sg))
    
    # Cox stratified by randomization stratification
    fit <- try(summary(coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR),
                             data = df_sg)), TRUE)
    if (!inherits(fit, "try-error")) {
      subgroup_hrs1[ss, gg] <- c(fit$conf.int[1, 1])
      subgroup_ubs1[ss, gg] <- c(fit$conf.int[1, 4])
      rm("fit")
    }
    
    # Non-stratified Cox
    fit <- try(summary(coxph(Surv(y.sim, event.sim) ~ treat.sim, data = df_sg)),
               TRUE)
    if (!inherits(fit, "try-error")) {
      subgroup_hrs2[ss, gg] <- c(fit$conf.int[1, 1])
      subgroup_ubs2[ss, gg] <- c(fit$conf.int[1, 4])
      rm("fit")
    }
    
    # Stratify by W (strong prognostic factor)
    fit <- try(summary(coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(w),
                             data = df_sg)), TRUE)
    if (!inherits(fit, "try-error")) {
      subgroup_hrs3[ss, gg] <- c(fit$conf.int[1, 1])
      subgroup_ubs3[ss, gg] <- c(fit$conf.int[1, 4])
      rm("fit")
    }
    
    # Stratify by biomarker cut
    n_opt <- with(df_sg, sum(cut_opt))
    m_opt <- nrow(df_sg) - n_opt
    if (n_opt >= 5 & m_opt >= 5) {
      fit <- try(summary(coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(cut_opt),
                               data = df_sg)), TRUE)
      if (!inherits(fit, "try-error")) {
        subgroup_hrs4[ss, gg] <- c(fit$conf.int[1, 1])
        subgroup_ubs4[ss, gg] <- c(fit$conf.int[1, 4])
        rm("fit")
      }
    }
    
    # Stratify by both W and rand strata
    fit <- try(summary(coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(w) +
                               strata(strata.simR), data = df_sg)), TRUE)
    if (!inherits(fit, "try-error")) {
      subgroup_hrs5[ss, gg] <- c(fit$conf.int[1, 1])
      subgroup_ubs5[ss, gg] <- c(fit$conf.int[1, 4])
      rm("fit")
    }
    
    # Adjust for w and age as covariates
    fit <- try(summary(coxph(Surv(y.sim, event.sim) ~ treat.sim + w + age + 
                               strata(strata.simR), data = df_sg)), TRUE)
    if (!inherits(fit, "try-error")) {
      subgroup_hrs6[ss, gg] <- c(fit$conf.int[1, 1])
      subgroup_ubs6[ss, gg] <- c(fit$conf.int[1, 4])
      rm("fit")
    }
    
  }
  
  if (ss <= 10) {
    kmH.fit <- KM.plot.2sample.weighted(df = df_sim, tte.name = "y.sim", event.name = "event.sim",
                                        treat.name = "treat.sim", strata.name = "strata.simO", stop.onerror = FALSE,
                                        risk.set = TRUE, by.risk = 12, risk.cex = 0.8, risk_offset = 0.15, risk_delta = 0.075,
                                        Xlab = "Months", Ylab = "Overall Survival", details = FALSE, show.ticks = TRUE,
                                        prob.points = seq(0, 1, by = 0.1), col.1 = "black", col.2 = "blue", 
                                        ltys = c(1,1), lwds = c(2, 2), cox.cex = 0.8, show.logrank = TRUE, lr.digits = 3,
                                        lr.eps = 1e-04, put.legend.lr = "top", show.med = TRUE, med.cex = 0.8,
                                        ymed.offset = 0.3, show.cox = TRUE, cox.digits = 3, cox.eps = 1e-04,
                                        show_arm_legend = TRUE, arms = c("Treat", "Control"), arm.cex = 0.8)
  }
  
}

if(!is.null(outfile)){save(dgm, subgroup_ns, subgroups_name, subgroup_hrs1,
     subgroup_hrs2, subgroup_ubs1, subgroup_ubs2, subgroup_hrs3, subgroup_hrs4, subgroup_ubs3,
     subgroup_ubs4, subgroup_hrs5, subgroup_hrs6, subgroup_ubs5, subgroup_ubs6, file = outfile)
}

res$subgroups_hrs1 <- subgroup_hrs1 
res$subgroups_hrs2 <- subgroup_hrs2 
res$subgroups_hrs3 <- subgroup_hrs3 
res$subgroups_hrs4 <- subgroup_hrs4 
res$subgroups_hrs5 <- subgroup_hrs5 
res$subgroups_hrs6 <- subgroup_hrs6 


res$subgroups_ubs1 <- subgroup_ubs1 
res$subgroups_ubs2 <- subgroup_ubs2 
res$subgroups_ubs3 <- subgroup_ubs3 
res$subgroups_ubs4 <- subgroup_ubs4 
res$subgroups_ubs5 <- subgroup_ubs5 
res$subgroups_ubs6 <- subgroup_ubs6 

return(res)

}

