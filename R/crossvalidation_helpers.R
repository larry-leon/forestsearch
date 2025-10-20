#' ForestSearch K-Fold Cross-Validation Output Summary
#'
#' Summarizes cross-validation results for ForestSearch, including subgroup agreement and performance metrics.
#'
#' @param res List. Result object from ForestSearch cross-validation, must contain elements: \code{cv_args}, \code{sg_analysis}, \code{sg0.name}, \code{sg1.name}, \code{Kfolds}, \code{resCV}.
#' @param details Logical. Print details during execution (default: FALSE).
#' @param outall Logical. If TRUE, returns all summary tables; if FALSE, returns only metrics (default: FALSE).
#'
#' @return If \code{outall=FALSE}, a list with \code{sens_metrics_original} and \code{find_metrics}. If \code{outall=TRUE}, a list with summary tables and metrics.
#'
#' @importFrom weightedsurv df_counting
#' @importFrom stringr str_sub str_length
#' @export

forestsearch_KfoldOut <- function(res, details = FALSE,outall = FALSE){
  if (!requireNamespace("weightedsurv", quietly = TRUE)) {
    stop("Package 'weightedsurv' needed for this function to work. Please install it via install_github('larry-leon/weightedsurv').")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
  stop("Package 'stringr' required: If issues loading consider substituting
 'stringr::str_sub(x,start,end)' with 'substr(x,start,end)' and 'stringr::str_length(x)' with 'nchar(x)'.")
  }

  # List of required elements in res
  required_elements <- c(
    "cv_args", "sg_analysis", "sg0.name", "sg1.name", "Kfolds", "resCV"
  )
    # Check for existence in res
  missing_elements <- required_elements[!sapply(required_elements, function(x) !is.null(res[[x]]))]
  if (length(missing_elements) > 0) {
    stop("The following elements are missing in 'res': ", paste(missing_elements, collapse = ", "))
  }
    # For nested elements in cv_args
  cv_args_elements <- c("confounders.name", "outcome.name", "event.name", "treat.name")
  missing_cv_args <- cv_args_elements[!sapply(cv_args_elements, function(x) !is.null(res$cv_args[[x]]))]
  if (length(missing_cv_args) > 0) {
    stop("The following elements are missing in 'res$cv_args': ", paste(missing_cv_args, collapse = ", "))
  }
  confounders.name <- c(res$cv_args$confounders.name)
  outcome.name <- c(res$cv_args$outcome.name)
  event.name <- c(res$cv_args$event.name)
  treat.name <- c(res$cv_args$treat.name)
  sg_analysis <- c(res$sg_analysis)
  sg0.name <- c(res$sg0.name)
  sg1.name <- c(res$sg1.name)

  Kfolds <- res$Kfolds
  df_CV <- res$resCV

  ## outall=FALSE --> only output find_metrics and sens_metrics_original
  ## Extract sg1 and sg2
  sg1 <- sg2 <- rep(NA,Kfolds)
  for(ks in 1:Kfolds){
    ## First element since all identical for same cvindex (=ks)
    sg1[ks] <- subset(df_CV, cvindex==ks)$sg1[1]
    sg2[ks] <- subset(df_CV, cvindex==ks)$sg2[1]
  }
  SGs_found <- cbind(sg1,sg2)

  CV_summary <- CV_sgs(sg1 = sg1, sg2 = sg2, confs = confounders.name, sg_analysis = sg_analysis)

  if(details){
    cat("Any found",c(mean(CV_summary$any_found)),"\n")
    cat("Exact match",c(mean(CV_summary$exact_match)),"\n")
    cat("At least 1 match",c(mean(CV_summary$one_match)),"\n")
    cat("Cov 1 any",c(mean(CV_summary$cov1_any)),"\n")
    cat("Cov 2 any",c(mean(CV_summary$cov2_any)),"\n")
    cat("Cov 1 and 2 any",c(mean(ifelse(CV_summary$cov1_any & CV_summary$cov2_any,1,0))),"\n")
    cat("Cov 1 exact",c(mean(CV_summary$cov1_exact)),"\n")
    cat("Cov 2 exact",c(mean(CV_summary$cov2_exact)),"\n")
  }

  find_metrics <- c(mean(CV_summary$any_found),mean(CV_summary$exact_match),mean(CV_summary$one_match),
                    mean(CV_summary$cov1_any),mean(CV_summary$cov2_any),mean(ifelse(CV_summary$cov1_any & CV_summary$cov2_any,1,0)),
                    mean(CV_summary$cov1_exact),mean(CV_summary$cov2_exact))

  names(find_metrics) <- c("Any","Exact","At least 1","Cov1","Cov2","Cov 1 & 2","Cov1 exact","Cov2 exact")

  # Propn agreement in H and H^c
  # between sample estimate (original) and cross-validation
  # Possible all folds return no subgroup
  n_sgfound <- length(unique(df_CV$treat.recommend))
  if(n_sgfound==2){
    tabit <- with(df_CV,table(treat.recommend,treat.recommend.original))
    sensH <- tabit[1,1]/sum(tabit[,1])
    sensHc <- tabit[2,2]/sum(tabit[,2])
    ppvH <- tabit[1,1]/sum(tabit[1,])
    ppvHc <- tabit[2,2]/sum(tabit[2,])
  }
  if(n_sgfound==1){
    tabit <- with(df_CV,table(treat.recommend,treat.recommend.original))
    sensH <- tabit[1,1]/sum(tabit[,1])
    sensHc <- NA
    ppvH <- tabit[1,1]/sum(tabit[1,])
    ppvHc <- NA
  }
  sens_metrics_original <- c(sensH,sensHc,ppvH,ppvHc)
  names(sens_metrics_original) <- c("sens_H","sens_Hc","ppv_H","ppv_Hc")

  if(details) cat("Agreement (sens, ppv) in H and Hc:",c(sensH, sensHc, ppvH, ppvHc),"\n")

  if(outall){
    itt_tab <- SG_tab_estimates(df = as.data.frame(df_CV), SG_flag = "ITT", draws = 0, details = FALSE,
                                outcome.name = outcome.name,
                                event.name = event.name,
                                treat.name = treat.name,
                                strata.name = NULL,
                                potentialOutcome.name = NULL,
                                est.scale = est.scale, sg0_name = "Questionable", sg1_name = "Recommend")
    if(n_sgfound == 2){
      SG_tab_Kfold <- SG_tab_estimates(df = as.data.frame(df_CV),SG_flag = "treat.recommend", sg1_name = sg1.name,sg0_name = sg0.name,
                                       outcome.name = outcome.name, event.name = event.name, treat.name = treat.name, strata.name =NULL, draws = 0, details = FALSE)
    }
    if(n_sgfound == 1){
      SG_tab_Kfold <- itt_tab
    }
    # Data analysis subgroups
    SG_tab_original <- SG_tab_estimates(df=as.data.frame(df_CV), SG_flag="treat.recommend.original", sg1_name = sg1.name, sg0_name = sg0.name,
                                        outcome.name = outcome.name,event.name = event.name,treat.name = treat.name,draws = 0)

    # Combine tables starting with ITT

    estimates_toget <-c("Subgroup","n","n1","m1","m0","RMST","HR (95% CI)")

    if(n_sgfound==2){
      temp1 <- itt_tab[estimates_toget]
      temp2a <- SG_tab_original[1,estimates_toget]
      temp2b <- SG_tab_original[2,estimates_toget]
      temp3a <- SG_tab_Kfold[1,estimates_toget]
      temp3b <- SG_tab_Kfold[2,estimates_toget]
      tab_all <- rbind(temp1,temp2a,temp3a,temp2b,temp3b)
      colnames(tab_all) <- c("Subgroup","n","n1","m1","m0","RMST","Hazard ratio")
      rownames(tab_all) <- c("Overall","FA_0","KfA_0","FA_1","KfA_1")
    }

    if(n_sgfound == 1){
      temp1 <- itt_tab[estimates_toget]
      temp2a <- SG_tab_original[1,estimates_toget]
      temp2b <- SG_tab_original[2,estimates_toget]
      tab_all <- rbind(temp1,temp2a,temp2b)
      colnames(tab_all) <- c("Subgroup","n","n1","n0","m1","m0","RMST","Hazard ratio")
      rownames(tab_all) <- c("Overall","FA_0","FA_1")
    }

    if(details) print(tab_all)

    out <- list(itt_tab=itt_tab,SG_tab_original=SG_tab_original,SG_tab_Kfold=SG_tab_Kfold,
                CV_summary=CV_summary,sens_metrics_original=sens_metrics_original,find_metrics=find_metrics,
                SGs_found=SGs_found,tab_all=tab_all)
  }
  if(!outall){
    out <- list(sens_metrics_original=sens_metrics_original,find_metrics=find_metrics)
  }
  return(out)
}

#' Cross-Validation Subgroup Match Summary
#'
#' Summarizes the match between cross-validation subgroups and analysis subgroups.
#'
#' @param sg1 Character vector. Subgroup 1 labels for each fold.
#' @param sg2 Character vector. Subgroup 2 labels for each fold.
#' @param confs Character vector. Confounder names.
#' @param sg_analysis Character vector. Subgroup analysis labels.
#'
#' @return List with indicators for any match, exact match, one match, and covariate-specific matches.
#'
#' @importFrom stringr str_sub str_length
#' @export

CV_sgs <- function(sg1, sg2, confs, sg_analysis){
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' required.")
  }

  any_found <- ifelse(!is.na(sg1) | !is.na(sg2),1,0)
  sg_depth <- length(sg_analysis)

  if(sg_depth==2){
    sg1a <- sg_analysis[1]
    sg2a <- sg_analysis[2]
    ## Exact match on both to analysis data
    exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a) & (sg1 == sg2a | sg2 == sg2a),1,0)
    exact_match[is.na(exact_match)] <- 0.0
    ## At least 1 exact match
    one_match <- ifelse((sg1 == sg1a | sg2 == sg1a) | (sg1 == sg2a | sg2 == sg2a),1,0)
    one_match[is.na(one_match)] <- 0.0
    ## Cov 1 exact
    cov1_match <- ifelse((sg1 == sg1a | sg2 == sg1a),1,0)
    cov1_match[is.na(cov1_match)] <- 0.0
    cov2_match <- ifelse((sg1 == sg2a | sg2 == sg2a),1,0)
    cov2_match[is.na(cov2_match)] <- 0.0
    ## Find confounder names involved in sg1a and sg2a
    ## Add { or !{ to names for matching (a bit tedious, but let's see)
    dda <- charmatch("{",sg1a, nomatch=0)
    ddb <- charmatch("!{",sg1a, nomatch=0)
    if(dda ==1) aa <- rep("{",length(confs))
    if(ddb ==1) aa <- rep("!{",length(confs))

    temp <- paste0(aa,confs)
    temp <-paste0(temp,"}")
    ## Is sg1a confounder (NOT necessarily same cut) involved in any
    loc_name <- charmatch(temp,sg1a)
    index_name <- which(loc_name==1)

    if(length(index_name)==0){
      temp <- paste0(aa,confs)
      loc_name <- charmatch(temp,sg1a)
      index_name <- which(loc_name==1)
    }

    # Find exact match
    if(length(index_name)>1){
      confs2 <- confs[index_name]
      lc <- stringr::str_length(confs2)
      if(dda==1) ctoget <- stringr::str_sub(sg1a,2,max(lc))
      if(ddb==1) ctoget <- stringr::str_sub(sg1a,3,max(lc))
      itoget <- which(confs == ctoget)
      cfs <- confs[itoget]
    }
    if(length(index_name)==1) cfs <- confs[which(loc_name==1)]

    bb1 <- grepl(cfs,sg1)
    bb2 <- grepl(cfs,sg2)
    cov1_any <- ifelse(bb1 | bb2, 1,0)
    rm("bb1","bb2","cfs")


    ## Second
    dda <- charmatch("{",sg2a, nomatch=0)
    ddb <- charmatch("!{",sg2a, nomatch=0)
    if(dda ==1) aa <- rep("{",length(confs))
    if(ddb ==1) aa <- rep("!{",length(confs))
    temp <- paste0(aa,confs)
    temp <-paste0(temp,"}")

    loc_name <- charmatch(temp,sg2a)
    index_name <- which(loc_name==1)

    if(length(index_name)==0){
      temp <- paste0(aa,confs)
      loc_name <- charmatch(temp,sg2a)
      index_name <- which(loc_name==1)
    }

    # Find exact match
    if(length(index_name)>1){
      confs2 <- confs[index_name]
      lc <- stringr::str_length(confs2)
      if(dda==1) ctoget <- stringr::str_sub(sg2a,2,max(lc))
      if(ddb==1) ctoget <- stringr::str_sub(sg2a,3,max(lc))
      itoget <- which(confs == ctoget)
      cfs <- confs[itoget]
    }
    if(length(index_name)==1) cfs <- confs[which(loc_name==1)]

    bb1 <- grepl(cfs,sg1)
    bb2 <- grepl(cfs,sg2)
    cov2_any <- ifelse(bb1 | bb2, 1,0)
  }
  if(sg_depth==1){
    sg1a <- sg_analysis[1]
    ## Exact match on both to analysis data
    exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a),1,0)
    exact_match[is.na(exact_match)] <- 0.0
    dda <- charmatch("{",sg1a, nomatch=0)
    ddb <- charmatch("!{",sg1a, nomatch=0)
    if(dda ==1) aa <- rep("{",length(confs))
    if(ddb ==1) aa <- rep("!{",length(confs))
    temp <- paste0(aa,confs)
    temp <-paste0(temp,"}")
    ## Is sg1a confounder (NOT necessarily same cut) involved in any
    loc_name <- charmatch(temp,sg1a)
    index_name <- which(loc_name==1)

    if(length(index_name)==0){
      temp <- paste0(aa,confs)
      loc_name <- charmatch(temp,sg1a)
      index_name <- which(loc_name==1)
    }

    # If names have numbers which may not be unique
    # E.g., "z1", "z11", this will match both
    # Find exact
    # Find exact match
    if(length(index_name)>1){
      confs2 <- confs[index_name]
      lc <- stringr::str_length(confs2)
      if(dda==1) ctoget <- stringr::str_sub(sg1a,2,max(lc))
      if(ddb==1) ctoget <- stringr::str_sub(sg1a,3,max(lc))
      itoget <- which(confs == ctoget)
      cfs <- confs[itoget]
    }
    if(length(index_name)==1) cfs <- confs[which(loc_name==1)]
    bb1 <- grepl(cfs,sg1)
    bb2 <- grepl(cfs,sg2)
    cov1_any <- ifelse(bb1 | bb2, 1,0)
    one_match <- exact_match
    cov2_any <- NA

    cov1_match <- exact_match
    cov2_match <- NA

  }
  return(list(any_found=any_found, exact_match=exact_match, one_match=one_match, cov1_any=cov1_any, cov2_any=cov2_any, cov1_exact=cov1_match, cov2_exact=cov2_match))
}

