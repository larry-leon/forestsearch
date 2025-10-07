#' GRF Subgroup Identification for Survival Data
#'
#' Identifies subgroups with differential treatment effect using generalized random forests (GRF) and policy trees.
#'
#' @param data Data frame containing the analysis data.
#' @param confounders.name Character vector of confounder variable names.
#' @param outcome.name Character. Name of outcome variable (e.g., time-to-event).
#' @param event.name Character. Name of event indicator variable (0/1).
#' @param id.name Character. Name of ID variable.
#' @param treat.name Character. Name of treatment group variable (0/1).
#' @param frac.tau Numeric. Fraction of tau for GRF horizon (default: 1.0).
#' @param n.min Integer. Minimum subgroup size (default: 60).
#' @param dmin.grf Numeric. Minimum difference in subgroup mean (default: 0.0).
#' @param RCT Logical. Is the data from a randomized controlled trial? (default: TRUE)
#' @param details Logical. Print details during execution (default: FALSE).
#' @param sg.criterion Character. Subgroup selection criterion ("mDiff" or "Nsg").
#' @param maxdepth Integer. Maximum tree depth (1, 2, or 3; default: 2).
#' @param seedit Integer. Random seed (default: 8316951).
#'
#' @return A list with GRF results, including subgroup definition, tree cuts, and policy tree objects.
#'
#' @importFrom grf causal_survival_forest
#' @importFrom policytree double_robust_scores policy_tree
#' @importFrom stats aggregate
#' @export

grf.subg.harm.survival <- function(data,confounders.name ,outcome.name ,event.name, id.name, treat.name, frac.tau = 1.0,n.min = 60, dmin.grf = 0.0,
RCT = TRUE, details = FALSE,sg.criterion = "mDiff", maxdepth = 2,seedit = 8316951){
if(maxdepth>3) stop("Max depth at most 3")

# original data.frame version
temp <- as.matrix(data[,confounders.name])
# Convert to numeric
X <- apply(temp,2,as.numeric)
Y <- data[,outcome.name]
W <- data[,treat.name]
D <- data[,event.name]

tau.rmst <- frac.tau * min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))

if(RCT) cs.forest <- try(suppressWarnings(grf::causal_survival_forest(X,Y,W,D,W.hat=0.5,horizon=tau.rmst,seed=seedit)),TRUE)
if(!RCT) cs.forest <- try(suppressWarnings(grf::causal_survival_forest(X,Y,W,D,horizon=tau.rmst,seed=seedit)),TRUE)

######################################################
# NOTE: GRF can assign all subjects to *a* treatment
# Restrict to trees with at least 1 split
# Equivalently, restrict to SG's with #|H| < N
# SG cannot include all subjects
# This is for the SG identification stage:
# If no SG found, then H^c=ITT (all subjects)
######################################################

n.max <- length(Y)

# Compute doubly robust scores.
dr.scores <- double_robust_scores(cs.forest)

# Fit trees of depths 1-3
tree1 <- policytree::policy_tree(X, dr.scores, depth = 1)
# Print and plot the tree - action 1 corresponds to control, and 2 treated.
# predicted=1 --> treat rec. = control
# predicted=2 --> treat rec. = treated
#data$predict.action<-predict(tree, X)
data$predict1.node <- predict(tree1,X,type="node.id")
values1 <- stats::aggregate(dr.scores, by = list(leaf.node = data$predict1.node),
                    FUN = function(x) c(mean = mean(x), size=length(x), se = sd(x) / sqrt(length(x))))
values1$diff <- values1$control[,"mean"]-values1$treated[,"mean"]
values1$Nsg <- values1$control[,"size"]
# Only subgroups with control.size > nmin.harm
values1 <- values1[which(values1$control[,"size"] >= n.min),]
values1$depth <- 1.0
if(maxdepth==1) values <- values1
tree2 <- NULL
if(maxdepth >= 2){
tree2 <- policytree::policy_tree(X, dr.scores, depth = 2)
data$predict2.node <- predict(tree2,X,type="node.id")
values2 <- stats::aggregate(dr.scores, by = list(leaf.node = data$predict2.node),
                     FUN = function(x) c(mean = mean(x), size=length(x), se = sd(x) / sqrt(length(x))))
values2$diff <- values2$control[,"mean"]-values2$treated[,"mean"]
values2$Nsg <- values2$control[,"size"]
values2 <- values2[which(values2$control[,"size"]>=n.min),]
values2$depth <- 2.0
}
if(maxdepth==2) values <- rbind(values1,values2)


tree3 <- NULL
# At most depth=3
if(maxdepth==3){
tree3 <- policytree::policy_tree(X, dr.scores, depth = 3)
data$predict3.node <- predict(tree3,X,type="node.id")
values3 <- stats::aggregate(dr.scores, by = list(leaf.node = data$predict3.node),
                     FUN = function(x) c(mean = mean(x), size=length(x), se = sd(x) / sqrt(length(x))))
values3$diff <- values3$control[,"mean"]-values3$treated[,"mean"]
values3$Nsg <- values3$control[,"size"]
values3 <- values3[which(values3$control[,"size"]>=n.min),]
values3$depth <- 3.0
values <- rbind(values1,values2,values3)
}

# Check that 'values' has required columns
required_cols <- c("diff", "depth")
if (!all(required_cols %in% names(values))) {
  stop("The 'values' data.frame must contain columns: ", paste(required_cols, collapse = ", "))
}
# Initialize max.diff as NULL
max.diff <- NULL
# Select max.diff based on criterion
if (identical(sg.criterion, "mDiff")) {
  if (nrow(values) == 0) stop("'values' is empty.")
  loc.max <- which.max(values$diff)
  max.diff <- values[loc.max, , drop = FALSE]
} else if (identical(sg.criterion, "Nsg")) {
  if (!"Nsg" %in% names(values)) stop("'values' must contain column 'Nsg' for Nsg criterion.")
  values_new <- values[values$diff >= dmin.grf, , drop = FALSE]
  if (nrow(values_new) == 0) stop("No rows in 'values' meet the diff >= dmin.grf criterion.")
  loc.max <- which.max(values_new$Nsg)
  max.diff <- values_new[loc.max, , drop = FALSE]
} else {
  stop("Unknown sg.criterion: ", sg.criterion)
}
# Check that max.diff has a 'depth' column
if (!"depth" %in% names(max.diff)) stop("'max.diff' must contain column 'depth'.")
# Assign predict.node and tree based on depth
if (max.diff$depth == 1) {
  data$predict.node <- data$predict1.node
  tree <- tree1
} else if (maxdepth >= 2 && max.diff$depth == 2) {
  data$predict.node <- data$predict2.node
  tree <- tree2
} else if (maxdepth >= 3 && max.diff$depth == 3) {
  data$predict.node <- data$predict3.node
  tree <- tree3
} else {
  stop("No matching depth found for assignment.")
}
if(details){
cat("tau, maxdepth=",c(tau.rmst,maxdepth),"\n")
temp <- values[,c("leaf.node","control","depth")]
print(round(temp,2))
temp <- max.diff[,c("leaf.node","control","depth")]
print(round(temp,2))
}
harm.est <- NULL
sg.harm.id <- NULL
# Restricting to SG less than N
# As GRF can identify ITT population as a
# "tree" meeting dmin criteria;
# However, this is not a subgroup

if(max.diff[,"diff"] >= dmin.grf && max.diff[,"Nsg"] < n.max){
harm.est<-max.diff
data$treat.recommend <- ifelse(data$predict.node == max.diff$leaf.node,0,1)
sg_node <- c(max.diff$leaf.node)
grf_names <- tree$columns
tnodes <- tree$nodes
# Find terminal leaf for identified subgroup
# Search though indices of tnodes
sg_cov<-NULL
sg_cut<-NULL
for(tt in 1:length(tnodes)){
  temp <- tnodes[[tt]]
  if(!temp$is_leaf){
    if(temp$left_child == sg_node | temp$right_child == sg_node){
      sg_cov <- temp$split_variable
      sg_cut <- temp$split_value
    }
  }
}
vmax <- paste0(grf_names[sg_cov]," <= ")
vmax <- paste0(vmax,sg_cut,sep="")

# Find all cuts for selected tree
Vsg_cuts <- NULL
Vsg_names <- NULL
for(tt in 1:length(tnodes)){
  temp <- tnodes[[tt]]
  if(!temp$is_leaf){
    sg_cov <- temp$split_variable
    sg_cut <- round(temp$split_value,2)
    vcut <- paste0(grf_names[sg_cov]," <= ")
    vcut <- paste0(vcut,sg_cut,sep="")
    Vsg_cuts <- c(Vsg_cuts,vcut)
    Vsg_names <- c(Vsg_names,grf_names[sg_cov])
  }
}

# Tree 1
grf_names <- tree1$columns
tnodes <- tree1$nodes
Vsg1_cuts <- NULL
Vsg1_names <- NULL
for(tt in 1:length(tnodes)){
  temp <- tnodes[[tt]]
  if(!temp$is_leaf){
    sg_cov <- temp$split_variable
    sg_cut <- round(temp$split_value,2)
    vcut <- paste0(grf_names[sg_cov]," <= ")
    vcut <- paste0(vcut,sg_cut,sep="")
    Vsg1_cuts <- c(Vsg1_cuts,vcut)
    Vsg1_names <- c(Vsg1_names,grf_names[sg_cov])
  }
}

# Trees 2 and 3 (depth >= 2)
if(maxdepth >=2){
grf_names <- tree2$columns
tnodes <- tree2$nodes
Vsg2_cuts <- NULL
Vsg2_names <- NULL
for(tt in 1:length(tnodes)){
  temp <- tnodes[[tt]]
  if(!temp$is_leaf){
    sg_cov <- temp$split_variable
    sg_cut <- round(temp$split_value,2)
    vcut <- paste0(grf_names[sg_cov]," <= ")
    vcut <- paste0(vcut,sg_cut,sep="")
    Vsg2_cuts <- c(Vsg2_cuts,vcut)
    Vsg2_names <- c(Vsg2_names,grf_names[sg_cov])
  }
}

Vsg3_cuts <- NULL
Vsg3_names <- NULL
if(maxdepth == 3){
grf_names <- tree3$columns
tnodes <- tree3$nodes
for(tt in 1:length(tnodes)){
  temp <- tnodes[[tt]]
  if(!temp$is_leaf){
    sg_cov <- temp$split_variable
    sg_cut <- round(temp$split_value,2)
    vcut <- paste0(grf_names[sg_cov]," <= ")
    vcut <- paste0(vcut,sg_cut,sep="")
    Vsg3_cuts <- c(Vsg3_cuts,vcut)
    Vsg3_names <- c(Vsg3_names,grf_names[sg_cov])
  }
 }
}
}

if(details && max.diff[,"diff"] >= dmin.grf & max.diff[,"Nsg"] < n.max){
cat("GRF subgroup found","\n")
cat("All splits","\n")
print(Vsg_cuts,2)
cat("Terminating node at max.diff (sg.harm.id)","\n")
print(vmax)
}

result <- list(data=data,
grf.gsub=harm.est,
sg.harm.id=vmax,
tree.cuts=Vsg_cuts,
tree.names=unique(Vsg_names),
tree1.cuts=Vsg1_cuts,
tree1.names=unique(Vsg1_names),
tree2.cuts=Vsg2_cuts,
tree2.names=unique(Vsg2_names),
tree3.cuts=Vsg3_cuts,
tree3.names=unique(Vsg3_names),
harm.any=values[which(values$diff>0),],
tree=tree,tau.rmst=tau.rmst,tree1=tree1,tree2=tree2,tree3=tree3)
}
if(max.diff[,"diff"]<dmin.grf | max.diff[,"Nsg"] >= n.max){
if(details){
cat("GRF subgroup NOT found","\n")
}
harm.any <- NULL
cand.sgs <- which(values$diff>0)
if(length(cand.sgs)>0) harm.any<-values[cand.sgs,]
result <- list(data=data,grf.gsub=harm.est,sg.harm.id=sg.harm.id,harm.any=harm.any,tree=tree,tau.rmst=tau.rmst,dmin.grf=dmin.grf,
frac.tau=frac.tau,maxdepth=maxdepth,n.min=n.min)
}
return(result)
}

#' GRF Subgroup Evaluation and Performance Metrics
#'
#' Evaluates the performance of GRF-identified subgroups, including hazard ratios, bias, and predictive values.
#'
#' @param df Data frame containing the analysis data.
#' @param grf.est List. Output from \code{grf.subg.harm.survival}.
#' @param dgm List. Data-generating mechanism (truth) for simulation.
#' @param cox.formula.sim Formula for unadjusted Cox model.
#' @param cox.formula.adj.sim Formula for adjusted Cox model.
#' @param analysis Character. Analysis label (default: "GRF").
#' @param frac.tau Numeric. Fraction of tau for GRF horizon (default: 1.0).
#'
#' @return Data frame with performance metrics (hazard ratios, bias, PPV, NPV, etc.).
#'
#' @importFrom survival coxph Surv
#' @export

grf.estimates.out <- function(df,grf.est=NULL,dgm=NULL,cox.formula.sim=NULL,cox.formula.adj.sim=NULL,analysis="GRF",frac.tau=1.0){
  p.cens<-mean(1-df[,event.name])

  Y<-df[,outcome.name]
  W<-df[,treat.name]
  D<-df[,event.name]
  taumax <- frac.tau * min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))

  if(!is.null(grf.est)) sg.harm.grf<-grf.est$sg.harm.id
  if(is.null(grf.est)) sg.harm.grf<-NULL

    fit <- summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
    hr.itt<-c(fit[1])
    l.itt<-c(fit[3])
    u.itt<-c(fit[4])
    rm("fit")

    fit <- summary(coxph(cox.formula.adj.sim,data=df,robust=FALSE))$conf.int
    hr.adj.itt<-c(fit[1,1])
    l.adj.itt<-c(fit[1,3])
    u.adj.itt<-c(fit[1,4])
    rm("fit")


  if(dgm$model=="null" & !is.null(dgm$grf.harm.true)) stop("For dgm model null grf.harm.true should be null")

  # dgm$grf.harm.true=NULL --> H does not exist

  # If H exists and found something
  if(!is.null(dgm$grf.harm.true) & !is.null(sg.harm.grf)){
    any.H<-1.0
    dfout<-grf.est$data
    # PPV for H membership: Hpred --> treat.reccomend=0 & flag.harm=1
    # Pred H and true H
    aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
    # Pred Hc and true H
    bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
    # Pred H and true Hc
    cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
    # Pred Hc and true Hc
    dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))

    size.H<-with(dfout,sum(treat.recommend==0))
    size.Hc<-with(dfout,sum(treat.recommend==1))


    # hr estimate with true H
    #hr.H.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=TRUE))$conf.int[1])
    #hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=TRUE))$conf.int[1])

    fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=FALSE))$conf.int
    hr.H.true<-c(fit[1])
    l.H.true<-c(fit[3])
    u.H.true<-c(fit[4])
    rm("fit")

    fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=FALSE))$conf.int
    hr.Hc.true<-c(fit[1])
    l.Hc.true<-c(fit[3])
    u.Hc.true<-c(fit[4])
    rm("fit")


    # hr estimate with estimated H
    #hr.H.hat<-c(summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=TRUE))$conf.int[1])
    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=TRUE))$conf.int[1])

    fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=FALSE))$conf.int
    hr.H.hat<-c(fit[1])
    l.H.hat<-c(fit[3])
    u.H.hat<-c(fit[4])
    rm("fit")

    fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=FALSE))$conf.int
    hr.Hc.hat<-c(fit[1])
    l.Hc.hat<-c(fit[3])
    u.Hc.hat<-c(fit[4])
    rm("fit")


    # Bias1 w.r.t hr.H.true
    b1.H<-c(hr.H.hat-hr.H.true)
    b2.H<-c(hr.H.hat-dgm$hr.H.true)

    b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
    b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)

    ppv<-aa/(aa+bb)
    npv<-dd/(cc+dd)
    specificity<-dd/(bb+dd)
    sensitivity<-aa/(aa+cc)

    #cat("PPV and NPV=",c(ppv,npv),"\n")
    #cat("Sensitivity and Specificity=",c(specificity,sensitivity),"\n")

    # Identify truth in dgm
    #found.1<-sum(grepl(dgm$grf.harm.true[1],sg.harm.grf))
    #found.2<-sum(grepl(dgm$grf.harm.true[2],sg.harm.grf))
    #found.both<-as.numeric(c(found.1+found.2)==2)

    found.1<-found.2<-found.both<-NA

    #found.al3<-as.numeric(length(sg.harm.grf)>=3)
    found.al3<-NA
  }

  # If H exists and did NOT find something
  if(!is.null(dgm$grf.harm.true) & is.null(sg.harm.grf)){
    any.H<-0
    size.H<-0
    size.Hc<-nrow(df)

    # hr estimate with true H
    #hr.H.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=TRUE))$conf.int[1])
    #hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=TRUE))$conf.int[1])

    fit <- summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=FALSE))$conf.int
    hr.H.true<-c(fit[1])
    l.H.true<-c(fit[3])
    u.H.true<-c(fit[4])
    rm("fit")

    fit <- summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=FALSE))$conf.int
    hr.Hc.true<-c(fit[1])
    l.Hc.true<-c(fit[3])
    u.Hc.true<-c(fit[4])
    rm("fit")


    # hr estimate with estimated H, in this case, does not exist
    #hr.H.hat<-NA
    # Define HC as ITT
    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=TRUE))$conf.int[1])

    # hr estimate with estimated H, in this case, does not exist
    hr.H.hat<-NA
    l.H.hat<-NA
    u.H.hat<-NA

    # Define HC as ITT
    # hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int[1])

    fit <- summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
    hr.Hc.hat<-c(fit[1])
    l.Hc.hat<-c(fit[3])
    u.Hc.hat<-c(fit[4])
    rm("fit")



    # Bias1 w.r.t hr.H.true
    b1.H<-NA
    b2.H<-NA

    b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
    b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)

    # no H found, treat.recommend==1 for all
    # Pred H and true H
    #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
    # NO H found
    aa<-0
    # Pred Hc and true H
    #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
    # All treat.recommend=1
    bb<-with(df,sum(flag.harm==1))
    # Pred H and true Hc
    #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
    cc<-0
    # Pred Hc and true Hc
    #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
    dd<-with(df,sum(flag.harm==0))

    ppv<-aa/(aa+bb)
    npv<-dd/(cc+dd)
    specificity<-dd/(bb+dd)
    sensitivity<-ifelse(aa==0,0,aa/(aa+cc))

    # Identify truth in dgm
    #found.1<-0
    #found.2<-0
    #found.both<-0

    found.1<-found.2<-found.both<-NA


    found.al3<-0

  }


  # If H does NOT exist and found something
  if(is.null(dgm$grf.harm.true) & !is.null(sg.harm.grf)){
    any.H<-1.0
    dfout<-grf.est$data
    # PPV for H membership: Hpred --> treat.reccomend=0 & flag.harm=1
    # Pred H and true H
    #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
    # flag.harm=0 for all
    aa<-0
    # Pred Hc and true H
    #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
    bb<-0
    # Pred H and true Hc
    #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
    cc<-with(dfout,sum(treat.recommend==0))
    # Pred Hc and true Hc
    #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
    dd<-with(dfout,sum(treat.recommend==1))

    size.H<-with(dfout,sum(treat.recommend==0))
    size.Hc<-with(dfout,sum(treat.recommend==1))

    # hr estimate with true H
    #hr.H.true<-NA
    # ITT
    #hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=df,robust=TRUE))$conf.int[1])

    # hr estimate with true H
    hr.H.true<-NA
    l.H.true<-NA
    u.H.true<-NA

    # no H --> ITT
    fit <- summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
    hr.Hc.true<-c(fit[1])
    l.Hc.true<-c(fit[3])
    u.Hc.true<-c(fit[4])
    rm("fit")


    # hr estimate with estimated H
    #hr.H.hat<-c(summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=TRUE))$conf.int[1])
    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=TRUE))$conf.int[1])

    # hr estimate with estimated H
    fit <- summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=FALSE))$conf.int
    hr.H.hat<-c(fit[1])
    l.H.hat<-c(fit[3])
    u.H.hat<-c(fit[4])
    rm("fit")

    fit <- summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=FALSE))$conf.int
    hr.Hc.hat<-c(fit[1])
    l.Hc.hat<-c(fit[3])
    u.Hc.hat<-c(fit[4])
    rm("fit")

    # Bias1 w.r.t hr.H.true
    b1.H<-NA
    b2.H<-NA

    b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
    b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)

    ppv<-NA
    npv<-dd/(cc+dd)
    specificity<-dd/(bb+dd)
    sensitivity<-aa/(aa+cc)

    # Identify truth in dgm
    # fs.harm.true will be null
    # Still want to track how often v1.1 and v4.1 are found
    #found.1<-sum(grepl(c("v1=1"),sg.harm.grf))
    #found.2<-sum(grepl(c("v4=1"),sg.harm.grf))
    #found.both<-as.numeric(c(found.1+found.2)==2)

    found.1<-found.2<-found.both<-NA

    #found.al3<-as.numeric(length(sg.harm.grf)>=3)
    found.al3<-NA
    }

  # If H does NOT exist and did NOT find something
  # Define hat(Hc)=ITT
  if(is.null(dgm$grf.harm.true) & is.null(sg.harm.grf)){
    any.H<-0
    # **flag.harm=0 for all**
    # Pred H and true H
    #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
    # flag.harm=0 for all
    aa<-0
    # Pred Hc and true H
    #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
    bb<-0
    # Pred H and true Hc
    #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
    cc<-0
    # Pred Hc and true Hc
    #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
    # Did NOT find anything --> treat.recommend=1 for all
    dd<-nrow(df)

    size.H<-0
    size.Hc<-nrow(df)

    # hr estimate with true H
    #hr.H.true<-NA
    # ITT
    #hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=df,robust=TRUE))$conf.int[1])

    # hr estimate with true H
    hr.H.true<-NA
    l.H.true<-NA
    u.H.true<-NA

    # ITT
    fit <- summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
    hr.Hc.true<-c(fit[1])
    l.Hc.true<-c(fit[3])
    u.Hc.true<-c(fit[4])
    rm("fit")


    # hr estimate with estimated H
    #hr.H.hat<-NA
    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=TRUE))$conf.int[1])

    # hr estimate with estimated H
    hr.H.hat<-NA
    l.H.hat<-NA
    u.H.hat<-NA

    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int[1])
    fit <- summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
    hr.Hc.hat<-c(fit[1])
    l.Hc.hat<-c(fit[3])
    u.Hc.hat<-c(fit[4])
    rm("fit")


    # Bias1 w.r.t hr.H.true
    b1.H<-NA
    b2.H<-NA

    b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
    b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)

    ppv<-NA
    npv<-dd/(cc+dd)
    specificity<-dd/(bb+dd)
    sensitivity<-NA

    # Identify truth in dgm
    # fs.harm.true will be null
    # Still want to track how often v1.1 and v4.1 are found

    # I don't think these are correct
    #found.1<-sum(grepl(c("v1=1"),sg.harm.grf))
    #found.2<-sum(grepl(c("v4=1"),sg.harm.grf))
    #found.both<-as.numeric(c(found.1+found.2)==2)
    found.1<-found.2<-found.both<-NA
    #found.al3<-as.numeric(length(sg.harm.grf)>=3)
    found.al3<-NA
  }

  df.res <- data.frame(any.H,size.H,size.Hc,ppv,npv,specificity,sensitivity,found.1,found.2,found.both,found.al3,hr.H.true,
                     hr.Hc.true,hr.H.hat,hr.Hc.hat,b1.H,b2.H,b1.Hc,b2.Hc,p.cens,analysis,taumax,hr.itt,l.itt,u.itt,
                     hr.adj.itt,l.adj.itt,u.adj.itt,
                     l.H.true,u.H.true,
                     l.Hc.true,u.Hc.true,
                     l.H.hat,u.H.hat,
                     l.Hc.hat,u.Hc.hat)
  return(df.res)
}


