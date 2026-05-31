suppressMessages({library(survival); library(data.table)})
R <- "fs_extract/R"
for (f in c(file.path(R, c("subgroup_consistency_helpers.R","dina.R","dina_df.R")),
            "dina_subgroup.R","forestsearch_helpers.R")) source(f, local=FALSE)
recompute <- function(labels, df){ m<-rep(TRUE,nrow(df)); for(lab in labels){
  cl<-gsub("^!?\\{(.*)\\}$","\\1",lab); op<-regmatches(cl,regexpr("(<=|>=|<|>)",cl))
  p<-strsplit(cl,"(<=|>=|<|>)")[[1]]; v<-trimws(p[1]); val<-as.numeric(trimws(p[2]))
  m<-m & switch(op,"<="=df[[v]]<=val,">="=df[[v]]>=val,"<"=df[[v]]<val,">"=df[[v]]>val)}; m }

## PART 1: planted data, real gaussian dina fit, scan floor to elicit depth-2
set.seed(11); n<-1500
df1<-data.frame(x1=runif(n),x2=runif(n),x3=runif(n),x4=runif(n),x5=runif(n),W=rbinom(n,1,.5))
df1$Y<-df1$W*(df1$x1+df1$x2)+rnorm(n); covs1<-paste0("x",1:5)
fit1<-dina(df=df1,outcome="Y",treatment="W",covariates=covs1,family="gaussian",seed=1L)
# scan m_diff to find a floor where maxSG selects a depth-2 conjunction
Fgrid<-seq(1.20,1.70,by=0.02); Fstar<-NA
for(F in Fgrid){ s<-dina_subgroup(fit1,df1,covs1,m_diff=F,n_min=60L,max_depth=2L,
       sg_focus="maxSG"); if(isTRUE(s$found)&&s$depth==2L){Fstar<-F;break} }
cat("======== PART 1  synthetic gaussian, depth-2-eliciting floor =========\n")
cat("floor where maxSG picks a conjunction: m_diff =", Fstar, "\n")
common1<-list(df=df1,df.predict=NULL,df.test=NULL,confounders.name=covs1,
  outcome.name="Y",event.name="W",treat.name="W",id.name="id",
  outcome_type="continuous",hr.threshold=Fstar,n.min=60L,sg_focus="maxSG",
  selection_rule="neighborhood",effect_neighborhood=0.10,dina_res=fit1,
  seedit=1L,details=FALSE)
p1<-do.call(.forestsearch_dina_select,c(common1,list(dina_args=list(family="gaussian",max_depth=1L))))
p2<-do.call(.forestsearch_dina_select,c(common1,list(dina_args=list(family="gaussian",max_depth=2L))))
cat(sprintf("max_depth=1: found=%s  sg.harm=%s  harm_n=%d\n",p1$found,
    paste(p1$sg.harm,collapse=" & "),sum(p1$df.est$treat.recommend==0L)))
cat(sprintf("max_depth=2: found=%s  sg.harm=%s  harm_n=%d\n",p2$found,
    paste(p2$sg.harm,collapse=" & "),sum(p2$df.est$treat.recommend==0L)))
cat("[A] depth-1 single label via real path: ",isTRUE(p1$found)&&length(p1$sg.harm)==1L,"\n")
cat("[B] depth-2 conjunction (2 labels) via real path: ",isTRUE(p2$found)&&length(p2$sg.harm)==2L,"\n")
cat("[C] conjunction uses only planted x1/x2: ",all(gsub("\\{| .*","",p2$sg.harm)%in%c("x1","x2")),"\n")
cat("[D] 2-label sg.harm AND-composes to stored treat.recommend==0: ",
    all(recompute(p2$sg.harm,df1)==(p2$df.est$treat.recommend==0L)),"\n")
cat("[E] grp.consistency shape intact & algorithm=='dina': ",
    all(c("out_sg","sg.harm","sg.harm.id","df_flag","algorithm")%in%names(p2$grp.consistency))
    && p2$grp.consistency$algorithm=="dina","\n")
cat("[F] depth-2 harm subgroup >= depth-1 under maxSG: ",
    sum(p2$df.est$treat.recommend==0L)>=sum(p1$df.est$treat.recommend==0L),"\n")
cat("[G] df_flag rows == nrow(df): ",nrow(p2$grp.consistency$df_flag)==nrow(df1),"\n")

## PART 2: gbsg cox real, scan floor to elicit depth-2 too
cat("\n======== PART 2  gbsg cox (real data) ================================\n")
g<-survival::gbsg; covsg<-c("age","size","nodes","pgr","er")
fitg<-dina(df=g,outcome="rfstime",treatment="hormon",covariates=covsg,status="status",family="cox",seed=1L)
Fg<-seq(quantile(as.numeric(coef(fitg)[1]+as.matrix(g[covsg])%*%coef(fitg)[-1]),0.5),
        max(as.numeric(coef(fitg)[1]+as.matrix(g[covsg])%*%coef(fitg)[-1]))-0.05,length.out=25)
Fgs<-NA
for(F in Fg){ s<-dina_subgroup(fitg,g,covsg,m_diff=F,n_min=60L,max_depth=2L,sg_focus="maxSG")
  if(isTRUE(s$found)&&s$depth==2L){Fgs<-F;break} }
cat("floor eliciting depth-2 on gbsg (maxSG): m_diff =",
    if(is.na(Fgs)) "none in range (depth-1 always optimal)" else round(Fgs,3),"\n")
useF<-if(is.na(Fgs)) as.numeric(quantile(as.numeric(coef(fitg)[1]+as.matrix(g[covsg])%*%coef(fitg)[-1]),0.7)) else Fgs
commong<-list(df=g,df.predict=NULL,df.test=NULL,confounders.name=covsg,outcome.name="rfstime",
  event.name="status",treat.name="hormon",id.name="pid",outcome_type="survival",
  hr.threshold=exp(useF),n.min=60L,sg_focus="maxSG",selection_rule="neighborhood",
  effect_neighborhood=0.10,dina_res=fitg,seedit=1L,details=FALSE)
g2<-do.call(.forestsearch_dina_select,c(commong,list(dina_args=list(family="cox",max_depth=2L))))
cat(sprintf("max_depth=2: found=%s  depth=%d  sg.harm=%s  harm_n=%d\n",g2$found,
    length(g2$sg.harm),paste(g2$sg.harm,collapse=" & "),
    if(isTRUE(g2$found)) sum(g2$df.est$treat.recommend==0L) else 0))
cat("[H] gbsg real cox path runs, found, composes: ",
    isTRUE(g2$found)&&all(recompute(g2$sg.harm,g)==(g2$df.est$treat.recommend==0L)),"\n")
