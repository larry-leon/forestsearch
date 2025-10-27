
sg_UB <- function(x,z,analysis){
  rej <- 100*round(mean(x<=1, na.rm=TRUE),3)
  rej <- paste0(rej,"%")
  medx <- median(x,na.rm=TRUE)
  xlow <- quantile(x,c(0.025),na.rm=TRUE)
  xup <- quantile(x,c(0.975),na.rm=TRUE)
  # IQR
  xiqr <- c(quantile(x,c(0.75),na.rm=TRUE) - quantile(x,c(0.25),na.rm=TRUE))
  bb <- c(z,NA,medx,xlow,xup,xiqr)
  names(bb) <- c("N","UB<1","est","low","hi","se")
  bb <- as.data.frame(t(bb))
  bb$Subgroup <- c(analysis)  
  bb[,c("UB<1")] <- c(rej)
  bb<- bb[,c("Subgroup","N","UB<1","est","low","hi","se")]
  return(bb)
}


getSG_dfUB <- function(x,y,z){
  sg_tabs <- NULL
  sg_labels <- colnames(x)
  for(gg in 1:length(sg_labels)){
    ## bb is stratified analysis (1)
    # % HR<1 
    bb <- sg_UB(x[,gg],z[gg],"Stratified")
    cc <- sg_UB(y[,gg],z[gg],"Un-Stratified")
    # Add header 
    aa <- bb 
    # Change "Analysis" to subgroup name
    aa$Subgroup <- sg_labels[gg]
    aa[,c("est","low","hi","se")] <- NA
    aa[,c("N","UB<1")] <- ""
    sg_tab <- rbind(aa,bb,cc)
    sg_tabs <- rbind(sg_tabs,sg_tab)
  }
  return(sg_tabs)
}


sg_hr <- function(x,z,ub,analysis,alpha=0.01,est.threshold=0.80){
  rej <- 100*round(mean(ub<=1, na.rm=TRUE),3)
  rej <- paste0(rej,"%")
  
  est1 <- 100*round(mean(x<=est.threshold, na.rm=TRUE),2)
  est1 <- paste0(est1,"%")
  label_threshold <- paste("HR<",est.threshold)
  
  medx <- median(x,na.rm=TRUE)
  if(alpha==0.025){
    xlow <- quantile(x,c(0.025),na.rm=TRUE)
    xup <- quantile(x,c(0.975),na.rm=TRUE)
  }
  if(alpha==0.01){
    xlow <- quantile(x,c(0.01),na.rm=TRUE)
    xup <- quantile(x,c(0.99),na.rm=TRUE)
  }
  
  # IQR
  xiqr <- c(quantile(x,c(0.75),na.rm=TRUE) - quantile(x,c(0.25),na.rm=TRUE))
  #se_emp <- sqrt(var(x,na.rm=TRUE))
  bb <- c(z,NA,medx,xlow,xup,xiqr)
  names(bb) <- c("N","UB<1","est","low","hi","se")
  bb <- as.data.frame(t(bb))
  bb$Subgroup <- c(analysis)  
  bb[,c("UB<1")] <- c(rej)
  bb[,c(label_threshold)] <- c(est1)
  bb<- bb[,c("Subgroup","N",label_threshold,"UB<1","est","low","hi","se")]
  return(bb)
}


sg_ub2 <- function(x,z,ub,analysis,alpha=0.01,est.threshold=0.80){
  rej <- 100*round(mean(ub<=1, na.rm=TRUE),3)
  rej <- paste0(rej,"%")
  
  est1 <- 100*round(mean(x<=est.threshold, na.rm=TRUE),2)
  est1 <- paste0(est1,"%")
  label_threshold <- paste("HR<",est.threshold)
  
  # x=ub
  medx <- median(ub,na.rm=TRUE)
  if(alpha==0.025){
    xlow <- quantile(ub,c(0.025),na.rm=TRUE)
    xup <- quantile(ub,c(0.975),na.rm=TRUE)
  }
  if(alpha==0.01){
    xlow <- quantile(ub,c(0.01),na.rm=TRUE)
    xup <- quantile(ub,c(0.99),na.rm=TRUE)
  }
  
  # IQR
  xiqr <- c(quantile(ub,c(0.75),na.rm=TRUE) - quantile(ub,c(0.25),na.rm=TRUE))
  #se_emp <- sqrt(var(ub,na.rm=TRUE))
  bb <- c(z,NA,medx,xlow,xup,xiqr)
  names(bb) <- c("N","UB<1","est","low","hi","se")
  bb <- as.data.frame(t(bb))
  bb$Subgroup <- c(analysis)  
  bb[,c("UB<1")] <- c(rej)
  bb[,c(label_threshold)] <- c(est1)
  bb<- bb[,c("Subgroup","N",label_threshold,"UB<1","est","low","hi","se")]
  return(bb)
}




getSG_dfhrTWO <- function(x,y,z,ubx,uby,analysisx="x",analysisy="y"){
  sg_tabs <- NULL
  sg_labels <- colnames(x)
  for(gg in 1:length(sg_labels)){
    ## bb is stratified analysis (1)
    # % HR<1 
    bb <- sg_hr(x[,gg],z[gg],ubx[,gg],analysisx)
    cc <- sg_hr(y[,gg],z[gg],uby[,gg],analysisy)
    # Add header 
    aa <- bb 
    # Change "Analysis" to subgroup name
    aa$Subgroup <- sg_labels[gg]
    aa[,c("est","low","hi","se")] <- NA
    aa[,c("N","HR<1","UB<1")] <- ""
    sg_tab <- rbind(aa,bb,cc)
    sg_tabs <- rbind(sg_tabs,sg_tab)
  }
  return(sg_tabs)
}


getSG_dfhrSIX <- function(z,x1,ubx1,x2,ubx2,x3,ubx3,
                          x4,ubx4,x5,ubx5,x6,ubx6,analysisx1="x1",analysisx2="x2",
                          analysisx3="x3",analysisx4="x4",
                          analysisx5="x5",analysisx6="x6",which_sgs = c("All Patients"),alpha=0.01,est.threshold=0.80){
  sg_tabs <- NULL
  
  label_threshold <- paste("HR<",est.threshold)
  
  sg_labels <- colnames(x1)
  
  for(sgALL in 1:length(which_sgs)){
    ## bb is stratified analysis (1)
    # % HR<1 
    sg <- sg_labels %in% which_sgs[sgALL]  
    
    bb <- sg_hr(x1[,sg],z[sg],ubx1[,sg],analysisx1,alpha=alpha,est.threshold=est.threshold)
    cc <- sg_hr(x2[,sg],z[sg],ubx2[,sg],analysisx2,alpha=alpha,est.threshold=est.threshold)
    dd <- sg_hr(x3[,sg],z[sg],ubx3[,sg],analysisx3,alpha=alpha,est.threshold=est.threshold)
    ee <- sg_hr(x4[,sg],z[sg],ubx4[,sg],analysisx4,alpha=alpha,est.threshold=est.threshold)
    ff <- sg_hr(x5[,sg],z[sg],ubx5[,sg],analysisx5,alpha=alpha,est.threshold=est.threshold)
    gg <- sg_hr(x6[,sg],z[sg],ubx6[,sg],analysisx6,alpha=alpha,est.threshold=est.threshold)
    # Add header 
    aa <- bb 
    # Change "Analysis" to subgroup name
    aa$Subgroup <- sg_labels[sg]
    aa[,c("est","low","hi","se")] <- NA
    aa[,c("N",label_threshold,"UB<1")] <- ""
    sg_tab <- rbind(aa,bb,cc,dd,ee,ff,gg)
    sg_tabs <- rbind(sg_tabs,sg_tab)
  }
  return(sg_tabs)
}


getSG_dfubSIX <- function(z,x1,ubx1,x2,ubx2,x3,ubx3,
                          x4,ubx4,x5,ubx5,x6,ubx6,analysisx1="x1",analysisx2="x2",
                          analysisx3="x3",analysisx4="x4",
                          analysisx5="x5",analysisx6="x6",which_sgs = c("All Patients"),alpha=0.01,est.threshold=0.80){
  sg_tabs <- NULL
  
  label_threshold <- paste("HR<",est.threshold)
  
  sg_labels <- colnames(x1)
  
  for(sgALL in 1:length(which_sgs)){
    ## bb is stratified analysis (1)
    # % HR<1 
    sg <- sg_labels %in% which_sgs[sgALL]  
    
    bb <- sg_ub2(x1[,sg],z[sg],ubx1[,sg],analysisx1,alpha=alpha,est.threshold=est.threshold)
    cc <- sg_ub2(x2[,sg],z[sg],ubx2[,sg],analysisx2,alpha=alpha,est.threshold=est.threshold)
    dd <- sg_ub2(x3[,sg],z[sg],ubx3[,sg],analysisx3,alpha=alpha,est.threshold=est.threshold)
    ee <- sg_ub2(x4[,sg],z[sg],ubx4[,sg],analysisx4,alpha=alpha,est.threshold=est.threshold)
    ff <- sg_ub2(x5[,sg],z[sg],ubx5[,sg],analysisx5,alpha=alpha,est.threshold=est.threshold)
    gg <- sg_ub2(x6[,sg],z[sg],ubx6[,sg],analysisx6,alpha=alpha,est.threshold=est.threshold)
    # Add header 
    aa <- bb 
    # Change "Analysis" to subgroup name
    aa$Subgroup <- sg_labels[sg]
    aa[,c("est","low","hi","se")] <- NA
    aa[,c("N",label_threshold,"UB<1")] <- ""
    sg_tab <- rbind(aa,bb,cc,dd,ee,ff,gg)
    sg_tabs <- rbind(sg_tabs,sg_tab)
  }
  return(sg_tabs)
}



getSG_dfhrTHREE <- function(z,x1,ubx1,x2,ubx2,x3,ubx3,
                            x4,ubx4,x5,ubx5,x6,ubx6,analysisx1="x1",analysisx2="x2",
                            analysisx3="x3",which_sgs = c("All Patients"),alpha=0.01,est.threshold=0.80){
  sg_tabs <- NULL
  
  label_threshold <- paste("HR<",est.threshold)
  
  sg_labels <- colnames(x1)
  
  for(sgALL in 1:length(which_sgs)){
    ## bb is stratified analysis (1)
    # % HR<1 
    sg <- sg_labels %in% which_sgs[sgALL]  
    
    bb <- sg_hr(x1[,sg],z[sg],ubx1[,sg],analysisx1,alpha=alpha,est.threshold=est.threshold)
    cc <- sg_hr(x2[,sg],z[sg],ubx2[,sg],analysisx2,alpha=alpha,est.threshold=est.threshold)
    dd <- sg_hr(x3[,sg],z[sg],ubx3[,sg],analysisx3,alpha=alpha,est.threshold=est.threshold)
    # Add header 
    aa <- bb 
    # Change "Analysis" to subgroup name
    aa$Subgroup <- sg_labels[sg]
    aa[,c("est","low","hi","se")] <- NA
    aa[,c("N",label_threshold,"UB<1")] <- ""
    sg_tab <- rbind(aa,bb,cc,dd)
    sg_tabs <- rbind(sg_tabs,sg_tab)
  }
  return(sg_tabs)
}


getSG_dfubTHREE <- function(z,x1,ubx1,x2,ubx2,x3,ubx3,
                            x4,ubx4,x5,ubx5,x6,ubx6,analysisx1="x1",analysisx2="x2",
                            analysisx3="x3", which_sgs = c("All Patients"),alpha=0.01,est.threshold=0.80){
  sg_tabs <- NULL
  
  label_threshold <- paste("HR<",est.threshold)
  
  sg_labels <- colnames(x1)
  
  for(sgALL in 1:length(which_sgs)){
    ## bb is stratified analysis (1)
    # % HR<1 
    sg <- sg_labels %in% which_sgs[sgALL]  
    
    bb <- sg_ub2(x1[,sg],z[sg],ubx1[,sg],analysisx1,alpha=alpha,est.threshold=est.threshold)
    cc <- sg_ub2(x2[,sg],z[sg],ubx2[,sg],analysisx2,alpha=alpha,est.threshold=est.threshold)
    dd <- sg_ub2(x3[,sg],z[sg],ubx3[,sg],analysisx3,alpha=alpha,est.threshold=est.threshold)
    # Add header 
    aa <- bb 
    # Change "Analysis" to subgroup name
    aa$Subgroup <- sg_labels[sg]
    aa[,c("est","low","hi","se")] <- NA
    aa[,c("N",label_threshold,"UB<1")] <- ""
    sg_tab <- rbind(aa,bb,cc,dd)
    sg_tabs <- rbind(sg_tabs,sg_tab)
  }
  return(sg_tabs)
}



getSG_dfhrONE <- function(z,x1,ubx1,x2,ubx2,x3,ubx3,
                            x4,ubx4,x5,ubx5,x6,ubx6,analysisx1="x1",
                          which_sgs = c("All Patients"),alpha=0.01,est.threshold=0.80){
  sg_tabs <- NULL
  
  label_threshold <- paste("HR<",est.threshold)
  
  sg_labels <- colnames(x1)
  
  for(sgALL in 1:length(which_sgs)){
    ## bb is stratified analysis (1)
    # % HR<1 
    sg <- sg_labels %in% which_sgs[sgALL]  
    
    bb <- sg_hr(x1[,sg],z[sg],ubx1[,sg],analysisx1,alpha=alpha,est.threshold=est.threshold)
    # Add header 
    aa <- bb 
    # Change "Analysis" to subgroup name
    aa$Subgroup <- sg_labels[sg]
    aa[,c("est","low","hi","se")] <- NA
    aa[,c("N",label_threshold,"UB<1")] <- ""
    sg_tab <- rbind(aa,bb)
    sg_tabs <- rbind(sg_tabs,sg_tab)
  }
  return(sg_tabs)
}


getSG_dfubONE <- function(z,x1,ubx1,x2,ubx2,x3,ubx3,
                            x4,ubx4,x5,ubx5,x6,ubx6,analysisx1="x1", which_sgs = c("All Patients"),alpha=0.01,est.threshold=0.80){
  sg_tabs <- NULL
  
  label_threshold <- paste("HR<",est.threshold)
  
  sg_labels <- colnames(x1)
  
  for(sgALL in 1:length(which_sgs)){
    ## bb is stratified analysis (1)
    # % HR<1 
    sg <- sg_labels %in% which_sgs[sgALL]  
    
    bb <- sg_ub2(x1[,sg],z[sg],ubx1[,sg],analysisx1,alpha=alpha,est.threshold=est.threshold)
    # Add header 
    aa <- bb 
    # Change "Analysis" to subgroup name
    aa$Subgroup <- sg_labels[sg]
    aa[,c("est","low","hi","se")] <- NA
    aa[,c("N",label_threshold,"UB<1")] <- ""
    sg_tab <- rbind(aa,bb)
    sg_tabs <- rbind(sg_tabs,sg_tab)
  }
  return(sg_tabs)
}
