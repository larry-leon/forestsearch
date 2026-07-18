#MONET1 study code

library(survival)

#caculate the upper bound of the hazard ratio
#D: dataset; the first three columns are aval1(survival time), CNSR(censoring indicator) and treatment(treatment indicator) 
#and the rest are the subgroup indicators
#B: the number of bootstrap; r: the tuning parameter; 
#n: the sample size of the dataset; level: the level of the upper bound
BC=function(D,B,r,n,level)
{ 
  X_bar=0
  Y_bar=0
  cc=dim(D)[2]-3
  for(l in 1:cc) #calculate the observed negative log-hazard ratio for each subgroup
  {
    X=D[D[,3+l]==1,]
    Y=D[D[,3+l]==0,]
    X_bar[l]=-coxph(Surv(aval1,CNSR)~treatment,data=X)$coef #change to the negative log-hazard ratio to make it consistent to the framework of the algorithm
    Y_bar[l]=-coxph(Surv(aval1,CNSR)~treatment,data=Y)$coef
  }
  MM=max(X_bar,Y_bar) #the observed negative log-hazard ratio of the best selected subgroup
  c_n=0
  d_n=0
  for(l in 1:cc) #calculate the correction term with r as the tuning
  {
    c_n[l]=(n^(1/2)-n^r)*(MM-X_bar[l])/sqrt(n)
    d_n[l]=(n^(1/2)-n^r)*(MM-Y_bar[l])/sqrt(n)
  }
  BS=0
  for(i in 1:B) #bootstrap
  {
    index=sample(1:n,n,replace=TRUE) #pair bootstrap the sample
    DD=D[index,]
    X_b=0
    Y_b=0
    for(l in 1:cc) #calculate the observed negative log hazard ratio for the bootstrap sample
    {
      XB=DD[DD[,3+l]==1,]
      YB=DD[DD[,3+l]==0,]
      X_b[l]=-coxph(Surv(aval1,CNSR)~treatment,data=XB)$coef 
      Y_b[l]=-coxph(Surv(aval1,CNSR)~treatment,data=YB)$coef
    }
    BS[i]=sqrt(n)*(max(X_b+c_n,Y_b+d_n)-MM) #calculate the quantity
  }
  return(exp(quantile(BS,1-level)/sqrt(n)-max(X_bar,Y_bar))) #the upper bound of the hazard ratio of the best subgroup
}

#caculate the bias-reduced estimate for the best subgroup
#the inputs are defined the same as those in CB
CB=function(D,B,r,n)
{
  X_bar=0
  Y_bar=0
  cc=dim(D)[2]-3
  for(l in 1:cc)
  {
    X=D[D[,3+l]==1,]
    Y=D[D[,3+l]==0,]
    X_bar[l]=-coxph(Surv(aval1,CNSR)~treatment,data=X)$coef 
    Y_bar[l]=-coxph(Surv(aval1,CNSR)~treatment,data=Y)$coef
  }
  MM=max(X_bar,Y_bar)
  c_n=0
  d_n=0
  for(l in 1:cc)
  {
    c_n[l]=(n^(1/2)-n^r)*(MM-X_bar[l])/sqrt(n)
    d_n[l]=(n^(1/2)-n^r)*(MM-Y_bar[l])/sqrt(n)
  }
  BS=0
  for(i in 1:B)
  {
    index=sample(1:n,n,replace=TRUE)
    DD=D[index,]
    X_b=0
    Y_b=0
    for(l in 1:cc)
    {
      XB=DD[DD[,3+l]==1,]
      YB=DD[DD[,3+l]==0,]
      X_b[l]=-coxph(Surv(aval1,CNSR)~treatment,data=XB)$coef 
      Y_b[l]=-coxph(Surv(aval1,CNSR)~treatment,data=YB)$coef
    }
    BS[i]=(max(X_b+c_n,Y_b+d_n)-MM)
  }
  return(exp(mean(BS)-max(X_bar,Y_bar))) #the bias-reduced for the hazard ratio of the best subgroup
}

#read the data
D=read.csv('syntheticdata.csv')

#basic setup
B=2000
r=0.03
n=1090
level=0.05
bc=0
cb=0

#calculate the upper bound and bias-reduced estimate for the best subgroup with different number of subgroup candidates
for(i in 1:8)
{ 
  D2=D[,1:(3+i)] #consider 2i subgroups
  bc[i]=BC(D2,B,r,n,level) #calculate the upper bound
  cb[i]=CB(D2,B,r,n) #calculate the bias-reduced estimate
}

