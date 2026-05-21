# Changing the path
# Will need to do this for every session
.libPaths("C:/Users/leolarr2/OneDrive - Merck Sharp & Dohme, Corp/documents/MyRlibrary")

rm(list=ls())


library(speff2trial)
library(Hmisc)
library(ggplot2)


# df_analysis24 dataset
load("outdata_analysis/df_week24.Rdata")


# Plot longitudinal
p<-ggplot(data=df_long,aes(x=ady,y=chg,group=subjid))
p<-p+xlab("study day")+ylab("Change-from-baseline")
p<-p+geom_point(position=position_jitter(height=0.03,width=0))+stat_smooth(aes(group=1))+facet_grid(. ~ trtp)

plot(p)


# Rename as df0

dfa<-df_analysis24

dfa<-within(dfa,{
offtreat_day<-ifelse(offtrt_day<=168,offtrt_day,168)
})

# First try speff2trial approach

table(dfa$ontrtfl_24)
table(dfa$miss_24,df0$trtp)
# Missingness matches explorations doc


bl_summary<-summary(trtp ~ y_0+sex+age+chg_4_impute+chg_8_impute+last_treat8+last_treat12+last_treat16+last_treat20+miss_24+chg_24+Hchg_24,method="reverse", data=dfa, test=TRUE)

print(bl_summary,digits=2,prtest=c('P'))


mod1<-speff(Hchg_24 ~ y_0+age+female,trt.id="treat_45mg",data=dfa)
print(summary(mod1))

mod2<-speff(Hchg_24 ~ y_0+age+female+y_4_impute, postrandom=c("y_4_impute"), trt.id="treat_45mg",data=dfa)
print(summary(mod2))

mod3<-speff(Hchg_24 ~ y_0+age+female+y_4_impute, postrandom=c("y_4_impute"), trt.id="treat_45mg",data=dfa,force.in=c(1,2,3))
print(summary(mod3))


mod4<-speff(Hchg_24 ~ y_0+age+female+y_8_impute, 
postrandom=c("y_8_impute"), 
trt.id="treat_45mg",data=dfa)
print(summary(mod4))

mod5<-speff(Hchg_24 ~ y_0+age+female+y_8_impute, 
postrandom=c("y_8_impute"), 
trt.id="treat_45mg",data=dfa,force.in=c(1,2,3))
print(summary(mod5))

#mod4<-speff(chg_24 ~ y_0+age+female+y_8_impute+last_treat16, postrandom=c("y_8_impute","last_treat16"), trt.id="treat_45mg",data=df0)
#mod5<-speff(chg_24 ~ y_0+age+female+y_8_impute+last_treat8, postrandom=c("y_8_impute","last_treat8"), trt.id="treat_45mg",data=df0)
#mod6<-speff(Hchg_24 ~ y_0+age+female+y_8_impute+, postrandom=c("y_8_impute","last_treat8"), trt.id="treat_45mg",data=df0)


library(knitr)
library(gee)
#library(tab)
library("kableExtra")
require(reshape2)
require(reshape)
library("Matching")
#library("rms") # For MAR1 robust variance estimation



ep.name<-"Hchg_24" 
or.names<-c("y_0","age","female","y_8_impute")
ps.names<-c("y_0","y_8_impute")
miss.name<-c("miss_H24")
treat.name<-c("treat_45mg")
direction.superiority<-c("negative") # Negative values in favor of experimental

dr.fit<-get.dr.analyses(df=dfa,ep.name=ep.name,or.names=or.names,ps.names=ps.names,miss.name=miss.name,
                treat.name=treat.name,direction.superiority=direction.superiority,details=FALSE)

print(dr.fit$dr1.out[c("m0","m1","d","sed","dL","dU","pval.1sided")])
print(dr.fit$dr2.out[c("m0","m1","d","sed","dL","dU","pval.1sided")])

#ExactMatch.names<-c("female")
#allow.ties<-TRUE
#M<-1
#min.cell<-1
#match.fit<-get.matching.analysis(df=df0,ep.name=ep.name,treat.name=treat.name,miss.name=miss.name,match.names=ExactMatch.names,min.cell=min.cell,allow.ties=allow.ties,ps.names=ps.names,details=FALSE)


