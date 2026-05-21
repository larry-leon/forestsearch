# Changing the path
# Will need to do this for every session
.libPaths("C:/Users/leolarr2/OneDrive - Merck Sharp & Dohme, Corp/documents/MyRlibrary")

df_long <- read.table("outdata/eff_cough.csv",header=TRUE, sep=",")

names(df_long) <- tolower(names(df_long))

#View(df_long)

head(df_long)

df_long <- within(df_long,{
analysis_day<-ady
treat_45mg<-ifelse(trtp=="Placebo",0,1)
flag_baseline<-ifelse(avisitn==0,1,0)
y_0<-base
chg<-aval-base
}
)

# Sort by subjid and analysis_day
# Note: can also use data.table
df_long<-df_long[order(df_long$subjid,df_long$analysis_day),]

df_baseline <- subset(df_long,flag_baseline==1)
table(df_baseline$treat_45mg)

df_baseline$base<-df_baseline
# Only keep subjid and base2=aval

df_bl<-df_baseline[,c("subjid","y_0","age","sex","trtp","treat_45mg","offtrt_day")]

# Week 4
# Note: chg_ontrt_4 --> CBL at week 4 *while on treatment*
# If missing then CBL "off treatment"
df_week4 <- subset(df_long,avisitn==4)
# 12 subjects missing "week 4"
df_week4<-within(df_week4,{
y_4<-aval  
chg_ontrt_4<-chg_ontrt  
})

df_week4<-df_week4[,c("subjid","y_4","chg_ontrt_4")]

temp<-merge(df_bl,df_week4,by="subjid",all.x=TRUE)

# Any off treatment befor Week 4?
check<-subset(temp,offtrt_day<=7*4)
cat("# of subects off treatment by day 28=",c(nrow(check)),"\n")
#subj_look<-subset(df_long,subjid=="402459")

#Baseline and week 4
df_0_4<-temp

# Week 8
df_week8 <- subset(df_long,avisitn==8)
df_week8<-within(df_week8,{
  y_8<-aval  
  chg_ontrt_8<-chg_ontrt  
})
df_week8<-df_week8[,c("subjid","y_8","chg_ontrt_8")]

temp<-merge(df_bl,df_week8,by="subjid",all.x=TRUE)

# Any off treatment befor Week 4?
check<-subset(temp,offtrt_day<=7*8)
cat("# of subects off treatment by day 56=",c(nrow(check)),"\n")

# Merge week8 with already merged baseline and week4

df_0_4_8<-merge(df_0_4,df_week8,by="subjid",all.x=TRUE)

# Missing y_4 AND y_8
miss_both<-subset(df_0_4_8,is.na(y_4) & is.na(y_8))
# 4 subjects missing both y_4 AND y_8

# Impute missing y_4 by y_0
# Impute missing y_8 by y_4(impute)

df_0_4_8<-within(df_0_4_8,{
miss_4_8<-ifelse(is.na(y_4) & is.na(y_8),1,0)
y_4_impute<-ifelse(!is.na(y_4),y_4,y_0)
y_8_impute<-ifelse(!is.na(y_8),y_8,y_4_impute)
chg_4_impute<-y_4_impute-y_0
chg_8_impute<-y_8_impute-y_0
})

# Week 24
df_week24 <- subset(df_long,avisitn==24)
df_week24<-within(df_week24,{
  y_24<-aval  
  chg_ontrt_24<-chg_ontrt  
  ontrtfl_24<-ontrtfl
})
df_week24<-df_week24[,c("subjid","y_24","chg_ontrt_24","ontrtfl_24")]

# Merge with df_0_4_8

df_analysis24<-merge(df_0_4_8,df_week24,by="subjid",all.x=TRUE)

# Note: offtrt_day refers to last treatment day
# create flags for last treatment prior to week 4, week 8, 12, 16, 20, 24

# Define missingness for "general reason"
# and "hypothetical" setting to missing if ontrtfl=N at week 24
# Define Hy_24 as "hypothetical"
df_analysis24<-within(df_analysis24,{
Hy_24<-ifelse(is.na(chg_ontrt_24) & !is.na(y_24),NA,y_24)
Hchg_24<-Hy_24-y_0
chg_24<-y_24-y_0
}
)

df_analysis24<-within(df_analysis24,{
last_treat8<-ifelse(offtrt_day<=8*7,1,0)
last_treat12<-ifelse(offtrt_day<=12*7,1,0)
last_treat16<-ifelse(offtrt_day<=16*7,1,0)
last_treat20<-ifelse(offtrt_day<=20*7,1,0)
miss_offtreat<-ifelse(is.na(Hy_24) & !is.na(y_24),1,0)
miss_24<-ifelse(is.na(y_24),1,0)
miss_H24<-ifelse(is.na(Hy_24),1,0)
female<-ifelse(sex=="F",1,0)
})

# Set to missing for off treatment
check<-subset(df_analysis24,is.na(Hy_24) & !is.na(y_24))
# 5 subjects



save(df_analysis24,df_long,file="outdata_analysis/df_week24.Rdata")


