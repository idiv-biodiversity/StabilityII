rm(list=ls()) 

library(lmerTest)
library(nlme)
library(dplyr)
library(lm.beta)

# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,Study_length,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,LeafN, LeafP,Plot_Biomassxbar, Plot_Biomasssd,Plot_Asynchrony,Plot_TempStab)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony


#Drop studies without monocultures

study_div<-summarize(group_by(stab_444,Site),min_spp=min(SppN),max_spp=max(SppN))

study_div$keep<-ifelse(study_div$min_spp>1,NA,1)

study_div<-select(study_div,Site,keep)

stab_g<-merge(stab_444,study_div,by.y="Site")
stab_g<-filter(stab_g,keep==1)

###########################################
# Step 1: lm regression / loop through all#
###########################################


n<-length(unique(stab_g$Site))

outt=c();

for(i in 1:n){
  
  test=subset(stab_g, stab_g$Site==(unique(stab_g$Site))[i])  
  Site<-as.character(unique(test$Site))
  Study_length<-max(test$Study_length)

biom_lm<-lm(log(Plot_Biomassxbar)~log(SppN),data=test)

beta_b<-lm.beta(biom_lm)

biom_beta<-coef(beta_b)[2]

biom_slope<-summary(biom_lm)$coefficient[2,1]

biom_var<-vcov(biom_lm)[2,2]

sd_lm<-lm(log(Plot_Biomasssd)~log(SppN),data=test)

beta_sd<-lm.beta(sd_lm)

sd_beta<-coef(beta_sd)[2]

sd_slope<-summary(sd_lm)$coefficient[2,1]
sd_var<-vcov(sd_lm)[2,2]


study<-cbind.data.frame(Site,Study_length,biom_beta,sd_beta,biom_slope,biom_var,sd_slope,sd_var)

outt[[i]]<-rbind.data.frame(study)

}

jjj<-do.call(rbind,outt)  
jjj<-data.frame(jjj)

write.table(jjj,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Div_ts.csv",sep=",",row.names=F)


######################################
## Meta-analysis #####################
######################################

require(metaphor)

div_ts<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Div_ts.csv",sep=",",header=T)


biom_re<-rma.uni(yi=biom_slope,vi=biom_var,method="REML",knha=TRUE,data=div_ts)

sd_re<-rma.uni(yi=sd_slope,vi=sd_var,method="REML",knha=TRUE,data=div_ts)


