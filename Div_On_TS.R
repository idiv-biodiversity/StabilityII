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



##################################
# Gross Figures 1 & 2 ############
##################################

#find studies without monocultures


study_div<-summarize(group_by(stab_444,Site),min_spp=min(SppN),max_spp=max(SppN))

study_div$keep<-ifelse(study_div$min_spp>1,NA,1)

study_div<-select(study_div,Site,keep)

stab_g<-merge(stab_444,study_div,by.y="Site")
stab_g<-filter(stab_g,keep==1)

###########################################
# Step 1: lm regression / loop through all#
###########################################

library(QuantPsyc)

n<-length(unique(stab_g$Site))

outt=c();
for(i in 1:n){
  
  test=subset(stab_g, stab_g$Site==(unique(stab_g$Site))[i])  
  Site<-as.character(unique(test$Site))
  Study_length<-unique(test$Study_length)

biom_lm<-lm(log(Plot_Biomassxbar)~log(SppN),data=test)

lm.beta(lm1)

lm.D9.beta <- lm.beta(lm.D9)
coef(lm.D9.beta)

sd_lm<-lm(log(Plot_Biomasssd)~log(SppN),data=test)

sd_beta<-summary(sd_lm)$coefficient[2,1]

study<-cbind.data.frame(Site,Study_length,biom_beta,sd_beta)

outt[[i]]<-rbind.data.frame(study)

}

jjj<-do.call(rbind,outt)  
jjj<-data.frame(jjj)
