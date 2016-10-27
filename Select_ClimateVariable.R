rm(list=ls()) 

library(dplyr)

library(ggplot2)

library(lme4)
library(AICcmodavg)
#library(lars)

# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony,meanPrecip,annualTemp,meanPET,CV_Temp,CV_Precip)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)


stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony

cc<-lmeControl(opt="optim")


stab_clim<-summarize(group_by(stab_444,Site),meanTS=mean(Plot_TempStab),meanPrecip=mean(meanPrecip),meanTemp=mean(annualTemp),meanPET=mean(meanPET),
                    CV_Temp=mean(CV_Temp),CV_Precip=mean(CV_Precip))

stab_clim$TS_lg2<-log(stab_clim$meanTS,base=2)

############################
## div + climate ###########
############################

Cand.set <- list( )
Cand.set[[1]]<-lmer(TS_lg2~meanPrecip+(1+lg2SppN|Site),data=stab_444)
Cand.set[[2]]<-lmer(TS_lg2~annualTemp+(1+lg2SppN|Site),data=stab_444)
Cand.set[[3]]<-lmer(TS_lg2~meanPET+(1+lg2SppN|Site),data=stab_444)
Cand.set[[4]]<-lmer(TS_lg2~CV_Temp+(1+lg2SppN|Site),data=stab_444)
Cand.set[[5]]<-lmer(TS_lg2~CV_Precip+(1+lg2SppN|Site),data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

# select mean PET :  delta AICc = 3.68 

#########################

Cand.set <- list( )
Cand.set[[1]]<-lm(TS_lg2~meanPrecip,data=stab_clim)
Cand.set[[2]]<-lm(TS_lg2~meanTemp,data=stab_clim)
Cand.set[[3]]<-lm(TS_lg2~meanPET,data=stab_clim)
Cand.set[[4]]<-lm(TS_lg2~CV_Temp,data=stab_clim)
Cand.set[[5]]<-lm(TS_lg2~CV_Precip,data=stab_clim)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

# select CV Temp :  delta AICc = 3.68 



##########################

require(MuMIn)

r.squaredGLMM(Cand.set[[5]])
