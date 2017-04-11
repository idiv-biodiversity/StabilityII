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

stab_4$PlotAsynchrony_s<-stab_4$Plot_Asynchrony*-1

stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony



cc<-lmeControl(opt="optim")


stab_clim<-summarize(group_by(stab_444,Site),meanTS=mean(Plot_TempStab),meanPrecip=mean(meanPrecip),meanTemp=mean(annualTemp),meanPET=mean(meanPET),
                    CV_Temp=mean(CV_Temp),CV_Precip=mean(CV_Precip))

stab_clim$TS_lg2<-log(stab_clim$meanTS,base=2)


stab_clim$meanPrecips<-scale(stab_clim$meanPrecip,scale=T,center=T)
stab_clim$meanTemps<-scale(stab_clim$meanTemp,scale=T,center=T)
stab_clim$meanPETs<-scale(stab_clim$meanPET,scale=T,center=T)
stab_clim$CV_Temps<-scale(stab_clim$CV_Temp,scale=T,center=T)
stab_clim$CV_Precips<-scale(stab_clim$CV_Precip,scale=T,center=T)


############################
## div + climate ###########
############################


#########################

Cand.set <- list( )
Cand.set[[1]]<-lm(TS_lg2~meanPrecips,data=stab_clim)
Cand.set[[2]]<-lm(TS_lg2~meanTemps,data=stab_clim)
Cand.set[[3]]<-lm(TS_lg2~meanPETs,data=stab_clim)
Cand.set[[4]]<-lm(TS_lg2~CV_Temps,data=stab_clim)
Cand.set[[5]]<-lm(TS_lg2~CV_Precips,data=stab_clim)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

# select CV Temp :  delta AICc = 1.13

######################

require(MuMIn)

options(na.action = 'na.fail')
all_clim<-lm(TS_lg2~meanPrecips+meanTemps+meanPETs+CV_Temps+CV_Precips,data=stab_clim)

vif(all_clim)

dd<-dredge(all_clim)  # only 

importance(dd)

# choose CV_Precips based on relative importance (ie sum of Akaike mdoels)


