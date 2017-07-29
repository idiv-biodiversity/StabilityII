rm(list=ls()) 

library(dplyr)

library(ggplot2)

library(lme4)
library(AICcmodavg)
require(MuMIn)

require(car)

# Data

stab<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_VI.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,Study_length,UniqueID,SppN,eMPD,eMNTD,ePSE,sMPD,sMNTD,FDis4,FRic4,PCAdim1_4trts,SLA, LDMC, LeafN, LeafP,
               Plot_TempStab,Plot_Biomassxbar, Plot_Biomasssd,Plot_Asynchrony, Gross_synchrony, Loreau_synchrony,annualTemp,meanPrecip,meanPET,CV_Temp,CV_Precip)


# for plots with ONLY 1 spp, we assume that a species #is perfectly synchronized with itself

stab_4$Plot_Asynchrony<-ifelse(is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) 
stab_4$Gross_synchrony<-ifelse(is.na(stab_4$Gross_synchrony)==TRUE,1,stab_4$Gross_synchrony) 
stab_4$Loreau_synchrony<-ifelse(is.na(stab_4$Loreau_synchrony)==TRUE,1,stab_4$Loreau_synchrony) 

# convert synchrony metrics to different scale

stab_4$PlotAsynchrony_s<-stab_4$Plot_Asynchrony*-1

stab_4$GrossAsynchrony_s<-stab_4$Gross_synchrony*-1


# further adjustments

stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)


# Filter out NAs for Asynchrony and  FRic4
stab_444<-filter(stab_4, is.na(PlotAsynchrony_s)==FALSE)
stab_444<-filter(stab_444, is.na(FRic4)==FALSE)

### summarize


stab_s<-summarize(group_by(stab_444, Site),MAT=mean(annualTemp),MAP=mean(meanPrecip),PET=mean(meanPET),CV_Temp=mean(CV_Temp),CV_Precip=mean(CV_Precip),
                       TS=mean(Plot_TempStab), Async=mean(GrossAsynchrony_s), meanBiomass=mean(Plot_Biomassxbar), sdd=sd(Plot_Biomasssd))

stab_s$TS_lg2<-log(stab_s$TS,base=2)

############################
#Predictor(s) of Stability #
############################

options(na.action = 'na.fail')

all_clim<-lm(TS_lg2~MAT+MAP+PET+CV_Temp+CV_Precip,data=stab_s)

vif(all_clim)

dd<-dredge(all_clim)  # only 

importance(dd)

# choose CV_Precips (.42 importance)

#############################
# Predictors of Asynchrony ##
#############################

options(na.action = 'na.fail')

all_clim<-lm(Async~MAT+MAP+PET+CV_Temp+CV_Precip,data=stab_s)

vif(all_clim)

dd<-dredge(all_clim)  # only 

importance(dd)  # 

#choose CV_Precips (.89 importance)

#############################
# Predictors of mean/sd of ##
# biomass                 ###  
#############################


# biomass

options(na.action = 'na.fail')

all_clim<-lm(meanBiomass~MAT+MAP+PET+CV_Temp+CV_Precip,data=stab_s)

vif(all_clim)

dd<-dredge(all_clim)  # only 

importance(dd)  # 

#choose CV_Precips (.36 importance), temp also high

#sd

options(na.action = 'na.fail')

all_clim<-lm(sdd~MAT+MAP+PET+CV_Temp+CV_Precip,data=stab_s)

vif(all_clim)

dd<-dredge(all_clim)  # only 

importance(dd)  # 

#choose CV_Precips (.96 importance)