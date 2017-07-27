rm(list=ls()) 

library(dplyr)

library(ggplot2)

library(lme4)
library(AICcmodavg)
require(MuMIn)

require(car)
#library(lars)

# Data
stab<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,
               Plot_Biomassxbar, Plot_Biomasssd,Plot_Asynchrony,meanPrecip,annualTemp,meanPET,CV_Temp,CV_Precip)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself

stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

stab_4$PlotAsynchrony_s<-stab_4$Plot_Asynchrony*-1

stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony


############################
#Predictor(s) of Stability #
############################


options(na.action = 'na.fail')
all_clim<-lme(TS_lg2~meanPrecip+annualTemp+CV_Temp+CV_Precip,random=~1+lg2SppN|Site,data=stab_444)

vif(all_clim)

dd<-dredge(all_clim)  # only 

importance(dd)

# choose CV_Precips (.06 importance)

#############################
# Predictors of Asynchrony ##
#############################

options(na.action = 'na.fail')
all_clim<-lme(Plot_Asynchrony~meanPrecip+annualTemp+CV_Temp+CV_Precip,random=~1+lg2SppN|Site,data=stab_444)

vif(all_clim)

dd<-dredge(all_clim)  # only 

importance(dd)  # all very low...


#############################
# Predictors of mean/sd of ##
# biomass                 ###  
#############################


options(na.action = 'na.fail')
all_clim<-lme(Plot_Biomassxbar~meanPET+meanPrecip+CV_Temp+CV_Precip,random=~1+lg2SppN|Site,data=stab_444)

vif(all_clim)

dd<-dredge(all_clim)  # only 

importance(dd)  # 
      # xbar biomass: annual temp, cv_temp, cv_precip, mean_precip
      # sd biomass: annual temp, cv_precip, cv_temp, mean precip



