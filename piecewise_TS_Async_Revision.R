rm(list=ls()) 
#require(devtools)
#install_github("jslefche/piecewiseSEM")

require(dplyr)
require(piecewiseSEM)
library(semPlot)
library(lmerTest)
library(nlme)
library(car) 

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

# Filter out NAs for all functional traits

stab_555<-filter(stab_444, is.na(LDMC)==FALSE)
stab_555<-filter(stab_555, is.na(LeafN)==FALSE)
stab_555<-filter(stab_555, is.na(LeafP)==FALSE)


####################
# Prelim       #####
####################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")

####################################
# EXTENDED SEM including mean + SD #
####################################


######################
# piecewise SEM  #####
# with Asynchrony  #
######################

##################
# FDis_eMNTD #####
##################

modList2=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(PlotAsynchrony_s~lg2SppN+FDis4+eMNTD+meanPrecip+CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Biomassxbar~PlotAsynchrony_s+PCAdim1_4trts+lg2SppN+eMNTD+FDis4+meanPrecip, random=~1+lg2SppN|Site,control=cc, data=stab_444),
  
  lme(Plot_Biomasssd~PlotAsynchrony_s+PCAdim1_4trts+lg2SppN+CV_Precip, random=~1+lg2SppN|Site,control=cc, data=stab_444),
  
  lme(TS_lg2~Plot_Biomassxbar+Plot_Biomasssd,random=~1|Site, control=cc,data=stab_444)
)

lapply(modList2, plot)

# Explore distribution of residuals

lapply(modList2, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList2[3:6], vif)


sem.fit(modList2,stab_444,corr.errors=c("eMNTD~~FDis4","Plot_Biomasssd~~Plot_Biomassxbar","FDis4 ~~ PCAdim1_4trts",
                                        "TS_lg2~~FDis4","TS_lg2~~eMNTD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #initial model


emntdfdis.fit<-sem.fit(modList2,stab_444,
                       corr.errors=c("eMNTD~~FDis4","Plot_Biomasssd~~eMNTD","Plot_Biomasssd~~Plot_Biomassxbar","FDis4 ~~ PCAdim1_4trts",
                                     "TS_lg2~~FDis4","TS_lg2~~eMNTD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"),
                       conditional=T,
                       model.control = list(lmeControl(opt = "optim"))) 

# add correlated error for eMNTD~~ sd of biomass

emntdfdis.fit<-cbind(emntdfdis.fit$Fisher.C,emntdfdis.fit$AIC)
emntdfdis.fit$ModClass<-"FDis_eMNTD"

ts_emntd2<-sem.coefs(modList2,stab_444,standardize="scale",
                     corr.errors=c("eMNTD~~FDis4","Plot_Biomasssd~~eMNTD","Plot_Biomasssd~~Plot_Biomassxbar","FDis4 ~~ PCAdim1_4trts",
                                   "TS_lg2~~FDis4","TS_lg2~~eMNTD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"))
ts_emntd2$ModClass<-"FDis_eMNTD"

sem.plot(modList2, stab_444, standardize = "scale")


mf_ts_emntd<-sem.model.fits(modList2)
mf_ts_emntd$ResponseVars<-c("eMNTD","FDis4","Asynchrony","xbar_Biomass", "sd_Biomass","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FDis4, meanPrecip, cvPrecip","Asynchrony,lg2SppN,F-S,FD,PD, meanPrecip",
                        "Asynchrony,lg2SppN,F-S, cvPrecip","meanBiomass,sdBiomass")
mf_ts_emntd$ModClass<-"FDis_eMNTD"

write.table(ts_emntd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_sem_coefs_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_model_fits_July2017.csv",sep=",",row.names=F)
write.table(emntdfdis.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_semfit_July2017.csv",sep=",",row.names=F)

#######################
## FRic4 - eMNTD    ###
#######################


modList22=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(PlotAsynchrony_s~lg2SppN+FRic4+eMNTD+meanPrecip+CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Biomassxbar~PlotAsynchrony_s+PCAdim1_4trts+lg2SppN+eMNTD+FRic4+meanPrecip, random=~1+lg2SppN|Site,control=cc, data=stab_444),
  
  lme(Plot_Biomasssd~PlotAsynchrony_s+PCAdim1_4trts+lg2SppN+CV_Precip, random=~1+lg2SppN|Site,control=cc, data=stab_444),
  
  lme(TS_lg2~Plot_Biomassxbar+Plot_Biomasssd,random=~1|Site, control=cc,data=stab_444)
  
)

lapply(modList22, plot)

# Explore distribution of residuals

lapply(modList22, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList22[3:6], vif)


sem.fit(modList22,stab_444,corr.errors=c("eMNTD~~FRic4","Plot_Biomasssd~~Plot_Biomassxbar","FRic4 ~~ PCAdim1_4trts",
                                        "TS_lg2~~FRic4","TS_lg2~~eMNTD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #initial model


emntdfric.fit<-sem.fit(modList22,stab_444,
                       corr.errors=c("eMNTD~~FRic4","Plot_Biomasssd~~eMNTD","Plot_Biomasssd~~FRic4","Plot_Biomasssd~~Plot_Biomassxbar","FRic4 ~~ PCAdim1_4trts",
                                     "TS_lg2~~FRic4","TS_lg2~~eMNTD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"),
                       conditional=T,
                       model.control = list(lmeControl(opt = "optim"))) 


#add correlated errors for emntd and FRicto sd biomass

emntdfric.fit<-cbind(emntdfric.fit$Fisher.C,emntdfric.fit$AIC)
emntdfric.fit$ModClass<-"FRic_eMNTD"

ts_emntd2<-sem.coefs(modList22,stab_444,standardize="scale",
                     corr.errors=c("eMNTD~~FRic4","Plot_Biomasssd~~eMNTD","Plot_Biomasssd~~FRic4","Plot_Biomasssd~~Plot_Biomassxbar","FRic4 ~~ PCAdim1_4trts",
                                   "TS_lg2~~FRic4","TS_lg2~~eMNTD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"))
ts_emntd2$ModClass<-"FRic_eMNTD"

mf_ts_emntd<-sem.model.fits(modList22)
mf_ts_emntd$ResponseVars<-c("eMNTD","FRic4","Asynchrony","xbar_Biomass", "sd_Biomass","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FRic4, meanPrecip, cvPrecip","Asynchrony,lg2SppN,F-S,FRic4,PD, meanPrecip",
                        "Asynchrony,lg2SppN,F-S, cvPrecip","meanBiomass,sdBiomass")
mf_ts_emntd$ModClass<-"FRic_eMNTD"

write.table(ts_emntd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_sem_coefs_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_model_fits_July2017.csv",sep=",",row.names=F)
write.table(emntdfric.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_semfit_July2017.csv",sep=",",row.names=F)

##################
# FDis_MPD ######
##################

modList3=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(PlotAsynchrony_s~lg2SppN+FDis4+eMPD+meanPrecip+CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Biomassxbar~PlotAsynchrony_s+PCAdim1_4trts+lg2SppN+eMPD+FDis4+meanPrecip, random=~1+lg2SppN|Site,control=cc, data=stab_444),
  
  lme(Plot_Biomasssd~PlotAsynchrony_s+PCAdim1_4trts+lg2SppN+CV_Precip, random=~1+lg2SppN|Site,control=cc, data=stab_444),
  
  lme(TS_lg2~Plot_Biomassxbar+Plot_Biomasssd,random=~1|Site, control=cc,data=stab_444)
)

lapply(modList3, plot)

# Explore distribution of residuals

lapply(modList3, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList3[3:6], vif)


sem.fit(modList3,stab_444,corr.errors=c("eMPD~~FDis4","Plot_Biomasssd~~Plot_Biomassxbar","FDis4 ~~ PCAdim1_4trts",
                                        "TS_lg2~~FDis4","TS_lg2~~eMPD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #initial model


empdfdis.fit<-sem.fit(modList3,stab_444,
                       corr.errors=c("eMPD~~FDis4","Plot_Biomasssd~~Plot_Biomassxbar","FDis4 ~~ PCAdim1_4trts",
                                     "TS_lg2~~FDis4","TS_lg2~~eMPD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"),
                       conditional=T,
                       model.control = list(lmeControl(opt = "optim"))) 

# 
empdfdis.fit<-cbind(empdfdis.fit$Fisher.C,empdfdis.fit$AIC)
empdfdis.fit$ModClass<-"FDis_eMPD"

ts_empd2<-sem.coefs(modList3,stab_444,standardize="scale",
                     corr.errors=c("eMPD~~FDis4","Plot_Biomasssd~~Plot_Biomassxbar","FDis4 ~~ PCAdim1_4trts",
                                   "TS_lg2~~FDis4","TS_lg2~~eMPD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"))
ts_empd2$ModClass<-"FDis_eMPD"

sem.plot(modList3, stab_444, standardize = "scale")


mf_ts_empd<-sem.model.fits(modList3)
mf_ts_empd$ResponseVars<-c("eMPD","FDis4","Asynchrony","xbar_Biomass", "sd_Biomass","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FDis4, meanPrecip, cvPrecip","Asynchrony,lg2SppN,F-S,FD,PD, meanPrecip",
                        "Asynchrony,lg2SppN,F-S, cvPrecip","meanBiomass,sdBiomass")
mf_ts_empd$ModClass<-"FDis_eMPD"

write.table(ts_empd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_sem_coefs_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_empd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_modelfits_July2017.csv",sep=",",row.names=F)
write.table(empdfdis.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_semfits_July2017.csv",sep=",",row.names=F)


##################
# FRic_MPD #######
##################


modList44=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(PlotAsynchrony_s~lg2SppN+FRic4+eMPD+meanPrecip+CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Biomassxbar~PlotAsynchrony_s+PCAdim1_4trts+lg2SppN+eMPD+FRic4+meanPrecip, random=~1+lg2SppN|Site,control=cc, data=stab_444),
  
  lme(Plot_Biomasssd~PlotAsynchrony_s+PCAdim1_4trts+lg2SppN+CV_Precip, random=~1+lg2SppN|Site,control=cc, data=stab_444),
  
  lme(TS_lg2~Plot_Biomassxbar+Plot_Biomasssd,random=~1|Site, control=cc,data=stab_444)
  
)

lapply(modList44, plot)

# Explore distribution of residuals

lapply(modList44, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList44[3:6], vif)


sem.fit(modList44,stab_444,corr.errors=c("eMPD~~FRic4","Plot_Biomasssd~~Plot_Biomassxbar","FRic4 ~~ PCAdim1_4trts",
                                         "TS_lg2~~FRic4","TS_lg2~~eMPD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #initial model


empdfric.fit<-sem.fit(modList44,stab_444,
                       corr.errors=c("eMPD~~FRic4","Plot_Biomasssd~~FRic4","Plot_Biomasssd~~Plot_Biomassxbar","FRic4 ~~ PCAdim1_4trts",
                                     "TS_lg2~~FRic4","TS_lg2~~eMPD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"),
                       conditional=T,
                       model.control = list(lmeControl(opt = "optim"))) 


#add correlated errors for emntd and FRicto sd biomass

empdfric.fit<-cbind(empdfric.fit$Fisher.C,empdfric.fit$AIC)
empdfric.fit$ModClass<-"FRic_eMPD"

ts_empd2<-sem.coefs(modList44,stab_444,standardize="scale",
                     corr.errors=c("eMPD~~FRic4","Plot_Biomasssd~~FRic4","Plot_Biomasssd~~Plot_Biomassxbar","FRic4 ~~ PCAdim1_4trts",
                                   "TS_lg2~~FRic4","TS_lg2~~eMPD","TS_lg2~~PlotAsynchrony_s","TS_lg2 ~~ lg2SppN"))
ts_empd2$ModClass<-"FRic_eMPD"

mf_ts_empd2<-sem.model.fits(modList44)
mf_ts_empd2$ResponseVars<-c("eMPD","FRic4","Asynchrony","xbar_Biomass", "sd_Biomass","Temp_Stability")
mf_ts_empd2$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FRic4, meanPrecip, cvPrecip","Asynchrony,lg2SppN,F-S,FRic4,PD, meanPrecip",
                        "Asynchrony,lg2SppN,F-S, cvPrecip","meanBiomass,sdBiomass")
mf_ts_empd2$ModClass<-"FRic_eMPD"

write.table(ts_empd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_sem_coefs_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_empd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_model_fits_July2017.csv",sep=",",row.names=F)
write.table(empdfric.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_semfit_July2017.csv",sep=",",row.names=F)

