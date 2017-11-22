###############################################
# sensitivity analysis of long-term data sets #
# using basic model in 'TS_alternatemodels.R' #
###############################################
#update 21.11.2017 add direct paths of FD+PD to stability
#
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


####################
# Data subsets  ####
####################

stab_666<-filter(stab_444, Study_length>4)

######################
# piecewise SEM  #####
# with synchrony  ####
######################

# adjustments
#1)add meanPrecip, with path to asynchrony + stability
#2) add path from CV_Precip to asynchrony

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")



##################
# FDis_eMNTD #####
##################

modList2=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(GrossAsynchrony_s~lg2SppN+FDis4+eMNTD+CV_Precip+meanPrecip,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+eMNTD+FDis4+lg2SppN+CV_Precip+meanPrecip,random=~1+lg2SppN|Site, control=cc,data=stab_666)
)


lapply(modList2, plot)

# Explore distribution of residuals

lapply(modList2, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList2[3:4], vif)

sem.fit(modList2,stab_666,corr.errors=c("eMNTD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

emntdfdis.fit<-sem.fit(modList2,stab_666,
                       corr.errors=c("eMNTD~~FDis4","GrossAsynchrony_s~~PCAdim1_4trts",
                                     "eMNTD~~PCAdim1_4trts","eMNTD~~meanPrecip"),conditional=T,
                       model.control = list(lmeControl(opt = "optim")))  


emntdfdis.fit<-cbind(emntdfdis.fit$Fisher.C,emntdfdis.fit$AIC)
emntdfdis.fit$ModClass<-"FDis_eMNTD_LongTerm"

ts_emntd2<-sem.coefs(modList2,stab_666,standardize="scale",corr.errors=c("eMNTD ~~ FDis4","GrossAsynchrony_s ~~ PCAdim1_4trts",
                                                                         "eMNTD ~~ PCAdim1_4trts","eMNTD ~~ meanPrecip"))
ts_emntd2$ModClass<-"FDis_eMNTD_LongTerm"

mf_ts_emntd<-sem.model.fits(modList2)
mf_ts_emntd$ResponseVars<-c("eMNTD","FDis4","Asynchrony","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FDis4, meanPrecip, CV_Precip","eMNTD,FDis4,eMNTD,Asynchrony,lg2SppN,meanPrecip,CV_Precip")
mf_ts_emntd$ModClass<-"FDis_eMNTD_LongTerm"

write.table(ts_emntd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_sem_coefs_LongTerm_SIMPLE_Nov2017.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_model_fits_LongTerm_SIMPLE_Nov2017.csv",sep=",",row.names=F)
write.table(emntdfdis.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_semfit_LongTerm_SIMPLE_Nov2017.csv",sep=",",row.names=F)

#######################
## FRic4 - eMNTD    ###
#######################
#start here

modList22=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(GrossAsynchrony_s~lg2SppN+FRic4+eMNTD+CV_Precip+meanPrecip,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+eMNTD+FRic4+lg2SppN+CV_Precip+meanPrecip,random=~1+lg2SppN|Site, control=cc,data=stab_666)
)

lapply(modList22, plot)

# Explore distribution of residuals

lapply(modList22, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList22[3:4], vif)

sem.fit(modList22,stab_666,corr.errors=c("eMNTD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

emntdfric.fit<-sem.fit(modList22,stab_666,corr.errors=c("FRic4~~eMNTD","eMNTD ~~ meanPrecip", " eMNTD ~~ PCAdim1_4trts",
                                                        "GrossAsynchrony_s ~~ PCAdim1_4trts"),conditional=T,
                       model.control = list(lmeControl(opt = "optim")))  # add eMNTD as predictor of TS

emntdfric.fit<-cbind(emntdfric.fit$Fisher.C,emntdfric.fit$AIC)
emntdfric.fit$ModClass<-"FRic_eMNTD_LongTerm"

ts_emntd2<-sem.coefs(modList22,stab_666,standardize="scale",corr.errors=c("FRic4~~eMNTD","eMNTD ~~ meanPrecip", " eMNTD ~~ PCAdim1_4trts",
                                                                          "GrossAsynchrony_s ~~ PCAdim1_4trts"))
ts_emntd2$ModClass<-"FRic_eMNTD_LongTerm"

mf_ts_emntd<-sem.model.fits(modList22)
mf_ts_emntd$ResponseVars<-c("eMNTD","FRic4","Asynchrony","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FRic4,meanPrecip,CV_Precip","Asynchrony,FRic4, eMNTD,lg2SppN,meanPrecip, CV_Precip")
mf_ts_emntd$ModClass<-"FRic_eMNTD_LongTerm"

write.table(ts_emntd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_sem_coefs_LongTerm_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_model_fits_LongTerm_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(emntdfric.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_semfit_LongTerm_SIMPLE_July2017.csv",sep=",",row.names=F)

##################
# FDis_eMPD ######
##################

modList3=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(GrossAsynchrony_s~lg2SppN+FDis4+eMPD+meanPrecip+CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(TS_lg2~GrossAsynchrony_s+FDis4+PCAdim1_4trts+lg2SppN+meanPrecip+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_666)
)


lapply(modList3, plot)

# Explore distribution of residuals

lapply(modList3, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList3[3:4], vif)


sem.fit(modList3,stab_666,corr.errors=c("eMPD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

empdfdis.fit<-sem.fit(modList3,stab_666,corr.errors=c("eMPD~~FDis4","eMPD ~~ PCAdim1_4trts","GrossAsynchrony_s ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  # 

empdfdis.fit<-cbind(empdfdis.fit$Fisher.C,empdfdis.fit$AIC)
empdfdis.fit$ModClass<-"FDis_eMPD_LongTerm"

ts_empd2<-sem.coefs(modList3,stab_666,standardize="scale",corr.errors=c("eMPD ~~ FDis4","eMPD ~~ PCAdim1_4trts","GrossAsynchrony_s ~~ PCAdim1_4trts"))
ts_empd2$ModClass<-"FDis_eMPD_LongTerm"


mf_ts_empd<-sem.model.fits(modList3)
mf_ts_empd$ResponseVars<-c("eMPD","FDis4","Asynchrony","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FDis4,meanPrecip,CV_Precip","Asynchrony,lg2SppN,FDis4,CV_Precip,meanPrecip,CV_Precip")
mf_ts_empd$ModClass<-"FDis_eMPD_LongTerm"

write.table(ts_empd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_sem_coefs_LongTerm_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_empd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_model_fits_LongTermSIMPLE_July2017.csv",sep=",",row.names=F)
write.table(empdfdis.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_semfit_LongTerm_SIMPLE_July2017.csv",sep=",",row.names=F)

#######################
## FRic - eMPD #######
#######################


modList33=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(GrossAsynchrony_s~lg2SppN+FRic4+eMPD+meanPrecip+CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_666),
  lme(TS_lg2~GrossAsynchrony_s+FRic4+PCAdim1_4trts+lg2SppN+meanPrecip+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_666)
)


lapply(modList33, plot)

# Explore distribution of residuals

lapply(modList33, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList33[3:4], vif)

sem.fit(modList33,stab_666,corr.errors=c("eMPD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

smpdfdis.fit<-sem.fit(modList33,stab_666,corr.errors=c("eMPD~~FRic4","eMPD ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  


smpdfdis.fit<-cbind(smpdfdis.fit$Fisher.C,smpdfdis.fit$AIC)
smpdfdis.fit$ModClass<-"FRic4_eMPD_LongTerm"

ts_smpd2<-sem.coefs(modList33,stab_666,standardize="scale",corr.errors=c("eMPD~~FRic4","eMPD ~~ PCAdim1_4trts"))
ts_smpd2$ModClass<-"FRic4_eMPD_LongTerm"


mf_ts_smpd<-sem.model.fits(modList33)
mf_ts_smpd$ResponseVars<-c("eMPD","FRic4","Asynchrony","Temp_Stability")
mf_ts_smpd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FRic4,meanPrecip,CV_Precip","Asynchrony,FRic4,F-S,lg2SppN,meanPrecip,CV_Precip")
mf_ts_smpd$ModClass<-"FRic4_eMPD_LongTerm"

write.table(ts_smpd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_sem_coefs_LongTerm_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_smpd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_model_fits_LongTerm_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(smpdfdis.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_semfit_LongTermSIMPLE_July2017.csv",sep=",",row.names=F)
