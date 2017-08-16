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

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


######################
# piecewise SEM  #####
# with synchrony  ####
######################

# adjustments
#1)add meanPrecip, with path to asynchrony + stability
#2) add path from CV_Precip to asynchrony

##################
# FDis_eMNTD #####
##################


modList2=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(GrossAsynchrony_s~lg2SppN+FDis4+eMNTD+CV_Precip+meanPrecip,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+lg2SppN+eMNTD+CV_Precip++meanPrecip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


lapply(modList2, plot)

# Explore distribution of residuals

lapply(modList2, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList2[3:4], vif)

sem.fit(modList2,stab_444,corr.errors=c("eMNTD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

emntdfdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("eMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))  


emntdfdis.fit<-cbind(emntdfdis.fit$Fisher.C,emntdfdis.fit$AIC)
emntdfdis.fit$ModClass<-"FDis_eMNTD"

ts_emntd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("eMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_emntd2$ModClass<-"FDis_eMNTD"

mf_ts_emntd<-sem.model.fits(modList2)
mf_ts_emntd$ResponseVars<-c("eMNTD","FDis4","Asynchrony","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FDis4, meanPrecip, CV_Precip","eMNTD,Asynchrony,lg2SppN,CV_Precip")
mf_ts_emntd$ModClass<-"FDis_eMNTD"

write.table(ts_emntd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_sem_coefs_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_model_fits_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(emntdfdis.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_semfit_SIMPLE_July2017.csv",sep=",",row.names=F)

#######################
## FRic4 - eMNTD    ###
#######################


modList22=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(GrossAsynchrony_s~lg2SppN+FRic4+eMNTD+CV_Precip+meanPrecip,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+lg2SppN+eMNTD+CV_Precip+meanPrecip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)

lapply(modList22, plot)

# Explore distribution of residuals

lapply(modList22, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList22[3:4], vif)

sem.fit(modList22,stab_444,corr.errors=c("eMNTD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

emntdfric.fit<-sem.fit(modList22,stab_444,corr.errors=c("FRic4~~eMNTD","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                       model.control = list(lmeControl(opt = "optim")))  # add eMNTD as predictor of TS

emntdfric.fit<-cbind(emntdfric.fit$Fisher.C,emntdfric.fit$AIC)
emntdfric.fit$ModClass<-"FRic_eMNTD"

ts_emntd2<-sem.coefs(modList22,stab_444,standardize="scale",corr.errors=c("eMNTD~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_emntd2$ModClass<-"FRic_eMNTD"

mf_ts_emntd<-sem.model.fits(modList22)
mf_ts_emntd$ResponseVars<-c("eMNTD","FRic4","Asynchrony","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FRic4,meanPrecip,CV_Precip","Asynchrony,eMNTD,lg2SppN,meanPrecip, CV_Precip")
mf_ts_emntd$ModClass<-"FRic_eMNTD"

write.table(ts_emntd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_sem_coefs_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_model_fits_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(emntdfric.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_semfit_SIMPLE_July2017.csv",sep=",",row.names=F)

##################
# FDis_eMPD ######
##################

modList3=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(GrossAsynchrony_s~lg2SppN+FDis4+eMPD+meanPrecip+CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+lg2SppN+meanPrecip+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


lapply(modList3, plot)

# Explore distribution of residuals

lapply(modList3, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList3[3:4], vif)


sem.fit(modList3,stab_444,corr.errors=c("eMPD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

empdfdis.fit<-sem.fit(modList3,stab_444,corr.errors=c("eMPD~~FDis4","FDis4 ~~ PCAdim1_4trts","eMPD ~~ PCAdim1_4trts"),conditional=T,
                       model.control = list(lmeControl(opt = "optim")))  # 

empdfdis.fit<-cbind(empdfdis.fit$Fisher.C,empdfdis.fit$AIC)
empdfdis.fit$ModClass<-"FDis_eMPD"

ts_empd2<-sem.coefs(modList3,stab_444,standardize="scale",corr.errors=c("eMPD~~FDis4","FDis4 ~~ PCAdim1_4trts","eMPD ~~ PCAdim1_4trts"))
ts_empd2$ModClass<-"FDis_eMPD"


mf_ts_empd<-sem.model.fits(modList3)
mf_ts_empd$ResponseVars<-c("eMPD","FDis4","Asynchrony","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FDis4,meanPrecip,CV_Precip","Asynchrony,lg2SppN,CV_Precip,meanPrecip,CV_Precip")
mf_ts_empd$ModClass<-"FDis_eMPD"

write.table(ts_empd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_sem_coefs_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_empd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_model_fits_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(empdfdis.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_semfit_SIMPLE_July2017.csv",sep=",",row.names=F)

#######################
## FRic - eMPD #######
#######################

modList33=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(GrossAsynchrony_s~lg2SppN+FRic4+eMPD+meanPrecip+CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+lg2SppN+meanPrecip+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


lapply(modList33, plot)

# Explore distribution of residuals

lapply(modList33, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList33[3:4], vif)


sem.fit(modList33,stab_444,corr.errors=c("eMPD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

smpdfdis.fit<-sem.fit(modList33,stab_444,corr.errors=c("eMPD~~FRic4","FRic4 ~~ PCAdim1_4trts"," eMPD ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  


smpdfdis.fit<-cbind(smpdfdis.fit$Fisher.C,smpdfdis.fit$AIC)
smpdfdis.fit$ModClass<-"FRic4_eMPD"

ts_smpd2<-sem.coefs(modList33,stab_444,standardize="scale",corr.errors=c("eMPD~~FRic4","FRic4 ~~ PCAdim1_4trts"," eMPD ~~ PCAdim1_4trts"))
ts_smpd2$ModClass<-"FRic4"


mf_ts_smpd<-sem.model.fits(modList33)
mf_ts_smpd$ResponseVars<-c("eMPD","FRic4","Asynchrony","Temp_Stability")
mf_ts_smpd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FRic4,meanPrecip,CV_Precip","Asynchrony,F-S,lg2SppN,meanPrecip,CV_Precip")
mf_ts_smpd$ModClass<-"FRic4_eMPD"

write.table(ts_smpd2,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_sem_coefs_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(mf_ts_smpd,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_model_fits_SIMPLE_July2017.csv",sep=",",row.names=F)
write.table(smpdfdis.fit,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_semfit_SIMPLE_July2017.csv",sep=",",row.names=F)


####################
# Merge models #####
####################

[[not yet]] # 15.08.2017

fdis_paths<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_sem_coefs_Oct2016.csv",sep=",",header=T)
fdis_R2.1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_model_fits_Oct2016.csv",sep=",",header=T)
fdis_semfit<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_semfit_Oct2016.csv",sep=",",header=T)

fdis_paths2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_sem_coefs_Oct2016.csv",sep=",",header=T)
fdis_R2.2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_model_fits_Oct2016.csv",sep=",",header=T)
fdis_semfit2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_semfit_Oct2016.csv",sep=",",header=T)

fdis_paths3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fdis_sem_coefs_Oct2016.csv",sep=",",header=T)
fdis_R2.3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fdis_model_fits_Oct2016.csv",sep=",",header=T)
fdis_semfit3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fdis_semfit_Oct2016.csv",sep=",",header=T)

fdis_paths4<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fdis_sem_coefs_Oct2016.csv",sep=",",header=T)
fdis_R2.4<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fdis_model_fits_Oct2016.csv",sep=",",header=T)
fdis_semfit4<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fdis_semfit_Oct2016.csv",sep=",",header=T)

#####

fric_paths1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_sem_coefs_Oct2016.csv",sep=",",header=T)
fric_R2.1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_model_fits_Oct2016.csv",sep=",",header=T)
fric_semfit<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_semfit_Oct2016.csv",sep=",",header=T)

fric_paths2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_sem_coefs_Oct2016.csv",sep=",",header=T)
fric_R2.2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_model_fits_Oct2016.csv",sep=",",header=T)
fric_semfit2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_semfit_Oct2016.csv",sep=",",header=T)

fric_paths3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fric_sem_coefs_Oct2016.csv",sep=",",header=T)
fric_R2.3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fric_model_fits_Oct2016.csv",sep=",",header=T)
fric_semfit3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fric_semfit_Oct2016.csv",sep=",",header=T)

fric_paths4<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fric_sem_coefs_Oct2016.csv",sep=",",header=T)
fric_R2.4<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fric_model_fits_Oct2016.csv",sep=",",header=T)
fric_semfit4<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fric_semfit_Oct2016.csv",sep=",",header=T)

#

sem_paths_sync<-rbind.data.frame(fdis_paths,fdis_paths2,fdis_paths3,fdis_paths4,fric_paths1,fric_paths2,fric_paths3,fric_paths4)
sem_paths_sync$X<-NULL

write.table(sem_paths_sync,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_sync_paths_Oct2016.csv",sep=",",row.names=F)

sem_R2_sync<-rbind.data.frame(fdis_R2.1,fdis_R2.2,fdis_R2.3,fdis_R2.4,fric_R2.1,fric_R2.2,fric_R2.3,fric_R2.4)
write.table(sem_R2_sync,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_sync_R2_Oct2016.csv",sep=",",row.names=F)

sem_fits_sync<-rbind.data.frame(fdis_semfit,fdis_semfit2,fdis_semfit3,fdis_semfit4,fric_semfit,fric_semfit2,fric_semfit3,fric_semfit4)
write.table(sem_fits_sync,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_sync_semfits_Oct2016.csv",sep=",",row.names=F)

