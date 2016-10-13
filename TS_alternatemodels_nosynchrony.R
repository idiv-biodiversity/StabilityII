rm(list=ls()) 
#require(devtools)
#install_github("jslefche/piecewiseSEM")

require(dplyr)
require(piecewiseSEM)
library(semPlot)
library(lmerTest)
library(nlme)

# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony,annualTemp,meanPrecip,meanPET,CV_Temp,CV_Precip)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)


stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


######################
# piecewise SEM  #####
# without synchrony  #
######################


modList=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  #lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("eMNTD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

fdis.fit1<-sem.fit(modList,stab_444,corr.errors=c("eMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) 

fdis.fit1<-cbind(fdis.fit1$Fisher.C,fdis.fit1$AIC)
fdis.fit1$ModClass<-"FDis4_eMNTD"


ts_emntd<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("eMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_emntd$ModClass<-"FDis4_eMNTD"


mf_ts_emntd<-sem.model.fits(modList)
mf_ts_emntd$ResponseVars<-c("eMNTD","FDis4","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","F-S,eMNTD,FDis4,lg2SppN,CV_Precip")
mf_ts_emntd$ModClass<-"FDis4_eMNTD"


write.table(ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_nosync_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_nosync__model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(fdis.fit1,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_nosync__semfit_Oct2016.csv",sep=",",row.names=F)

###############
## FDis -eMPD #
###############


modList=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  #lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FDis4+eMPD+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("eMPD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

fdis.fit2<-sem.fit(modList,stab_444,corr.errors=c("eMPD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
                   model.control = list(lmeControl(opt = "optim"))) 

fdis.fit2<-cbind(fdis.fit2$Fisher.C,fdis.fit2$AIC)
fdis.fit2$ModClass<-"FDis4_eMPD"


ts_empd<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("eMPD~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_empd$ModClass<-"FDis4_eMPD"

mf_ts_empd<-sem.model.fits(modList)
mf_ts_empd$ResponseVars<-c("eMPD","FDis4","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","F-S,eMPD,FDis4,lg2SppN,CV_Precip")
mf_ts_empd$ModClass<-"FDis4_eMPD"


write.table(ts_empd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_nosync_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_empd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_nosync__model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(fdis.fit2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_nosync__semfit_Oct2016.csv",sep=",",row.names=F)

##############
# FDis - ePSE#
##############

modList=list(
  lme(ePSE~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  #lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FDis4+ePSE+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("ePSE~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

fdis.fit3<-sem.fit(modList,stab_444,corr.errors=c("ePSE~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
                   model.control = list(lmeControl(opt = "optim"))) 

fdis.fit3<-cbind(fdis.fit3$Fisher.C,fdis.fit3$AIC)
fdis.fit3$ModClass<-"FDis4_ePSE"


ts_epse<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("ePSE~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_epse$ModClass<-"FDis4_ePSE"

mf_ts_epse<-sem.model.fits(modList)
mf_ts_epse$ResponseVars<-c("ePSE","FDis4","Temp_Stability")
mf_ts_epse$PredVars<-c("lg2SppN","lg2SppN","F-S,ePSE,FDis4,lg2SppN,CV_Precip")
mf_ts_epse$ModClass<-"FDis4_ePSE"


write.table(ts_epse,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_nosync_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_epse,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_nosync__model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(fdis.fit3,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_nosync__semfit_Oct2016.csv",sep=",",row.names=F)

#################
## FRic - eMNTD #
#################


modList=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  #lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FRic4+eMNTD+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("eMNTD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

fric.fit1<-sem.fit(modList,stab_444,corr.errors=c("eMNTD~~FRic4","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                   model.control = list(lmeControl(opt = "optim"))) 

fric.fit1<-cbind(fric.fit1$Fisher.C,fric.fit1$AIC)
fric.fit1$ModClass<-"FRic4_eMNTD"


ts_emntd<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("eMNTD~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_emntd$ModClass<-"FRic4_eMNTD"


mf_ts_emntd<-sem.model.fits(modList)
mf_ts_emntd$ResponseVars<-c("eMNTD","FRic4","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","F-S,eMNTD,FRic4,lg2SppN,CV_Precip")
mf_ts_emntd$ModClass<-"FRic4_eMNTD"


write.table(ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_nosync_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_nosync__model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(fric.fit1,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_nosync__semfit_Oct2016.csv",sep=",",row.names=F)

###############
## FDis -eMPD #
###############


modList=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  #lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FRic4+eMPD+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("eMPD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

fric.fit2<-sem.fit(modList,stab_444,corr.errors=c("eMPD~~FRic4","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                   model.control = list(lmeControl(opt = "optim"))) 

fric.fit2<-cbind(fric.fit2$Fisher.C,fric.fit2$AIC)
fric.fit2$ModClass<-"FRic4_eMPD"


ts_empd<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("eMPD~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_empd$ModClass<-"FRic4_eMPD"

mf_ts_empd<-sem.model.fits(modList)
mf_ts_empd$ResponseVars<-c("eMPD","FRic4","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","F-S,eMPD,FRic4,lg2SppN,CV_Precip")
mf_ts_empd$ModClass<-"FRic4_eMPD"


write.table(ts_empd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_nosync_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_empd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_nosync__model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(fric.fit2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_nosync__semfit_Oct2016.csv",sep=",",row.names=F)

##############
# FRic - ePSE#
##############

modList=list(
  lme(ePSE~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  #lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FRic4+ePSE+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("ePSE~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

fric.fit3<-sem.fit(modList,stab_444,corr.errors=c("ePSE~~FRic4","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                   model.control = list(lmeControl(opt = "optim"))) 

fric.fit3<-cbind(fric.fit3$Fisher.C,fric.fit3$AIC)
fric.fit3$ModClass<-"FRic4_ePSE"


ts_epse<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("ePSE~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_epse$ModClass<-"FRic4_ePSE"

mf_ts_epse<-sem.model.fits(modList)
mf_ts_epse$ResponseVars<-c("ePSE","FRic4","Temp_Stability")
mf_ts_epse$PredVars<-c("lg2SppN","lg2SppN","F-S,ePSE,FRic4,lg2SppN,CV_Precip")
mf_ts_epse$ModClass<-"FRic4_ePSE"


write.table(ts_epse,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_nosync_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_epse,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_nosync__model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(fric.fit3,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_nosync__semfit_Oct2016.csv",sep=",",row.names=F)


###############
# merge files #
###############


fdis_paths1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_nosync_sem_coefs_Oct2016.csv",sep=",",header=T)
fdis_R2.1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_nosync__model_fits_Oct2016.csv",sep=",",header=T)
fdis_semfit1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_nosync__semfit_Oct2016.csv",sep=",",header=T)


fdis_paths2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_nosync_sem_coefs_Oct2016.csv",sep=",",header=T)
fdis_R2.2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_nosync__model_fits_Oct2016.csv",sep=",",header=T)
fdis_semfit2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_nosync__semfit_Oct2016.csv",sep=",",header=T)

fdis_paths3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_nosync_sem_coefs_Oct2016.csv",sep=",",header=T)
fdis_R2.3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_nosync__model_fits_Oct2016.csv",sep=",",header=T)
fdis_semfit3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_nosync__semfit_Oct2016.csv",sep=",",header=T)

fric_paths1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_nosync_sem_coefs_Oct2016.csv",sep=",",header=T)
fric_R2.1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_nosync__model_fits_Oct2016.csv",sep=",",header=T)
fric_semfit1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_nosync__semfit_Oct2016.csv",sep=",",header=T)

fric_paths2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_nosync_sem_coefs_Oct2016.csv",sep=",",header=T)
fric_R2.2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_nosync__model_fits_Oct2016.csv",sep=",",header=T)
fric_semfit2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_nosync__semfit_Oct2016.csv",sep=",",header=T)


fric_paths3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_nosync_sem_coefs_Oct2016.csv",sep=",",header=T)
fric_R2.3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_nosync__model_fits_Oct2016.csv",sep=",",header=T)
fric_semfit3<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_nosync__semfit_Oct2016.csv",sep=",",header=T)


sem_paths_no_sync<-rbind.data.frame(fdis_paths1,fdis_paths2,fdis_paths3,fric_paths1,fric_paths2,fric_paths3)
sem_paths_no_sync$X<-NULL

write.table(sem_paths_no_sync,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_no_sync_paths_Oct2016.csv",sep=",",row.names=F)

sem_R2_no_sync<-rbind.data.frame(fdis_R2.1,fdis_R2.2,fdis_R2.3,fric_R2.1,fric_R2.2,fric_R2.3)
write.table(sem_R2_no_sync,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_no_sync_R2_Oct2016.csv",sep=",",row.names=F)

sem_fits_no_sync<-rbind.data.frame(fdis_semfit1,fdis_semfit2,fdis_semfit3,fric_semfit1,fric_semfit2,fric_semfit3)
write.table(sem_fits_no_sync,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_no_sync_semfits_Oct2016.csv",sep=",",row.names=F)
