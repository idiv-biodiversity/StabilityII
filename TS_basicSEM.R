###################################
# Basic piecewise SEM             #
# for each combination of         #
# phylogenetic and functional     #
# diversity metric:               #
# 1) FDis + eMNTD                 #
# 2) FRic + eMNTD                 #
# 3) FDis + eMPD                  #
# 4) FRic + eMPD                  #
###################################

require(dplyr)
require(piecewiseSEM)
library(semPlot)
library(lmerTest)
library(nlme)
library(car) 

# Prepare data for analysis
stab<-read.delim("data.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,Study_length,UniqueID,SppN,eMPD,eMNTD,ePSE,sMPD,sMNTD,FDis4,FRic4,PCAdim1_4trts,SLA, LDMC, LeafN, LeafP,
               Plot_TempStab,Plot_Biomassxbar, Plot_Biomasssd,Gross_synchrony,
               mArid,sdAridity)

# for plots with ONLY 1 spp, we assume that a species #is perfectly synchronized with itself

stab_4$Gross_synchrony<-ifelse(is.na(stab_4$Gross_synchrony)==TRUE,1,stab_4$Gross_synchrony) 

# convert synchrony metrics to different scale

stab_4$GrossAsynchrony_s<-stab_4$Gross_synchrony*-1

# further adjustments

stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

stab_4$lg2_mArid  <-log(stab_4$mArid,base=2)

# Filter out NAs for Asynchrony and  FRic4

stab_4<-filter(stab_4, is.na(PlotAsynchrony_s)==FALSE)
stab_444<-filter(stab_4, is.na(FRic4)==FALSE)

############################
# fit piecewise SEM models #
############################

# Control list set up for LMM in nlme 

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")

##################
# FDis_eMNTD #####
##################

modList2=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(GrossAsynchrony_s~lg2SppN+FDis4+eMNTD+sdAridity+lg2_mArid,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+lg2SppN+eMNTD+FDis4+sdAridity+lg2_mArid,random=~1+lg2SppN|Site, control=cc,data=stab_444)
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
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FDis4, lg2_mArid, Arid_sd","eMNTD, FDis4,Asynchrony,lg2SppN,lg2_mArid")
mf_ts_emntd$ModClass<-"FDis_eMNTD"

# write results to file
write.table(ts_emntd2,"TS_emntd_fdis_sem_coefs_SIMPLE.csv",sep=",",row.names=F) # sem path coefficients
write.table(mf_ts_emntd,"TS_emntd_fdis_model_fits_SIMPLE.csv",sep=",",row.names=F) # model fits (marginal + conditional R2)
write.table(emntdfdis.fit,"TS_emntd_fdis_semfit_SIMPLE.csv",sep=",",row.names=F) #  sem model fit

#######################
## FRic4 - eMNTD    ###
#######################

modList22=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(GrossAsynchrony_s~lg2SppN+FRic4+eMNTD+sdAridity+lg2_mArid,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+lg2SppN+eMNTD+FRic4+sdAridity+lg2_mArid,random=~1+lg2SppN|Site, control=cc,data=stab_444)
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
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FRic4,mArid,Arid_sd","Asynchrony,FRic4,eMNTD,lg2SppN,mArid,Arid_sd")
mf_ts_emntd$ModClass<-"FRic_eMNTD"

# write results to file
write.table(ts_emntd2,"TS_emntd_fric_sem_coefs_SIMPLE.csv",sep=",",row.names=F) # sem path coefficients
write.table(mf_ts_emntd,"TS_emntd_fric_model_fits_SIMPLE.csv",sep=",",row.names=F) # model fits (marginal + conditional R2)
write.table(emntdfric.fit,"TS_emntd_fric_semfit_SIMPLE.csv",sep=",",row.names=F) #  sem model fit

##################
# FDis_eMPD ######
##################

modList3=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(GrossAsynchrony_s~lg2SppN+FDis4+eMPD+sdAridity+lg2_mArid,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+lg2SppN+FDis4+eMPD+sdAridity+lg2_mArid,random=~1+lg2SppN|Site, control=cc,data=stab_444)
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
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FDis4,mArid,Arid_sd","Asynchrony,FDis, eMPD,lg2SppN,mArid,Arid_sd")
mf_ts_empd$ModClass<-"FDis_eMPD"

# write results to file
write.table(ts_empd2,"TS_empd_fdis_sem_coefs_SIMPLE.csv",sep=",",row.names=F) # sem path coefficients
write.table(mf_ts_empd,"TS_empd_fdis_model_fits_SIMPLE.csv",sep=",",row.names=F) # model fits (marginal + conditional R2)
write.table(empdfdis.fit,"TS_empd_fdis_semfit_SIMPLE.csv",sep=",",row.names=F) #  sem model fit

#######################
## FRic - eMPD #######
#######################

modList33=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(GrossAsynchrony_s~lg2SppN+FRic4+eMPD+sdAridity+lg2_mArid,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~GrossAsynchrony_s+PCAdim1_4trts+lg2SppN+FRic4+eMPD+sdAridity+lg2_mArid,random=~1+lg2SppN|Site, control=cc,data=stab_444)
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
mf_ts_smpd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FRic4,mArid,Arid_sd","Asynchrony,F-S,FRic, eMPD,lg2SppN,mArid,Arid_sd")
mf_ts_smpd$ModClass<-"FRic4_eMPD"

# write results to file
write.table(ts_smpd2,"TS_empd_fric_sem_coefs_SIMPLE.csv",sep=",",row.names=F) # sem path coefficients
write.table(mf_ts_smpd,"TS_empd_fric_model_fits_SIMPLE_April2018.csv",sep=",",row.names=F) # model fits (marginal + conditional R2)
write.table(smpdfdis.fit,"TS_empd_fric_semfit_SIMPLE_April2018.csv",sep=",",row.names=F) #  sem model fit
