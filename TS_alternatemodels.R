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
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,sMPD,sMNTD,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony,annualTemp,meanPrecip,meanPET,CV_Temp,CV_Precip)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)


stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony


####################
# Prelim       #####
####################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


a<-lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
a1<-lme(eMNTD~lg2SppN,random=~1|Site,control=cc,data=stab_444)

AIC(a,a1)

plot(a1)
qqnorm(a1)

aa<-lme(sMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
aa1<-lme(sMNTD~lg2SppN,random=~1|Site,control=cc,data=stab_444)

AIC(aa,aa1)

plot(aa)
qqnorm(aa)



b<-lme(FDis4~lg2SppN,random=~1|Site, control=cc,data=stab_444)
b1<-lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
AIC(b,b1)

plot(b1)
qqnorm(b1)


c<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_444)
c1<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)


AIC(c,c1)

plot(c1)
qqnorm(c1)

c<-lme(Plot_Asynchrony~FDis4,random=~1|Site,control=cc,data=stab_444)
c1<-lme(Plot_Asynchrony~FDis4,random=~1+lg2SppN|Site,control=cc,data=stab_444)
c2<-lme(Plot_Asynchrony~FDis4,random=~1+FDis4|Site/SppN,control=cc,data=stab_444)

AIC(c,c1,c2)

plot(c2)
qqnorm(c2)

d<-lme(Plot_Asynchrony~PCAdim1_4trts,random=~1|Site,control=cc,data=stab_444)
d1<-lme(Plot_Asynchrony~PCAdim1_4trts,random=~1+lg2SppN|Site,control=cc,data=stab_444)
d2<-lme(Plot_Asynchrony~PCAdim1_4trts,random=~1+PCAdim1_4trts|Site/SppN,control=cc,data=stab_444)


AIC(d,d1,d2)

plot(d1)
qqnorm(d1)


e<-lme(Plot_Asynchrony~eMNTD,random=~1|Site,control=cc,data=stab_444)
e1<-lme(Plot_Asynchrony~eMNTD,random=~1+lg2SppN|Site,control=cc,data=stab_444)
e2<-lme(Plot_Asynchrony~eMNTD,random=~1+eMNTD|Site/SppN,control=cc,data=stab_444)

AIC(e,e1,e2)

plot(e2)
qqnorm(e2)

f<-lme(Plot_Asynchrony~CV_Precip,random=~1|Site,control=cc,data=stab_444)
f1<-lme(Plot_Asynchrony~CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_444)

AIC(f,f1)

plot(f1)
qqnorm(f1)

g<-lme(TS_lg2~CV_Precip,random=~1|Site,control=cc,data=stab_444)
g1<-lme(TS_lg2~CV_Precip,random=~1+lg2SppN|Site,control=cc,data=stab_444)

AIC(g,g1)

plot(f1)
qqnorm(f1)



g<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1|Site,control=cc,data=stab_444)
g1<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+lg2SppN|Site,control=cc,data=stab_444)
g2<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+eMNTD|Site,control=cc,data=stab_444)
g3<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+FDis4|Site,control=cc,data=stab_444)
g4<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+lg2SppN+PCAdim1_4trts|Site,control=cc,data=stab_444)
g5<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+lg2SppN+eMNTD|Site,control=cc,data=stab_444)
g6<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+lg2SppN+FDis4|Site,control=cc,data=stab_444)
g7<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+PCAdim1_4trts+FDis4|Site,control=cc,data=stab_444)
g8<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+PCAdim1_4trts+eMNTD|Site,control=cc,data=stab_444)
g9<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+PCAdim1_4trts+FDis4|Site,control=cc,data=stab_444)
g10<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+PCAdim1_4trts+FDis4|Site,control=cc,data=stab_444)


AICc(g,g1,g2,g3,g4)

plot(e2)
qqnorm(e2)



f<-lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_444)
f1<-lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN,random=~1|Site, control=cc,data=stab_444)
f3<-lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN,random=~1|Site/SppN, control=cc,data=stab_444)



AIC(f,f1,f3)


######################
# piecewise SEM  #####
# with synchrony  #
######################


##################
# FDis_eMNTD #####
##################

modList2=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FDis4+PCAdim1_4trts+eMNTD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+eMNTD+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


lapply(modList2, plot)

# Explore distribution of residuals

lapply(modList2, function(i) hist(resid(i)))

# Look at variance inflation factors
lapply(modList2[3:4], vif)


sem.fit(modList2,stab_444,corr.errors=c("eMNTD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

emntdfdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("eMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))  # add eMNTD as predictor of TS

emntdfdis.fit<-cbind(emntdfdis.fit$Fisher.C,emntdfdis.fit$AIC)
emntdfdis.fit$ModClass<-"FDis_eMNTD"

ts_emntd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("eMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_emntd2$ModClass<-"FDis_eMNTD"


mf_ts_emntd<-sem.model.fits(modList2)
mf_ts_emntd$ResponseVars<-c("eMNTD","FDis4","Asynchrony","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FDis4,F-S","eMNTD,Asynchrony,lg2SppN,CV_Precip")
mf_ts_emntd$ModClass<-"FDis_eMNTD"

write.table(ts_emntd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(emntdfdis.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_semfit_Oct2016.csv",sep=",",row.names=F)

##################
# FDis_eMPD ######
##################


modList2=list(
  lme(sMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FDis4+PCAdim1_4trts+sMPD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("sMPD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

empdfdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("sMPD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
                       model.control = list(lmeControl(opt = "optim")))  # did not add eMPD as predictor of TS


empdfdis.fit<-cbind(empdfdis.fit$Fisher.C,empdfdis.fit$AIC)
empdfdis.fit$ModClass<-"FDis_eMPD"

ts_empd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("sMPD~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_empd2$ModClass<-"FDis_eMPD"


mf_ts_empd<-sem.model.fits(modList2)
mf_ts_empd$ResponseVars<-c("eMPD","FDis4","Asynchrony","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FDis4,F-S","Asynchrony,lg2SppN,CV_Precip")
mf_ts_empd$ModClass<-"FDis_eMPD"

write.table(ts_empd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_empd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(empdfdis.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_semfit_Oct2016.csv",sep=",",row.names=F)

#######################
## FDis - sMPD #######
#######################

modList2=list(
  lme(sMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FDis4+PCAdim1_4trts+sMPD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+lg2SppN+sMPD+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("sMPD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

smpdfdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("sMPD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  


smpdfdis.fit<-cbind(smpdfdis.fit$Fisher.C,smpdfdis.fit$AIC)
smpdfdis.fit$ModClass<-"FDis_sMPD"

ts_smpd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("sMPD~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_smpd2$ModClass<-"FDis_sMPD"


mf_ts_smpd<-sem.model.fits(modList2)
mf_ts_smpd$ResponseVars<-c("sMPD","FDis4","Asynchrony","Temp_Stability")
mf_ts_smpd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,sMPD,FDis4,F-S","Asynchrony,sMPD,lg2SppN,CV_Precip")
mf_ts_smpd$ModClass<-"FDis_sMPD"

write.table(ts_smpd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fdis_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_smpd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fdis_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(smpdfdis.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fdis_semfit_Oct2016.csv",sep=",",row.names=F)

#######################
# FDis -sMNTD #########
#######################

modList2=list(
  lme(sMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FDis4+PCAdim1_4trts+sMNTD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+sMNTD+PCAdim1_4trts+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("sMNTD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

smmntdfdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("sMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  # did not add eMPD as predictor of TS


smmntdfdis.fit<-cbind(smmntdfdis.fit$Fisher.C,smmntdfdis.fit$AIC)
smmntdfdis.fit$ModClass<-"FDis_sMNTD"

ts_smntd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("sMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_smntd2$ModClass<-"FDis_sMNTD"


mf_ts_smntd<-sem.model.fits(modList2)
mf_ts_smntd$ResponseVars<-c("sMNTD","FDis4","Asynchrony","Temp_Stability")
mf_ts_smntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,sMPD,FDis4,F-S","Asynchrony,sMNTD,lg2SppN,CV_Precip")
mf_ts_smntd$ModClass<-"FDis_sMNTD"

write.table(ts_smntd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fdis_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_smntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fdis_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(smmntdfdis.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fdis_semfit_Oct2016.csv",sep=",",row.names=F)


#######################
## FRic4 - eMNTD    ###
#######################



modList2=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FRic4+PCAdim1_4trts+eMNTD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+eMNTD+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("eMNTD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

emntdfric.fit<-sem.fit(modList2,stab_444,corr.errors=c("FRic4~~eMNTD","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                       model.control = list(lmeControl(opt = "optim")))  # add eMNTD as predictor of TS

emntdfric.fit<-cbind(emntdfric.fit$Fisher.C,emntdfric.fit$AIC)
emntdfric.fit$ModClass<-"FRic_eMNTD"

ts_emntd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("eMNTD~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_emntd2$ModClass<-"FRic_eMNTD"

mf_ts_emntd<-sem.model.fits(modList2)
mf_ts_emntd$ResponseVars<-c("eMNTD","FRic4","Asynchrony","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMNTD,FRic4,F-S","Asynchrony,eMNTD,lg2SppN,CV_Precip")
mf_ts_emntd$ModClass<-"FRic_eMNTD"

write.table(ts_emntd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(emntdfric.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_semfit_Oct2016.csv",sep=",",row.names=F)

##################
# FRic_eMPD ######
##################


modList2=list(
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FRic4+PCAdim1_4trts+eMPD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)

sem.fit(modList2,stab_444,corr.errors=c("eMPD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

empdfric.fit<-sem.fit(modList2,stab_444,corr.errors=c("eMPD~~FRic4","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  # did not add eMPD as predictor of TS


empdfric.fit<-cbind(empdfric.fit$Fisher.C,empdfric.fit$AIC)
empdfric.fit$ModClass<-"FRic_eMPD"

ts_empd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("eMPD~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_empd2$ModClass<-"FRic_eMPD"


mf_ts_empd<-sem.model.fits(modList2)
mf_ts_empd$ResponseVars<-c("eMPD","FRic4","Asynchrony","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FRic4,F-S","Asynchrony,lg2SppN,CV_Precip")
mf_ts_empd$ModClass<-"FRic_eMPD"

write.table(ts_empd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_empd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(empdfric.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_semfit_Oct2016.csv",sep=",",row.names=F)

#######################
## FRic - sMPD #######
#######################

modList2=list(
  lme(sMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FRic4+PCAdim1_4trts+sMPD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+sMPD+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("sMPD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

smpdfdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("sMPD~~FRic4","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  


smpdfdis.fit<-cbind(smpdfdis.fit$Fisher.C,smpdfdis.fit$AIC)
smpdfdis.fit$ModClass<-"FRic_sMPD"

ts_smpd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("sMPD~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_smpd2$ModClass<-"FRic_sMPD"


mf_ts_smpd<-sem.model.fits(modList2)
mf_ts_smpd$ResponseVars<-c("sMPD","FRic4","Asynchrony","Temp_Stability")
mf_ts_smpd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,sMPD,FRic4,F-S","Asynchrony,sMPD,lg2SppN,CV_Precip")
mf_ts_smpd$ModClass<-"FRic_sMPD"

write.table(ts_smpd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fric_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_smpd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fric_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(smpdfdis.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smpd_fric_semfit_Oct2016.csv",sep=",",row.names=F)

#######################
# FRic -sMNTD #########
#######################

modList2=list(
  lme(sMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FRic4+PCAdim1_4trts+sMNTD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+sMNTD+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("sMNTD~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

smmntdfdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("sMNTD~~FRic4","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                        model.control = list(lmeControl(opt = "optim")))  # 

smmntdfdis.fit<-cbind(smmntdfdis.fit$Fisher.C,smmntdfdis.fit$AIC)
smmntdfdis.fit$ModClass<-"FRic_sMNTD"

ts_smntd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("sMNTD~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_smntd2$ModClass<-"FRic_sMNTD"


mf_ts_smntd<-sem.model.fits(modList2)
mf_ts_smntd$ResponseVars<-c("sMNTD","FRic4","Asynchrony","Temp_Stability")
mf_ts_smntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,sMPD,FRic4,F-S","Asynchrony,sMNTD,lg2SppN,CV_Precip")
mf_ts_smntd$ModClass<-"FRic_sMNTD"

write.table(ts_smntd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fric_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_smntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fric_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(smmntdfdis.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_smntd_fric_semfit_Oct2016.csv",sep=",",row.names=F)


####################
# Merge models #####
####################

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



##################################
## JUNK: ePSE models #############
###################################



###############
## FD - ePSE ##
###############

modList2=list(
  lme(ePSE~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FDis4+PCAdim1_4trts+ePSE,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("ePSE~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

epsefdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("ePSE~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  # did not add ePSE as predictor of TS

epsefdis.fit<-cbind(epsefdis.fit$Fisher.C,epsefdis.fit$AIC)
epsefdis.fit$ModClass<-"FDis_ePSE"

ts_epse2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("ePSE~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_epse2$ModClass<-"FDis_ePSE"


mf_ts_epse<-sem.model.fits(modList2)
mf_ts_epse$ResponseVars<-c("ePSE","FDis4","Asynchrony","Temp_Stability")
mf_ts_epse$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,ePSE,FDis4,F-S","Asynchrony,lg2SppN,CV_Precip")
mf_ts_epse$ModClass<-"FDis_ePSE"

write.table(ts_epse2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_epse,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(epsefdis.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fdis_semfit_Oct2016.csv",sep=",",row.names=F)


################
## Fric - ePSE #
################

modList2=list(
  lme(ePSE~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FRic4+PCAdim1_4trts+ePSE,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+lg2SppN+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("ePSE~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

epsefric.fit<-sem.fit(modList2,stab_444,corr.errors=c("ePSE~~FRic4","FRic4 ~~ PCAdim1_4trts"),conditional=T,
                      model.control = list(lmeControl(opt = "optim")))  # did not add ePSE as predictor of TS

epsefric.fit<-cbind(epsefric.fit$Fisher.C,epsefric.fit$AIC)
epsefric.fit$ModClass<-"FRic_ePSE"

ts_epse2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("ePSE~~FRic4","FRic4 ~~ PCAdim1_4trts"))
ts_epse2$ModClass<-"FRic_ePSE"


mf_ts_epse<-sem.model.fits(modList2)
mf_ts_epse$ResponseVars<-c("ePSE","FRic4","Asynchrony","Temp_Stability")
mf_ts_epse$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,ePSE,FRic4,F-S","Asynchrony,lg2SppN,CV_Precip")
mf_ts_epse$ModClass<-"FRic_ePSE"

write.table(ts_epse2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_epse,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(epsefric.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_semfit_Oct2016.csv",sep=",",row.names=F)



