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


modList2=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FDis4+PCAdim1_4trts+eMNTD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+lg2SppN+eMNTD+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


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
  lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN+FDis4+PCAdim1_4trts+eMPD,random=~1+lg2SppN|Site,control=bb,data=stab_444),
  lme(TS_lg2~Plot_Asynchrony+PCAdim1_4trts+lg2SppN+eMPD+CV_Precip,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList2,stab_444,corr.errors=c("eMPD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

empdfdis.fit<-sem.fit(modList2,stab_444,corr.errors=c("eMPD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
                       model.control = list(lmeControl(opt = "optim")))  # add eMNTD as predictor of TS

empdfdis.fit<-cbind(empdfdis.fit$Fisher.C,empdfdis.fit$AIC)
empdfdis.fit$ModClass<-"FDis_eMPD"

ts_empd2<-sem.coefs(modList2,stab_444,standardize="scale",corr.errors=c("eMPD~~FDis4","FDis4 ~~ PCAdim1_4trts"))
ts_empd2$ModClass<-"FDis_eMPD"


mf_ts_empd<-sem.model.fits(modList2)
mf_ts_empd$ResponseVars<-c("eMPD","FDis4","Asynchrony","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN,eMPD,FDis4,F-S","eMPD,Asynchrony,lg2SppN,CV_Precip")
mf_ts_empd$ModClass<-"FDis_eMPD"

write.table(ts_empd2,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_sem_coefs_Oct2016.csv",sep=",",row.names=F)
write.table(mf_ts_empd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_model_fits_Oct2016.csv",sep=",",row.names=F)
write.table(empdfdis.fit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_semfit_Oct2016.csv",sep=",",row.names=F)

[STOP HERE]