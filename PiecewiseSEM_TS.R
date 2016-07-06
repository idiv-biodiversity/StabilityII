rm(list=ls()) 
#require(devtools)
#install_github("jslefche/piecewiseSEM")

require(piecewiseSEM)
library(semPlot)
library(lmerTest)
library(nlme)
library(dplyr)

# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
      #is perfectly synchronized with itself



stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

#stab_44<-filter(stab_4,SppN>1) eliminate monocultures

stab_44<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony

stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony


####################
# eMNTD + FDis #####
####################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


a<-lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
a1<-lme(eMNTD~lg2SppN,random=~1|Site,control=cc,data=stab_44)
AIC(a,a1)


re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$eMNTD)
qqline(stab_44$eMNTD)   ## use the log data


b<-lme(FDis4~lg2SppN,random=~1|Site, control=cc,data=stab_44)
b1<-lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
AIC(b,b1)


re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$FDis4)
qqline(stab_44$FDis4)   ## use the log data


c<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_44)
c1<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$Plot_Asynchrony)
qqline(stab_44$Plot_Asynchrony)   ## use the log data



d<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FDis4+eMNTD+Plot_Asynchrony,random=~1|Site, control=cc,data=stab_44)
d1<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FDis4+eMNTD+Plot_Asynchrony,random=~1+lg2SppN|Site, control=cc,data=stab_44)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$TS_lg2)
qqline(stab_44$TS_lg2)   


modList=list(
  lme(eMNTD~lg2SppN,random=~1|Site,control=cc,data=stab_44),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44),
  lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44),
  lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+Plot_Asynchrony+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_44)
)

[[STOP HERE]]
sem.fit(modList,stab_44,corr.errors=c("eMNTD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model



sem.fit(modList,stab_44,corr.errors=c("eMNTD~~FDis4","PCAdim1_4trts~~eMNTD","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model


sem.fit(modList,stab_44,corr.errors=c("eMNTD~~FDis4","PCAdim1_4trts~~eMNTD","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN",
                                      "Plot_Asynchrony~~FDis4","Plot_Asynchrony~~eMNTD"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))

ts_emntd<-sem.coefs(modList,stab_44,standardize="scale",corr.errors=c("eMNTD~~FDis4","PCAdim1_4trts~~eMNTD","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN",
                                                             "Plot_Asynchrony~~FDis4","Plot_Asynchrony~~eMNTD"))

sem.missing.paths(modList, stab_4)


mf_ts_emntd<-sem.model.fits(modList)
mf_ts_emntd$ResponseVars<-c("eMNTD","FDis4","Plot_Asynchrony","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN","F-S,eMNTD,Asychrony,FDis4,lg2SppN")

sem.plot(modList,stab_44,show.nonsig = FALSE,scaling=5)

resids.df1<-partial.resid(TS_lg2~lg2SppN,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df2<-partial.resid(TS_lg2~eMNTD,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df3<-partial.resid(TS_lg2~FDis4,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df4<-partial.resid(TS_lg2~Plot_Asynchrony,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df5<-partial.resid(TS_lg2~PCAdim1_4trts,modList,data=stab_44,list(lmeControl(opt="optim")))


write.table(ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_sem_coefs.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_model_fits.csv",sep=",",row.names=F)


###############
# eMPD + FDis #
###############

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


a<-lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
a1<-lme(eMPD~lg2SppN,random=~1|Site,control=cc,data=stab_44)
AIC(a,a1)


re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$eMPD)
qqline(stab_44$eMPD)   ## use the log data


b<-lme(FDis4~lg2SppN,random=~1|Site, control=cc,data=stab_44)
b1<-lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
AIC(b,b1)


re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$FDis4)
qqline(stab_44$FDis4)   ## use the log data


c<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_44)
c1<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$Plot_Asynchrony)
qqline(stab_44$Plot_Asynchrony)   ## use the log data

d<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FDis4+eMPD+Plot_Asynchrony,random=~1|Site, control=cc,data=stab_44)
d1<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FDis4+eMPD+Plot_Asynchrony,random=~1+lg2SppN|Site, control=cc,data=stab_44)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$TS_lg2)
qqline(stab_44$TS_lg2)   


modList=list(
  lme(eMPD~lg2SppN,random=~1|Site,control=cc,data=stab_44),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44),
  lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44),
  lme(TS_lg2~PCAdim1_4trts+FDis4+eMPD+Plot_Asynchrony+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_44)
)

sem.fit(modList,stab_44,corr.errors=c("eMPD~~FDis4","PCAdim1_4trts~~eMPD","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))

sem.fit(modList,stab_44,corr.errors=c("eMPD~~FDis4","PCAdim1_4trts~~eMPD","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN",
                                      "Plot_Asynchrony~~FDis4","Plot_Asynchrony~~eMPD"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))

ts_empd<-sem.coefs(modList,stab_44,standardize="scale",corr.errors=c("eMPD~~FDis4","PCAdim1_4trts~~eMPD","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN",
                                                                      "Plot_Asynchrony~~FDis4","Plot_Asynchrony~~eMPD"))

sem.missing.paths(modList, stab_4)


mf_ts_empd<-sem.model.fits(modList)
mf_ts_empd$ResponseVars<-c("eMPD","FDis4","Plot_Asynchrony","Temp_Stability")
mf_ts_empd$PredVars<-c("lg2SppN","lg2SppN","lg2SppN","F-S,eMPD,Asychrony,FDis4,lg2SppN")

sem.plot(modList,stab_44,show.nonsig = FALSE,scaling=5)

resids.df1<-partial.resid(TS_lg2~lg2SppN,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df2<-partial.resid(TS_lg2~eMPD,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df3<-partial.resid(TS_lg2~FDis4,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df4<-partial.resid(TS_lg2~Plot_Asynchrony,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df5<-partial.resid(TS_lg2~PCAdim1_4trts,modList,data=stab_44,list(lmeControl(opt="optim")))


write.table(ts_empd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_sem_coefs.csv",sep=",",row.names=F)
write.table(mf_ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_model_fits.csv",sep=",",row.names=F)

#################
## ePSE + FDis4 #
#################


bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


a<-lme(ePSE~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
a1<-lme(ePSE~lg2SppN,random=~1|Site,control=cc,data=stab_44)
AIC(a,a1)


re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$eMPD)
qqline(stab_44$eMPD)   ## use the log data


b<-lme(FDis4~lg2SppN,random=~1|Site, control=cc,data=stab_44)
b1<-lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
AIC(b,b1)


re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$FDis4)
qqline(stab_44$FDis4)   ## use the log data


c<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_44)
c1<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$Plot_Asynchrony)
qqline(stab_44$Plot_Asynchrony)   ## use the log data

d<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FDis4+ePSE+Plot_Asynchrony,random=~1|Site, control=cc,data=stab_44)
d1<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FDis4+ePSE+Plot_Asynchrony,random=~1+lg2SppN|Site, control=cc,data=stab_44)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_44$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$TS_lg2)
qqline(stab_44$TS_lg2)   


modList=list(
  lme(ePSE~lg2SppN,random=~1|Site,control=cc,data=stab_44),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44),
  lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_44),
  lme(TS_lg2~PCAdim1_4trts+FDis4+ePSE+Plot_Asynchrony+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_44)
)


sem.fit(modList,stab_44,corr.errors=c("ePSE~~FDis4","PCAdim1_4trts~~ePSE","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) # naive model

sem.fit(modList,stab_44,corr.errors=c("ePSE~~FDis4","PCAdim1_4trts~~ePSE","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN",
                                      "Plot_Asynchrony~~FDis4","Plot_Asynchrony~~ePSE"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))

ts_epse<-sem.coefs(modList,stab_44,standardize="scale",corr.errors=c("ePSE~~FDis4","PCAdim1_4trts~~ePSE","PCAdim1_4trts~~FDis4","PCAdim1_4trts~~lg2SppN",
                                                                     "Plot_Asynchrony~~FDis4","Plot_Asynchrony~~ePSE"))

sem.missing.paths(modList, stab_44)


mf_ts_epse<-sem.model.fits(modList)
mf_ts_epse$ResponseVars<-c("ePSE","FDis4","Plot_Asynchrony","Temp_Stability")
mf_ts_epse$PredVars<-c("lg2SppN","lg2SppN","lg2SppN","F-S,ePSE,Asychrony,FDis4,lg2SppN")

sem.plot(modList,stab_44,show.nonsig = FALSE,scaling=5)

resids.df1<-partial.resid(TS_lg2~lg2SppN,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df2<-partial.resid(TS_lg2~ePSE,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df3<-partial.resid(TS_lg2~FDis4,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df4<-partial.resid(TS_lg2~Plot_Asynchrony,modList,data=stab_44,list(lmeControl(opt="optim")))
resids.df5<-partial.resid(TS_lg2~PCAdim1_4trts,modList,data=stab_44,list(lmeControl(opt="optim")))


write.table(ts_epse,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empse_fdis_sem_coefs.csv",sep=",",row.names=F)
write.table(mf_ts_epse,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empse_fdis_model_fits.csv",sep=",",row.names=F)

################
# eMPD +FRic4  #
################


bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


a<-lme(eMPD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
a1<-lme(eMPD~lg2SppN,random=~1|Site,control=cc,data=stab_444)
AIC(a,a1)


re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$eMPD)
qqline(stab_444$eMPD)   ## use the log data


b<-lme(FRic4~lg2SppN,random=~1|Site, control=cc,data=stab_444)
b1<-lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
AIC(b,b1)


re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$FRic4)
qqline(stab_444$FRic4)   ## use the log data


c<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_444)
c1<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$Plot_Asynchrony)
qqline(stab_44$Plot_Asynchrony)   ## use the log data

d<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FRic4+eMPD+Plot_Asynchrony,random=~1|Site, control=cc,data=stab_444)
d1<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FRic4+eMPD+Plot_Asynchrony,random=~1+lg2SppN|Site, control=cc,data=stab_444)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$TS_lg2)
qqline(stab_444$TS_lg2)   


modList=list(
  lme(eMPD~lg2SppN,random=~1|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FRic4+eMPD+Plot_Asynchrony+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("eMPD~~FRic4","PCAdim1_4trts~~eMPD","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) # naive model

sem.fit(modList,stab_444,corr.errors=c("eMPD~~FRic4","PCAdim1_4trts~~eMPD","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN",
                                      "Plot_Asynchrony~~FRic4"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))

ts_empd_fric<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("eMPD~~FRic4","PCAdim1_4trts~~eMPD","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN",
                                                                      "Plot_Asynchrony~~FRic4"))

sem.missing.paths(modList, stab_444)


mf_empd_fric<-sem.model.fits(modList)
mf_empd_fric$ResponseVars<-c("eMPD","FRic4","Plot_Asynchrony","Temp_Stability")
mf_empd_fric$PredVars<-c("lg2SppN","lg2SppN","lg2SppN","F-S,eMPD,Asychrony,FRic4,lg2SppN")

sem.plot(modList,stab_444,show.nonsig = FALSE,scaling=5)

resids.df1<-partial.resid(TS_lg2~lg2SppN,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df2<-partial.resid(TS_lg2~ePSE,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df3<-partial.resid(TS_lg2~FRic4,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df4<-partial.resid(TS_lg2~Plot_Asynchrony,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df5<-partial.resid(TS_lg2~PCAdim1_4trts,modList,data=stab_444,list(lmeControl(opt="optim")))


write.table(ts_empd_fric,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_sem_coefs.csv",sep=",",row.names=F)
write.table(mf_empd_fric,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fric_model_fits.csv",sep=",",row.names=F)

##################
# eMNTD + FRic4  #
##################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


a<-lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
a1<-lme(eMNTD~lg2SppN,random=~1|Site,control=cc,data=stab_444)
AIC(a,a1)


re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$eMNTD)
qqline(stab_444$eMNTD)   ## use the log data


b<-lme(FRic4~lg2SppN,random=~1|Site, control=cc,data=stab_444)
b1<-lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
AIC(b,b1)

re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$FRic4)
qqline(stab_444$FRic4)   ## use the log data

c<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_444)
c1<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$Plot_Asynchrony)
qqline(stab_44$Plot_Asynchrony)   ## use the log data

d<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FRic4+eMNTD+Plot_Asynchrony,random=~1|Site, control=cc,data=stab_444)
d1<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FRic4+eMNTD+Plot_Asynchrony,random=~1+lg2SppN|Site, control=cc,data=stab_444)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$TS_lg2)
qqline(stab_444$TS_lg2)   


modList=list(
  lme(eMNTD~lg2SppN,random=~1|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FRic4+eMNTD+Plot_Asynchrony+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("eMNTD~~FRic4","PCAdim1_4trts~~eMNTD","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) # naive model

sem.fit(modList,stab_444,corr.errors=c("eMNTD~~FRic4","PCAdim1_4trts~~eMNTD","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN",
                                       "Plot_Asynchrony~~FRic4","Plot_Asynchrony~~eMNTD"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))

ts_emntd_fric<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("eMNTD~~FRic4","PCAdim1_4trts~~eMNTD","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN",
                                                                           "Plot_Asynchrony~~FRic4","Plot_Asynchrony~~eMNTD"))

sem.missing.paths(modList, stab_444)


mf_emntd_fric<-sem.model.fits(modList)
mf_emntd_fric$ResponseVars<-c("eMNTD","FRic4","Plot_Asynchrony","Temp_Stability")
mf_emntd_fric$PredVars<-c("lg2SppN","lg2SppN","lg2SppN","F-S,eMNTD,Asychrony,FRic4,lg2SppN")

sem.plot(modList,stab_444,show.nonsig = FALSE,scaling=5)

resids.df1<-partial.resid(TS_lg2~lg2SppN,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df2<-partial.resid(TS_lg2~eMNTD,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df3<-partial.resid(TS_lg2~FRic4,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df4<-partial.resid(TS_lg2~Plot_Asynchrony,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df5<-partial.resid(TS_lg2~PCAdim1_4trts,modList,data=stab_444,list(lmeControl(opt="optim")))


write.table(ts_emntd_fric,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_sem_coefs.csv",sep=",",row.names=F)
write.table(mf_emntd_fric,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fric_model_fits.csv",sep=",",row.names=F)

###################
## ePSE & FRic ####
###################


bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


a<-lme(ePSE~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
a1<-lme(ePSE~lg2SppN,random=~1|Site,control=cc,data=stab_444)
AIC(a,a1)


re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$eMNTD)
qqline(stab_444$eMNTD)   ## use the log data


b<-lme(FRic4~lg2SppN,random=~1|Site, control=cc,data=stab_444)
b1<-lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
AIC(b,b1)

re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$FRic4)
qqline(stab_444$FRic4)   ## use the log data

c<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_444)
c1<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_44$Plot_Asynchrony)
qqline(stab_44$Plot_Asynchrony)   ## use the log data

d<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FRic4+ePSE+Plot_Asynchrony,random=~1|Site, control=cc,data=stab_444)
d1<-lme(TS_lg2~lg2SppN+PCAdim1_4trts+FRic4+ePSE+Plot_Asynchrony,random=~1+lg2SppN|Site, control=cc,data=stab_444)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$TS_lg2)
qqline(stab_444$TS_lg2)   


modList=list(
  lme(ePSE~lg2SppN,random=~1|Site,control=cc,data=stab_444),
  lme(FRic4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FRic4+ePSE+Plot_Asynchrony+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("ePSE~~FRic4","PCAdim1_4trts~~ePSE","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) # naive model

sem.fit(modList,stab_444,corr.errors=c("ePSE~~FRic4","PCAdim1_4trts~~ePSE","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN",
                                       "Plot_Asynchrony~~FRic4","Plot_Asynchrony~~ePSE"),conditional=T,
        model.control = list(lmeControl(opt = "optim")))

ts_epse_fric<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("ePSE~~FRic4","PCAdim1_4trts~~ePSE","PCAdim1_4trts~~FRic4","PCAdim1_4trts~~lg2SppN",
                                                                           "Plot_Asynchrony~~FRic4","Plot_Asynchrony~~ePSE"))

sem.missing.paths(modList, stab_444)


mf_epse_fric<-sem.model.fits(modList)
mf_epse_fric$ResponseVars<-c("ePSE","FRic4","Plot_Asynchrony","Temp_Stability")
mf_epse_fric$PredVars<-c("lg2SppN","lg2SppN","lg2SppN","F-S,ePSE,Asychrony,FRic4,lg2SppN")

sem.plot(modList,stab_444,show.nonsig = FALSE,scaling=5)

resids.df1<-partial.resid(TS_lg2~lg2SppN,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df2<-partial.resid(TS_lg2~ePSE,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df3<-partial.resid(TS_lg2~FRic4,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df4<-partial.resid(TS_lg2~Plot_Asynchrony,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df5<-partial.resid(TS_lg2~PCAdim1_4trts,modList,data=stab_444,list(lmeControl(opt="optim")))


write.table(ts_epse_fric,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_sem_coefs.csv",sep=",",row.names=F)
write.table(mf_epse_fric,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_epse_fric_model_fits.csv",sep=",",row.names=F)

