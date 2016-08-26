rm(list=ls()) 

library(nlme)
library(dplyr)

##########################################
# why does stability increase with       #
#increasing dominance of resource        #
#acquisitive species?                    #
##########################################


# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Biomassxbar,Plot_Biomasssd,Plot_Asynchrony)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

stab_4$lg2Biomass <- log(stab_4$Plot_Biomassxbar,2)

stab_4$lg2SDBiomass <- log(stab_4$Plot_Biomasssd,2)

stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony


#################
# fit model #####
#################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


#Model 1: Biomass and SppN

a<-lme(lg2Biomass~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
a1<-lme(lg2Biomass~lg2SppN,random=~1|Site,control=cc,data=stab_444)
AIC(a,a1)


re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$lg2Biomass)
qqline(stab_444$lg2Biomass)   ## use the log data

# Plot Biomass increases with SppN (coefficient = 0.594344)


#Model 2: SD of biomass and SppN


b<-lme(lg2SDBiomass~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b1<-lme(lg2SDBiomass~lg2SppN,random=~1|Site,control=cc,data=stab_444)
AIC(b,b1)


re<-resid(b, type="normalized")
fi<-fitted(b)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$lg2SppN,y=re,xlab="SppN",ylab="residuals")

qqnorm(stab_444$lg2SDBiomass)
qqline(stab_444$lg2SDBiomass)   ## use the log data

# SD of plot biomass also increased with SppN (coefficient = 0.311863)


# Model 3:  increase of biomass with fast-slow

a<-lme(lg2Biomass~PCAdim1_4trts,random=~1+lg2SppN|Site,control=cc,data=stab_444)
a1<-lme(lg2Biomass~PCAdim1_4trts,random=~1|Site,control=cc,data=stab_444)
AIC(a,a1)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_444$PCAdim1_4trts,y=re,xlab="SppN",ylab="residuals")


# Plot Biomass increases with f-s axis (coefficient = 0.184024)


# Model 4:  increase of  SD of biomass with fast-slow

a<-lme(lg2SDBiomass~PCAdim1_4trts,random=~1+lg2SppN|Site,control=cc,data=stab_4)
a1<-lme(lg2SDBiomass~PCAdim1_4trts,random=~1|Site,control=cc,data=stab_4)
AIC(a,a1)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=stab_4$PCAdim1_4trts,y=re,xlab="SppN",ylab="residuals")


# Plot Biomass increases with f-s axis (coefficient = 0.060625, not significant)


