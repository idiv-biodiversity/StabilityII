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
stab<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,Study_length,UniqueID,SppN,eMPD,eMNTD,ePSE,sMPD,sMNTD,FDis4,FRic4,PCAdim1_4trts,SLA, LDMC, LeafN, LeafP,
               Plot_TempStab,Plot_Biomassxbar, Plot_Biomasssd, Plot_Asynchrony,annualTemp,meanPrecip,meanPET,CV_Temp,CV_Precip)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself

stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

stab_4$PlotAsynchrony_s<-stab_4$Plot_Asynchrony*-1


# Filter out NAs for Asynchrony and  FRic4
stab_444<-filter(stab_4, is.na(PlotAsynchrony_s)==FALSE)
stab_444<-filter(stab_444, is.na(FRic4)==FALSE)

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")

##############
##############

## SR on biomass

mod0<-lme(log(Plot_Biomassxbar)~log(SppN),random=~1|Site,control=cc,data=stab_444)
mod1<-lme(log(Plot_Biomassxbar)~log(SppN),random=~1+log(SppN)|Site,control=cc,data=stab_444)
AIC(mod0,mod1) #mod1

div_ranef<-ranef(mod1)
colnames(div_ranef)[2]<-"DivSlope_xbarBiomass"
div_ranef$Site<-rownames(div_ranef)
div_ranef<-select(div_ranef, Site,DivSlope_xbarBiomass)

#mod0<-lme(log(Plot_Biomassxbar)~log(SppN),random=~1|Site,control=cc,data=stab_444)
mod2<-lme(log(Plot_Biomasssd)~log(SppN),random=~1+log(SppN)|Site,control=cc,data=stab_444)

div_ranef2<-ranef(mod2)
colnames(div_ranef2)[2]<-"DivSlope_sdBiomass"
div_ranef2$Site<-rownames(div_ranef2)
div_ranef2<-select(div_ranef2, Site,DivSlope_sdBiomass)


SR_gross<-merge(div_ranef,div_ranef2,by.y="Site",all.x=TRUE)

SR_cv<- ggplot(SR_gross, aes(x=DivSlope_xbarBiomass,y=DivSlope_sdBiomass))+
         geom_point()+geom_abline(slope=1,intercept=0,linetype=2)+
         geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  labs(x="Slope of log(biomass) vs.  log(species richness)",y="Slope of log(SD) vs.  log(species richness)")+
  theme_bw()+theme(axis.title.y=element_text(colour=c("black"),face="bold",size=6),
                   axis.title.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=6),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### eMNTD on biomass

mod11<-lme(log(Plot_Biomassxbar)~log(eMNTD+0.001),random=~1+log(eMNTD+0.001)|Site,control=cc,data=stab_444)

div_ranef<-ranef(mod11)
colnames(div_ranef)[2]<-"PDSlope_xbarBiomass"
div_ranef$Site<-rownames(div_ranef)
div_ranef<-select(div_ranef, Site,PDSlope_xbarBiomass)

#mod0<-lme(log(Plot_Biomassxbar)~log(SppN),random=~1|Site,control=cc,data=stab_444)
mod22<-lme(log(Plot_Biomasssd)~log(eMNTD+0.001),random=~1+log(eMNTD+0.001)|Site,control=cc,data=stab_444)

div_ranef2<-ranef(mod22)
colnames(div_ranef2)[2]<-"PDSlope_sdBiomass"
div_ranef2$Site<-rownames(div_ranef2)
div_ranef2<-select(div_ranef2, Site,PDSlope_sdBiomass)


eMNTD_gross<-merge(div_ranef,div_ranef2,by.y="Site",all.x=TRUE)

PD_cv<- ggplot(eMNTD_gross, aes(x=PDSlope_xbarBiomass,y=PDSlope_sdBiomass))+
  geom_point()+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  labs(x="Slope of log(biomass) vs.  log(eMNTD)",y="Slope of log(SD) vs.  log(eMNTD)")+
  theme_bw()+theme(axis.title.y=element_text(colour=c("black"),face="bold",size=6),
                   axis.title.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=6),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## FD on biomass

mod33<-lme(log(Plot_Biomassxbar)~log(FDis4+0.001),random=~1+log(FDis4+0.001)|Site,control=cc,data=stab_444)

div_ranef<-ranef(mod33)
colnames(div_ranef)[2]<-"FDSlope_xbarBiomass"
div_ranef$Site<-rownames(div_ranef)
div_ranef<-select(div_ranef, Site,FDSlope_xbarBiomass)

#mod0<-lme(log(Plot_Biomassxbar)~log(SppN),random=~1|Site,control=cc,data=stab_444)
mod44<-lme(log(Plot_Biomasssd)~log(FDis4+0.001),random=~1+log(FDis4+0.001)|Site,control=cc,data=stab_444)

div_ranef2<-ranef(mod44)
colnames(div_ranef2)[2]<-"FDSlope_sdBiomass"
div_ranef2$Site<-rownames(div_ranef2)
div_ranef2<-select(div_ranef2, Site,FDSlope_sdBiomass)

FD_gross<-merge(div_ranef,div_ranef2,by.y="Site",all.x=TRUE)

FD_cv<- ggplot(FD_gross, aes(x=FDSlope_xbarBiomass,y=FDSlope_sdBiomass))+
  geom_point()+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  labs(x="Slope of log(biomass) vs.  log(FD)",y="Slope of log(SD) vs.  log(FD)")+
  theme_bw()+theme(axis.title.y=element_text(colour=c("black"),face="bold",size=6),
                   axis.title.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=6),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## F-S on biomass


mod55<-lme(log(Plot_Biomassxbar)~log(PCAdim1_4trts+3.62),random=~1+log(PCAdim1_4trts+3.62)|Site,control=cc,data=stab_444)

div_ranef<-ranef(mod55)
colnames(div_ranef)[2]<-"FSSlope_xbarBiomass"
div_ranef$Site<-rownames(div_ranef)
div_ranef<-select(div_ranef, Site,FSSlope_xbarBiomass)

#mod0<-lme(log(Plot_Biomassxbar)~log(SppN),random=~1|Site,control=cc,data=stab_444)
mod66<-lme(log(Plot_Biomasssd)~log(PCAdim1_4trts+3.62),random=~1+log(PCAdim1_4trts+3.62)|Site,control=cc,data=stab_444)

div_ranef2<-ranef(mod66)
colnames(div_ranef2)[2]<-"FSSlope_sdBiomass"
div_ranef2$Site<-rownames(div_ranef2)
div_ranef2<-select(div_ranef2, Site,FSSlope_sdBiomass)

FS_gross<-merge(div_ranef,div_ranef2,by.y="Site",all.x=TRUE)

FS_cv<- ggplot(FS_gross, aes(x=FSSlope_xbarBiomass,y=FSSlope_sdBiomass))+
  geom_point()+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  labs(x="Slope of log(biomass) vs.  log(F-S)",y="Slope of log(SD) vs.  log(F-S)")+
  theme_bw()+theme(axis.title.y=element_text(colour=c("black"),face="bold",size=6),
                   axis.title.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=6),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##########
# merge ##
##########

require(cowplot)


png(filename="/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/BDonXbarSD_GrossFig2.png", 
    units="in", 
    width=7, 
    height=7, 
    pointsize=2, 
    res=600)


plot_grid(SR_cv, PD_cv,FD_cv, FS_cv,ncol=2, labels = c("a) SR","b) PD","c) FD", "d) F-S"),label_size=6,align = "hv")


dev.off()
