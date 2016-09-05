rm(list=ls()) 

library(lmerTest)
library(nlme)
library(dplyr)
library(lm.beta)

# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,Study_length,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,LeafN, LeafP,SLA, LDMC,
               Plot_Biomassxbar, Plot_Biomasssd,Plot_Asynchrony,Plot_TempStab)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg<-log(stab_4$Plot_TempStab)

stab_4$lgSppN <- log(stab_4$SppN)

stab_4$lg_meanbiom <- log(stab_4$Plot_Biomassxbar)
stab_4$lg_sdbiom <- log(stab_4$Plot_Biomasssd)

stab_4$lg_FDis<-log(stab_4$FDis4+.01)
stab_4$lg_FRic<-log(stab_4$FRic4+.01)

stab_4$lg_eMPD<-log(stab_4$eMPD+.01)
stab_4$lg_eMNTD<-log(stab_4$eMNTD+.01)
stab_4$lg_ePSE<-log(stab_4$ePSE+.01)

stab_4$lg_FS<-log(stab_4$PCAdim1_4trts+4.275)
stab_4$lg_P<-log(stab_4$LeafP)
stab_4$lg_N<-log(stab_4$LeafN)
stab_4$lg_SLA<-log(stab_4$SLA)
stab_4$lg_LDMC<-log(stab_4$LDMC)


stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony


#Drop studies without monocultures

study_div<-summarize(group_by(stab_444,Site),min_spp=min(SppN),max_spp=max(SppN))

study_div$keep<-ifelse(study_div$min_spp>1,NA,1)

study_div<-select(study_div,Site,keep)

stab_g<-merge(stab_444,study_div,by.y="Site")
stab_g<-filter(stab_g,keep==1)

###########################################
# Step 1: lm regression / loop through all#
###########################################


n<-length(unique(stab_g$Site))

outt=c();

for(i in 1:n){
  
  test=subset(stab_g, stab_g$Site==(unique(stab_g$Site))[i])  
  Site<-as.character(unique(test$Site))
  Study_length<-max(test$Study_length)

  ##SPP
biom_lm<-lm(lg_meanbiom~lgSppN,data=test)

biom_slope<-summary(biom_lm)$coefficient[2,1]

biom_var<-vcov(biom_lm)[2,2]

sd_lm<-lm(lg_sdbiom~lgSppN,data=test)

sd_slope<-summary(sd_lm)$coefficient[2,1]
sd_var<-vcov(sd_lm)[2,2]

##FD
fd_biom_lm<-lm(lg_meanbiom~lg_FDis,data=test)

fd_biom_slope<-summary(fd_biom_lm)$coefficient[2,1]

fd_biom_var<-vcov(fd_biom_lm)[2,2]

fd_sd_lm<-lm(lg_sdbiom~lg_FDis,data=test)

fd_sd_slope<-summary(fd_sd_lm)$coefficient[2,1]
fd_sd_var<-vcov(fd_sd_lm)[2,2]

#FRic
fr_biom_lm<-lm(lg_meanbiom~lg_FRic,data=test)

fr_biom_slope<-summary(fr_biom_lm)$coefficient[2,1]

fr_biom_var<-vcov(fr_biom_lm)[2,2]

fr_sd_lm<-lm(lg_sdbiom~lg_FRic,data=test)

fr_sd_slope<-summary(fr_sd_lm)$coefficient[2,1]
fr_sd_var<-vcov(fr_sd_lm)[2,2]


### PD

#ePSE
pse_biom_lm<-lm(lg_meanbiom~lg_ePSE,data=test)

pse_biom_slope<-summary(pse_biom_lm)$coefficient[2,1]

pse_biom_var<-vcov(pse_biom_lm)[2,2]

pse_sd_lm<-lm(lg_sdbiom~lg_ePSE,data=test)

pse_sd_slope<-summary(pse_sd_lm)$coefficient[2,1]
pse_sd_var<-vcov(pse_sd_lm)[2,2]

#eMNTD
mntd_biom_lm<-lm(lg_meanbiom~lg_eMNTD,data=test)

mntd_biom_slope<-summary(mntd_biom_lm)$coefficient[2,1]

mntd_biom_var<-vcov(mntd_biom_lm)[2,2]

mntd_sd_lm<-lm(lg_sdbiom~lg_eMNTD,data=test)

mntd_sd_slope<-summary(mntd_sd_lm)$coefficient[2,1]
mntd_sd_var<-vcov(mntd_sd_lm)[2,2]

#eMPD
mpd_biom_lm<-lm(lg_meanbiom~lg_eMPD,data=test)

mpd_biom_slope<-summary(mpd_biom_lm)$coefficient[2,1]

mpd_biom_var<-vcov(mpd_biom_lm)[2,2]

mpd_sd_lm<-lm(lg_sdbiom~lg_eMPD,data=test)

mpd_sd_slope<-summary(mpd_sd_lm)$coefficient[2,1]
mpd_sd_var<-vcov(mpd_sd_lm)[2,2]

##CWMs

fs_biom_lm<-lm(lg_meanbiom~lg_FS,data=test)

fs_biom_slope<-summary(fs_biom_lm)$coefficient[2,1]

fs_biom_var<-vcov(fs_biom_lm)[2,2]

fs_sd_lm<-lm(lg_sdbiom~lg_FS,data=test)

fs_sd_slope<-summary(fs_sd_lm)$coefficient[2,1]
fs_sd_var<-vcov(fs_sd_lm)[2,2]

##leafN

##leaf P

##SLA

##LDMC


####
study<-cbind.data.frame(Site,Study_length,biom_slope,biom_var,sd_slope,sd_var)

outt[[i]]<-rbind.data.frame(study)

}

jjj<-do.call(rbind,outt)  
jjj<-data.frame(jjj)

write.table(jjj,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Div_ts.csv",sep=",",row.names=F)


######################################
## LME           #####################
######################################


#########
# SR ####
#########

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")

#biomass
a<-lme(lg2_meanbiom~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_g)
a1<-lme(lg2_meanbiom~lg2SppN,random=~1|Site,control=cc,data=stab_g)
AIC(a,a1)

qqnorm(a)
plot(a)

sr_biom<-intervals(a)$fixed[2,]


#sd
b<-lme(lg2_sdbiom~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_g)
b1<-lme(lg2_sdbiom~lg2SppN,random=~1|Site,control=cc,data=stab_g)
AIC(b,b1)

qqnorm(b)
plot(b)

sr_biomsd<-intervals(b)$fixed[2,]


#########
# FD ####
#########

#biomass

c<-lme(lg2_meanbiom~lg2_FDis,random=~1+lg2SppN|Site,control=cc,data=stab_g)
c1<-lme(lg2_meanbiom~lg2_FDis,random=~1|Site,control=cc,data=stab_g)
AIC(c,c1)

qqnorm(c)
plot(c)

fd_biom<-intervals(c)$fixed[2,]


#sd
b<-lme(lg2_sdbiom~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_g)
b1<-lme(lg2_sdbiom~lg2SppN,random=~1|Site,control=cc,data=stab_g)
AIC(b,b1)

qqnorm(b)
plot(b)

sr_biomsd<-intervals(b)$fixed[2,]


#########
# PD ####
#########

#eMPD
c<-lme(lg2_meanbiom~lg2_FDis,random=~1+lg2SppN|Site,control=cc,data=stab_g)
c1<-lme(lg2_meanbiom~lg2_FDis,random=~1|Site,control=cc,data=stab_g)
AIC(c,c1)

qqnorm(c)
plot(c)

fd_biom<-intervals(c)$fixed[2,]


#sd
b<-lme(lg2_sdbiom~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_g)
b1<-lme(lg2_sdbiom~lg2SppN,random=~1|Site,control=cc,data=stab_g)
AIC(b,b1)

qqnorm(b)
plot(b)

sr_biomsd<-intervals(b)$fixed[2,]

#eMNTD

#ePSE

##########
# F-S ####
##########

#F-S

#leafN

#leafP

#SLA

#LDMC