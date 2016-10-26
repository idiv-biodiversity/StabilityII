rm(list=ls()) 

library(lmerTest)
library(nlme)
library(dplyr)
library(lm.beta)

# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,Study_length,SppN,eMPD,eMNTD,ePSE,sMPD,sMNTD,FDis4,FRic4,PCAdim1_4trts,LeafN, LeafP,SLA, LDMC,
               Plot_Biomassxbar, Plot_Biomasssd,Plot_Asynchrony,Plot_TempStab)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) 

# for monocultures, we assume that a species is perfectly synchronized with itself


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

stab_4$lg_FS<-log(stab_4$PCAdim1_4trts+3.62)
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

biom_b<-lm.beta(biom_lm)
biom_beta<-coef(biom_b)[2]

sd_lm<-lm(lg_sdbiom~lgSppN,data=test)

sd_slope<-summary(sd_lm)$coefficient[2,1]
sd_var<-vcov(sd_lm)[2,2]

sd_b<-lm.beta(sd_lm)
sd_beta<-coef(sd_b)[2]

##FD
fd_biom_lm<-lm(lg_meanbiom~lg_FDis,data=test)

fd_biom_slope<-summary(fd_biom_lm)$coefficient[2,1]

fd_biom_var<-vcov(fd_biom_lm)[2,2]


fd_b<-lm.beta(fd_biom_lm)
fd_beta<-coef(fd_b)[2]


fd_sd_lm<-lm(lg_sdbiom~lg_FDis,data=test)

fd_sd_slope<-summary(fd_sd_lm)$coefficient[2,1]
fd_sd_var<-vcov(fd_sd_lm)[2,2]

fd_sd<-lm.beta(fd_sd_lm)
fd_sd_beta<-coef(fd_sd)[2]


#FRic
fr_biom_lm<-lm(lg_meanbiom~lg_FRic,data=test)

fr_biom_slope<-summary(fr_biom_lm)$coefficient[2,1]

fr_biom_var<-vcov(fr_biom_lm)[2,2]

fr_b<-lm.beta(fr_biom_lm)
fr_beta<-coef(fr_b)[2]


fr_sd_lm<-lm(lg_sdbiom~lg_FRic,data=test)

fr_sd_slope<-summary(fr_sd_lm)$coefficient[2,1]
fr_sd_var<-vcov(fr_sd_lm)[2,2]

fr_sd<-lm.beta(fr_sd_lm)
fr_sd_beta<-coef(fr_sd)[2]


### PD

#ePSE
pse_biom_lm<-lm(lg_meanbiom~lg_ePSE,data=test)

pse_biom_slope<-summary(pse_biom_lm)$coefficient[2,1]

pse_biom_var<-vcov(pse_biom_lm)[2,2]

pse_b<-lm.beta(pse_biom_lm)
pse_beta<-coef(pse_b)[2]

pse_sd_lm<-lm(lg_sdbiom~lg_ePSE,data=test)

pse_sd_slope<-summary(pse_sd_lm)$coefficient[2,1]
pse_sd_var<-vcov(pse_sd_lm)[2,2]

pse_sd<-lm.beta(pse_sd_lm)
pse_sd_beta<-coef(pse_sd)[2]


#eMNTD
mntd_biom_lm<-lm(lg_meanbiom~lg_eMNTD,data=test)

mntd_biom_slope<-summary(mntd_biom_lm)$coefficient[2,1]

mntd_biom_var<-vcov(mntd_biom_lm)[2,2]

mntd_sd_lm<-lm(lg_sdbiom~lg_eMNTD,data=test)

mntd_sd_slope<-summary(mntd_sd_lm)$coefficient[2,1]
mntd_sd_var<-vcov(mntd_sd_lm)[2,2]

mntd_b<-lm.beta(mntd_biom_lm)
mntd_beta<-coef(mntd_b)[2]

mntd_sd<-lm.beta(mntd_sd_lm)
mntd_sd_beta<-coef(mntd_sd)[2]


#eMPD
mpd_biom_lm<-lm(lg_meanbiom~lg_eMPD,data=test)

mpd_biom_slope<-summary(mpd_biom_lm)$coefficient[2,1]

mpd_biom_var<-vcov(mpd_biom_lm)[2,2]

mpd_sd_lm<-lm(lg_sdbiom~lg_eMPD,data=test)

mpd_sd_slope<-summary(mpd_sd_lm)$coefficient[2,1]
mpd_sd_var<-vcov(mpd_sd_lm)[2,2]

mpd_b<-lm.beta(mpd_biom_lm)
mpd_beta<-coef(mpd_b)[2]

mpd_sd<-lm.beta(mpd_sd_lm)
mpd_sd_beta<-coef(mpd_sd)[2]


##CWMs

fs_biom_lm<-lm(lg_meanbiom~lg_FS,data=test)

fs_biom_slope<-summary(fs_biom_lm)$coefficient[2,1]

fs_biom_var<-vcov(fs_biom_lm)[2,2]

fs_sd_lm<-lm(lg_sdbiom~lg_FS,data=test)

fs_sd_slope<-summary(fs_sd_lm)$coefficient[2,1]
fs_sd_var<-vcov(fs_sd_lm)[2,2]

fs_b<-lm.beta(fs_biom_lm)
fs_beta<-coef(fs_b)[2]

fs_sd<-lm.beta(fs_sd_lm)
fs_sd_beta<-coef(fs_sd)[2]



##leafN


n_biom_lm<-lm(lg_meanbiom~lg_N,data=test)

n_biom_slope<-summary(n_biom_lm)$coefficient[2,1]

n_biom_var<-vcov(n_biom_lm)[2,2]

n_sd_lm<-lm(lg_sdbiom~lg_N,data=test)

n_sd_slope<-summary(n_sd_lm)$coefficient[2,1]
n_sd_var<-vcov(n_sd_lm)[2,2]

n_b<-lm.beta(n_biom_lm)
n_beta<-coef(n_b)[2]

n_sd<-lm.beta(n_sd_lm)
n_sd_beta<-coef(n_sd)[2]


##leaf P

p_biom_lm<-lm(lg_meanbiom~lg_P,data=test)

p_biom_slope<-summary(p_biom_lm)$coefficient[2,1]

p_biom_var<-vcov(p_biom_lm)[2,2]

p_sd_lm<-lm(lg_sdbiom~lg_P,data=test)

p_sd_slope<-summary(p_sd_lm)$coefficient[2,1]
p_sd_var<-vcov(p_sd_lm)[2,2]


p_b<-lm.beta(p_biom_lm)
p_beta<-coef(p_b)[2]

p_sd<-lm.beta(p_sd_lm)
p_sd_beta<-coef(p_sd)[2]


##SLA

sla_biom_lm<-lm(lg_meanbiom~lg_SLA,data=test)

sla_biom_slope<-summary(sla_biom_lm)$coefficient[2,1]

sla_biom_var<-vcov(sla_biom_lm)[2,2]

sla_sd_lm<-lm(lg_sdbiom~lg_SLA,data=test)

sla_sd_slope<-summary(sla_sd_lm)$coefficient[2,1]
sla_sd_var<-vcov(sla_sd_lm)[2,2]

sla_b<-lm.beta(sla_biom_lm)
sla_beta<-coef(sla_b)[2]

sla_sd<-lm.beta(sla_sd_lm)
sla_sd_beta<-coef(sla_sd)[2]


##LDMC

ldmc_biom_lm<-lm(lg_meanbiom~lg_LDMC,data=test)

ldmc_biom_slope<-summary(ldmc_biom_lm)$coefficient[2,1]

ldmc_biom_var<-vcov(ldmc_biom_lm)[2,2]

ldmc_sd_lm<-lm(lg_sdbiom~lg_LDMC,data=test)

ldmc_sd_slope<-summary(ldmc_sd_lm)$coefficient[2,1]
ldmc_sd_var<-vcov(ldmc_sd_lm)[2,2]

ldmc_b<-lm.beta(ldmc_biom_lm)
ldmc_beta<-coef(ldmc_b)[2]

ldmc_sd<-lm.beta(ldmc_sd_lm)
ldmc_sd_beta<-coef(ldmc_sd)[2]



####
study<-cbind.data.frame(Site,Study_length,biom_slope,biom_var,sd_slope,sd_var,
                        fd_biom_slope,fd_biom_var,fd_sd_slope,fd_sd_var,
                        fr_biom_slope,fr_biom_var,fr_sd_slope,fr_sd_var,
                        mpd_biom_slope,mpd_biom_var,mpd_sd_slope,mpd_sd_var,
                        pse_biom_slope,pse_biom_var,pse_sd_slope,pse_sd_var,
                        mntd_biom_slope,mntd_biom_var,mntd_sd_slope,mntd_sd_var,
                        fs_biom_slope,fs_biom_var,fs_sd_slope,fs_sd_var,
                        n_biom_slope,n_biom_var,n_sd_slope,n_sd_var,
                        p_biom_slope,p_biom_var,p_sd_slope,p_sd_var,
                        sla_biom_slope,sla_biom_var,sla_sd_slope,sla_sd_var,
                        ldmc_biom_slope,ldmc_biom_var,ldmc_sd_slope,ldmc_sd_var,
                        
                        biom_beta,sd_beta,
                        fd_beta, fd_sd_beta,
                        fr_beta, fr_sd_beta,
                        
                        mpd_beta, mpd_sd_beta,
                        pse_beta, pse_sd_beta,
                        mntd_beta, mntd_sd_beta,
                        
                        fs_beta, fs_sd_beta,
                        n_beta, n_sd_beta,
                        p_beta, p_sd_beta,
                        sla_beta, sla_sd_beta,
                        ldmc_beta, ldmc_sd_beta
                        )

outt[[i]]<-rbind.data.frame(study)

}

jjj<-do.call(rbind,outt)  
jjj<-data.frame(jjj)

write.table(jjj,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Div_ts.csv",sep=",",row.names=F)

