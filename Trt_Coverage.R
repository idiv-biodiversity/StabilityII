require(reshape2)
require(dplyr)

# Data

traits<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/Data/Traits/StabII_FStraits.csv",sep=",",header=T)


stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-unique(select(stab,Site,Study))

##########


togg<-merge(traits,stab_4,by.y.=c("Study","Site"))

spp_site<-summarise(group_by(togg, Site), SppN=length(unique(SppCode6)))

sla<-subset(togg,Trait=="SLA")
sla<-filter(sla,is.na(Value)==FALSE)
sla_site<-summarise(group_by(sla, Site), SLA_spp=length(unique(SppCode6)))

ldmc<-subset(togg,Trait=="LDMC")
ldmc<-filter(ldmc,is.na(Value)==FALSE)
ldmc_site<-summarise(group_by(ldmc, Site), LDMC_spp=length(unique(SppCode6)))

leafN<-subset(togg,Trait=="LeafN")
leafN<-filter(leafN,is.na(Value)==FALSE)
leafN_site<-summarise(group_by(leafN, Site), leafN_spp=length(unique(SppCode6)))

leafP<-subset(togg,Trait=="LeafP")
leafP<-filter(leafP,is.na(Value)==FALSE)
leafP_site<-summarise(group_by(leafP, Site), leafP_spp=length(unique(SppCode6)))


togg2<-merge(spp_site,sla_site,by.y="Site",all.x=TRUE)
togg22<-merge(togg2,ldmc_site,by.y="Site",all.x=TRUE)
togg222<-merge(togg22,leafN_site,by.y="Site",all.x=TRUE)
togg222<-merge(togg222,leafP_site,by.y="Site",all.x=TRUE)

togg222$SLA_spp<-togg222$SLA_spp/togg222$SppN
togg222$LDMC_spp<-togg222$LDMC_spp/togg222$SppN
togg222$leafN_spp<-togg222$leafN_spp/togg222$SppN
togg222$leafP_spp<-togg222$leafP_spp/togg222$SppN

trt_cov<-summarize(group_by(togg222),avgSLA=mean(SLA_spp),avgLDMC=mean(LDMC_spp),avgLeafN=mean(leafN_spp),avgLeafP=mean(leafP_spp))


