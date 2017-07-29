#################################
# Calculate Loreau Synchrony  ##
#################################

rm(list=ls()) 
require(reshape2)
require(dplyr)
#require(plyr)
require(tidyr)
require(stringr)
require(codyn)

length_unique<-function(x){length(unique(x))}


# Data

all<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/Data/Community/StabilityRedux_spp_plot_biomass_24062016.csv",sep=",",header=T)
all<-all[order(all$Study,all$Site,all$BlockUnique,all$PlotUnique,all$SppN,all$SppCode6),]
all$UniqueID<-paste(all$Study,"_",all$Site,"_",all$Block,"_",all$Plot)
all$UniqueID<-str_replace_all(all$UniqueID," ","")

##############
# set up data#
##############


# quality check, species numbers per site

spp_site<-data.frame(tapply(all$SppCode6,all$Site,length_unique))
spp_site$Site<-rownames(spp_site)
colnames(spp_site)[1]<-"Spp_Num"

## plots for which asynchrony should be calculated

stab<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stabb<-select(stab, Site,UniqueID, Plot_Asynchrony)
stabb<-filter(stabb, is.na(Plot_Asynchrony)==FALSE)
stabb<-select(stabb,-Plot_Asynchrony)

all_sync<-select(all, Study, Site, UniqueID, SppN, SppCode6,Year, Spp_Biomass)

all_syncc<-merge(stabb, all_sync, by.y=c("Site", "UniqueID"))

all_syncc$Year<-as.numeric(all_syncc$Year)

all_syncc$Site<-droplevels(all_syncc$Site)

all_syncc<-arrange(all_syncc,Site, UniqueID, Year)

all_syncc_NoLanta<-filter(all_syncc,Site!="Lanta_drt")

### fix Lanta_drt

Lantaa<-filter(all_syncc, Site=="Lanta_drt")

des<-unique(select(Lantaa, Study, Site, UniqueID, SppN,SppCode6))

dess<-des[rep(seq_len(nrow(des)), 3), ]

dess$Year<-seq(from=1, to=3, by=)
dess$Year<-rep(1:3, each =384)

comb<-merge(dess, Lantaa, by.y=c("Study","Site","UniqueID","SppN","SppCode6","Year"),all.x=TRUE)

comb$Spp_Biomass<-ifelse(is.na(comb$Spp_Biomass)==TRUE,0,comb$Spp_Biomass)

# add Lanta back

all_syncc_tog<-rbind.data.frame(all_syncc_NoLanta,comb)



####################
# Loreau synchrony #
####################

nn<-length_unique(all_syncc_tog$Site)
outt=c();

for(i in 1:nn){
  
  test=subset(all_syncc_tog, all_syncc_tog$Site==(unique(all_syncc_tog$Site))[i])  


  # output
cat("progress", i, sep=' ','\n')

loreau<-synchrony(df=test, time.var="Year", species.var="SppCode6", abundance.var="Spp_Biomass",replicate.var="UniqueID")
colnames(loreau)[2]<-"Loreau_synchrony"
  
gross<-synchrony(df=test, time.var="Year", species.var="SppCode6", abundance.var="Spp_Biomass",replicate.var="UniqueID",metric="Gross")
colnames(gross)[2]<-"Gross_synchrony"

togg<-merge(loreau,gross, by.y="UniqueID")

outt[[i]]<-rbind.data.frame(togg)
}

synzz<-do.call(rbind.data.frame,outt)


# merge with big database

stab2<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)


stabb_out<-merge(stab2, synzz, by.y=c("UniqueID"),all.x=TRUE)


write.table(stabb_out,"/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_VI.csv",sep=",",row.names=F)


