require(dplyr)


sla<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/Data/Traits/StabII_SLA_spplevel.csv",sep=",",header=T)
sla$Trait<-"SLA"
sla$Value<-sla$SLA

sla<-select(sla,-SLA)

leafn<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/Data/Traits/StabII_LeafN_spplevel.csv",sep=",",header=T)
leafn$Trait<-"LeafN"
leafn$Value<-leafn$LeafN
leafn<-select(leafn,-LeafN)

ldmc<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/Data/Traits/StabII_LDMC_spplevel.csv",sep=",",header=T)
ldmc$Trait<-"LDMC"
ldmc$Value<-ldmc$LDMC
ldmc<-select(ldmc,-LDMC)

leafp<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/Data/Traits/StabII_leafP_spplevel.csv",sep=",",header=T)
leafp$Trait<-"LeafP"
leafp$Value<-leafp$LeafP
leafp<-select(leafp,-LeafP)

traits_spp<-rbind.data.frame(sla,leafn,ldmc,leafp)

spp_site<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/Data/Traits/Site_Study_SppCode6.txt",sep="\t",header=T)
spp_site<-unique(select(spp_site,Site,SppCode6))

stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab<-unique(select(stab,Site))

spp_site2<-merge(stab,spp_site,by.y="Site",all.x=TRUE)
spp_site2<-unique(select(spp_site2,SppCode6))


traits_spp<-merge(spp_site2,traits_spp,,by.y="SppCode6",all.x=TRUE)
traits_spp<-na.omit(traits_spp)


#filter out individual studies
traits_f<-filter(traits_spp,Source!="MEND")
traits_f<-filter(traits_f,Source!="Catford_CC")
traits_f<-filter(traits_f,Source!="Grime_2007")
traits_f<-filter(traits_f,Source!="Wacker")
traits_f<-filter(traits_f,Source!="Jena_2006")
traits_f<-filter(traits_f,Source!="Jena_Mono_2011")
traits_f<-filter(traits_f,Source!="Cedar_Creek")

traits_ff<-unique(select(traits_f,AccSpeciesID,Source))

write.table(traits_ff,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Trait_datasources.csv",sep=",",row.names=F)
