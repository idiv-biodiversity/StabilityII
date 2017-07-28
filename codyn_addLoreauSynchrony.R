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

# quality check, species numbers per site

spp_site<-data.frame(tapply(all$SppCode6,all$Site,length_unique))
spp_site$Site<-rownames(spp_site)
colnames(spp_site)[1]<-"Spp_Num"


########################
# use codyn ############
########################

all_sync<-select(all, Study, Site, UniqueID, SppN, Year, SppCode6, Spp_Biomass)



[[stop here]]
spp_site2<-unique(select(all_sync, Study, Site, UniqueID, SppCode6))




all_syncc<-merge(spp_site2, all_sync, by.y=c("Study","Site", "UniqueID","SppCode6"),all.x=TRUE)


all_Loreau<-synchrony(df=all_sync, time.var="Year", species.var="SppCode6", abundance.var="Spp_Biomass",replicate.var="UniqueID")


Error in check_multispp(df, species.var, replicate.var) : 
  One or more replicates consists of only a single species;
please remove these replicates prior to calculations 

synch_onerep(df, time.var, species.var, abundance.var, metric = "Loreau")

check_multispp(alls, species.var="SppCode6", replicate.var="UniqueID")

alls2<-filter(alls, Spp_Biomass>0)

Spp1 <-summarize(group_by(alls2, Study, Site, UniqueID,Year), realSppN=length_unique(SppCode6))
Spp1<-filter(Spp1,realSppN>1)

alls3<-merge(Spp1,alls2,by.y=c("Study","Site","UniqueID","Year"))

belg<-filter(alls3, Study=="BigBio")

codyn_Loreau_sync<-synchrony(df=belg, time.var="Year", species.var="SppCode6", abundance.var="Spp_Biomass",replicate.var="UniqueID")

codyn_Gross_sync<-synchrony(df=alls3, time.var="Year", species.var="SppCode6", abundance.var="Spp_Biomass",replicate.var="UniqueID",metric = "Gross")


### or 'synchrony'

install.packages("synchrony")


require(synchrony)

#######################

c<-merge(fun,jjj,by.y=c("UniqueID"),all.x=TRUE)
colnames(c)[10]<-"Plot_Biomassxbar"
colnames(c)[11]<-"Plot_Biomasssd"
colnames(c)[13]<-"Plot_Asynchrony"


write.table(c,"/home/dylan/Dropbox/leipzigPhyTrt/Data/Community/StabRedux_Plot_Level_Stability_24062016.csv",sep=",",row.names=F)



#####################################
## community-level              #####
## asynchrony sensu Gross et al 2013#
# ave corr of spp biomass with  #####
## biomass of all other spp, for ####
## each spp. ########################
#####################################

require(plyr)
alls<-subset(all,SppN>1)
alls$PlotBiomass_minus1Spp<-alls$Plot_Biomass-alls$Spp_Biomass

nn<-length_unique(alls$UniqueID)
outt=c();
for(i in 1:nn){
  
  test=subset(alls, alls$UniqueID==(unique(alls$UniqueID))[i])  
  
  func1 <- function(test)
  {
    return(data.frame(corr=cor(test$Spp_Biomass,test$PlotBiomass_minus1Spp)))
  }
  
  out1<-ddply(test,.(SppCode6),func1)
  
  dd<-subset(test,Spp_Biomass>0)
  realized_spp<- length_unique(dd$SppCode6)
  r_xbar<-mean(out1$corr,na.rm=T)
  UniqueID<-unique(test$UniqueID)
  
  outt[[i]]<-cbind.data.frame(UniqueID,r_xbar,realized_spp)
  
}

jjj<-do.call(rbind,outt)  
jjj <- jjj[!(is.nan(jjj$r_xbar)),]
colnames(jjj)[2]<-"Community_Asynchrony"
colnames(jjj)[3]<-"Realized_Spp"

