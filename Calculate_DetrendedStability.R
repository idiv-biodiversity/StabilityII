############################################
# Calculate detrended ecosystem stability  #
# following Tilman et al. 2006             #
############################################
# 1. for each plot, fit biomass ~ year 
# 2. extract model residuals
# 4. calculate detrended stability as mean biomass/ detrended SD
# 3. calculate detrended SD as the sd of model residuals

require(reshape2)
require(dplyr)
require(tidyr)
require(stringr)
require(mgcv)

len_un<-function(x){length(unique(x))}

# prepare data
all<-read.delim("data.csv",sep=",",header=T)
all<-all[order(all$Study,all$Site,all$BlockUnique,all$PlotUnique,all$SppN,all$SppCode6),]
all$UniqueID<-paste(all$Study,"_",all$Site,"_",all$Block,"_",all$Plot)
all$UniqueID<-str_replace_all(all$UniqueID," ","")

all<-summarize(group_by(all, Study, Site, BlockUnique, PlotUnique, UniqueID, SppN, Year), Plot_Biomass=sum(Spp_Biomass))
all<-ungroup(all)

####################
# detrended SD     #
####################

un_un<-unique(all$UniqueID) # vector of unique plot ids for loop
IDn<-seq(from=1, to=length(un_un))

unn<-cbind.data.frame(un_un,IDn)
colnames(unn)[1]<-"UniqueID"

all<-left_join(all, unn,by="UniqueID") # add unique plot id to year level data 

nn<-len_un(all$IDn)
outt=c();

for(i in 1:nn){
  
  test<-subset(all, all$IDn==(unique(all$IDn))[i])  
  cat("progress", i, sep=' ','\n')
  
  modd<-lm(Plot_Biomass~Year, data=test) 
   
  test$resid<-residuals(modd, type="response")
  
  detrend_sd<-sd(test$resid)
  
  UniqueID<-unique(test$UniqueID)
  
  outt[[i]]<-cbind.data.frame(UniqueID, detrend_sd)
}

synzz<-do.call(rbind.data.frame,outt)

# merge with big data

stab<-read.delim("data_II.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stabb_out<-left_join(stab, synzz, by=c("UniqueID"))

# calculate de-treneded stability

stabb_out$Plot_TempStab_Detrend<-stabb_out$Plot_Biomassxbar/stabb_out$detrend_sd

write.table(stabb_out,"data_II.csv",sep=",",row.names=F)
