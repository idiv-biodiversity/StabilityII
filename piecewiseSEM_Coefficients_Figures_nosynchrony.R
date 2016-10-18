rm(list=ls()) 
#####################
require(ggplot2)
require(cowplot)
require(dplyr)
require(tidyr)

#####################

nosync_p<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_no_sync_paths_Oct2016.csv",sep=",",header=TRUE)
nosync1<-filter(nosync_p,response=="TS_lg2")

nosync1<-nosync1 %>%
  separate(ModClass, c("FDclass", "PDclass"), "_")
nosync1$FDclass2<-"FD"

sync1$predictor<-as.character(sync1$predictor)
sync1$predictor<-ifelse(sync1$predictor=="PCAdim1_4trts","Fast-Slow Axis",sync1$predictor)
sync1$predictor<-ifelse(sync1$predictor=="lg2SppN","Species Richness",sync1$predictor)
sync1$predictor<-ifelse(sync1$predictor=="FDis4" |sync1$predictor=="FRic4","FD",sync1$predictor)
sync1$predictor<-ifelse(sync1$predictor=="eMNTD" |sync1$predictor=="ePSE"|sync1$predictor=="eMPD","PD",sync1$predictor)

sync1$predictor<-as.factor(sync1$predictor)

sync1$response<-"Synchrony"


sync1$lowCI<-sync1$estimate-(sync1$std.error*1.96)
sync1$highCI<-sync1$estimate+(sync1$std.error*1.96)

sync1$PathSig<-ifelse(sync1$p.value<0.05,1,0)
sync1$PathSig<-as.factor(sync1$PathSig)