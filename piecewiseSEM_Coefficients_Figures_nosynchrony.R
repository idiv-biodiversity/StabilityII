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

nosync1$predictor<-as.character(nosync1$predictor)
nosync1$predictor<-ifelse(nosync1$predictor=="PCAdim1_4trts","Fast-Slow Axis",nosync1$predictor)
nosync1$predictor<-ifelse(nosync1$predictor=="lg2SppN","Species Richness",nosync1$predictor)
nosync1$predictor<-ifelse(nosync1$predictor=="FDis4" |nosync1$predictor=="FRic4","FD",nosync1$predictor)
nosync1$predictor<-ifelse(nosync1$predictor=="eMNTD" |nosync1$predictor=="ePSE"|nosync1$predictor=="eMPD","PD",nosync1$predictor)

nosync1$predictor<-as.factor(nosync1$predictor)

nosync1$response<-"Synchrony"


nosync1$lowCI<-nosync1$estimate-(nosync1$std.error*1.96)
nosync1$highCI<-nosync1$estimate+(nosync1$std.error*1.96)

nosync1$PathSig<-ifelse(nosync1$p.value<0.05,1,0)
nosync1$PathSig<-as.factor(nosync1$PathSig)


[[stop here]]
########
##SppN #
########

sppn2<-filter(sync2,predictor=="Species Richness")

p11<-ggplot(data=sppn2,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = expression(bold(paste("Ecosystem stability ( ",mu," / ",sigma," )")))) + ggtitle("Plant species richness")+
  
  facet_grid(~FDclass)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8,face="bold"),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())

fss2<-filter(sync2,predictor=="Fast-Slow Axis")

p22<-ggplot(data=fss2,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y=expression(bold(paste("Ecosystem stability ( ",mu," / ",sigma," )"))))+ ggtitle("Fast-slow axis")+
  
  facet_grid(~FDclass)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.title.y=element_blank(),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())


fdd2<-filter(sync2,predictor=="CV_Precip")


p33<-ggplot(data=fdd2,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y=expression(bold(paste("Ecosystem stability ( ",mu," / ",sigma," )"))))+ ggtitle("Precipitation seasonality")+
  
  facet_grid(~FDclass)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8,face="bold"),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())


pdd2<-filter(sync2,predictor=="PD")

empty<-filter(sync2,predictor=="Fast-Slow Axis")
empty$estimate<-0
empty$predictor<-"PD"
empty$lowCI<-0
empty$highCI<-0
empty$PathSig<-0
empty<-filter(empty,Model!="FDis_eMNTD")
empty<-filter(empty,Model!="FRic_eMNTD")
pdd2<-rbind.data.frame(pdd2,empty)

p44<-ggplot(data=pdd2,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="white"))+
  
  labs(x = "", y=expression(bold(paste("Ecosystem stability ( ",mu," / ",sigma," )")))) + ggtitle("Phylogenetic diversity")+
  
  facet_grid(~FDclass)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_blank(),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())

fdd5<-filter(sync2,predictor=="eta")


p55<-ggplot(data=fdd5,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y=expression(bold(paste("Ecosystem stability ( ",mu," / ",sigma," )"))))+ ggtitle(expression(eta))+
  
  facet_grid(~FDclass)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_blank(),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())


#########

png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_TS_nosync_NEW.png", 
    type="cairo",
    units="in", 
    width=8, 
    height=8, 
    pointsize=2, 
    res=200)

plot_grid(p11,p22,p44,p33,p55,ncol=3,nrow=2,labels=c("(a)","(b)","(c)","(d)","(e)"),label_size=7,align="hv")


dev.off()

