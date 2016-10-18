rm(list=ls()) 
#####################
require(ggplot2)
require(cowplot)
require(dplyr)
require(tidyr)

sync_p<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_sync_paths_Oct2016.csv",sep=",",header=TRUE)

# predictors of synchrony

sync1<-filter(sync_p,response=="Plot_Asynchrony")
sync1<-sync1 %>%
  separate(ModClass, c("FDclass", "PDclass"), "_")
sync1$FDclass2<-"FD"

sync1$predictor<-as.character(sync1$predictor)
sync1$predictor<-ifelse(sync1$predictor=="PCAdim1_4trts","Fast-Slow Axis",sync1$predictor)
sync1$predictor<-ifelse(sync1$predictor=="lg2SppN","Species Richness",sync1$predictor)
sync1$predictor<-ifelse(sync1$predictor=="FDis4" |sync1$predictor=="FRic4","FD",sync1$predictor)
sync1$predictor<-ifelse(sync1$predictor=="eMNTD" |sync1$predictor=="ePSE"|sync1$predictor=="eMPD","PD",sync1$predictor)

sync1$predictor<-as.factor(sync1$predictor)

sync1$response<-"Synchrony"

sync1$Model<-paste(sync1$FDclass,sync1$PDclass,sep="_")

sync1$lowCI<-sync1$estimate-(sync1$std.error*1.96)
sync1$highCI<-sync1$estimate+(sync1$std.error*1.96)

sync1$PathSig<-ifelse(sync1$p.value<0.05,1,0)
sync1$PathSig<-as.factor(sync1$PathSig)




##SppN

sppn<-filter(sync1,predictor=="Species Richness")

p1<-ggplot(data=sppn,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = expression(eta)) + ggtitle("Plant species richness")+
  
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

fss<-filter(sync1,predictor=="Fast-Slow Axis")

p2<-ggplot(data=fss,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = expression(eta)) + ggtitle("Fast-slow axis")+
  
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


fdd<-filter(sync1,predictor=="FD")


p3<-ggplot(data=fdd,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = expression(eta)) + ggtitle("Functional diversity")+
  
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


pdd<-filter(sync1,predictor=="PD")


p4<-ggplot(data=pdd,aes(x=PDclass,y=estimate,group=FDclass,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-1,1))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = expression(eta)) + ggtitle("Phylogenetic diversity")+
  
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

#########

png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_Sync_NEW.png", 
    type="cairo",
    units="in", 
    width=8, 
    height=8, 
    pointsize=2, 
    res=200)

plot_grid(p1,p2,p3,p4,ncol=2,nrow=2,labels=c("(a)","(b)","(c)","(d)"),label_size=7,align="hv")


dev.off()

#########

# predictors of TS

sync2<-filter(sync_p,response=="TS_lg2")
sync2<-sync2 %>%
  separate(ModClass, c("FDclass", "PDclass"), "_")
sync2$FDclass2<-"FD"

sync2$predictor<-as.character(sync2$predictor)
sync2$predictor<-ifelse(sync2$predictor=="PCAdim1_4trts","Fast-Slow Axis",sync2$predictor)
sync2$predictor<-ifelse(sync2$predictor=="lg2SppN","Species Richness",sync2$predictor)
sync2$predictor<-ifelse(sync2$predictor=="Plot_Asynchrony" ,"eta",sync2$predictor)
sync2$predictor<-ifelse(sync2$predictor=="eMNTD" |sync2$predictor=="ePSE"|sync2$predictor=="eMPD","PD",sync2$predictor)

sync2$predictor<-as.factor(sync2$predictor)

sync2$response<-"Stability"

sync2$Model<-paste(sync2$FDclass,sync2$PDclass,sep="_")

sync2$lowCI<-sync2$estimate-(sync2$std.error*1.96)
sync2$highCI<-sync2$estimate+(sync2$std.error*1.96)

sync2$PathSig<-ifelse(sync2$p.value<0.05,1,0)
sync2$PathSig<-as.factor(sync2$PathSig)


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

png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_TS_NEW.png", 
    type="cairo",
    units="in", 
    width=8, 
    height=8, 
    pointsize=2, 
    res=200)

plot_grid(p11,p22,p44,p33,p55,ncol=3,nrow=2,labels=c("(a)","(b)","(c)","(d)","(e)"),label_size=7,align="hv")


dev.off()

