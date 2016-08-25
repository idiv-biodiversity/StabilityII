###########################
## Resilience #############
###########################

# path coefficients across climate events

rm(list=ls()) 
#####################
require(ggplot2)
require(cowplot)
require(dplyr)

#####################


########################
# Model 1: FDis, eMNTD #
########################

dry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_emntd_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
dry1<-filter(dry,response=="lg2Rsl12")
dry1$PD_met<-"eMNTD"
dry1$FD_met<-"FDisp"
dry1$X<-NULL

exdry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXDRY_emntd_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
exdry1<-filter(exdry,response=="lg2Rsl12")
exdry1$PD_met<-"eMNTD"
exdry1$FD_met<-"FDisp"
exdry1$X<-NULL

wet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_WET_emntd_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
wet1<-filter(wet,response=="lg2Rsl12")
wet1$PD_met<-"eMNTD"
wet1$FD_met<-"FDisp"
wet1$X<-NULL

exwet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXWET_emntd_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
exwet1<-filter(exwet,response=="lg2Rsl12")
exwet1$PD_met<-"eMNTD"
exwet1$FD_met<-"FDisp"
exwet1$X<-NULL

paths<-rbind.data.frame(dry1,exdry1, wet1,exwet1)

paths$lowCI<-paths$estimate-(paths$std.error*1.96)
paths$highCI<-paths$estimate+(paths$std.error*1.96)

paths$predictor<-as.character(paths$predictor)
paths$predictor<-ifelse(paths$predictor=="lg2SppN","Species Richness",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="Ye_PCAcwm4trts","Fast-Slow (CWM)",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="Ye_eMNTD","Phylogenetic Diversity",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="Ye_FDis4","Functional Diversity",paths$predictor)
paths$predictor<-as.factor(paths$predictor)

paths$PathSig<-ifelse(paths$p.value<0.05,1,0)
paths$PathSig<-as.factor(paths$PathSig)

##Figure


p1<-ggplot(data=paths,aes(x=Climate_Bin,y=estimate,group=predictor,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-0.45,.65),breaks=c(-0.4,-.3,-.2,0,-.1,0,.1,.2,.3,0.4,0.5,0.6))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = "Standardized path coefficient (95 % CI)") + ggtitle("Resilience (SPEI12)")+
  
  facet_wrap(~predictor,ncol=2,nrow=2)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8,face="bold"),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())



png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_Rsl_FDisEMNTD.png", 
    type="cairo",
    units="in", 
    width=5 ,
    height=5, 
    pointsize=2, 
    res=200)

p1

dev.off()


########################
# Model 2: FDis, eMPD #
########################


dry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_empd_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
dry1<-filter(dry,response=="lg2Rsl12")
dry1$PD_met<-"eMPD"
dry1$FD_met<-"FDisp"
dry1$X<-NULL

exdry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXDRY_empd_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
exdry1<-filter(exdry,response=="lg2Rsl12")
exdry1$PD_met<-"eMPD"
exdry1$FD_met<-"FDisp"
exdry1$X<-NULL

wet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_WET_empd_fdis_sem_coefs_V1_NEW.csv",sep=",",header=TRUE)
wet1<-filter(wet,response=="lg2Rsl12")

# no models fit data with FDis included as exogenous variable

wet2<-cbind.data.frame("lg2Rsl12","Ye_FDis4",0,0,1,"","Wet")
colnames(wet2)[1]<-"response"
colnames(wet2)[2]<-"predictor"
colnames(wet2)[3]<-"estimate"
colnames(wet2)[4]<-"std.error"
colnames(wet2)[5]<-"p.value"
colnames(wet2)[6]<-"X"
colnames(wet2)[7]<-"Climate_Bin"

wet1<-rbind.data.frame(wet1,wet2)

wet1$PD_met<-"eMPD"
wet1$FD_met<-"FDisp"
wet1$X<-NULL

# separate models fit including either eMPD or FDis; path coefficients of both included for visualisation

exwet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXWET_empd_fdis_sem_coefs_V1_NEW.csv",sep=",",header=TRUE)
exwet1<-filter(exwet,response=="lg2Rsl12")
exwet1$PD_met<-"eMPD"
exwet1$FD_met<-"FDisp"
exwet1$X<-NULL


exwet2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXWET_empd_fdis_sem_coefs_V2_NEW.csv",sep=",",header=TRUE)
exwet2<-filter(exwet2,response=="lg2Rsl12")
exwet2$PD_met<-"eMPD"
exwet2$FD_met<-"FDisp"
exwet2$X<-NULL

exwet22<-exwet2[1,]

exwet11<-rbind.data.frame(exwet1,exwet22)


paths2<-rbind.data.frame(dry1,exdry1, wet1,exwet11)

paths2$lowCI<-paths2$estimate-(paths2$std.error*1.96)
paths2$highCI<-paths2$estimate+(paths2$std.error*1.96)

paths2$predictor<-as.character(paths2$predictor)
paths2$predictor<-ifelse(paths2$predictor=="lg2SppN","Species Richness",paths2$predictor)
paths2$predictor<-ifelse(paths2$predictor=="Ye_PCAcwm4trts","Fast-Slow (CWM)",paths2$predictor)
paths2$predictor<-ifelse(paths2$predictor=="Ye_eMPD","Phylogenetic Diversity",paths2$predictor)
paths2$predictor<-ifelse(paths2$predictor=="Ye_FDis4","Functional Diversity",paths2$predictor)
paths2$predictor<-as.factor(paths2$predictor)

paths2$PathSig<-ifelse(paths2$p.value<0.05,1,0)
paths2$PathSig<-as.factor(paths2$PathSig)

##Figure


p2<-ggplot(data=paths2,aes(x=Climate_Bin,y=estimate,group=predictor,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-0.45,.68),breaks=c(-0.4,-.3,-.2,0,-.1,0,.1,.2,.3,0.4,0.5,0.6))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = "Standardized path coefficient (95 % CI)") + ggtitle("Resilience (SPEI12)")+
  
  facet_wrap(~predictor,ncol=2,nrow=2)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8,face="bold"),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())

png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_Rsl_FDisEMPD.png", 
    type="cairo",
    units="in", 
    width=5 ,
    height=5, 
    pointsize=2, 
    res=200)

p2

dev.off()



############
## Model 3 #
############

dry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_epse_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
dry1<-filter(dry,response=="lg2Rsl12")
dry1$PD_met<-"ePSE"
dry1$FD_met<-"FDisp"
dry1$X<-NULL

exdry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXDRY_epse_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
exdry1<-filter(exdry,response=="lg2Rsl12")
exdry1$PD_met<-"ePSE"
exdry1$FD_met<-"FDisp"
exdry1$X<-NULL

wet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_WET_epse_fdis_sem_coefs_V1_NEW.csv",sep=",",header=TRUE)
wet1<-filter(wet,response=="lg2Rsl12")
wet1$PD_met<-"ePSE"
wet1$FD_met<-"FDisp"
wet1$X<-NULL

wet2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_WET_epse_fdis_sem_coefs_V2_NEW.csv",sep=",",header=TRUE)
wet2<-filter(wet2,response=="lg2Rsl12")
wet2$PD_met<-"ePSE"
wet2$FD_met<-"FDisp"
wet2$X<-NULL

wet22<-wet2[2,]

wet11<-rbind.data.frame(wet1,wet22)


exwet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXWET_epse_fdis_sem_coefs_V1_NEW.csv",sep=",",header=TRUE)
eexwet1<-filter(exwet,response=="lg2Rsl12")
exwet1$PD_met<-"ePSE"
exwet1$FD_met<-"FDisp"
exwet1$X<-NULL

exwet2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXWET_epse_fdis_sem_coefs_V2_NEW.csv",sep=",",header=TRUE)
exwet2<-filter(exwet2,response=="lg2Rsl12")
exwet2$PD_met<-"ePSE"
exwet2$FD_met<-"FDisp"
exwet2$X<-NULL

exwet22<-exwet2[1,]

exwet11<-rbind.data.frame(exwet1,exwet22)

paths3<-rbind.data.frame(dry1,exdry1, wet11,exwet11)

paths3$lowCI<-paths3$estimate-(paths3$std.error*1.96)
paths3$highCI<-paths3$estimate+(paths3$std.error*1.96)

paths3$predictor<-as.character(paths3$predictor)
paths3$predictor<-ifelse(paths3$predictor=="lg2SppN","Species Richness",paths3$predictor)
paths3$predictor<-ifelse(paths3$predictor=="Ye_PCAcwm4trts","Fast-Slow (CWM)",paths3$predictor)
paths3$predictor<-ifelse(paths3$predictor=="Ye_ePSE","Phylogenetic Diversity",paths3$predictor)
paths3$predictor<-ifelse(paths3$predictor=="Ye_eMPD","Phylogenetic Diversity",paths3$predictor)
paths3$predictor<-ifelse(paths3$predictor=="Ye_FDis4","Functional Diversity",paths3$predictor)
paths3$predictor<-as.factor(paths3$predictor)

paths3$PathSig<-ifelse(paths3$p.value<0.05,1,0)
paths3$PathSig<-as.factor(paths3$PathSig)

##Figure

p3<-ggplot(data=paths3,aes(x=Climate_Bin,y=estimate,group=predictor,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-0.45,.65),breaks=c(-0.4,-.3,-.2,0,-.1,0,.1,.2,.3,0.4,0.5,0.6))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = "Standardized path coefficient (95 % CI)") + ggtitle("Resilience (SPEI12)")+
  
  facet_wrap(~predictor,ncol=2,nrow=2)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8,face="bold"),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())

png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_Rsl_FDisEPSE.png", 
    type="cairo",
    units="in", 
    width=5 ,
    height=5, 
    pointsize=2, 
    res=200)

p3

dev.off()

##############
# Model 4  ###
##############


dry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_epse_fric_sem_coefs.csv",sep=",",header=TRUE)
dry1<-filter(dry,response=="lg2Rsl12")
dry1$PD_met<-"ePSE"
dry1$FD_met<-"FRic"
dry1$X<-NULL

exdry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXDRY_epse_fric_sem_coefs.csv",sep=",",header=TRUE)
exdry1<-filter(exdry,response=="lg2Rsl12")
exdry1$PD_met<-"ePSE"
exdry1$FD_met<-"FRic"
exdry1$X<-NULL


# include path coefficients for FRic and ePSE from different models for visualisation 

wet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_WET_epse_fric_sem_coefs_V1.csv",sep=",",header=TRUE)
wet1<-filter(wet,response=="lg2Rsl12")

wet2<-cbind.data.frame("lg2Rsl12","Ye_FRic4",0,0,1,"","Wet")
colnames(wet2)[1]<-"response"
colnames(wet2)[2]<-"predictor"
colnames(wet2)[3]<-"estimate"
colnames(wet2)[4]<-"std.error"
colnames(wet2)[5]<-"p.value"
colnames(wet2)[6]<-"X"
colnames(wet2)[7]<-"Climate_Bin"

wet1<-rbind.data.frame(wet1,wet2)
wet1$PD_met<-"ePSE"
wet1$FD_met<-"FRic"
wet1$X<-NULL



exwet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXWET_epse_fric_sem_coefs.csv",sep=",",header=TRUE)
exwet1<-filter(exwet,response=="lg2Rsl12")
exwet1$PD_met<-"ePSE"
exwet1$FD_met<-"FRic"
exwet1$X<-NULL

paths4<-rbind.data.frame(dry1,exdry1, wet1,exwet1)

paths4$lowCI<-paths4$estimate-(paths4$std.error*1.96)
paths4$highCI<-paths4$estimate+(paths4$std.error*1.96)

paths4$predictor<-as.character(paths4$predictor)
paths4$predictor<-ifelse(paths4$predictor=="lg2SppN","Species Richness",paths4$predictor)
paths4$predictor<-ifelse(paths4$predictor=="Ye_PCAcwm4trts","Fast-Slow (CWM)",paths4$predictor)
paths4$predictor<-ifelse(paths4$predictor=="Ye_ePSE","Phylogenetic Diversity",paths4$predictor)
paths4$predictor<-ifelse(paths4$predictor=="Ye_FRic4","Functional Diversity",paths4$predictor)
paths4$predictor<-ifelse(paths4$predictor=="Ye_FDis4","Functional Diversity",paths4$predictor)

paths4$predictor<-as.factor(paths4$predictor)

paths4$PathSig<-ifelse(paths4$p.value<0.05,1,0)
paths4$PathSig<-as.factor(paths4$PathSig)

##Figure

p4<-ggplot(data=paths4,aes(x=Climate_Bin,y=estimate,group=predictor,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-0.45,.65),breaks=c(-0.4,-.3,-.2,0,-.1,0,.1,.2,.3,0.4,0.5,0.6))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = "Standardized path coefficient (95 % CI)") + ggtitle("Resilience (SPEI12)")+
  
  facet_wrap(~predictor,ncol=2,nrow=2)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8,face="bold"),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())

png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_Rsl_FRicEPSE.png", 
    type="cairo",
    units="in", 
    width=5 ,
    height=5, 
    pointsize=2, 
    res=200)

p4

dev.off()


###################
### Model 5 #######
###################


dry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_emntd_fric_sem_coefs.csv",sep=",",header=TRUE)
dry1<-filter(dry,response=="lg2Rsl12")
dry1$PD_met<-"eMNTD"
dry1$FD_met<-"FRic"
dry1$X<-NULL

exdry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXDRY_emntd_fric_sem_coefs.csv",sep=",",header=TRUE)
exdry1<-filter(exdry,response=="lg2Rsl12")
exdry1$PD_met<-"eMNTD"
exdry1$FD_met<-"FRic"
exdry1$X<-NULL

#  use path coefficients from separate models for visualization purposes
wet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_WET_emntd_fric_sem_coefs.csv",sep=",",header=TRUE)
wet1<-filter(wet,response=="lg2Rsl12")
wet1$PD_met<-"eMNTD"
wet1$FD_met<-"FRic"
wet1$X<-NULL


exwet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_EXWET_emntd_fric_sem_coefs.csv",sep=",",header=TRUE)
exwet1<-filter(exwet,response=="lg2Rsl12")
exwet1$PD_met<-"eMNTD"
exwet1$FD_met<-"FRic"
exwet1$X<-NULL

paths5<-rbind.data.frame(dry1,exdry1, wet1,exwet1)

paths5$lowCI<-paths5$estimate-(paths5$std.error*1.96)
paths5$highCI<-paths5$estimate+(paths5$std.error*1.96)

paths5$predictor<-as.character(paths5$predictor)
paths5$predictor<-ifelse(paths5$predictor=="lg2SppN","Species Richness",paths5$predictor)
paths5$predictor<-ifelse(paths5$predictor=="Ye_PCAcwm4trts","Fast-Slow (CWM)",paths5$predictor)
paths5$predictor<-ifelse(paths5$predictor=="Ye_eMNTD","Phylogenetic Diversity",paths5$predictor)
paths5$predictor<-ifelse(paths5$predictor=="Ye_FRic4","Functional Diversity",paths5$predictor)
#paths5$predictor<-ifelse(paths5$predictor=="Yn_FDis4","Functional Diversity",paths5$predictor)

paths5$predictor<-as.factor(paths5$predictor)

paths5$PathSig<-ifelse(paths5$p.value<0.05,1,0)
paths5$PathSig<-as.factor(paths5$PathSig)

##Figure

p5<-ggplot(data=paths5,aes(x=Climate_Bin,y=estimate,group=predictor,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-0.45,.6),breaks=c(-0.4,-.3,-.2,0,-.1,0,.1,.2,.3,0.4,0.5,0.6))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = "Standardized path coefficient (95 % CI)") + ggtitle("Resilience (SPEI12)")+
  
  facet_wrap(~predictor,ncol=2,nrow=2)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8,face="bold"),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())

png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_Rsl_FRicEMNTD.png", 
    type="cairo",
    units="in", 
    width=5 ,
    height=5, 
    pointsize=2, 
    res=200)

p5

dev.off()


##################
# Model 6  #######
##################

[Stop here]

dry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_DRY_empd_fric_sem_coefs.csv",sep=",",header=TRUE)
dry1<-filter(dry,response=="lg2Rst12")
dry1$PD_met<-"eMPD"
dry1$FD_met<-"FRic"
dry1$X<-NULL

exdry<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_EXDRY_empd_fric_sem_coefs.csv",sep=",",header=TRUE)
exdry1<-filter(exdry,response=="lg2Rst12")
exdry1$PD_met<-"eMPD"
exdry1$FD_met<-"FRic"
exdry1$X<-NULL



#  use path coefficients from separate models for visualization purposes
wet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_empd_fric_sem_coefs_V1.csv",sep=",",header=TRUE)
wet1<-filter(wet,response=="lg2Rst12")
wet1$PD_met<-"eMPD"
wet1$FD_met<-"FRic"
wet1$X<-NULL


wet2<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_empd_fric_sem_coefs_V2.csv",sep=",",header=TRUE)
wet22<-filter(wet2,response=="lg2Rst12")
wet22$PD_met<-"eMPD"
wet22$FD_met<-"FRic"
wet22$X<-NULL

wet22<-wet22[3,]

wet11<-rbind.data.frame(wet1,wet22)

exwet<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_EXWET_empd_fric_sem_coefs.csv",sep=",",header=TRUE)
exwet1<-filter(exwet,response=="lg2Rst12")
exwet1$PD_met<-"eMPD"
exwet1$FD_met<-"FRic"
exwet1$X<-NULL

paths6<-rbind.data.frame(dry1,exdry1, wet11,exwet1)

paths6$lowCI<-paths6$estimate-(paths6$std.error*1.96)
paths6$highCI<-paths6$estimate+(paths6$std.error*1.96)

paths6$predictor<-as.character(paths6$predictor)
paths6$predictor<-ifelse(paths6$predictor=="lg2SppN","Species Richness",paths6$predictor)
paths6$predictor<-ifelse(paths6$predictor=="Yn_PCAcwm4trts","Fast-Slow (CWM)",paths6$predictor)
paths6$predictor<-ifelse(paths6$predictor=="Yn_eMPD","Phylogenetic Diversity",paths6$predictor)
paths6$predictor<-ifelse(paths6$predictor=="Yn_FRic4","Functional Diversity",paths6$predictor)
#paths5$predictor<-ifelse(paths5$predictor=="Yn_FDis4","Functional Diversity",paths5$predictor)

paths6$predictor<-as.factor(paths6$predictor)

paths6$PathSig<-ifelse(paths6$p.value<0.05,1,0)
paths6$PathSig<-as.factor(paths6$PathSig)

##Figure

p6<-ggplot(data=paths5,aes(x=Climate_Bin,y=estimate,group=predictor,colour=PathSig))+
  geom_hline(yintercept=0,linetype=3,color="gray40")+
  geom_point(size=2,shape=19)+
  geom_errorbar(width=0.2,aes(ymin=lowCI,ymax=highCI))+
  scale_y_continuous(limits=c(-0.45,.5),breaks=c(-0.4,-.3,-.2,0,-.1,0,.1,.2,.3,0.4))+
  
  scale_colour_manual(name="",values = c("1"= "#de2d26","0" ="#3182bd"))+
  
  labs(x = "", y = "Standardized path coefficient (95 % CI)") + ggtitle("Resistance (SPEI12)")+
  
  facet_wrap(~predictor,ncol=2,nrow=2)+
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8,face="bold"),
        strip.text.x=element_text(size=6,face="bold"),
        legend.position="none",
        panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank())

png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/PathCoefficients_Rst_FRicEMPD.png", 
    type="cairo",
    units="in", 
    width=5 ,
    height=5, 
    pointsize=2, 
    res=200)

p6

dev.off()


