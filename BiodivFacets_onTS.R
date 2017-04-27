rm(list=ls()) 

library(lmerTest)
library(nlme)
library(dplyr)
library(AICcmodavg)
library(MuMIn)
library(ggplot2)

ICClme <- function(out) {
  varests <- as.numeric(VarCorr(out)[1:2])
  return(paste("ICC =", varests[1]/sum(varests)))
}


# Data
stab<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,eMPD,eMNTD,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

stab_4$PlotAsynchrony_s<-stab_4$Plot_Asynchrony*-1

stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony

cc<-lmeControl(opt="optim")

###########################################
# Effects on TS of biodiversity facets   ##
###########################################


###########
## SR  ####
###########

Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~lg2SppN,random=~1|Site,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~lg2SppN,random=~1|Site/SppN,control=cc,data=stab_444)

Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot(a)
qqnorm(a)


big<-lme(TS_lg2~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=~1+lg2SppN|Site,control=cc,data=stab_444,method="ML")

anova(big,small)

#final
final<-lme(TS_lg2~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)

r.squaredGLMM(final)

intervals(final, which="fixed")

VarCorr(final)

ICClme(final)

###############
# predictions #
###############
spp_ts<-select(stab_444,Site,UniqueID,SppN,TS_lg2,lg2SppN)

spp_ts$pred<-predict(final,spp_ts,re.form=~(~1+lg2SppN|Site))

spp_ts$TS<-2^(spp_ts$TS_lg2)
spp_ts$pred_t<-2^(spp_ts$pred)

newdat <- expand.grid(lg2SppN=seq(from=0,to=6,by=0.1))
newdat$pred <- predict(final, newdat, level = 0)

Designmat <- model.matrix(formula(final)[-2], newdat)
predvar <- diag(Designmat %*% vcov(final) %*% t(Designmat)) 
newdat$SE <- sqrt(predvar) 


newdat$p_lCI<-newdat$pred-(newdat$SE*1.96)
newdat$p_uCI<-newdat$pred+(newdat$SE*1.96)

newdat$pred_t<-2^(newdat$pred)
newdat$SppN<-2^(newdat$lg2SppN)
newdat$pred_uCIt<-2^(newdat$p_uCI)
newdat$pred_lCIt<-2^(newdat$p_lCI)


c<-ggplot(data=spp_ts,aes(x=SppN,y=pred_t))+

    geom_smooth(data=spp_ts,aes(y=pred_t,x=SppN,group=Site),method="lm",formula=y~x,size=0.5,color="gray80",se=FALSE)+
  geom_smooth(data=newdat,aes(y=pred_t,x=SppN),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+


  labs(x="Plant species richness",y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous(trans="log2",breaks=c(1,2,4,8,16,32,60)) +   scale_y_continuous(trans="log2",lim=c(0.8,16),breaks=c(1,2,4,8,16))

SppN<-c+ theme(axis.title.x=element_text(colour="black",face="bold",size=10),
                axis.title.y=element_text(colour="black",face="bold",size=10,vjust=1),
                axis.text.y=element_text(colour="black",face="bold",size=8),
                axis.text.x=element_text(colour="black",face="bold",size=8),
                plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
               panel.background = element_rect(fill = "white"))


###################
#Plot_Asynchrony  #
###################

Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~PlotAsynchrony_s,random=~1+PlotAsynchrony_s|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~PlotAsynchrony_s,random=~1+PlotAsynchrony_s|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~PlotAsynchrony_s,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~PlotAsynchrony_s,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~PlotAsynchrony_s,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~PlotAsynchrony_s,random=~1+lg2SppN*PlotAsynchrony_s|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~PlotAsynchrony_s,random=list(~1+lg2SppN+PlotAsynchrony_s|Site),control=cc,data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot(Cand.set[[6]])
qqnorm(Cand.set[[6]])

### LRT

big<-lme(TS_lg2~PlotAsynchrony_s,random=list(~1+lg2SppN*PlotAsynchrony_s|Site),control=cc,data=stab_444,method="ML")


small<-lme(TS_lg2~1,random=list(~1+lg2SppN*PlotAsynchrony_s|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~PlotAsynchrony_s,random=list(~1+lg2SppN*PlotAsynchrony_s|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

VarCorr(final)

ICClme(final)

#predictions

sync_ts<-select(stab_444,Site,UniqueID,PlotAsynchrony_s,TS_lg2,lg2SppN)

sync_ts$pred<-predict(final,sync_ts,re.form=~(~1+lg2SppN*PlotAsynchrony_s|Site))

sync_ts$TS<-2^(sync_ts$TS_lg2)
sync_ts$pred_t<-2^(sync_ts$pred)

newdat <- expand.grid(PlotAsynchrony_s=seq(from=-1,to=1,by=0.01))
newdat$pred <- predict(final, newdat, level = 0)

Designmat <- model.matrix(formula(final)[-2], newdat)
predvar <- diag(Designmat %*% vcov(final) %*% t(Designmat)) 
newdat$SE <- sqrt(predvar) 


newdat$p_lCI<-newdat$pred-(newdat$SE*1.96)
newdat$p_uCI<-newdat$pred+(newdat$SE*1.96)

newdat$pred_t<-2^(newdat$pred)
newdat$pred_uCIt<-2^(newdat$p_uCI)
newdat$pred_lCIt<-2^(newdat$p_lCI)


c<-ggplot(data=sync_ts,aes(x=PlotAsynchrony_s,y=pred_t))+
  
  geom_smooth(data=sync_ts,aes(y=pred_t,x=PlotAsynchrony_s,group=Site),method="lm",formula=y~x,size=0.5,color="gray80",se=FALSE)+
  geom_smooth(data=newdat,aes(y=pred_t,x=PlotAsynchrony_s),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+
  
  
  labs(x=expression(bold(paste("Species asynchrony (",-eta," )"))),y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous() + scale_y_continuous(trans="log2",lim=c(0.8,16),breaks=c(1,2,4,8,16))


Sync<-c+ theme(axis.title.x=element_text(colour="black",face="bold",size=10),
               axis.title.y=element_blank(),
               axis.text.y=element_text(colour="black",face="bold",size=8),
               axis.text.x=element_text(colour="black",face="bold",size=8),
               plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
               panel.background = element_rect(fill = "white"))

##############################
# Div and Sync on Stability  #
##############################

require(cowplot)

png(filename="/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/StabilityII/Div_Sync_onTS_April2017.png", 
    type="cairo",
    units="in", 
    width=7, 
    height=4, 
    pointsize=2, 
    res=200)


cairo_ps(filename="/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/StabilityII/Div_Sync_onTSDiv_Sync_onTS_April2017.eps",
         family="sans",
         height=4,width=7,
         bg="white")


plot_grid(SppN, Sync, labels=c("(a)","(b)"),label_size = 8,align ="hv")

dev.off()


##########
# PD #####
##########

#eMNTD


Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~eMNTD,random=~1+eMNTD|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~eMNTD,random=~1+eMNTD|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~eMNTD,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~eMNTD,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~eMNTD,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~eMNTD,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~eMNTD,random=list(~1+lg2SppN+eMNTD|Site),control=cc,data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table


plot(Cand.set[[6]])
qqnorm(Cand.set[[6]])



### LRT

big<-lme(TS_lg2~eMNTD*CV_Precip,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444,method="ML")

small<-lme(TS_lg2~1,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~eMNTD,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444)
r.squaredGLMM(final)

ICClme(final)
VarCorr(final)

################
# predictions  #
################

emntd_ts<-select(stab_444,Site,UniqueID,lg2SppN,TS_lg2,eMNTD)

emntd_ts$pred<-predict(final,emntd_ts,re.form=~(~1+lg2SppN*eMNTD|Site))

emntd_ts$TS<-2^(emntd_ts$TS_lg2)
emntd_ts$pred_t<-2^(emntd_ts$pred)

newdat <- expand.grid(eMNTD=seq(from=0,to=20,by=0.1))
newdat$pred <- predict(final, newdat, level = 0)

Designmat <- model.matrix(formula(final)[-2], newdat)
predvar <- diag(Designmat %*% vcov(final) %*% t(Designmat)) 
newdat$SE <- sqrt(predvar) 


newdat$p_lCI<-newdat$pred-(newdat$SE*1.96)
newdat$p_uCI<-newdat$pred+(newdat$SE*1.96)

newdat$pred_t<-2^(newdat$pred)
newdat$pred_uCIt<-2^(newdat$p_uCI)
newdat$pred_lCIt<-2^(newdat$p_lCI)


d<-ggplot(data=emntd_ts,aes(x=eMNTD,y=pred_t))+
  
  geom_smooth(data=emntd_ts,aes(y=pred_t,x=eMNTD,group=Site),method="lm",formula=y~x,size=0.5,color="gray80",se=FALSE)+
  geom_smooth(data=newdat,aes(y=pred_t,x=eMNTD),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+
  
  
  labs(x="Phylogenetic diversity (MNTD)",y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous() + scale_y_continuous(trans="log2",lim=c(0.8,16),breaks=c(1,2,4,8,16))

emNTD<-d+ theme(axis.title.x=element_text(colour="black",face="bold",size=10),
               axis.title.y=element_text(colour="black",face="bold",size=10,vjust=1),
               axis.text.y=element_text(colour="black",face="bold",size=8),
               axis.text.x=element_text(colour="black",face="bold",size=8),
               plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
               panel.background = element_rect(fill = "white"))



#eMPD


Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~eMPD,random=~1+eMPD|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~eMPD,random=~1+eMPD|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~eMPD,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~eMPD,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~eMPD,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~eMPD,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~eMPD,random=list(~1+lg2SppN+eMPD|Site),control=cc,data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table


plot(Cand.set[[7]])
qqnorm(Cand.set[[7]])

### LRT

big<-lme(TS_lg2~eMPD,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~eMPD,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444)
r.squaredGLMM(final)



################
# predictions  #
################

empd_ts<-select(stab_444,Site,UniqueID,lg2SppN,TS_lg2,eMPD)

empd_ts$pred<-predict(final,empd_ts,re.form=~(~1+lg2SppN*eMPD|Site))

empd_ts$TS<-2^(empd_ts$TS_lg2)
empd_ts$pred_t<-2^(empd_ts$pred)

newdat <- expand.grid(eMPD=seq(from=0,to=13.2,by=0.1))
newdat$pred <- predict(final, newdat, level = 0)

Designmat <- model.matrix(formula(final)[-2], newdat)
predvar <- diag(Designmat %*% vcov(final) %*% t(Designmat)) 
newdat$SE <- sqrt(predvar) 


newdat$p_lCI<-newdat$pred-(newdat$SE*1.96)
newdat$p_uCI<-newdat$pred+(newdat$SE*1.96)

newdat$pred_t<-2^(newdat$pred)
newdat$pred_uCIt<-2^(newdat$p_uCI)
newdat$pred_lCIt<-2^(newdat$p_lCI)


e<-ggplot(data=empd_ts,aes(x=eMPD,y=pred_t))+
  
  geom_smooth(data=empd_ts,aes(y=pred_t,x=eMPD,group=Site),method="lm",formula=y~x,size=0.5,color="gray80",se=FALSE)+
  geom_smooth(data=newdat,aes(y=pred_t,x=eMPD),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+
  
  
  labs(x="Phylogenetic diversity (eMPD)",y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous() + scale_y_continuous(trans="log2",lim=c(0.8,16),breaks=c(1,2,4,8,16))

eMPD<-e+ theme(axis.title.x=element_text(colour="black",face="bold",size=8),
                axis.title.y=element_text(colour="black",face="bold",size=8,vjust=1),
                axis.text.y=element_text(colour="black",face="bold",size=8),
                axis.text.x=element_text(colour="black",face="bold",size=8),
                plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
                panel.background = element_rect(fill = "white"))



#sMPD

#sMNTD


#FDis

Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~FDis4,random=~1+FDis4|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~FDis4,random=~1+FDis4|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~FDis4,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~FDis4,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~FDis4,random=~1+FDis4|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~FDis4,random=~1+lg2SppN*FDis4|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~FDis4,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot(Cand.set[[7]])
qqnorm(Cand.set[[7]])

### LRT

big<-lme(TS_lg2~FDis4,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~FDis4,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

ICClme(final)
VarCorr(final)

################
# predictions  #
################


fdis_ts<-select(stab_444,Site,UniqueID,lg2SppN,TS_lg2,FDis4)

fdis_ts$pred<-predict(final,fdis_ts,re.form=~(~1+lg2SppN+FDis4|Site))

fdis_ts$TS<-2^(fdis_ts$TS_lg2)
fdis_ts$pred_t<-2^(fdis_ts$pred)

newdat <- expand.grid(FDis4=seq(from=0,to=2.1,by=0.1))
newdat$pred <- predict(final, newdat, level = 0)

Designmat <- model.matrix(formula(final)[-2], newdat)
predvar <- diag(Designmat %*% vcov(final) %*% t(Designmat)) 
newdat$SE <- sqrt(predvar) 


newdat$p_lCI<-newdat$pred-(newdat$SE*1.96)
newdat$p_uCI<-newdat$pred+(newdat$SE*1.96)

newdat$pred_t<-2^(newdat$pred)
newdat$pred_uCIt<-2^(newdat$p_uCI)
newdat$pred_lCIt<-2^(newdat$p_lCI)


g<-ggplot(data=fdis_ts,aes(x=FDis4,y=pred_t))+
  
  geom_smooth(data=fdis_ts,aes(y=pred_t,x=FDis4,group=Site),method="lm",formula=y~x,size=0.5,color="gray80",se=FALSE)+
  geom_smooth(data=newdat,aes(y=pred_t,x=FDis4),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+
  
  
  labs(x="Functional diversity (FD)",y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous() + scale_y_continuous(trans="log2",lim=c(0.8,16),breaks=c(1,2,4,8,16))

FDis<-g+ theme(axis.title.x=element_text(colour="black",face="bold",size=10),
                axis.title.y=element_blank(),
                axis.text.y=element_text(colour="black",face="bold",size=8),
                axis.text.x=element_text(colour="black",face="bold",size=8),
                plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
                panel.background = element_rect(fill = "white"))


#FRic

Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~FRic4,random=~1+FRic4|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~FRic4,random=~1+FRic4|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~FRic4,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~FRic4,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~FRic4,random=~1+FRic4|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~FRic4,random=~1+lg2SppN*FRic4|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~FRic4,random=list(~1+lg2SppN+FRic4|Site),control=cc,data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot(Cand.set[[6]])
qqnorm(Cand.set[[6]])


### LRT

big<-lme(TS_lg2~FRic4,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~FRic4,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

ICClme(final)

VarCorr(final)

################
# predictions  #
################


fric_ts<-select(stab_444,Site,UniqueID,lg2SppN,TS_lg2,FRic4)

fric_ts$pred<-predict(final,fric_ts,re.form=~(~1+lg2SppN*FRic4|Site))

fric_ts$TS<-2^(fric_ts$TS_lg2)
fric_ts$pred_t<-2^(fric_ts$pred)

newdat <- expand.grid(FRic4=seq(from=0,to=10.3,by=0.1))
newdat$pred <- predict(final, newdat, level = 0)

Designmat <- model.matrix(formula(final)[-2], newdat)
predvar <- diag(Designmat %*% vcov(final) %*% t(Designmat)) 
newdat$SE <- sqrt(predvar) 


newdat$p_lCI<-newdat$pred-(newdat$SE*1.96)
newdat$p_uCI<-newdat$pred+(newdat$SE*1.96)

newdat$pred_t<-2^(newdat$pred)
newdat$pred_uCIt<-2^(newdat$p_uCI)
newdat$pred_lCIt<-2^(newdat$p_lCI)


f<-ggplot(data=fric_ts,aes(x=FRic4,y=pred_t))+
  
  geom_smooth(data=fric_ts,aes(y=pred_t,x=FRic4,group=Site),method="lm",formula=y~x,size=0.5,color="gray80",se=FALSE)+
  geom_smooth(data=newdat,aes(y=pred_t,x=FRic4),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+
  
  
  labs(x="Functional diversity (FRic)",y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous() + scale_y_continuous(trans="log2",lim=c(0.8,16),breaks=c(1,2,4,8,16))

FRic<-f+ theme(axis.title.x=element_text(colour="black",face="bold",size=8),
               axis.title.y=element_text(colour="black",face="bold",size=8,vjust=1),
               axis.text.y=element_text(colour="black",face="bold",size=8),
               axis.text.x=element_text(colour="black",face="bold",size=8),
               plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
               panel.background = element_rect(fill = "white"))


#PCAdim1_4trts

Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~PCAdim1_4trts,random=~1+PCAdim1_4trts|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~PCAdim1_4trts,random=~1+PCAdim1_4trts|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~PCAdim1_4trts,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~PCAdim1_4trts,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~PCAdim1_4trts,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~PCAdim1_4trts,random=~1+lg2SppN*PCAdim1_4trts|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot(Cand.set[[7]])
qqnorm(Cand.set[[7]])


### LRT

big<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)


ICClme(final)
VarCorr(final)
################
# predictions  #
################


fs_ts<-select(stab_444,Site,UniqueID,lg2SppN,TS_lg2,PCAdim1_4trts)

fs_ts$pred<-predict(final,fs_ts,re.form=~(~1+lg2SppN+PCAdim1_4trts|Site))

fs_ts$TS<-2^(fs_ts$TS_lg2)
fs_ts$pred_t<-2^(fs_ts$pred)

newdat <- expand.grid(PCAdim1_4trts=seq(from=-3.62,to=4.47,by=0.1))
newdat$pred <- predict(final, newdat, level = 0)

Designmat <- model.matrix(formula(final)[-2], newdat)
predvar <- diag(Designmat %*% vcov(final) %*% t(Designmat)) 
newdat$SE <- sqrt(predvar) 


newdat$p_lCI<-newdat$pred-(newdat$SE*1.96)
newdat$p_uCI<-newdat$pred+(newdat$SE*1.96)

newdat$pred_t<-2^(newdat$pred)
newdat$pred_uCIt<-2^(newdat$p_uCI)
newdat$pred_lCIt<-2^(newdat$p_lCI)


gg<-ggplot(data=fs_ts,aes(x=PCAdim1_4trts,y=pred_t))+
  
  geom_smooth(data=fs_ts,aes(y=pred_t,x=PCAdim1_4trts,group=Site),method="lm",formula=y~x,size=0.5,color="gray80",se=FALSE)+
  #geom_smooth(data=newdat,aes(y=pred_t,x=PCAdim1_4trts),method="lm",formula=y~x,size=1,color="gray35",se=FALSE)+
  #geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+
  
  
  labs(x="Fast-slow spectrum ",y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous() + scale_y_continuous(trans="log2",lim=c(0.8,16),breaks=c(1,2,4,8,16))

FStrt<-gg+ theme(axis.title.x=element_text(colour="black",face="bold",size=10),
               axis.title.y=element_text(colour="black",face="bold",size=10),
               axis.text.y=element_text(colour="black",face="bold",size=8),
               axis.text.x=element_text(colour="black",face="bold",size=8),
               plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
               panel.background = element_rect(fill = "white"))


#################
################

cairo_ps(filename="/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/eMNTD_FD_CWM_onTS_April2017.eps",
                   family="sans",
                   height=6,width=6,
                bg="white")

png(filename="/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/eMNTD_FD_CWM_onTS_April2017.png", 
    units="in", 
    width=6, 
    height=6, 
    pointsize=2, 
    res=200)

plot_grid(emNTD,FDis,FStrt, labels=c("(a)","(b)","(c)"), label_size = 7, ncol=2,align="hv")

dev.off()

