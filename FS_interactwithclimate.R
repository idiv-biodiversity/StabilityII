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

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,eMPD,eMNTD,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony,CV_Precip,CV_Temp,annualTemp, meanPrecip)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself

stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)


stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony

cc<-lmeControl(opt="optim")

#PCAdim1_4trts

Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~PCAdim1_4trts*meanPrecip,random=~1+PCAdim1_4trts|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~PCAdim1_4trts*meanPrecip,random=~1+PCAdim1_4trts|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~PCAdim1_4trts*meanPrecip,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~PCAdim1_4trts*meanPrecip,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~PCAdim1_4trts*meanPrecip,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~PCAdim1_4trts*meanPrecip,random=~1+lg2SppN*PCAdim1_4trts|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~PCAdim1_4trts*meanPrecip,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot(Cand.set[[7]])
qqnorm(Cand.set[[7]])


### LRT

big<-lme(TS_lg2~PCAdim1_4trts*CV_Temp,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~PCAdim1_4trts+CV_Temp,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)


ICClme(final)

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

FStrt<-gg+ theme(axis.title.x=element_text(colour="black",face="bold",size=8),
                 axis.title.y=element_text(colour="black",face="bold",size=8),
                 axis.text.y=element_text(colour="black",face="bold",size=8),
                 axis.text.x=element_text(colour="black",face="bold",size=8),
                 plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
                 panel.background = element_rect(fill = "white"))

