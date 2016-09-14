rm(list=ls()) 

library(lmerTest)
library(nlme)
library(dplyr)
library(AICcmodavg)
library(MuMIn)
library(ggplot2)

# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)


stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony


###########################################
# Effects on TS of biodiversity facets   ##
###########################################

cc<-lmeControl(opt="optim")

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

    geom_smooth(data=spp_ts,aes(y=pred_t,x=SppN,group=Site),method="lm",formula=y~x,size=0.5,color="gray50",se=FALSE)+
  geom_smooth(data=newdat,aes(y=pred_t,x=SppN),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+


  labs(x="Plant species richness",y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous(trans="log2",breaks=c(1,2,4,8,16,32,60)) + scale_y_continuous(trans="log2",breaks=c(0.5,1,2,4,8,16))

SppN<-c+ theme(axis.title.x=element_text(colour="black",face="bold",size=8),
                axis.title.y=element_text(colour="black",face="bold",size=8,vjust=1),
                axis.text.y=element_text(colour="black",face="bold",size=8),
                axis.text.x=element_text(colour="black",face="bold",size=8),
                plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
               panel.background = element_rect(fill = "white"))



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

big<-lme(TS_lg2~eMNTD,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~eMNTD,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444)
r.squaredGLMM(final)


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
  
  geom_smooth(data=emntd_ts,aes(y=pred_t,x=eMNTD,group=Site),method="lm",formula=y~x,size=0.5,color="gray50",se=FALSE)+
  geom_smooth(data=newdat,aes(y=pred_t,x=eMNTD),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  geom_ribbon(data=newdat,aes(ymin=pred_lCIt,ymax=pred_uCIt),fill="gray50",colour="transparent",alpha=0.4)+
  
  
  labs(x="Phylogenetic diversity (eMNTD)",y=expression(bold(paste("Ecosystem stability ( ", mu," / ",sigma," )")))) +
  scale_x_continuous() + scale_y_continuous(trans="log2",breaks=c(0.5,1,2,4,8,16))

emNTD<-d+ theme(axis.title.x=element_text(colour="black",face="bold",size=8),
               axis.title.y=element_text(colour="black",face="bold",size=8,vjust=1),
               axis.text.y=element_text(colour="black",face="bold",size=8),
               axis.text.x=element_text(colour="black",face="bold",size=8),
               plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),panel.border=element_rect(fill=NA,colour="black"),
               panel.background = element_rect(fill = "white"))

[[STOP HERE]]

#eMPD

b<-lme(TS_lg2~eMPD,random=~1+eMPD|Site,control=cc,data=stab_444)
b1<-lme(TS_lg2~eMPD,random=~1+eMPD|Site/SppN,control=cc,data=stab_444)
b2<-lme(TS_lg2~eMPD,random=~1|Site,control=cc,data=stab_444)
b3<-lme(TS_lg2~eMPD,random=~1|Site/SppN,control=cc,data=stab_444)
b4<-lme(TS_lg2~eMPD,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b5<-lme(TS_lg2~eMPD,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444)
b6<-lme(TS_lg2~eMPD,random=list(~1+lg2SppN+eMPD|Site),control=cc,data=stab_444)

AICc(b,b1,b2,b3,b4,b5,b6)

plot(b6)
qqnorm(b6)

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

#ePSE

b<-lme(TS_lg2~ePSE,random=~1+ePSE|Site,control=cc,data=stab_444)
b1<-lme(TS_lg2~ePSE,random=~1+ePSE|Site/SppN,control=cc,data=stab_444)
b2<-lme(TS_lg2~ePSE,random=~1|Site,control=cc,data=stab_444)
b3<-lme(TS_lg2~ePSE,random=~1|Site/SppN,control=cc,data=stab_444)
b4<-lme(TS_lg2~ePSE,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b5<-lme(TS_lg2~ePSE,random=~1+lg2SppN*ePSE|Site,control=cc,data=stab_444)
b6<-lme(TS_lg2~ePSE,random=list(~1+lg2SppN+ePSE|Site),control=cc,data=stab_444)

AICc(b,b1,b2,b3,b4,b5,b6)

plot(b6)
qqnorm(b6)

### LRT

big<-lme(TS_lg2~ePSE,random=list(~1+lg2SppN+ePSE|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN+ePSE|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~ePSE,random=list(~1+lg2SppN+ePSE|Site),control=cc,data=stab_444)

r.squaredGLMM(final)



################
# predictions  #
################

#FDis

b<-lme(TS_lg2~FDis4,random=~1+FDis4|Site,control=cc,data=stab_444)
b1<-lme(TS_lg2~FDis4,random=~1+FDis4|Site/SppN,control=cc,data=stab_444)
b2<-lme(TS_lg2~FDis4,random=~1|Site,control=cc,data=stab_444)
b3<-lme(TS_lg2~FDis4,random=~1|Site/SppN,control=cc,data=stab_444)
b4<-lme(TS_lg2~FDis4,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b5<-lme(TS_lg2~FDis4,random=~1+lg2SppN*FDis4|Site,control=cc,data=stab_444)
b6<-lme(TS_lg2~FDis4,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444)

AICc(b,b1,b2,b3,b4,b5,b6)

plot(b5)
qqnorm(b5)

### LRT

big<-lme(TS_lg2~FDis4,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~FDis4,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)


################
# predictions  #
################


#FRic

b<-lme(TS_lg2~FRic4,random=~1+FRic4|Site,control=cc,data=stab_444)
b1<-lme(TS_lg2~FRic4,random=~1+FRic4|Site/SppN,control=cc,data=stab_444)
b2<-lme(TS_lg2~FRic4,random=~1|Site,control=cc,data=stab_444)
b3<-lme(TS_lg2~FRic4,random=~1|Site/SppN,control=cc,data=stab_444)
b4<-lme(TS_lg2~FRic4,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b5<-lme(TS_lg2~FRic4,random=~1+lg2SppN*FRic4|Site,control=cc,data=stab_444)
b6<-lme(TS_lg2~FRic4,random=list(~1+lg2SppN+FRic4|Site),control=cc,data=stab_444)

AICc(b,b1,b2,b3,b4,b5,b6)

plot(b4)
qqnorm(b4)

### LRT

big<-lme(TS_lg2~FRic4,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~FRic4,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)


################
# predictions  #
################

#PCAdim1_4trts
b<-lme(TS_lg2~PCAdim1_4trts,random=~1+PCAdim1_4trts|Site,control=cc,data=stab_444)
b1<-lme(TS_lg2~PCAdim1_4trts,random=~1+PCAdim1_4trts|Site/lg2SppN,control=cc,data=stab_444)
b2<-lme(TS_lg2~PCAdim1_4trts,random=~1|Site,control=cc,data=stab_444)
b3<-lme(TS_lg2~PCAdim1_4trts,random=~1|Site/SppN,control=cc,data=stab_444)
b4<-lme(TS_lg2~PCAdim1_4trts,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b5<-lme(TS_lg2~PCAdim1_4trts,random=~1+lg2SppN*PCAdim1_4trts|Site,control=cc,data=stab_444)
b6<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444)

AICc(b,b1,b2,b3,b4,b5,b6)

plot(b6)
qqnorm(b6)


### LRT

big<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN*PCAdim1_4trts|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

################
# predictions  #
################


#Plot_Asynchrony

b<-lme(TS_lg2~Plot_Asynchrony,random=~1+Plot_Asynchrony|Site,control=cc,data=stab_444)
b1<-lme(TS_lg2~Plot_Asynchrony,random=~1+Plot_Asynchrony|Site/SppN,control=cc,data=stab_444)
b2<-lme(TS_lg2~Plot_Asynchrony,random=~1|Site,control=cc,data=stab_444)
b3<-lme(TS_lg2~Plot_Asynchrony,random=~1|Site/SppN,control=cc,data=stab_444)
b4<-lme(TS_lg2~Plot_Asynchrony,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b5<-lme(TS_lg2~Plot_Asynchrony,random=~1+lg2SppN*Plot_Asynchrony|Site,control=cc,data=stab_444)
b6<-lme(TS_lg2~Plot_Asynchrony,random=list(~1+lg2SppN+Plot_Asynchrony|Site),control=cc,data=stab_444)

AICc(b,b1,b2,b3,b4,b5,b6)

plot(b4)
qqnorm(b4)

### LRT

big<-lme(TS_lg2~Plot_Asynchrony,random=list(~1+lg2SppN*Plot_Asynchrony|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN*Plot_Asynchrony|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(TS_lg2~Plot_Asynchrony,random=list(~1+lg2SppN*Plot_Asynchrony|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

