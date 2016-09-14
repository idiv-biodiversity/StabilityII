rm(list=ls()) 

library(lmerTest)
library(nlme)
library(dplyr)
library(AICcmodavg)
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


###############################################
# Effects on Synchrony of biodiversity facets #
###############################################

cc<-lmeControl(opt="optim")

###########
## SR  ####
###########


a<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
a1<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_444)
a2<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site/lg2SppN,control=cc,data=stab_444)

AICc(a,a1,a2)

plot(a2)
qqnorm(a2)


big<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site/SppN,control=cc,data=stab_444,method="ML")
small<-lme(Plot_Asynchrony~1,random=~1|Site/SppN,control=cc,data=stab_444,method="ML")

anova(big,small)

#final
final<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site/SppN,control=cc,data=stab_444)

r.squaredGLMM(final)



###############
# predictions #
###############
spp_ts<-select(stab_444,Site,UniqueID,SppN,TS_lg2,lg2SppN)

spp_ts$pred<-predict(a,spp_ts,re.form=~(~1+lg2SppN|Site))

Designmat <- model.matrix(formula(a1)[-2], spp_ts)
predvar <- diag(Designmat %*% vcov(a1) %*% t(Designmat)) 
spp_ts$SE <- sqrt(predvar) 
spp_ts$SE2 <- sqrt(predvar+a1$sigma^2)

c<-ggplot(data=spp_ts,aes(x=lg2SppN,y=TS_lg2))+
  geom_point(colour="white",fill="black",shape= 21,cex=3) + 
  
  geom_smooth(data=spp_ts,aes(y=pred,x=lg2SppN),method="lm",formula=y~x,size=1,color="black",se=FALSE)+
  labs(x="Time since abandoment (years)",y=expression(paste(bold("SES"["FRic"]))))+
  scale_y_continuous(limits=c(-5.0,5),breaks=c(-5,-4,-3,-2,-1,0 ,1,2,3,4,5))



##########
# PD #####
##########

#eMNTD

b<-lme(Plot_Asynchrony~eMNTD,random=~1+eMNTD|Site,control=cc,data=stab_444)
b1<-lme(Plot_Asynchrony~eMNTD,random=~1+eMNTD|Site/SppN,control=cc,data=stab_444)
b2<-lme(Plot_Asynchrony~eMNTD,random=~1|Site,control=cc,data=stab_444)
b3<-lme(Plot_Asynchrony~eMNTD,random=~1|Site/SppN,control=cc,data=stab_444)
b4<-lme(Plot_Asynchrony~eMNTD,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b5<-lme(Plot_Asynchrony~eMNTD,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444)
b6<-lme(Plot_Asynchrony~eMNTD,random=list(~1+lg2SppN+eMNTD|Site),control=cc,data=stab_444)


AICc(b,b1,b2,b3,b4,b5,b6)

plot(b5)
qqnorm(b5)



### LRT

big<-lme(Plot_Asynchrony~eMNTD,random=~1+eMNTD|Site/SppN,control=cc,data=stab_444,method="ML")
small<-lme(Plot_Asynchrony~1,random=~1+eMNTD|Site/SppN,control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(Plot_Asynchrony~eMNTD,random=~1+eMNTD|Site/SppN,control=cc,data=stab_444,method="REML")
r.squaredGLMM(final)


################
# predictions  #
################

#eMPD

b<-lme(Plot_Asynchrony~eMPD,random=~1+eMPD|Site,control=cc,data=stab_444)
b1<-lme(Plot_Asynchrony~eMPD,random=~1+eMPD|Site/SppN,control=cc,data=stab_444)
b2<-lme(Plot_Asynchrony~eMPD,random=~1|Site,control=cc,data=stab_444)
b3<-lme(Plot_Asynchrony~eMPD,random=~1|Site/SppN,control=cc,data=stab_444)
b4<-lme(Plot_Asynchrony~eMPD,random=~1+lg2SppN|Site,control=cc,data=stab_444)
b5<-lme(Plot_Asynchrony~eMPD,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444)
b6<-lme(Plot_Asynchrony~eMPD,random=list(~1+lg2SppN+eMPD|Site),control=cc,data=stab_444)

AICc(b,b1,b2,b3,b4,b5,b6)

plot(b6)
qqnorm(b6)

### LRT

big<-lme(Plot_Asynchrony~eMPD,random=~1+eMPD|Site/SppN,control=cc,data=stab_444,method="ML")
small<-lme(Plot_Asynchrony~1,random=~1+eMPD|Site/SppN,control=cc,data=stab_444,method="ML")

anova(big,small)


#final 

final<-lme(Plot_Asynchrony~eMPD,random=~1+eMPD|Site/SppN,control=cc,data=stab_444,method="REML")
r.squaredGLMM(final)


[[STOP HERE]]
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

