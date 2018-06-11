#######################################################################
# Effects on ecosystem stability of individual biodiversity facets   ##
#######################################################################
# 1. Species diversity
# 2. Species asynchrony
# 3. Phylogenetic diversity (MNTD)
# 4. Phylogenetic diversity (MPD)
# 5. Functional diversity (FDis)
# 6. Functional diversity (FRic)
# 7. Functional composition


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

# Prepare data for analysis

stab<-read.delim("data.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # remove site where we didn't have good trait coverage

stab_4<-select(stab,Site,Study_length,UniqueID,SppN,eMPD,eMNTD,FDis4,FRic4,PCAdim1_4trts,
               Plot_TempStab,Plot_Biomassxbar, Plot_Biomasssd,Gross_synchrony)

# for plots with ONLY 1 spp, we assume that a species is perfectly synchronized with itself (ie asynchrony = 1)

stab_4$Gross_synchrony<-ifelse(is.na(stab_4$Gross_synchrony)==TRUE,1,stab_4$Gross_synchrony) 

# convert synchrony metrics to different scale

stab_4$PlotAsynchrony_s<-stab_4$Plot_Asynchrony*-1

stab_4$GrossAsynchrony_s<-stab_4$Gross_synchrony*-1

# further adjustments

stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

# Filter out NAs for Asynchrony and FRic4
stab_444<-filter(stab_4, is.na(PlotAsynchrony_s)==FALSE)
stab_444<-filter(stab_444, is.na(FRic4)==FALSE)

#######################################################################
# Effects on ecosystem stability of individual biodiversity facets   ##
#######################################################################

# Control list set up for LMM in nlme 

cc<-lmeControl(opt="optim")

#########################
## Species Richness  ####
#########################

## Test multiple random effect structures; compare using AICc
Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~lg2SppN,random=~1|Site,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~lg2SppN,random=~1|Site/SppN,control=cc,data=stab_444)

Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot([[best_model]])
qqnorm([[best_model]])

# Log ratio test for significance of main effect

big<-lme(TS_lg2~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=~1+lg2SppN|Site,control=cc,data=stab_444,method="ML")

anova(big,small)

#final model

final<-lme(TS_lg2~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)

r.squaredGLMM(final)

intervals(final, which="fixed")

VarCorr(final)

ICClme(final)

#######################
# Species Asynchrony  #
#######################

## Test multiple random effect structures; compare using AICc 
Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~GrossAsynchrony_s,random=~1+GrossAsynchrony_s|Site,control=cc,data=stab_444)
Cand.set[[2]]<-lme(TS_lg2~GrossAsynchrony_s,random=~1+GrossAsynchrony_s|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~GrossAsynchrony_s,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~GrossAsynchrony_s,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~GrossAsynchrony_s,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~GrossAsynchrony_s,random=~1+lg2SppN*GrossAsynchrony_s|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~GrossAsynchrony_s,random=list(~1+lg2SppN+GrossAsynchrony_s|Site),control=cc,data=stab_444)


Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot(Cand.set[[best_model]])
qqnorm(Cand.set[[best_model]])

### Log ratio test for significance of main effect

big<-lme(TS_lg2~GrossAsynchrony_s,random=~1+lg2SppN*GrossAsynchrony_s|Site,control=cc,data=stab_444,method="ML")

small<-lme(TS_lg2~1,random=~1+lg2SppN*GrossAsynchrony_s|Site,control=cc,data=stab_444,method="ML")

anova(big,small)

#final model summary

final<-lme(TS_lg2~GrossAsynchrony_s,random=list(~1+lg2SppN*GrossAsynchrony_s|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

VarCorr(final)

ICClme(final)

##########################
# Phylogenetic diversity #
##########################

#eMNTD

## Test multiple random effect structures; compare using AICc 
Cand.set <- list( )
Cand.set[[1]]<-lme(TS_lg2~eMNTD,random=~1+eMNTD|Site,control=cc,data=cd)
Cand.set[[2]]<-lme(TS_lg2~eMNTD,random=~1+eMNTD|Site/SppN,control=cc,data=stab_444)
Cand.set[[3]]<-lme(TS_lg2~eMNTD,random=~1|Site,control=cc,data=stab_444)
Cand.set[[4]]<-lme(TS_lg2~eMNTD,random=~1|Site/SppN,control=cc,data=stab_444)
Cand.set[[5]]<-lme(TS_lg2~eMNTD,random=~1+lg2SppN|Site,control=cc,data=stab_444)
Cand.set[[6]]<-lme(TS_lg2~eMNTD,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444)
Cand.set[[7]]<-lme(TS_lg2~eMNTD,random=list(~1+lg2SppN+eMNTD|Site),control=cc,data=stab_444)

Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,second.ord = T)
res.table

plot(Cand.set[[best_model]])
qqnorm(Cand.set[[best_model]])

### Log ratio test for significance of main effect

big<-lme(TS_lg2~eMNTD*CV_Precip,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444,method="ML")

small<-lme(TS_lg2~1,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444,method="ML")

anova(big,small)


#final model summary

final<-lme(TS_lg2~eMNTD,random=~1+lg2SppN*eMNTD|Site,control=cc,data=stab_444)
r.squaredGLMM(final)

ICClme(final)
VarCorr(final)

#eMPD

## Test multiple random effect structures; compare using AICc 
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

plot(Cand.set[[best_model]])
qqnorm(Cand.set[[best_model]])

### Log ratio test for significance of main effect

big<-lme(TS_lg2~eMPD,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444,method="ML")

anova(big,small)

#final model

final<-lme(TS_lg2~eMPD,random=~1+lg2SppN*eMPD|Site,control=cc,data=stab_444)
r.squaredGLMM(final)

ICClme(final)
VarCorr(final)

########################
# Functional diversity # 
########################
#FDis (Functional Dispersion)

## Test multiple random effect structures; compare using AICc 
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

plot(Cand.set[[best_model]])
qqnorm(Cand.set[[best_model]])

### Log ratio test for significance of main effect

big<-lme(TS_lg2~FDis4,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="ML")

anova(big,small)


#final model summary

final<-lme(TS_lg2~FDis4,random=list(~1+lg2SppN+FDis4|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

ICClme(final)
VarCorr(final)

#FRic (functional richness)

## Test multiple random effect structures; compare using AICc 
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

plot(Cand.set[[best_model]])
qqnorm(Cand.set[[best_model]])

### Log ratio test for significance of main effect

big<-lme(TS_lg2~FRic4,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="ML")

anova(big,small)

#final model summary  

final<-lme(TS_lg2~FRic4,random=list(~1+lg2SppN*FRic4|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

ICClme(final)

VarCorr(final)


# PCAdim1_4trts (functional composition (CWM Fast-Slow))

## Test multiple random effect structures; compare using AICc 
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

plot(Cand.set[[best_model]])
qqnorm(Cand.set[[best_model]])


### Log ratio test for significance of main effect

big<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="ML")
small<-lme(TS_lg2~1,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="ML")

anova(big,small)

#final 

final<-lme(TS_lg2~PCAdim1_4trts,random=list(~1+lg2SppN+PCAdim1_4trts|Site),control=cc,data=stab_444,method="REML")

r.squaredGLMM(final)

ICClme(final)
VarCorr(final)