rm(list=ls()) 
#require(devtools)
#install_github("jslefche/piecewiseSEM")

require(piecewiseSEM)
library(semPlot)
library(lmerTest)
library(nlme)
library(dplyr)

# Data
rs_rl<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/StabII_ResistResil_FD_PD_PCA_CWM_29062016.csv",sep=",",header=T)

rs_rl$PlotUnique2<-paste(rs_rl$Study,rs_rl$Study,rs_rl$BlockUnique,rs_rl$PlotUnique,sep="_")
rs_rl$PlotUnique2<-as.factor(rs_rl$PlotUnique2)

rs_rl<-filter(rs_rl,Climate_Bin=="SPEI12")
rs_rl$Climate_Bin<-droplevels(rs_rl$Climate_Bin)

##########
## DRY ###
##########
rs_rll<-filter(rs_rl,Value=="Dry")
rs_rll$Value<-droplevels(rs_rll$Value)
rs_12<-select(rs_rll,Climate_Bin, Value,Dir, Int,ExpYear,Site,BlockUnique,PlotUnique,PlotUnique2,SppN,Yn_FRic4,Yn_FDis4,Yn_eMPD, Yn_eMNTD, Yn_ePSE,Yn_PCAcwm4trts,
              Yn_PlotBiomass,Ye_PlotBiomass,Rst_Plot)

#rs_12<-filter(rs_12,Int_SPEI12==1)
#rs_12$Bin_12<-droplevels(rs_12$Bin_12)
#rs_12<-subset(rs_12, !is.na(Rst_SPEI12))
#rs_12$FDis4_12_Yn<-ifelse(rs_12$SppN==1 & is.na(rs_12$FDis4_12_Yn)==TRUE,0,rs_12$FDis4_12_Yn)
#rs_12$eMPD_12_Yn<-ifelse(rs_12$SppN==1 & is.na(rs_12$eMPD_12_Yn)==TRUE,0,rs_12$eMPD_12_Yn)
#rs_12$eMNTD_12_Yn<-ifelse(rs_12$SppN==1 & is.na(rs_12$eMNTD_12_Yn)==TRUE,0,rs_12$eMNTD_12_Yn)
#rs_12$ePSE_12_Yn<-ifelse(rs_12$SppN==1 & is.na(rs_12$ePSE_12_Yn)==TRUE,0,rs_12$ePSE_12_Yn)


rs_12$PlotUnique2<-droplevels(rs_12$PlotUnique2)
rs_12$PlotUnique<-droplevels(rs_12$PlotUnique)
rs_12$Site<-droplevels(rs_12$Site)

rs_12<-arrange(rs_12,Site,PlotUnique2,ExpYear)

#Transformations

rs_12$Rst_Plot<-rs_122$Rst_Plot+.001

               
rs_12$SppN<-as.numeric(rs_12$SppN)
rs_12$lg2SppN <- log(rs_12$SppN,2)
rs_12$lg2Rst12<-log(rs_12$Rst_Plot,2)

rs_12$Dir<-as.factor(rs_12$Dir)
rs_12$Int<-as.factor(rs_12$Int)

rs_122<-subset(rs_12, !is.na(FRic4)) # for analysis with FRic4

####################
# fit SEMs #########
####################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpY_Ye |Site/PlotUnique2)
x2=corCompSymm(form=~ExpY_Ye |Site/PlotUnique2)


a<-lme(eMNTD_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
a1<-lme(eMNTD_12_Yn~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
AIC(a,a1)

b<-lme(FDis4_12_Yn~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
b1<-lme(FDis4_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(b,b1)

c<-lme(PCAdim1_4trts_12_Yn~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
c1<-lme(PCAdim1_4trts_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(c,c1)

d<-lme(lg2Rst12~lg2SppN+PCAdim1_4trts_12_Yn+FDis4_12_Yn+eMNTD_12_Yn,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
d1<-lme(lg2Rst12~lg2SppN+PCAdim1_4trts_12_Yn+FDis4_12_Yn+eMNTD_12_Yn,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(d,d1)

e<- lme(lg2Rst12~FDis4_12_Yn+PCAdim1_4trts_12_Yn+lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
e1<- lme(lg2Rst12~FDis4_12_Yn+PCAdim1_4trts_12_Yn+lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)



Rst_ModList=list(
  lme(eMNTD_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(FDis4_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(PCAdim1_4trts_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rst12~FDis4_12_Yn+PCAdim1_4trts_12_Yn+lg2SppN,random=~~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)

sem.fit(Rst_ModList,rs_122,corr.errors=c("eMNTD_12_Yn~~FDis4_12_Yn","PCAdim1_4trts_12_Yn~~eMNTD_12_Yn","PCAdim1_4trts_12_Yn~~FDis4_12_Yn"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

sem.fisher.c(Rst_ModList, data=rs_122,corr.errors=c("eMNTD_12_Yn~~FDis4_12_Yn","PCAdim1_4trts_12_Yn~~eMNTD_12_Yn","PCAdim1_4trts_12_Yn~~FDis4_12_Yn"),
             model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

sem.coefs(Rst_ModList,rs_122,standardize="scale")
d_pc<- sem.coefs(Rst_ModList,rs_122,standardize="scale")
#sem.coefs(Rst_ModList,rs_122,standardize="scale",corr.errors=c("eMNTD_12_Yn~~FDis4_12_Yn","PCAdim1_4trts_12_Yn~~eMNTD_12_Yn","PCAdim1_4trts_12_Yn~~FDis4_12_Yn"))



sem.missing.paths(Rst_ModList, rs_122,corr.errors=c("eMNTD_12_Yn~~FDis4_12_Yn","PCAdim1_4trts_12_Yn~~eMNTD_12_Yn","PCAdim1_4trts_12_Yn~~FDis4_12_Yn"),conditional=T)

#sem.model.fits(list(a,b,c,d),aicc=TRUE)

a<-lme(eMNTD_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
b<-lme(FDis4_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
c<-lme(PCAdim1_4trts_12_Yn~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
d<-lme(lg2Rst12~PCAdim1_4trts_12_Yn+lg2SppN+FDis4_12_Yn,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)

require(MuMIn)
r.squaredGLMM(d)

sem.plot(Rst_ModList,rs_122,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~lg2SppN,Rst_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

resids.df2<-partial.resid(lg2Rst12~PCAdim1_4trts_12_Yn,Rst_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

resids.df4<-partial.resid(lg2Rst12~FDis4_12_Yn,Rst_ModList,data=rs_122,
                         model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))


##################
## Extreme DRY ###
##################
#######################
# WET #################
#######################


#######################
# EXTREME WET #########
#######################



N