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

rs_rl$PlotUnique2<-paste(rs_rl$Study,rs_rl$Site,rs_rl$Block,rs_rl$Plot,sep="_")
rs_rl$PlotUnique2<-as.factor(rs_rl$PlotUnique2)

rs_rl<-filter(rs_rl,Climate_Bin=="SPEI12")
rs_rl$Climate_Bin<-droplevels(rs_rl$Climate_Bin)

##########
## WET ###
##########
rs_rll<-filter(rs_rl,Value=="Wet")
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
rs_12$Site<-droplevels(rs_12$Site)

rs_12<-arrange(rs_12,Site,PlotUnique2,ExpYear)

rs_12<-subset(rs_12, !is.na(Yn_eMPD))

#Transformations

rs_12$Rst_Plot<-rs_12$Rst_Plot+.001


rs_12$SppN<-as.numeric(rs_12$SppN)
rs_12$lg2SppN <- log(rs_12$SppN,2)
rs_12$lg2Rst12<-log(rs_12$Rst_Plot,2)

rs_12$Dir<-as.factor(rs_12$Dir)
rs_12$Int<-as.factor(rs_12$Int)

rs_122<-subset(rs_12, !is.na(Yn_FRic4)) # for analysis with FRic4

####################
# fit SEMs #########
####################
summarize(group_by(rs_12,Site),n_yrs=length(unique(ExpYear)))

#yes, temporal correlations are needed

####################
# Model 1         ##
# eMNTD + FDis #####
####################


bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpYear |Site/PlotUnique2)
x2=corCompSymm(form=~ExpYear |Site/PlotUnique2)


a<-lme(Yn_eMNTD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
a1<-lme(Yn_eMNTD~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
AIC(a,a1)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


b<-lme(Yn_FDis4~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
b1<-lme(Yn_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
b2<-lme(Yn_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x2,control=bb,data=rs_12)


re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


c<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
c1<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


d<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FDis4+Yn_eMNTD,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
d1<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FDis4+Yn_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 


Rst_Wet_ModList=list(
  lme(Yn_eMNTD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12),
  lme(Yn_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12),
  lme(lg2Rst12~Yn_PCAcwm4trts+lg2SppN+Yn_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
  
)

sem.fit(Rst_Wet_ModList,rs_12,corr.errors=c("Yn_eMNTD~~Yn_FDis4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))


# 

sem.fit(Rst_Wet_ModList,rs_12,corr.errors=c("Yn_eMNTD~~Yn_FDis4","Yn_FDis4~~Yn_PCAcwm4trts","Yn_eMNTD ~~ Yn_PCAcwm4trts"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))



#sem.coefs(Rst_Wet_ModList,rs_12,standardize="scale",corr.errors=c("Yn_eMNTD~~Yn_FDis4","Yn_FDis4~~Yn_PCAcwm4trts"))

wet_rst_fdis_emntd_pc<- sem.coefs(Rst_Wet_ModList,rs_12,standardize="scale")
wet_rst_fdis_emntd_pc$Climate_Bin<-"Wet"


wet_rst_fdis_emntd_modfit<-sem.model.fits(Rst_Wet_ModList)

wet_rst_fdis_emntd_modfit$ResponseVars<-c("eMNTD","FDis4","Resistance")
wet_rst_fdis_emntd_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,eMNTD,lg2SppN")
wet_rst_fdis_emntd_modfit$Climate_Bin<-"Wet"

#sem.plot(Rst_Dry_ModList,rs_12,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~lg2SppN,Rst_Wet_ModList,data=rs_12,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(wet_rst_fdis_emntd_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_emntd_fdis_sem_coefs.csv",sep=",",row.names=F)
write.table(wet_rst_fdis_emntd_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_emntd_fdis_model_fits.csv",sep=",",row.names=F)


#################
# Model 2      ##
# eMPD & FDis  ##
#################


bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpYear |Site/PlotUnique2)
x2=corCompSymm(form=~ExpYear |Site/PlotUnique2)


a<-lme(Yn_eMPD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
a1<-lme(Yn_eMPD~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
AIC(a,a1)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")



b<-lme(Yn_FDis4~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
b1<-lme(Yn_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(b,b1)

re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


c<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
c1<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


d<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FDis4+Yn_eMNTD,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
d1<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FDis4+Yn_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 


Rst_Wet_ModList=list(
  lme(Yn_eMPD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12),
  lme(Yn_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12),
  lme(lg2Rst12~Yn_PCAcwm4trts+lg2SppN+Yn_eMPD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
  
)

sem.fit(Rst_Wet_ModList,rs_12,corr.errors=c("Yn_eMPD~~Yn_FDis4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model

sem.fit(Rst_Wet_ModList,rs_12,corr.errors=c("Yn_eMPD~~Yn_FDis4","Yn_FDis4 ~~ Yn_PCAcwm4trts"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))


# no further changes
wet_rst_fdis_empd_pc<- sem.coefs(Rst_Wet_ModList,rs_12,standardize="scale")
wet_rst_fdis_empd_pc$Climate_Bin<-"Wet"


wet_rst_fdis_empd_modfit<-sem.model.fits(Rst_Wet_ModList)

wet_rst_fdis_empd_modfit$ResponseVars<-c("Yn_eMPD","FDis4","Resistance")
wet_rst_fdis_empd_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,Yn_eMPD,lg2SppN")
wet_rst_fdis_empd_modfit$Climate_Bin<-"Wet"

sem.plot(Rst_Dry_ModList,rs_12,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~Yn_PCAcwm4trts,Rst_Wet_ModList,data=rs_12,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(wet_rst_fdis_empd_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_empd_fdis_sem_coefs.csv",sep=",",row.names=F)
write.table(wet_rst_fdis_empd_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_empd_fdis_model_fits.csv",sep=",",row.names=F)


#################
## Model 3     ##
## FDis $ ePSE ##
#################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpYear |Site/PlotUnique2)
x2=corCompSymm(form=~ExpYear |Site/PlotUnique2)


a<-lme(Yn_ePSE~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
a1<-lme(Yn_ePSE~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
AIC(a,a1)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


b<-lme(Yn_FDis4~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
b1<-lme(Yn_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(b,b1)

re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


c<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
c1<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


d<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FDis4+Yn_eMNTD,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
d1<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FDis4+Yn_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 


Rst_Wet_ModList=list(
  lme(Yn_ePSE~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12),
  lme(Yn_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12),
  lme(lg2Rst12~Yn_PCAcwm4trts+lg2SppN+Yn_ePSE,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
  
)


sem.fit(Rst_Wet_ModList,rs_12,corr.errors=c("Yn_ePSE~~Yn_FDis4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model

sem.fit(Rst_Wet_ModList,rs_12,corr.errors=c("Yn_ePSE~~Yn_FDis4","Yn_FDis4 ~~ Yn_PCAcwm4trts","Yn_ePSE ~~ Yn_PCAcwm4trts"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))


# no further changes

wet_rst_fdis_epse_pc<- sem.coefs(Rst_Wet_ModList,rs_12,standardize="scale")
wet_rst_fdis_epse_pc$Climate_Bin<-"Wet"

wet_rst_fdis_epse_modfit<-sem.model.fits(Rst_Wet_ModList)

wet_rst_fdis_epse_modfit$ResponseVars<-c("Yn_ePSE","FDis4","Resistance")
wet_rst_fdis_epse_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,Yn_ePSE,lg2SppN")
wet_rst_fdis_epse_modfit$Climate_Bin<-"Wet"

sem.plot(Rst_Wet_ModList,rs_12,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~Yn_PCAcwm4trts,Rst_Wet_ModList,data=rs_12,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(wet_rst_fdis_epse_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_epse_fdis_sem_coefs.csv",sep=",",row.names=F)
write.table(wet_rst_fdis_epse_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_epse_fdis_model_fits.csv",sep=",",row.names=F)

############
# FRic4 ####
############

####################
# Model 4         ##
# eMNTD + FRic4 ####
####################


bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpYear |Site/PlotUnique2)
x2=corCompSymm(form=~ExpYear |Site/PlotUnique2)


a<-lme(Yn_eMNTD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
a1<-lme(Yn_eMNTD~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
AIC(a,a1)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")

b<-lme(Yn_FRic4~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
b1<-lme(Yn_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(b,b1)

re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")


c<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
c1<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


d<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FRic4+Yn_eMNTD,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
d1<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FRic4+Yn_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 


Rst_Wet_ModList=list(
  lme(Yn_eMNTD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Yn_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rst12~Yn_PCAcwm4trts+lg2SppN+Yn_FRic4,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)

sem.fit(Rst_Wet_ModList,rs_122,corr.errors=c("Yn_eMNTD~~Yn_FRic4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model

sem.fit(Rst_Wet_ModList,rs_122,corr.errors=c("Yn_eMNTD~~Yn_FRic4","Yn_eMNTD ~~ Yn_PCAcwm4trts","Yn_FRic4 ~~ Yn_PCAcwm4trts"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))


# no further changes

#sem.coefs(Rst_Dry_ModList,rs_122,standardize="scale",corr.errors=c("Yn_eMNTD~~Yn_FRic4"))

sem.coefs(Rst_Wet_ModList,rs_122,standardize="scale")

wet_rst_fric_emntd_pc<- sem.coefs(Rst_Wet_ModList,rs_122,standardize="scale")
wet_rst_fric_emntd_pc$Climate_Bin<-"Wet"

wet_rst_fric_emntd_modfit<-sem.model.fits(Rst_Wet_ModList)

wet_rst_fric_emntd_modfit$ResponseVars<-c("eMNTD","FRic4","Resistance")
wet_rst_fric_emntd_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,Yn_FRic4,lg2SppN")
wet_rst_fric_emntd_modfit$Climate_Bin<-"Wet"

sem.plot(Rst_Dry_ModList,rs_122,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~Yn_PCAcwm4trts,Rst_Wet_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(wet_rst_fric_emntd_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_emntd_fric_sem_coefs_V2.csv",sep=",",row.names=F)
write.table(wet_rst_fric_emntd_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_WET_emntd_fric_model_fits_V2.csv",sep=",",row.names=F)

#################
# Model 5 #######
# eMPD & FRIC  ##
#################


bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpYear |Site/PlotUnique2)
x2=corCompSymm(form=~ExpYear |Site/PlotUnique2)


a<-lme(Yn_eMPD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
a1<-lme(Yn_eMPD~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
AIC(a,a1)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")


b<-lme(Yn_FRic4~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
b1<-lme(Yn_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(b,b1)

re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")


c<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
c1<-lme(Yn_PCAcwm4trts~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(c,c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")


d<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FRic4+Yn_eMNTD,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
d1<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FRic4+Yn_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 


Rst_Wet_ModList=list(
  lme(Yn_eMPD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Yn_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rst12~Yn_FRic4+Yn_PCAcwm4trts+lg2SppN+Yn_eMPD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)

sem.fit(Rst_Wet_ModList,rs_122,corr.errors=c("Yn_eMPD~~Yn_FRic4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model

sem.fit(Rst_Wet_ModList,rs_122,corr.errors=c("Yn_eMPD~~Yn_FRic4","Yn_eMPD ~~ Yn_PCAcwm4trts","Yn_FRic4 ~~ Yn_PCAcwm4trts"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))


# no further changes

dry_rst_fric_empd_pc<- sem.coefs(Rst_Dry_ModList,rs_122,standardize="scale")
dry_rst_fric_empd_pc$Climate_Bin<-"Dry"

dry_rst_fric_empd_modfit<-sem.model.fits(Rst_Dry_ModList)
dry_rst_fric_empd_modfit$ResponseVars<-c("Yn_eMPD","FRic4","Resistance")
dry_rst_fric_empd_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,Yn_eMPD,Yn_FRic4,lg2SppN")
dry_rst_fric_empd_modfit$Climate_Bin<-"Dry"

sem.plot(Rst_Dry_ModList,rs_122,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~lg2SppN,Rst_Dry_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(dry_rst_fric_empd_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_DRY_empd_fric_sem_coefs.csv",sep=",",row.names=F)
write.table(dry_rst_fric_empd_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_DRY_empd_fric_model_fits.csv",sep=",",row.names=F)


#################
## Model 6    ###
## FDis $ ePSE ##
#################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpYear |Site/PlotUnique2)
x2=corCompSymm(form=~ExpYear |Site/PlotUnique2)


a<-lme(Yn_ePSE~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
a1<-lme(Yn_ePSE~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
AIC(a,a1)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")


b<-lme(Yn_FRic4~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
b1<-lme(Yn_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(b,b1)

re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


d<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FRic4+Yn_eMNTD,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
d1<-lme(lg2Rst12~lg2SppN+Yn_PCAcwm4trts+Yn_FRic4+Yn_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(d,d1)

re<-resid(d1, type="normalized")
fi<-fitted(d1)  
plot(x=fi,y=re,xlab="fitted values",ylab="residuals") 


Rst_Dry_ModList=list(
  lme(Yn_ePSE~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Yn_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rst12~Yn_FRic4+Yn_PCAcwm4trts+lg2SppN+Yn_ePSE,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)


sem.fit(Rst_Dry_ModList,rs_122,corr.errors=c("Yn_ePSE~~Yn_FRic4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model


dry_rst_fric_epse_pc<- sem.coefs(Rst_Dry_ModList,rs_122,standardize="scale")
dry_rst_fric_epse_pc$Climate_Bin<-"Dry"

dry_rst_fric_epse_modfit<-sem.model.fits(Rst_Dry_ModList)
dry_rst_fric_epse_modfit$ResponseVars<-c("Yn_eMPD","FRic4","Resistance")
dry_rst_fric_epse_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,Yn_ePSE,Yn_FRic4,lg2SppN")
dry_rst_fric_epse_modfit$Climate_Bin<-"Dry"

sem.plot(Rst_Dry_ModList,rs_12,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~lg2SppN,Rst_Dry_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(dry_rst_fric_epse_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_DRY_epse_fric_sem_coefs.csv",sep=",",row.names=F)
write.table(dry_rst_fric_epse_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rst_DRY_epse_fric_model_fits.csv",sep=",",row.names=F)

