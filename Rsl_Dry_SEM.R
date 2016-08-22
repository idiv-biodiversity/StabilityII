rm(list=ls()) 
#require(devtools)
#install_github("jslefche/piecewiseSEM")

require(piecewiseSEM)
library(semPlot)
#library(lmerTest)
library(nlme)
library(dplyr)

# Data
rs_rl<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/StabII_ResistResil_FD_PD_PCA_CWM_29062016.csv",sep=",",header=T)

rs_rl$PlotUnique2<-paste(rs_rl$Study,rs_rl$Study,rs_rl$BlockUnique,rs_rl$PlotUnique,sep="_")
rs_rl$PlotUnique2<-as.factor(rs_rl$PlotUnique2)

rs_rl<-filter(rs_rl,Climate_Bin=="SPEI12")
rs_rll<-filter(rs_rl,Value=="Dry")

##################
## Dry ###########
##################

rs_rll$Value<-droplevels(rs_rll$Value)
rs_12<-select(rs_rll,Climate_Bin, Value,Dir, Int,ExpYear,Site,BlockUnique,PlotUnique,PlotUnique2,SppN,Ye_FRic4,Ye_FDis4,Ye_eMPD, Ye_eMNTD, Ye_ePSE,Ye_PCAcwm4trts,
              Ye_PlotBiomass,Ye1_PlotBiomass,Rsl_Plot)


rs_12$PlotUnique2<-droplevels(rs_12$PlotUnique2)
rs_12$Site<-droplevels(rs_12$Site)

rs_12<-arrange(rs_12,Site,PlotUnique2,ExpYear)

rs_12<-subset(rs_12, !is.na(Rsl_Plot))
rs_12<-subset(rs_12, !is.na(Ye_eMPD))

#Transformations

rs_12$Rsl_Plot<-rs_12$Rsl_Plot+.001


rs_12$SppN<-as.numeric(rs_12$SppN)
rs_12$lg2SppN <- log(rs_12$SppN,2)
rs_12$lg2Rsl12<-log(rs_12$Rsl_Plot,2)

rs_12$Dir<-as.factor(rs_12$Dir)
rs_12$Int<-as.factor(rs_12$Int)

rs_122<-subset(rs_12, !is.na(Ye_FRic4)) # for analysis with FRic4

##################
# Preliminaries ##
##################


summarize(group_by(rs_12,Site),n_yrs=length(unique(ExpYear)))

#yes, temporal correlations are needed

######################
# preliminaries FDis #
######################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpYear |Site/PlotUnique2)
x2=corCompSymm(form=~ExpYear |Site/PlotUnique2)

#PD
a<-lme(Ye_eMNTD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
a1<-lme(Ye_eMNTD~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
AIC(a,a1)

plot(a)
qqnorm(a)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")

a<-lme(Ye_eMPD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
a1<-lme(Ye_eMPD~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
AIC(a,a1)

plot(a)
qqnorm(a)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")

a<-lme(Ye_ePSE~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
a1<-lme(Ye_ePSE~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_12)
AIC(a,a1)

plot(a)
qqnorm(a)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")



#FD

b<-lme(Ye_FDis4~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
b1<-lme(Ye_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)

AIC(b,b1)

plot(b1)
qqnorm(b1)


re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")

#F-S

c<-lme(Ye_PCAcwm4trts~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
c1<-lme(Ye_PCAcwm4trts~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(c,c1)

plot(c1)
qqnorm(c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=rs_12$lg2SppN,y=re,xlab="SppN",ylab="residuals")


# Rsl

d<-lme(lg2Rsl12~lg2SppN+Ye_PCAcwm4trts+Ye_FDis4+Ye_eMNTD,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
d1<-lme(lg2Rsl12~lg2SppN+Ye_PCAcwm4trts+Ye_FDis4+Ye_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_12)
AIC(d,d1)

plot(d)
qqnorm(d)

#######################
# preliminaries: FRic #
#######################



#PD
a<-lme(Ye_eMNTD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
a1<-lme(Ye_eMNTD~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
AIC(a,a1)

plot(a)
qqnorm(a)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")


re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")

a<-lme(Ye_eMPD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
a1<-lme(Ye_eMPD~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
AIC(a,a1)

plot(a)
qqnorm(a)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")

a<-lme(Ye_ePSE~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
a1<-lme(Ye_ePSE~lg2SppN,random=~1|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122)
AIC(a,a1)

plot(a)
qqnorm(a)

re<-resid(a, type="normalized")
fi<-fitted(a)  
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")


#FD

b<-lme(Ye_FRic4~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
b1<-lme(Ye_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)

AIC(b,b1)

plot(b1)
qqnorm(b1)


re<-resid(b1, type="normalized")
fi<-fitted(b1)  
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")

#F-S

c<-lme(Ye_PCAcwm4trts~lg2SppN,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
c1<-lme(Ye_PCAcwm4trts~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(c,c1)

plot(c1)
qqnorm(c1)

re<-resid(c1, type="normalized")
fi<-fitted(c1)  
plot(x=rs_122$lg2SppN,y=re,xlab="SppN",ylab="residuals")


# Rsl

d<-lme(lg2Rsl12~lg2SppN+Ye_PCAcwm4trts+Ye_FRic4+Ye_eMNTD,random=~1|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
d1<-lme(lg2Rsl12~lg2SppN+Ye_PCAcwm4trts+Ye_FRic4+Ye_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
AIC(d,d1)

plot(d)
qqnorm(d)


####################
# fit SEMs #########
####################
bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


x1=corAR1(form=~ExpYear |Site/PlotUnique2)
x2=corCompSymm(form=~ExpYear |Site/PlotUnique2)

####################
# Model 1         ##
# eMNTD + FDis #####
####################

Rsl_Dry_ModList=list(
  lme(Ye_eMNTD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Ye_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rsl12~Ye_PCAcwm4trts+lg2SppN+Ye_FDis4+Ye_eMNTD,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)

sem.fit(Rsl_Dry_ModList,rs_122,corr.errors=c("Ye_eMNTD~~Ye_FDis4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model


sem.fit(Rsl_Dry_ModList,rs_122,corr.errors=c("Ye_eMNTD~~Ye_FDis4","Ye_FDis4 ~~ Ye_PCAcwm4trts"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))


dry_rsl_fdis_emntd_pc<- sem.coefs(Rsl_Dry_ModList,rs_122,standardize="scale")
dry_rsl_fdis_emntd_pc$Climate_Bin<-"Dry"


dry_rsl_fdis_emntd_modfit<-sem.model.fits(Rsl_Dry_ModList)

dry_rsl_fdis_emntd_modfit$ResponseVars<-c("eMNTD","FDis4","Resistance")
dry_rsl_fdis_emntd_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,eMNTD,FDis4,lg2SppN")
dry_rsl_fdis_emntd_modfit$Climate_Bin<-"Dry"

#sem.plot(Rst_Dry_ModList,rs_12,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rsl12~Ye_PCAcwm4trts,Rsl_Dry_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(dry_rsl_fdis_emntd_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_emntd_fdis_sem_coefs_NEW.csv",sep=",",row.names=F)
write.table(dry_rsl_fdis_emntd_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_emntd_fdis_model_fits_NEW.csv",sep=",",row.names=F)


#################
# Model 2      ##
# eMPD & FDis  ##
#################



Rsl_Dry_ModList=list(
  lme(Ye_eMPD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Ye_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rsl12~Ye_eMPD+Ye_PCAcwm4trts+lg2SppN+Ye_FDis4,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)

sem.fit(Rsl_Dry_ModList,rs_122,corr.errors=c("Ye_eMPD~~Ye_FDis4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model

sem.fit(Rsl_Dry_ModList,rs_122,corr.errors=c("Ye_eMPD~~Ye_FDis4","Ye_FDis4 ~~ Ye_PCAcwm4trts"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))


# no further changes
dry_rsl_fdis_empd_pc<- sem.coefs(Rsl_Dry_ModList,rs_122,standardize="scale")
dry_rsl_fdis_empd_pc$Climate_Bin<-"Dry"


dry_rsl_fdis_empd_modfit<-sem.model.fits(Rsl_Dry_ModList)

dry_rsl_fdis_empd_modfit$ResponseVars<-c("Ye_eMPD","FDis4","Resilience")
dry_rsl_fdis_empd_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,Ye_eMPD,FDis4,lg2SppN")
dry_rsl_fdis_empd_modfit$Climate_Bin<-"Dry"

sem.plot(Rst_Dry_ModList,rs_12,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~Yn_PCAcwm4trts,Rst_Wet_ModList,data=rs_12,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(dry_rsl_fdis_empd_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_empd_fdis_sem_coefs_NEW.csv",sep=",",row.names=F)
write.table(dry_rsl_fdis_empd_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_empd_fdis_model_fits_NEW.csv",sep=",",row.names=F)

#################
## Model 3     ##
## FDis $ ePSE ##
#################


Rsl_Dry_ModList=list(
  lme(Ye_ePSE~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Ye_FDis4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rsl12~Ye_FDis4+Ye_PCAcwm4trts+lg2SppN+Ye_ePSE,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)


sem.fit(Rsl_Dry_ModList,rs_122,corr.errors=c("Ye_ePSE~~Ye_FDis4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model


sem.fit(Rsl_Dry_ModList,rs_122,corr.errors=c("Ye_ePSE~~Ye_FDis4","Ye_FDis4 ~~ Ye_PCAcwm4trts"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))



# no further changes

dry_rsl_fdis_epse_pc<- sem.coefs(Rsl_Dry_ModList,rs_122,standardize="scale")
dry_rsl_fdis_epse_pc$Climate_Bin<-"Dry"

dry_rsl_fdis_epse_modfit<-sem.model.fits(Rsl_Dry_ModList)

dry_rsl_fdis_epse_modfit$ResponseVars<-c("Ye_ePSE","FDis4","Resilience")
dry_rsl_fdis_epse_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,FDis4,Yn_ePSE,lg2SppN")
dry_rsl_fdis_epse_modfit$Climate_Bin<-"Dry"


resids.df1<-partial.resid(lg2Rsl12~Ye_PCAcwm4trts,Rsl_Dry_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(dry_rsl_fdis_epse_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_epse_fdis_sem_coefs_NEW.csv",sep=",",row.names=F)
write.table(dry_rsl_fdis_epse_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_epse_fdis_model_fits_NEW.csv",sep=",",row.names=F)

############
# FRic4 ####
############

####################
# Model 4         ##
# eMNTD + FRic4 ####
####################


Rsl_Dry_ModList=list(
  lme(Ye_eMNTD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Ye_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rsl12~Ye_PCAcwm4trts+lg2SppN+Ye_eMNTD+Ye_FRic4,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)

sem.fit(Rsl_Dry_ModList,rs_122,corr.errors=c("Ye_eMNTD~~Ye_FRic4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model


# no further changes

#sem.coefs(Rst_Dry_ModList,rs_122,standardize="scale",corr.errors=c("Yn_eMNTD~~Yn_FRic4"))

sem.coefs(Rsl_Dry_ModList,rs_122,standardize="scale")

dry_rsl_fric_emntd_pc<- sem.coefs(Rsl_Dry_ModList,rs_122,standardize="scale")
dry_rsl_fric_emntd_pc$Climate_Bin<-"Dry"

dry_rsl_fric_emntd_modfit<-sem.model.fits(Rsl_Dry_ModList)

dry_rsl_fric_emntd_modfit$ResponseVars<-c("eMNTD","FRic4","Resilience")
dry_rsl_fric_emntd_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,eMNTD,FRic4,lg2SppN")
dry_rsl_fric_emntd_modfit$Climate_Bin<-"Dry"

sem.plot(Rst_Dry_ModList,rs_122,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~Yn_PCAcwm4trts,Rst_Wet_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(dry_rsl_fric_emntd_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_emntd_fric_sem_coefs.csv",sep=",",row.names=F)
write.table(dry_rsl_fric_emntd_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_emntd_fric_model_fits.csv",sep=",",row.names=F)

#################
# Model 5 #######
# eMPD & FRIC  ##
#################



Rsl_Wet_ModList=list(
  lme(Ye_eMPD~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Ye_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rsl12~Ye_PCAcwm4trts+lg2SppN+Ye_eMPD+Ye_FRic4,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)

sem.fit(Rsl_Wet_ModList,rs_122,corr.errors=c("Ye_eMPD~~Ye_FRic4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model

# no further changes

dry_rsl_fric_empd_pc<- sem.coefs(Rsl_Wet_ModList,rs_122,standardize="scale")
dry_rsl_fric_empd_pc$Climate_Bin<-"Dry"

dry_rsl_fric_empd_modfit<-sem.model.fits(Rsl_Wet_ModList)
dry_rsl_fric_empd_modfit$ResponseVars<-c("Yn_eMPD","FRic4","Resilience")
dry_rsl_fric_empd_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,Ye_eMPD,FRic4,lg2SppN")
dry_rsl_fric_empd_modfit$Climate_Bin<-"Dry"

sem.plot(Rst_Dry_ModList,rs_122,show.nonsig = FALSE,scaling=20)

resids.df1<-partial.resid(lg2Rst12~Yn_PCAcwm4trts,Rst_Wet_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(dry_rsl_fric_empd_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_empd_fric_sem_coefs.csv",sep=",",row.names=F)
write.table(dry_rsl_fric_empd_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_empd_fric_model_fits.csv",sep=",",row.names=F)


#################
## Model 6    ###
## FDis $ ePSE ##
#################

Rsl_Dry_ModList=list(
  lme(Ye_ePSE~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2,correlation=x1,control=bb,data=rs_122),
  lme(Ye_FRic4~lg2SppN,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122),
  lme(lg2Rsl12~Ye_PCAcwm4trts+lg2SppN+Ye_ePSE+Ye_FRic4,random=~1+lg2SppN|Site/PlotUnique2, correlation=x1,control=bb,data=rs_122)
  
)


sem.fit(Rsl_Dry_ModList,rs_122,corr.errors=c("Ye_ePSE~~Ye_FRic4"),conditional=T,
        model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

#Naive Model

dry_rsl_fric_epse_pc<- sem.coefs(Rsl_Dry_ModList,rs_122,standardize="scale")
dry_rsl_fric_epse_pc$Climate_Bin<-"Dry"

dry_rsl_fric_epse_modfit<-sem.model.fits(Rsl_Dry_ModList)
dry_rsl_fric_epse_modfit$ResponseVars<-c("Ye_ePSE","FRic4","Resilience")
dry_rsl_fric_epse_modfit$PredVars<-c("lg2SppN","lg2SppN","F-S,Ye_ePSE,FRic4,lg2SppN")
dry_rsl_fric_epse_modfit$Climate_Bin<-"Dry"


resids.df1<-partial.resid(lg2Rst12~Yn_ePSE,Rst_Wet_ModList,data=rs_122,
                          model.control = list(lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")))

write.table(dry_rsl_fric_epse_pc,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_epse_fric_sem_coefs.csv",sep=",",row.names=F)
write.table(dry_rsl_fric_epse_modfit,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Rsl_DRY_epse_fric_model_fits.csv",sep=",",row.names=F)

#DONE
