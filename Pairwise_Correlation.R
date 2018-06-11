###############################################
# Correlations among facets of biodiversity   #
# estimated using multilevel meta-analytical  #
# regression models                           #
###############################################

require(dplyr)
require(psych)
require(reshape2)
require(metafor)
require(corrplot)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

# Data
stab<-read.delim("data.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # eliminate site where we didn't have good trait coverage

stab_4<-select(stab,Site,Study_length,UniqueID,SppN,eMPD,eMNTD,ePSE,sMPD,sMNTD,FDis4,FRic4,PCAdim1_4trts,SLA, LDMC, LeafN, LeafP,
               Plot_TempStab,Plot_Biomassxbar, Plot_Biomasssd,Gross_synchrony)

# for plots with ONLY 1 spp, we assume that a species #is perfectly synchronized with itself

stab_4$Gross_synchrony<-ifelse(is.na(stab_4$Gross_synchrony)==TRUE,1,stab_4$Gross_synchrony) 

# convert synchrony metrics to different scale

stab_4$GrossAsynchrony_s<-stab_4$Gross_synchrony*-1

# further adjustments

stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

# Filter out NAs for Asynchrony and  FRic4
stab_444<-filter(stab_4, is.na(PlotAsynchrony_s)==FALSE)
stab_444<-filter(stab_444, is.na(FRic4)==FALSE)

#################################
# calculate correlations among  #
# all predictors of TS ##########
#################################

stab_corr<-select(stab_444,Site,SppN, eMNTD,eMPD,FDis4,FRic4,PCAdim1_4trts,GrossAsynchrony_s,Plot_Biomassxbar,Plot_Biomasssd)

n<-length(unique(stab_corr$Site))

outt=c();

for(i in 1:n){
  
  test=subset(stab_corr, stab_corr$Site==(unique(stab_corr$Site))[i]) 
  
  Site<-as.character(unique(test$Site))
  Plot_n<-dim(test)[1]
  test<-select(test,-Site)
  
  p_out<-corr.test(test,method="pearson")
  p_outt<-as.matrix(p_out$r)
  p_outt<-data.frame(Var1=rownames(p_outt)[row(p_outt)[upper.tri(p_outt)]], 
                     Var2=colnames(p_outt)[col(p_outt)[upper.tri(p_outt)]], 
             corr=p_outt[upper.tri(p_outt)])
  
  p_outt$Site<-Site
  p_outt$Plot_n<-Plot_n

    outt[[i]]<-rbind.data.frame(p_outt)
  
}

jjj<-do.call(rbind,outt)  
jjj<-data.frame(jjj)

jjj<-filter(jjj, is.na(corr)==FALSE)

#Calculate effect size (raw correlation coefficient)

effect_size<-select(jjj,Site, Plot_n,Var1,Var2, corr)
effect_size<-escalc(measure="COR",ri=corr,ni=Plot_n,data=effect_size)

######### 

effect_size$Combn<-paste(effect_size$Var1,effect_size$Var2,sep="_")

n<-length(unique(effect_size$Combn))

outt=c();

for(i in 1:n){
  
  test=subset(effect_size, effect_size$Combn==(unique(effect_size$Combn))[i]) 

  Combn<-as.character(unique(test$Combn))
  Var1<-as.character(unique(test$Var1))
  Var2<-as.character(unique(test$Var2))
  
    Mod1<-rma.uni(yi,vi,measure="GEN",test="knha",method="REML",data=test)

  outt_p<-cbind.data.frame(Var1, Var2,Combn,Mod1$b,Mod1$ci.lb,Mod1$ci.ub)
  outt[[i]]<-rbind.data.frame(outt_p)
  
}

jjj<-do.call(rbind,outt)  
jjj<-data.frame(jjj)

colnames(jjj)[4]<-"r"
colnames(jjj)[5]<-"lower95"
colnames(jjj)[6]<-"upper95"

jjj<-arrange(jjj,Var1,Var2)


write.table(jjj,"Div_Corr_Effsizes.csv",sep=",",row.names=F)