require(dplyr)
require(psych)
require(dcast2)
require(metafor)
require(corrplot)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))


# Data
stab<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,sMPD, sMNTD,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony,CV_Temp)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself

stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)

#stab_44<-filter(stab_4,SppN>1) eliminate monocultures

stab_44<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony

stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony

#################################
# calculate correlations among  #
# all predictors of TS ##########
#################################

stab_corr<-select(stab_444,Site,SppN, eMNTD,eMPD,FDis4,FRic4,PCAdim1_4trts,Plot_Asynchrony)

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

jjj<-na.omit(jjj)


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


write.table(jjj,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Div_Corr_Effsizes.csv",sep=",",row.names=F)

############################
# make correlation matrix  #
############################

require(reshape2)
require(viridis)

jjj<-read.delim("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Div_Corr_Effsizes.csv",sep=",",header=T)

jjj<-filter(jjj,Var1=="SppN"|Var1=="FDis4" | Var1=="FRic4"|Var1=="eMNTD"|Var1=="PCAdim1_4trts" )
jjj<-filter(jjj,Var2=="SppN"|Var2=="FDis4" | Var2=="FRic4"|Var2=="eMNTD"|Var2=="PCAdim1_4trts")

corr_mat<-dcast(jjj,Var1~Var2,value.var="r",mean)

corr_mat<-arrange(corr_mat,-eMNTD)
rownames(corr_mat)<-corr_mat$Var1

corr_mat$Var1<-as.character(corr_mat$Var1)

corr_mat$Var1<-ifelse(corr_mat$Var1=="PCAdim1_4trts","Fast-slow",corr_mat$Var1)
#corr_mat$Var1<-ifelse(corr_mat$Var1=="Plot_Asynchrony","Synchrony",corr_mat$Var1)

#corr_mat<-select(corr_mat,-Plot_Asynchrony)

colnames(corr_mat)[5]<-"Fast-slow"
#colnames(corr_mat)[6]<-"Synchrony"

corr_mat$Var1<-NULL

corr_mat[is.nan(corr_mat)] <- 0
corr_mat<-as.matrix(corr_mat)

col<- colorRampPalette(c("red", "white", "blue"))(256)

col2<-plasma(256)

##################

png(filename="/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Div_Corr.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=6 , 
    pointsize=2, 
    res=200)


corrplot(corr_mat, method="ellipse",type="upper",col=col,is.corr=TRUE,diag=TRUE,bg="white",tl.pos=TRUE,tl.cex=6,tl.col="black",tl.srt=0,cl.cex=6)

dev.off()

cairo_ps("/homes/dc78cahe/Dropbox (iDiv)/Research_projects/leipzigPhyTrt/StabilityII_data/Community_Level/Div_Corr.eps",
         family="sans",
         height=6,width=6,
         bg="white")

corrplot(corr_mat, method="ellipse",type="upper",col=col,is.corr=TRUE,diag=TRUE,bg="white",tl.pos=TRUE,tl.col="black",tl.srt=0)

dev.off()