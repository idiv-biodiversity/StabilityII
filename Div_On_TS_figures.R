###############################
# Temp. stability: div effects#
#    on biomass (xbar & sd)####
###############################
rm(list=ls()) 

require(dplyr)
require(ggplot2)

########
# data #
########
div_ts<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Div_ts.csv",sep=",",header=T)



######## using unstandardized regression coefficients
  
## SppN
g1 <- ggplot(div_ts,aes(x=biom_slope,y=sd_slope))+ 
    geom_point(pch=21,alpha=0.5,fill="black",size=div_ts$Study_length)+
    geom_abline(intercept = 0, slope = 1,colour="black",linetype=3)+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    labs(x=expression(paste(beta [mu])),y =expression(paste(beta [sd])))
  
  tempstab1<-g1+ theme(axis.title.x=element_text(colour="black",face="bold",size=15,vjust=-1),
                      axis.title.y=element_text(colour="black",face="bold",size=15,vjust=2),
                      axis.text.y=element_text(colour="black",face="bold",size=12),
                      axis.text.x=element_text(colour="black",face="bold",size=12),
                      strip.text.x=element_text(colour="black",face="bold",size=12), 
                      panel.background =element_rect(fill="transparent",colour="black"),
                      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),panel.border=element_rect(fill=NA,colour="black"))
  

### FD
  
  g2 <- ggplot(div_ts,aes(x=fd_biom_slope,y=fd_sd_slope))+ 
    geom_point(pch=21,alpha=0.5,fill="black",size=div_ts$Study_length)+
    geom_abline(intercept = 0, slope = 1,colour="black",linetype=3)+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    labs(x=expression(paste(beta [mu])),y =expression(paste(beta [sd])))
  
  tempstab2<-g2+ theme(axis.title.x=element_text(colour="black",face="bold",size=15,vjust=-1),
                       axis.title.y=element_text(colour="black",face="bold",size=15,vjust=2),
                       axis.text.y=element_text(colour="black",face="bold",size=12),
                       axis.text.x=element_text(colour="black",face="bold",size=12),
                       strip.text.x=element_text(colour="black",face="bold",size=12), 
                       panel.background =element_rect(fill="transparent",colour="black"),
                       plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),panel.border=element_rect(fill=NA,colour="black"))
  
  
  
  png(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Div_onTS_standard.png", 
      type="cairo",
      units="in", 
      width=4, 
      height=4, 
      pointsize=2, 
      res=200)
  
  
  cairo_ps(filename="/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Div_onTS_standard.eps",
           family="sans",
           height=4,width=4,
           bg="white")
  
  
  
  tempstab1
  
  dev.off()
  
  
