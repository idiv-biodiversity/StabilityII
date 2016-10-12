rm(list=ls()) 
#require(devtools)
#install_github("jslefche/piecewiseSEM")

require(dplyr)
require(piecewiseSEM)
library(semPlot)
library(lmerTest)
library(nlme)

# Data
stab<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/Stab_Stability_FD_PD_CWM_PlotYearAverages_V.csv",sep=",",header=T)

stab<-filter(stab,Site!="BIODEPTH_GR")  # should get rid of site where we didn't have good trait coverage

stab_4<-select(stab,Site,UniqueID,SppN,eMPD,eMNTD,ePSE,FDis4,FRic4,PCAdim1_4trts,Plot_TempStab,Plot_Asynchrony,annualTemp,meanPrecip,meanPET,CV_Temp,CV_Precip)

stab_4$Plot_Asynchrony<-ifelse(stab_4$SppN==1 & is.na(stab_4$Plot_Asynchrony)==TRUE,1,stab_4$Plot_Asynchrony) # for monocultures, we assume that a species
#is perfectly synchronized with itself


stab_4$SppN<-as.numeric(stab_4$SppN)

stab_4$TS_lg2<-log(stab_4$Plot_TempStab,base=2)

stab_4$lg2SppN <- log(stab_4$SppN,2)


stab_444<-stab_4[!is.na(stab_4$Plot_Asynchrony),]  # no NAs for Plot Asynchrony
stab_444<-stab_444[!is.na(stab_444$FRic4),]  # no NAs for Plot Asynchrony


####################
# Prelim       #####
####################

bb<-lmeControl(msMaxIter=0,msVerbose = TRUE,opt="optim",maxIter=100,optimMEthod="L-BFGS-B")  ######## "msMaxIter=0" is important in here!!!
cc<-lmeControl(opt="optim")


a<-lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
a1<-lme(eMNTD~lg2SppN,random=~1|Site,control=cc,data=stab_444)

AIC(a,a1)

plot(a1)
qqnorm(a1)


b<-lme(FDis4~lg2SppN,random=~1|Site, control=cc,data=stab_444)
b1<-lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)
AIC(b,b1)

plot(b1)
qqnorm(b1)


c<-lme(Plot_Asynchrony~lg2SppN,random=~1|Site,control=cc,data=stab_444)
c1<-lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444)


AIC(c,c1)

plot(c1)
qqnorm(c1)

c<-lme(Plot_Asynchrony~FDis4,random=~1|Site,control=cc,data=stab_444)
c1<-lme(Plot_Asynchrony~FDis4,random=~1+lg2SppN|Site,control=cc,data=stab_444)
c2<-lme(Plot_Asynchrony~FDis4,random=~1+FDis4|Site/SppN,control=cc,data=stab_444)

AIC(c,c1,c2)

plot(c2)
qqnorm(c2)

d<-lme(Plot_Asynchrony~PCAdim1_4trts,random=~1|Site,control=cc,data=stab_444)
d1<-lme(Plot_Asynchrony~PCAdim1_4trts,random=~1+lg2SppN|Site,control=cc,data=stab_444)
d2<-lme(Plot_Asynchrony~PCAdim1_4trts,random=~1+PCAdim1_4trts|Site/SppN,control=cc,data=stab_444)


AIC(d,d1,d2)

plot(d1)
qqnorm(d1)


e<-lme(Plot_Asynchrony~eMNTD,random=~1|Site,control=cc,data=stab_444)
e1<-lme(Plot_Asynchrony~eMNTD,random=~1+lg2SppN|Site,control=cc,data=stab_444)
e2<-lme(Plot_Asynchrony~eMNTD,random=~1+eMNTD|Site/SppN,control=cc,data=stab_444)

AIC(e,e1,e2)

plot(e2)
qqnorm(e2)

g<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1|Site,control=cc,data=stab_444)
g1<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+lg2SppN|Site,control=cc,data=stab_444)
g2<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+eMNTD|Site,control=cc,data=stab_444)
g3<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+FDis4|Site,control=cc,data=stab_444)
g4<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+lg2SppN+PCAdim1_4trts|Site,control=cc,data=stab_444)
g5<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+lg2SppN+eMNTD|Site,control=cc,data=stab_444)
g6<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+lg2SppN+FDis4|Site,control=cc,data=stab_444)
g7<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+PCAdim1_4trts+FDis4|Site,control=cc,data=stab_444)
g8<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+PCAdim1_4trts+eMNTD|Site,control=cc,data=stab_444)
g9<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+PCAdim1_4trts+FDis4|Site,control=cc,data=stab_444)
g10<-lme(Plot_Asynchrony~eMNTD+lg2SppN+FDis4+PCAdim1_4trts,random=~1+PCAdim1_4trts+FDis4|Site,control=cc,data=stab_444)


AICc(g,g1,g2,g3,g4)

plot(e2)
qqnorm(e2)



f<-lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_444)
f1<-lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN,random=~1|Site, control=cc,data=stab_444)
f3<-lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN,random=~1|Site/SppN, control=cc,data=stab_444)



AIC(f,f1,f3)

######################
# piecewise SEM  #####
# without synchrony  #
######################


modList=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  #lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(modList,stab_444,corr.errors=c("eMNTD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

sem.fit(modList,stab_444,corr.errors=c("eMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) 

ts_emntd<-sem.coefs(modList,stab_444,standardize="scale",corr.errors=c("eMNTD~~FDis4","FDis4 ~~ PCAdim1_4trts"))


mf_ts_emntd<-sem.model.fits(modList)
mf_ts_emntd$ResponseVars<-c("eMNTD","FDis4","Temp_Stability")
mf_ts_emntd$PredVars<-c("lg2SppN","lg2SppN","F-S,eMNTD,FDis4,lg2SppN")


resids.df1<-partial.resid(TS_lg2~lg2SppN,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df2<-partial.resid(TS_lg2~eMNTD,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df3<-partial.resid(TS_lg2~FDis4,modList,data=stab_444,list(lmeControl(opt="optim")))
resids.df5<-partial.resid(TS_lg2~PCAdim1_4trts,modList,data=stab_444,list(lmeControl(opt="optim")))


write.table(ts_emntd,"/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_emntd_fdis_sem_coefs_noSYNC.csv",sep=",",row.names=F)

#############
# Visualize #
#############

#devtools::install_github('rich-iannone/DiagrammeR')
require(DiagrammeR)
require(stringr)
require(semPlot)


# Describe edges

paths <- ts_emntd %>%
  select(response,predictor, estimate,p.value)

#paths<-paths[1:8,]

paths$response<-as.character(paths$response)
paths$predictor<-as.character(paths$predictor)


paths$arrowtail<-grepl("~~ ", paths$response)
paths$arrowtail<-ifelse(paths$arrowtail==FALSE,"none","normal")
paths$response<-str_replace_all(paths$response, "~~ ", "")
paths$predictor<-str_replace_all(paths$predictor, "~~ ", "")


paths$response<-ifelse(paths$response=="PCAdim1_4trts","F-S",paths$response)
paths$response<-ifelse(paths$response=="eMNTD","PD",paths$response)
paths$response<-ifelse(paths$response=="FDis4","FD",paths$response)
#paths$response<-ifelse(paths$response=="Plot_Asynchrony","&#951;",paths$response)
paths$response<-ifelse(paths$response=="TS_lg2","1/CV",paths$response)
paths$predictor<-ifelse(paths$predictor=="lg2SppN","SR",paths$predictor)
#paths$predictor<-ifelse(paths$predictor=="Plot_Asynchrony","&#951;",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="PCAdim1_4trts","F-S",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="FDis4","FD",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="eMNTD","PD",paths$predictor)


#paths$predictor<-as.character(paths$predictor)
#paths$predictor<-str_replace_all(paths$predictor, "~~ ", "")

paths$color<-ifelse(paths$p.value<0.05,"black","gray50")

paths$style<-ifelse(paths$estimate>0,"solid","dashed")
paths$weight<-abs(paths$estimate)

colnames(paths)[1]<-"to"
colnames(paths)[2]<-"from"

paths$rel<-"leading_to"

#s<-round(paths$values,2)

paths$labels<-round(paths$estimate,2)
#paths<-select(paths,from, to, rel, weight,style,labels,color)

paths$labels<-as.character(paths$labels)

paths$dir<-ifelse(paths$arrowtail=="none","forward","both")

###Create nodes
nn1<-data.frame(unique(paths$from))
colnames(nn1)[1]<-"nodes"
nn2<-data.frame(unique(paths$to))
colnames(nn2)[1]<-"nodes"

nodes<-rbind.data.frame(nn1,nn2)
nodes<-unique(nodes)

nodess<-create_nodes(nodes = nodes$nodes,
                     label = as.character(nodes$nodes),
                     type = "lower",
                     style = "empty",
                     color = "black",
                     shape = c("rectangle"))

##Create edgges

edgess <-
  create_edges(from = paths$from,
               to = paths$to,
               rel = paths$rel,
               color = paths$color,
               style = paths$style,
               #label= paths$labels,
               arrowtail=paths$arrowtail,
               arrowhead="normal", dir=paths$dir,
               penwidth = paths$weight)

# create graph

my_graph <- create_graph(
  nodes_df = nodess, 
  edges_df = edgess, 
  edge_attrs= "fontsize = 6",
  graph_attrs = c("layout = circo"),
  directed=TRUE)

graph <-
  my_graph %>%
  #rescale_edge_attrs(
  #"weight", "gray80", "gray20", "color") %>%
  rescale_edge_attrs(
    "penwidth", 0.5, 4, "penwidth")

render_graph(graph)

######################
# hypothetical model:
# sans synchrony ####
#####################


hyp=list(
  lme(eMNTD~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(FDis4~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  #lme(Plot_Asynchrony~lg2SppN,random=~1+lg2SppN|Site,control=cc,data=stab_444),
  lme(TS_lg2~PCAdim1_4trts+FDis4+eMNTD+lg2SppN,random=~1+lg2SppN|Site, control=cc,data=stab_444)
)


sem.fit(hyp,stab_444,corr.errors=c("eMNTD~~FDis4"),conditional=T,
        model.control = list(lmeControl(opt = "optim"))) #naive model

hyp<-sem.coefs(hyp,stab_444,standardize="scale",corr.errors=c("eMNTD~~FDis4"))
hyp$estimate<-0.5

paths <- hyp %>%
  select(response,predictor, estimate,p.value)

#paths<-paths[1:8,]

paths$response<-as.character(paths$response)
paths$predictor<-as.character(paths$predictor)


paths$arrowtail<-grepl("~~ ", paths$response)
paths$arrowtail<-ifelse(paths$arrowtail==FALSE,"none","normal")
paths$response<-str_replace_all(paths$response, "~~ ", "")
paths$predictor<-str_replace_all(paths$predictor, "~~ ", "")


paths$response<-ifelse(paths$response=="PCAdim1_4trts","F-S",paths$response)
paths$response<-ifelse(paths$response=="eMNTD","PD",paths$response)
paths$response<-ifelse(paths$response=="FDis4","FD",paths$response)
#paths$response<-ifelse(paths$response=="Plot_Asynchrony","&#951;",paths$response)
paths$response<-ifelse(paths$response=="TS_lg2","1/CV",paths$response)
paths$predictor<-ifelse(paths$predictor=="lg2SppN","SR",paths$predictor)
#paths$predictor<-ifelse(paths$predictor=="Plot_Asynchrony","&#951;",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="PCAdim1_4trts","F-S",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="FDis4","FD",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="eMNTD","PD",paths$predictor)


#paths$predictor<-as.character(paths$predictor)
#paths$predictor<-str_replace_all(paths$predictor, "~~ ", "")

#paths$color<-ifelse(paths$p.value<0.05,"black","gray50")

#paths$style<-ifelse(paths$estimate>0,"solid","dashed")
paths$weight<-paths$estimate

colnames(paths)[1]<-"to"
colnames(paths)[2]<-"from"

paths$rel<-"leading_to"
paths$dir<-ifelse(paths$arrowtail=="none","forward","both")

###Create nodes
nn1<-data.frame(unique(paths$from))
colnames(nn1)[1]<-"nodes"
nn2<-data.frame(unique(paths$to))
colnames(nn2)[1]<-"nodes"

nodes<-rbind.data.frame(nn1,nn2)
nodes<-unique(nodes)

nodess<-create_nodes(nodes = nodes$nodes,
                     label = as.character(nodes$nodes),
                     type = "lower",
                     style = "empty",
                     color = "black",
                     shape = c("rectangle"))

##Create edgges

edgess <-
  create_edges(from = paths$from,
               to = paths$to,
               rel = paths$rel,
               color = "black",
               style ="solid",
               #label= paths$labels,
               arrowtail=paths$arrowtail,
               arrowhead="normal", dir=paths$dir,
               penwidth = paths$weight)

# create graph

my_graph <- create_graph(
  nodes_df = nodess, 
  edges_df = edgess, 
  edge_attrs= "fontsize = 6",
  graph_attrs = c("layout = circo"),
  directed=TRUE)

render_graph(my_graph)

