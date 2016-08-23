rm(list=ls()) 
#####################
require(dplyr)

devtools::install_github('rich-iannone/DiagrammeR')
require(DiagrammeR)
require(stringr)
require(semPlot)
#require(lavaan)

mod1<-read.delim("/home/dylan/Dropbox/leipzigPhyTrt/StabilityII_data/Community_Level/TS_empd_fdis_sem_coefs_NEW.csv",sep=",",header=TRUE)
#mod11<-filter(mod1,response=="TS_lg2")
#mod11$Model<-"eMPD"
#mod11$ModelGroup<-"FD"



# Describe edges

paths <- mod1 %>%
  select(response,predictor, estimate,p.value)

#paths<-paths[1:8,]

paths$response<-as.character(paths$response)
paths$predictor<-as.character(paths$predictor)


paths$arrowtail<-grepl("~~ ", paths$response)
paths$arrowtail<-ifelse(paths$arrowtail==FALSE,"none","normal")
paths$response<-str_replace_all(paths$response, "~~ ", "")
paths$predictor<-str_replace_all(paths$predictor, "~~ ", "")


paths$response<-ifelse(paths$response=="PCAdim1_4trts","F-S",paths$response)
paths$response<-ifelse(paths$response=="eMPD","PD",paths$response)
paths$response<-ifelse(paths$response=="FDis4","FD",paths$response)
paths$response<-ifelse(paths$response=="Plot_Asynchrony","&#951;",paths$response)
paths$response<-ifelse(paths$response=="TS_lg2","1/CV",paths$response)
paths$predictor<-ifelse(paths$predictor=="lg2SppN","SR",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="Plot_Asynchrony","&#951;",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="PCAdim1_4trts","F-S",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="FDis4","FD",paths$predictor)
paths$predictor<-ifelse(paths$predictor=="eMPD","PD",paths$predictor)


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
               label= paths$labels,
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
    "penwidth", 0.5, 2, "penwidth")

render_graph(graph)


