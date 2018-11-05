####### clean data #########
all<- read.csv("/Users/fanouyang/Desktop/sna_2018fall/all_data.csv")
all<-all %>% select(vert1_id,vert2_id)  

library(reshape2)
all_matrix<-all %>% dcast(vert1_id~vert2_id)
write.csv(all_matrix,file = "/Users/fanouyang/Desktop/sna_2018fall/all_matrix.csv") 
# manipulate on local and read in again, this is the final correct matrix
all_matrix<- read.csv("/Users/fanouyang/Desktop/sna_2018fall/all_matrix.csv",row.names=1)

all_edge <- melt(all_matrix)
write.csv(all_edge,file = "/Users/fanouyang/Desktop/sna_2018fall/all_edge.csv") 

# manipulate on local and read in again, this is the final correct edge list
all_edge<- read.csv("/Users/fanouyang/Desktop/sna_2018fall/all_edge.csv")

####### use package sna, network ##########
library(sna)
# plot the network, not pretty
gplot(all_matrix,gmode="digraph", displaylabels=TRUE,label.cex=0.8,vertex.col="darkolivegreen")

# convert matrix to network format for further analysis
overallnet=network(all_matrix)
overallnet

##################################
###### node level analysis###########
##################################
#calculate in/out degree
id<-degree(overallnet,gmode="digraph",cmode="indegree")
od<-degree(overallnet,gmode="digraph",cmode="outdegree")
id
od
id+od
#betweenness
bet1=betweenness(overallnet,rescale=T)
bet2=betweenness(overallnet)
bet1
bet2
#closeness
clo1=closeness(overallnet,rescale=T)
clo2=closeness(overallnet)
clo1
clo2
#Eigenvector
eig1=evcent(overallnet,rescale=T)
eig2=evcent(overallnet)
eig1
eig2

#plot node size and color based on in/out degree
gplot(overallnet, vertex.cex=(id+od)^0.5/2, gmode="graph",
      boxed.labels=FALSE,label.cex=0.7, label.pos=5, label.col="grey17",
      vertex.col=rgb((id+od)/max(id+od),0,(id+od)/max(id+od)),edge.col="grey17",
      label=network.vertex.names(overallnet),edge.lwd=all_matrix/2,mode = "fruchtermanreingold")

#sociomatrix of the overall network
plot.sociomatrix(overallnet, diaglab = FALSE, 
                 main = "The overall interaction", cex.lab = 0.4, asp = 0.5)


### function that takes a list of network objects and outputs classic ###
### centrality scores by nodes and ranks them ###
node.centrality<-function(list.net){
  lapply(list.net,
         function(z){
           x<-z
           central.nodes<-cbind(degree(x,cmode="indegree"), degree(x,cmode="outdegree"),evcent(x,rescale=T),betweenness(x,rescale=T),
                                closeness(x,cmode="directed",rescale=T))
           colnames(central.nodes)<-c("idegree","odegree","eigen","betweenness","closeness")
           rownames(central.nodes)<-z%v%"vertex.names"
           
           o1<-order(central.nodes[,1],decreasing =TRUE)
           o2<-order(central.nodes[,2],decreasing =TRUE)
           o3<-order(central.nodes[,3],decreasing =TRUE)
           o4<-order(central.nodes[,4],decreasing =TRUE)
           o5<-order(central.nodes[,5],decreasing =TRUE)
           list(ranking=central.nodes,order=cbind(o1,o2,o3,o4,o5))
         })
}
### Compute Centrality scores
nc<-node.centrality(list(overallnet))
nc


###########################################
######## graph level analysis  ######
###########################################
centralization(overallnet,degree)
centralization(overallnet,degree,cmode="outdegree")
centralization(overallnet,degree,cmode="indegree")
centralization(overallnet, betweenness)
centralization(overallnet, closeness)
centralization(overallnet, evcent)

network.size(overallnet)
gden(overallnet,mode="graph") #density
degree(overallnet)
mean(degree(overallnet)) #average degree
sum(id)/20 #average in degree
sum(od)/20 #average out degree

#transitivity(overallnet)
dyad.census(overallnet)
network.dyadcount(overallnet, na.omit = F)
network.edgecount(overallnet, na.omit = F)
grecip(overallnet, measure = "edgewise")
grecip(overallnet, measure = "dyadic")
grecip(overallnet, measure = "dyadic.nonnull")
gtrans(overallnet) #transitivity

hierarchy(overallnet, measure = "reciprocity")
hierarchy(overallnet, measure = "krackhardt")

# Returns the number of components within network
components(overallnet,connected="weak")
#returns the Krackhardt connectedness
connectedness(overallnet, g=NULL)

# geodist, how to explain? find the number and lengths of geodesics between all nodes
geodist(overallnet, inf.replace=Inf, count.paths=TRUE, predecessors=FALSE,
        ignore.eval=TRUE)
geo=geodist(overallnet)
max(geo$gdist)  #diameter

## Function for average path length
averagePathLength<-function(net){
  if(!is.network(net)){stop("Not a Network")}
  gd<-geodist(net)
  if(net%n%"directed"){
    return((1/choose(network.size(net),2))*sum(gd$gdist))
  }
  (1/(2*choose(network.size(net),2)))*sum(gd$gdist)
}

## Function for diameter
diameter<-function(net){
  gd<-geodist(net)
  max(gd$gdist)
}

## Compute Average Path Length
averagePathLength(overallnet)

## Compute Diameter
diameter(overallnet)

########### you can use tnet package to calculate those metrics ##########

########## igraph format analysis ##########
##### this is a good example http://kateto.net/networks-r-igraph #######
library(igraph)

# final edge list 
all_edge=read.csv("/Users/fanouyang/Desktop/sna_2018fall/all_edge.csv") 

# create igraph from the edge list
all_igraph <- graph_from_data_frame(d=all_edge, directed=T)
# cheeck node and edge
E(all_igraph)     
V(all_igraph)
# plot not pretty,  need further revise it later
plot(all_igraph, edge.arrow.size=.1)
## node-level analysis
igraph::degree(all_igraph)
igraph::degree(all_igraph,mode="out")
igraph::betweenness(all_igraph)
igraph::closeness(all_igraph, mode="in")
igraph::closeness(all_igraph, mode="out")
igraph::closeness(all_igraph, mode="all")
igraph::eigen_centrality(all_igraph)
## network level analysis
igraph::diameter(all_igraph, directed = TRUE)
dyad_census(all_igraph)
distances(all_igraph)
shortest_paths(all_igraph, 5)
centr_degree(all_igraph)$centralization
centr_clo(all_igraph, mode="all")$centralization
centr_eigen(all_igraph, directed=FALSE)$centralization

#### use visNetwork package########
# plot you should customize it later
library(visNetwork)
visIgraph(all_igraph, idToLabel = TRUE, layout = "layout_nicely",
          physics = FALSE, smooth = T)



