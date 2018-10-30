overallSNA<- as.matrix(read.csv("/Users/lt/Desktop/manuscripts/03 deep SNA/analysis process/OVERALL_onemode/OVERALL01.csv", row.names=1))
gplot(overallSNA,gmode="digraph", displaylabels=TRUE,label.cex=0.8,vertex.col="darkolivegreen")
overallnet=network(overallSNA)
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

#plot
plot(id, od, type="n", xlab="Indegree", ylab="Outdegree")
abline(0, 1, lty=3,col = "gray60")
text(jitter(id), jitter(od), network.vertex.names(overallnet), cex=0.7,
     col="black") 

plot(clo~bet)

#plot betweenness and closeness
plot(bet, clo, type="n", xlab="Betweenness", ylab="Closeness")
abline(0, 1, lty=3,col = "gray60")
text(jitter(bet), jitter(clo), network.vertex.names(overallnet), cex=0.7,
     col="black")

par(mfrow=c(1,1)) # Setup a 2 panel plot (for later)

#plot node size and color based on indegree
gplot(overallnet, vertex.cex=(id)^0.5/2, gmode="graph",
      boxed.labels=FALSE,label.cex=0.6, label.pos=5, label.col="grey17",
      vertex.col=rgb(id/max(id),0,id/max(id)),edge.col="grey17",
      label=network.vertex.names(overallnet),edge.lwd=overallSNA/2,mode = "fruchtermanreingold")

#plot node size and color based on outdegree
gplot(overallnet, vertex.cex=(od)^0.5/2, gmode="graph",
      boxed.labels=FALSE,label.cex=0.6, label.pos=5, label.col="grey17",
      vertex.col=rgb(od/max(od),0,od/max(od)),edge.col="grey17",
      label=network.vertex.names(overallnet),edge.lwd=overallSNA/2,mode = "fruchtermanreingold")

#plot node size and color based on in/out degree
gplot(overallnet, vertex.cex=(id+od)^0.5/2, gmode="graph",
      boxed.labels=FALSE,label.cex=0.7, label.pos=5, label.col="grey17",
      vertex.col=rgb((id+od)/max(id+od),0,(id+od)/max(id+od)),edge.col="grey17",
      label=network.vertex.names(overallnet),edge.lwd=overallSNA/2,mode = "fruchtermanreingold")

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
