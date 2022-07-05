# clean up
rm(list = ls())

# load the igraph library
library('igraph')

# load data
load('LAB 2/lab2_networks.Rdata')

# show graphs
advice
friend
report

# create the matrices
adviceY <- get.adjacency(advice, sparse = FALSE)
diag(adviceY) <- NA
friendY <- get.adjacency(friend, sparse = FALSE)
diag(friendY) <- NA
reportY <- get.adjacency(report, sparse = FALSE)
diag(reportY) <- NA


V(advice)
V(friend)
V(report)

E(advice)
E(friend)
E(report)

vcount(advice)
vcount(friend)
vcount(report)

ecount(advice)
ecount(friend)
ecount(report)

is.directed(advice)
is.directed(friend)
is.directed(report)

# graphical representaton of the network
par(mfrow = c(1,3), mar = c(0,0,2,0))
plot(advice, vertex.size = 20,  edge.arrow.size = 0.3,
     vertex.label.cex = 1.5, main = "Advice")
box()
plot(friend, vertex.size = 20,  edge.arrow.size = 0.3,
     vertex.label.cex = 1.5, main = "Friendship")
box()
plot(report, vertex.size = 20,  edge.arrow.size = 0.3,
     vertex.label.cex = 1.5, main = "Report")
box()

# density calculated manually
ecount(advice)/(vcount(advice)*(vcount(advice)-1))
ecount(friend)/(vcount(friend)*(vcount(friend)-1))
ecount(report)/(vcount(report)*(vcount(report)-1))

# same density using igraph function density
graph.density(advice)
graph.density(friend)
graph.density(report)

# density on undirected graphs
adviceUn = as.undirected(advice, mode = 'mutual')
friendUn = as.undirected(friend, mode = 'mutual')
reportUn = as.undirected(report, mode = 'mutual')

# density calculated manually
(ecount(adviceUn)*2)/(vcount(adviceUn)*(vcount(adviceUn)-1))
(ecount(friendUn)*2)/(vcount(friendUn)*(vcount(friendUn)-1))
(ecount(reportUn)*2)/(vcount(reportUn)*(vcount(reportUn)-1))

# same density using igraph function density
graph.density(adviceUn)
graph.density(friendUn)
graph.density(reportUn)


# Is there a tendency for the nodes in a directed network to return relations?
# Let's compute dyads
dyad.census(advice)
dyad.census(friend)
dyad.census(report)

unlist(dyad.census(advice))
unlist(dyad.census(friend))
unlist(dyad.census(report))

# there is a different counting between census and matrix calc
sum(adviceY * t(adviceY), na.rm = TRUE)
sum(friendY * t(friendY), na.rm = TRUE)
sum(reportY * t(reportY), na.rm = TRUE)

sum(adviceY * t(adviceY), na.rm = TRUE)/2
sum(friendY * t(friendY), na.rm = TRUE)/2
sum(reportY * t(reportY), na.rm = TRUE)/2

# reciprocity
sum(adviceY * t(adviceY), na.rm = TRUE)/sum(adviceY, na.rm = TRUE)
sum(friendY * t(friendY), na.rm = TRUE)/sum(friendY, na.rm = TRUE)
sum(reportY * t(reportY), na.rm = TRUE)/sum(reportY, na.rm = TRUE)

reciprocity(advice)
reciprocity(friend)
reciprocity(report)

# transitivity
rowNames <- c('empty',
  'A->B, C',
  'A<->B, C',
  'A<-B->C',
  'A->B<-C',
  'A->B->C',
  'A<->B<-C',
  'A<->B->C',
  'A->B<-C, A->C',
  'A<-B<-C, A->C',
  'A<->B<->C',
  'A<-B->C, A<->C',
  'A->B<-C, A<->C',
  'A->B->C, A<->C',
  'A->B<->C, A<->C',
  'A<->B<->C, A<->C')


adviceTr <- triad.census(advice)
friendTr <- triad.census(friend)
reportTr <- triad.census(report)
df <- data.frame(rowNames, row.names = 1, adviceTr, friendTr, reportTr)
colnames(df) <- c('Advice', 'Friend', 'Report')
df

transitivity(advice)
transitivity(friend)
transitivity(report)

transitivity(as.undirected(advice, mode = 'collapse'))
transitivity(as.undirected(friend, mode = 'collapse'))
transitivity(as.undirected(report, mode = 'collapse'))

# normalize in order to compare
# use log-odds ratio

# density
graph.density(advice)/(1-graph.density(advice))
graph.density(friend)/(1-graph.density(friend))
graph.density(report)/(1-graph.density(report))

# transitivity
transitivity(advice)/(1-transitivity(advice))
transitivity(friend)/(1-transitivity(friend))
transitivity(report)/(1-transitivity(report))

# assortativity
attrs = read.table('LAB 2/attributes.txt',  sep = '', head = TRUE)

V(advice)$Dept <- attrs$Dept
V(advice)$Age <- attrs$Age

V(friend)$Dept <- attrs$Dept
V(friend)$Age <- attrs$Age

V(report)$Dept <- attrs$Dept
V(report)$Age <- attrs$Age

par(mfrow = c(1,1))

plot(advice, vertex.size = 20,  edge.arrow.size = 0.4,
     vertex.label.cex = 1.5, main = "Advice", 
     vertex.color = attrs$Dept)

plot(friend, vertex.size = 20,  edge.arrow.size = 0.4,
     vertex.label.cex = 1.5, main = "Advice", 
     vertex.color = attrs$Dept)

plot(report, vertex.size = 20,  edge.arrow.size = 0.4,
     vertex.label.cex = 1.5, main = "Advice", 
     vertex.color = attrs$Dept)

assortativity(advice, V(advice)$Dept)
assortativity(friend, V(friend)$Dept)
assortativity(report, V(report)$Dept)

assortativity(advice, V(advice)$Age)
assortativity(friend, V(friend)$Age)
assortativity(report, V(report)$Age)

