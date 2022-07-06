#--------------------------------------------#
#--------------------------------------------#
#                  Lab 3                     #
#     Descriptive analysis of networks:      #
#          characterizing the nodes          #
#--------------------------------------------#
#--------------------------------------------#

rm(list = ls())

# ----------
# 1. Set-up 
# ----------
library(igraph)



# --------------------
# 2. load the data
# -------------------- 
load("LAB 3/HighTechNetworks.Rdata")
ls()
# Relations between employees of a high-tech firm,
# three types of relations are recorded for each node
# 1) to whom they reported -- report.net
# 2) with whom they were friends -- friendship.net
# 3) and to whom they go to for advice -- advice.net


# let us trasform the three networks in undirected
advice.net.sy = as.undirected(advice.net, mode = "collapse")
friend.net.sy = as.undirected(friend.net, mode = "collapse")
report.net.sy = as.undirected(report.net, mode = "collapse")


# how do the networks look like? 
par(mfrow = c(1,3), mar = c(0,0,2,0))
plot(advice.net.sy, vertex.size = 20,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Advice network")
plot(friend.net.sy, vertex.size = 20,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Friendship network")
plot(report.net.sy, vertex.size = 20,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Report network")


# let us consider the adjacency matrices
advice.Y.sy = get.adjacency(advice.net.sy, sparse = FALSE)
diag(advice.Y.sy) = NA
friend.Y.sy = get.adjacency(friend.net.sy, sparse = FALSE)
diag(friend.Y.sy) = NA
report.Y.sy = get.adjacency(report.net.sy, sparse = FALSE)
diag(report.Y.sy) = NA
# ATT: row/column ordering of the advice net is different from 
# that of the friendship and report net


# ---------------------
# CENTRALITY MEASURES 
# ---------------------

# 1) DEGREE CENTRALITY
# ---------------------
# A node is central if it is connected to many other nodes
# zeta_i^d = sum_j y_ij

zid.advice = rowSums(advice.Y.sy, na.rm = T)
zid.advice
# but also
colSums(advice.Y.sy, na.rm = T)
zid.friend = rowSums(friend.Y.sy, na.rm = T)
zid.friend
zid.report = rowSums(report.Y.sy, na.rm = T)
zid.report

# to fix ordering issues
ord1 = order(as.numeric(names(zid.friend)))
ord2 = order(as.numeric(names(zid.report)))


zid = cbind(advice = zid.advice, friend = zid.friend[ord1], report= zid.report[ord2])
# let us consider the summary of the degree centralities
summary(zid)

# As expected nodes in the advice network have higher degree


# let us graphically represent the degree distribution 
par(mfrow = c(3,1), mar = c(2.5,4,3.5,1.5))
hist(zid.advice, breaks = 10, xlim = c(0, 21))
hist(zid.friend, breaks = 10, xlim = c(0,21))
hist(zid.report, breaks = 10, xlim = c(0,21))


# to simplify comparison... 
# standardized degree centrality
# -------------------------------
n = vcount(advice.net)
zid.st.advice = zid.advice/(n-1)
zid.st.friend = zid.friend/(n-1)
zid.st.report = zid.report/(n-1)

zid.st = cbind(advice = zid.st.advice, friend = zid.st.friend[ord1], report = zid.st.report[ord2])

# let us consider the summary of the standardized degree centralities
summary(zid.st)

# what can we say about nodes?
# -------------------------------
# ADVICE NET: 
# most of the nodes have from moderately high to high centrality
# FRIEND NET: 
# most of the nodes have from moderately to low centrality, the maximum is however close to 1
# REPORT NET: 
# almost all nodes have very low centrality


# let us graphically represent standardized centrality
par(mfrow = c(1,1))
plot(advice.net.sy, vertex.size = zid.st.advice*20,  edge.arrow.size = 0.5,
     vertex.label.cex = 1, main = "Advice network")

# what can we say about centrality of the nodes?

# color in blue the top 3 nodes
ord = order(zid.st.advice, decreasing = T)
V(advice.net.sy)$color = "orange"
V(advice.net.sy)$color[ord[1:3]] = "lightblue"

plot(advice.net.sy, vertex.size = zid.st.advice*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1, main = "Advice network")

# let us do the same for the other networks
ord = order(zid.st.friend, decreasing = T)
V(friend.net.sy)$color = "orange"
V(friend.net.sy)$color[ord[1:3]] = "lightblue"
plot(friend.net.sy, vertex.size = zid.st.friend*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1, main = "Friendship")

# what can we say about nodes? 


# let us do the same for the report to network
ord = order(zid.st.report, decreasing = T)
V(report.net.sy)$color = "orange"
V(report.net.sy)$color[ord[1:3]] = "lightblue"
plot(report.net.sy, vertex.size = zid.st.report*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Report")
# what can we say now? 


# let us put all plots in the same window for comparison
par(mfrow = c(1,3))
plot(advice.net.sy, vertex.size = zid.st.advice*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1, main = "Advice network")
plot(friend.net.sy, vertex.size = zid.st.friend*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1, main = "Friendship")
plot(report.net.sy, vertex.size = zid.st.report*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1, main = "Report")

# what can we say now? 
# --------------------------------------------------
# ** ADVICE NET
# most of the nodes have quite a high degree. There are not big differences between them
# ** FRIENDSHIP NET
# nodes 17 and 11 present the highest number of connections
# ** REPORT NET
# node 14 presents the highest number of connections -- 
# this is clearly related to the number of subordinates



# more easily...
# ---------------

# degree centrality
degree(advice.net.sy)
zid.advice
degree(advice.net.sy, normalized = T)
zid.st.advice



# 2) CLOSENESS CENTRALITY
# ------------------------
# A node is central if it is "close" to many other nodes
# z_i^c = 1/sum_j d_ij

?distances

D.advice = distances(advice.net.sy)
D.advice 
zic.advice = 1/rowSums(D.advice)

D.friend = distances(friend.net.sy)
zic.friend = 1/rowSums(D.friend)

D.report = distances(report.net.sy)
zic.report = 1/rowSums(D.report)


zic = cbind(advice = zic.advice, friend = zic.friend[ord1], report = zic.report[ord2])
zic 
# NOTE --  closeness centralities of network nodes are not that different


# to simplify comparison ....
# standardized closeness centrality
# ----------------------------------
zic.st.advice = zic.advice*(n-1)
zic.st.friend = zic.friend*(n-1)
zic.st.report = zic.report*(n-1)
zic.st = cbind(zic.st.advice, zic.st.friend[ord1], zic.st.report[ord2])
summary(zic.st)


# let us graphically represent standardized centrality
par(mfrow = c(1,3), mar = c(0,0,2,0))
# color in blue the top 3 nodes

ord = order(zic.st.advice, decreasing = T)
V(advice.net.sy)$color = "orange"
V(advice.net.sy)$color[ord[1:3]] = "lightblue"
plot(advice.net.sy, vertex.size = zic.st.advice*30, edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Advice network")

ord = order(zic.st.friend, decreasing = T)
V(friend.net.sy)$color = "orange"
V(friend.net.sy)$color[ord[1:3]] = "lightblue"
plot(friend.net.sy, vertex.size = zic.st.friend*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Friendship")

ord = order(zic.st.report, decreasing = T)
V(report.net.sy)$color = "orange"
V(report.net.sy)$color[ord[1:3]] = "lightblue"
plot(report.net.sy, vertex.size = zic.st.report*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Report")


# what can we say about the nodes?
# ----------------------------------
# ** ADVICE NET
# most of the nodes have quite a high closeness, nodes 15, 2 and 3 are the most central ones
# ** FRIENDSHIP NET
# Nodes 17, 11 and 2 present the highest closenness. Information rapidly passes through them 
# ** REPORT NET
# Node 7 has the highest closeness, followed by 14 and 21


# more easily...
# -----------------
# closeness centrality
closeness(advice.net.sy)
zic.advice

closeness(advice.net.sy, normalized = T)
zic.st.advice



# 3) BETWEENNESS CENTRALITY
# --------------------------
# a node is central if it is located between many nodes

?betweenness
zib.advice = betweenness(advice.net.sy, directed = F, normalized = F)
zib.friend = betweenness(friend.net.sy, directed = F, normalized = F)
zib.report = betweenness(report.net.sy, directed = F, normalized = F)

zib = cbind(advice = zib.advice, friend = zib.friend[ord1], report = zib.report[ord2])
summary(zib)


# standardized betweennes centrality
zib.st.advice = zib.advice*2/((n-1)*(n-2))
# or simply... 
zib.st.advice = betweenness(advice.net.sy, directed = F, normalized = T)
zib.st.friend = betweenness(friend.net.sy, directed = F, normalized = T)
zib.st.report = betweenness(report.net.sy, directed = F, normalized = T)

zib.st = cbind(advice = zib.st.advice, friend = zib.st.friend[ord1], report = zib.st.report[ord2])
summary(zib.st)


# let us graphically represent standardized centrality
par(mfrow = c(1,3), mar = c(0,0,2,0))
# color in blue the top 3 nodes

# advice network
ord = order(zib.st.advice, decreasing = T)
V(advice.net.sy)$color = "orange"
V(advice.net.sy)$color[ord[1:3]] = "lightblue"
plot(advice.net.sy, vertex.size = zib.st.advice*50,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Advice network")

# friendship network
ord = order(zib.st.friend, decreasing = T)
V(friend.net.sy)$color = "orange"
V(friend.net.sy)$color[ord[1:3]] = "lightblue"
plot(friend.net.sy, vertex.size = zib.st.friend*50,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Friendship")

# report to network
ord = order(zib.st.report, decreasing = T)
V(report.net.sy)$color = "orange"
V(report.net.sy)$color[ord[1:3]] = "lightblue"
plot(report.net.sy, vertex.size = zib.st.report*50,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Report")


# what can we say about the nodes?
# ----------------------------------
# ** ADVICE AND FRIENDSHIP NET
# All nodes present in both networks low betwennes centrality
# ** REPORT NET
# Information in the network passes mainly through node 7, followed by 14, and 21



# 4) EIGENVECTOR CENTRALITY
# --------------------------
# a node is central if it is connected to other central nodes
?eigen_centrality

zie.adivice = eigen_centrality(advice.net.sy, scale = F)$vector
zie.friend = eigen_centrality(friend.net.sy, scale = F)$vector
zie.report = eigen_centrality(report.net.sy, scale = F)$vector

zie = cbind(advice = zie.adivice, friend = zie.friend[ord1], report = zie.report[ord2])
summary(zie)

# standardized eigen vector centrality
zie.st.advice = eigen_centrality(advice.net.sy, scale = T)$vector
zie.st.friend = eigen_centrality(friend.net.sy, scale = T)$vector
zie.st.report = eigen_centrality(report.net.sy, scale = T)$vector

zie.st = cbind(advice = zie.st.advice, friend = zie.st.friend[ord1], report = zie.st.report[ord2])
summary(zie.st)

# let us graphically represent standardized centrality
par(mfrow = c(1,3), mar = c(0,0,2,0))
# color in blue the top 3 nodes

ord = order(zie.st.advice, decreasing = T)
V(advice.net.sy)$color = "orange"
V(advice.net.sy)$color[ord[1:3]] = "lightblue"

plot(advice.net.sy, vertex.size = zie.st.advice*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Advice network")


ord = order(zie.st.friend, decreasing = T)
V(friend.net.sy)$color = "orange"
V(friend.net.sy)$color[ord[1:3]] = "lightblue"
plot(friend.net.sy, vertex.size = zie.st.friend*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Friendship")

ord = order(zie.st.report, decreasing = T)
V(report.net.sy)$color = "orange"
V(report.net.sy)$color[ord[1:3]] = "lightblue"
plot(report.net.sy, vertex.size = zie.st.report*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5, main = "Report")

# what can we say about the nodes?
# ----------------------------------
# ** ADVICE NET
# most of the nodes present quite high eigenvector centrality
# ** FRIENDSHIP NET
# Nodes 17, 11, and 19 seem to the the most central ones
# High number of connections with highly connected nodes
# ** REPORT NET
# Node 14 is the most central one
# Connected with many nodes and also with 7 to which all report


# --------------------------
# NETWORK CENTRALIZATION
# --------------------------

# 1) Degree centralization
# -------------------------
d.centralization = function(c){
  n = length(c)
  max.cd = max(c)
  num = sum(max(c) - c)
  den = (n-1)*(n-2)
  return(num/den)
}


# 2) Closeness centralization
# ----------------------------
c.centralization = function(c){
  n = length(c)
  max.cd = max(c)
  num = sum(max(c) - c)
  den = (n-2)/(2*n-3)
  return(num/den)
}


# 3) betweenness centralization
# ---------------------------
b.centralization = function(c){
  n = length(c)
  max.cd = max(c)
  num = sum(max(c) - c)
  den = (n-1)^2*(n-2)/2
  return(num/den)
}

centr.d = apply(zid, 2, d.centralization)
centr.c = apply(zic, 2, c.centralization)
centr.b = apply(zib, 2, b.centralization)
tab = rbind(centr.d, centr.c, centr.b)
colnames(tab) = c("advice", "friend", "report")
tab

# more simply ... 
centr_degree(advice.net.sy, loops = F)$centralization
centr_clo((advice.net.sy))$centralization
centr_betw(advice.net.sy, directed = F)$centralization


# comparing networks
rbind(centr.d, centr.c, centr.b)
# The friendship network is more centralized
# than the others with respect to degree and closeness
# few nodes have high degree 

# the report network turns to be the more centralized wrt to the betweenness



# how to deal with the direction of the ties? 

# in- and out degree 
# how many nodes do point to a given node? 
degree(advice.net, mode = "in")
# how many nodes are pointed by a given node? 
degree(advice.net, mode = "out")

# closeness centrality
# how easily can a node be reached from other nodes? 
closeness(advice.net, mode = "in")
# how easily can a node reach the others? 
closeness(advice.net, mode = "out")

# betweenness centrality
betweenness(advice.net, directed = T)

# eigenvector centrality
eigen_centrality(advice.net, directed = T)$vector
