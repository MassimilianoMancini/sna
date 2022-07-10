#--------------------------------------------#
#                  Lab 2                     #
#     Descriptive analysis of networks:      #
#           characterizing the network       #
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
load("lab2_networks.Rdata")

# Relations between employees of a high-tech firm,
# three types of relations are recorded for each node
# 1) report -- to whom they reported
# 2) friend -- with whom they are friends
# 3) advice -- to whom they go to for advice 


# Which are the properties of the observed networks?
# --------------------------------------------------
advice 
friend
report

# network vertices
V(advice)
V(friend)
V(report)

# vertex names
V(advice)$name
V(friend)$name
V(report)$name

# network edges 
E(advice)
E(friend)
E(report)

# how many nodes/vertices?
vcount(advice)
vcount(friend)
vcount(report)

# how many edges?
ecount(advice)
ecount(friend)
ecount(report)

# are the networks directed?
is.directed(advice)
is.directed(friend)
is.directed(report)

# graphical representaton of the network
par(mfrow = c(1,3), mar = c(0,0,2,0))
plot(advice, vertex.size = 20,  edge.arrow.size = 0.3,
     vertex.label.cex = 1.5, main = "Advice")
plot(friend, vertex.size = 20,  edge.arrow.size = 0.3,
     vertex.label.cex = 1.5, main = "Friendship")
plot(report, vertex.size = 20,  edge.arrow.size = 0.3,
     vertex.label.cex = 1.5, main = "Report")



# Q1: How connected is the network?
# ----------------------------------- 

# DENSITY
# --------
# let us start from the adjacency matrix
advice.Y = get.adjacency(advice, sparse = F)
diag(advice.Y) = NA
friend.Y = get.adjacency(friend, sparse = F)
diag(friend.Y) = NA
report.Y = get.adjacency(report, sparse = F)
diag(report.Y) = NA

n = ncol(advice.Y)
advice.rho = sum(advice.Y, na.rm = T)/(n*(n-1))
friend.rho = sum(friend.Y, na.rm = T)/(n*(n-1))
report.rho = sum(report.Y, na.rm = T)/(n*(n-1))

advice.rho
friend.rho
report.rho

# or more simply
graph.density(advice)
graph.density(friend)
graph.density(report)

# what can we say about the networks?

# what about undirected networks?
# let us transform networks in undirected ones -- only mutual relations
advice.sy = as.undirected(advice, mode = "mutual")
friend.sy = as.undirected(friend, mode = "mutual")
report.sy = as.undirected(report, mode = "mutual")
ecount(advice)
ecount(advice.sy)

advice.Y.sy = get.adjacency(advice.sy, sparse = F)
diag(advice.Y.sy) = NA
friend.Y.sy = get.adjacency(friend.sy, sparse = F)
diag(friend.Y.sy) = NA
report.Y.sy = get.adjacency(report.sy, sparse = F)
diag(report.Y.sy) = NA

# compute the density by hand
n = ncol(advice.Y.sy)
# upper/lower triangular elements only
advice.rho = sum(advice.Y.sy[lower.tri(advice.Y.sy)])/(n*(n-1)/2)
advice.rho

# but also 
sum(advice.Y.sy, na.rm = T)/(n*(n-1))

# or more simply
graph.density(advice.sy)
graph.density(friend.sy)
graph.density(report.sy)




# Q2: is there a tendency for the nodes in a directed network to return relations?
# ---------------------------------------------------------------------------------
# DYADS
# -------
# subgraph on two vertices
# (A,B) -- three possible relations between them
# null -- (0,0)
# asymmetric -- (0,1), (1,0)
# mutual -- (1,1)

# dyad census
dyad.census(advice)
unlist(dyad.census(advice))
unlist(dyad.census(friend))
unlist(dyad.census(report))

# for mutual relations, each pair is counted only ones in the dyad.census function
# Indeed, from the adjacency matrix...
sum(advice.Y * t(advice.Y), na.rm = T)
sum(advice.Y * t(advice.Y), na.rm=T)/2


# RECIPROCITY 
# ------------
# R = n. reciprocated ties / n. of observed ties

# from the adjacency matrix
num = sum(advice.Y * t(advice.Y), na.rm = T)
num
# ATT diagonal elements need to be zero
diag(advice.Y) = 0
num = sum(advice.Y * t(advice.Y), na.rm = T)
den = sum(advice.Y, na.rm = T)
R.advice = num/den
R.advice 

diag(friend.Y) = 0
num = sum(friend.Y * t(friend.Y))
den = sum(friend.Y)
R.friend = num/den
R.friend 

diag(report.Y) = 0
num = sum(report.Y * t(report.Y))
den = sum(report.Y)
R.report = num/den
R.report 

# or more simply
reciprocity(advice)
reciprocity(friend)
reciprocity(report)


# Q2: what can we say about the networks?
# in the advice network, reciprocal relations are slightly more frequent than in the friendship network
# reciprocity does not play a role in the report-to network
# this is due to the nature of the relation -- corporate hierarchy



# Q3: is there a tendency of the nodes in the network to cluster together?
# --------------------------------------------------------------------------

# TRIADS
# -------

# subgraph on three vertices
# (A,B,C) -- 8 possibile states for undirected networks
# (A,B,C) -- 16 possible states for directed networks

?triad.census

# let us focus on directed networks first
# Triad types (Davis & Leinhardt):
# 003  A, B, C, empty triad.
# 012  A->B, C 
# 102  A<->B, C  
# 021D A<-B->C 
# 021U A->B<-C 
# 021C A->B->C
# 111D A<->B<-C
# 111U A<->B->C
# 030T A->B<-C, A->C
# 030C A<-B<-C, A->C.
# 201  A<->B<->C.
# 120D A<-B->C, A<->C.
# 120U A->B<-C, A<->C.
# 120C A->B->C, A<->C.
# 210  A->B<->C, A<->C.
# 300  A<->B<->C, A<->C, completely connected.

# Let us create a vector of labels with all possible relations
census.labels = c('empty',
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

advice.tri = triad.census(advice)
friend.tri = triad.census(friend)
report.tri = triad.census(report)
data.frame(census.labels, advice.tri, friend.tri, report.tri)

# Q3: how do the three networks differ on this network statistic? 


# transitivity coefficient
?transitivity
# ATT: directed networks are transformed in undirected ones with mode = "collapse"

tr.advice = transitivity(advice)
tr.advice
 
advice.collapse = as.undirected(advice, mode = "collapse")
transitivity(advice.collapse)

tr.friend = transitivity(friend)
tr.report = transitivity(report)


c(tr.advice, tr.friend, tr.report)

# how do the three networks differ wrt this network statistic? 

# there is a higher tencency to create small clusters in the advice network than in the friendship one
# as for the reciprocity, transitivity does not play a role in the report net


# How to normalize the transitivity index?  
# let us compare it with the density
dens.advice = graph.density(advice)
dens.friend = graph.density(friend)
dens.report = graph.density(report)


# compute the log-odds ratio
# advice network
odd.dens.advice = dens.advice/(1-dens.advice)
odd.dens.advice

odd.tr.advice = tr.advice/(1-tr.advice)
odd.tr.advice

# odds ratio
odd.tr.advice/odd.dens.advice

# the chance of observing a relation between nodes sharing a common relation
# is more than three times higher than that of observing a relation between two randomly selected nodes

# all networks

# transitivity -- success and failure probabilities
tr.s = c(tr.advice, tr.friend, tr.report)
tr.f = 1-c(tr.advice, tr.friend, tr.report)

# density -- success and failure probabilities
de.s = c(dens.advice, dens.friend, dens.report)
de.f = 1 - c(dens.advice, dens.friend, dens.report)

odd.tr = tr.s / tr.f
odd.de = de.s / de.f
odd.tr/odd.de




# Q4: are highly connected nodes similar to each other?
# ------------------------------------------------------

# ASSORTATIVE MIXING
# -------------------
# Let us read the file with nodal attributes
attr = read.table("attributes.txt", sep = "", head = TRUE)
head(attr)

str(attr)
table(attr$Dept)
summary(attr$Age)

# add vertex attributes
V(advice)$Dept = attr$Dept
V(advice)$Age = attr$Age

V(friend)$Dept = attr$Dept
V(friend)$Age = attr$Age

V(report)$Dept = attr$Dept
V(report)$Age = attr$Age

# let us focus on the categorical variable department  -- 4 categories
par(mfrow = c(1,1))

plot(advice, vertex.size = 20,  edge.arrow.size = 0.4,
     vertex.label.cex = 1.5, main = "Advice", 
     vertex.color = attr$Dept)

plot(friend, vertex.size = 20,  edge.arrow.size = 0.4,
     vertex.label.cex = 1.5, main = "Advice", 
     vertex.color = attr$Dept)

plot(report, vertex.size = 20,  edge.arrow.size = 0.4,
     vertex.label.cex = 1.5, main = "Advice", 
     vertex.color = attr$Dept)


# qualititative attributes
# ---------------------------
# start from modularity
# Q = fraction of observed ties within the same category - 
#     expected fraction of ties within the same category in a random graph
# we always have Q < 1
# negative values --  disassortative mixing
# positive values -- assortative mixing

# then normalize it -- assortative coefficient
# -1 <= r <= 1
# r = -1?
# r = +1?
assortativity(advice, V(advice)$Dept)
assortativity(friend, V(friend)$Dept)
assortativity(report, V(report)$Dept)

# Q6: what can we say about the networks?

# scalar attributes
assortativity(advice, V(advice)$Age)
assortativity(friend, V(friend)$Age)
assortativity(report, V(report)$Age)

# what can we say about the networks?