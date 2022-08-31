# Project by Massimiliano Mancini

rm(list = ls())

library(igraph)
library(ergm)
library(intergraph)

setwd("~/unifi/sna/project")

Y <- as.matrix(read.table('HN34_10.DAT'))
attrs <- read.table('CBE10.DAT')

g <- graph_from_adjacency_matrix(Y, mode = 'directed')

V(g)$gender <- attrs$V1
V(g)$delinq <- attrs$V2
V(g)$friend <- attrs$V3

net = network(Y, directed = T)

# add network and edge attributes
net %v% "gender" <- attrs$V1
net %v% "delinq" <- attrs$V2
net %v% "friend" <- attrs$V3

# let us give a look at the properties
net

rho <- graph.density(g)
oddsRho <- rho/(1-rho)

cat('Density: ', rho)
cat('Density odds: ', oddsRho)

dyad.census(g)

rec <- reciprocity(g)
oddsRec <- rec/(1-rec)

cat('Reciprocity: ', rec)
cat ('Reciprocity odds: ', oddsRec)


censusLabels = c('empty',
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

data.frame(censusLabels, triad.census(g))

tran <- transitivity(g, type = 'global')
oddsTran <- tran/(1-tran)

standardizeTrans <- log(oddsTran) - log(oddsRho) 

cat('Transitivity: ', tran)
cat('Transitivity odds: ', oddsTran)
cat('Standardized transitivity: ', standardizeTrans)

fine = 500

genderColor <- ifelse(V(g)$gender == 1, 'pink', 'lightblue')

pal = colorRampPalette(c('white','red'))
delinqColor = pal(fine)[as.numeric(cut(V(g)$delinq,breaks = fine))]

pal = colorRampPalette(c('blue','white'))
friendColor = pal(fine)[as.numeric(cut(V(g)$friend,breaks = fine))]

set.seed(6)
png(filename = 'gender.png', width = 1024, height = 1024)

plot(g, 
     vertex.size = 10, 
     vertex.color = genderColor, 
     main = 'Gender')
dev.off()

png(filename = 'dt.png', width = 1024, height = 1024)

set.seed(6)
plot(g, 
     vertex.size = 10, 
     vertex.color = delinqColor, 
     main = 'Delinquent tendency')
dev.off()

set.seed(6)
png(filename = 'fi.png', width = 1024, height = 1024)
plot(g, 
     vertex.size = 10, 
     vertex.color = friendColor, 
     main = 'Friendship importance')
dev.off()


assortativity(g, V(g)$gender)
assortativity(g, V(g)$delinq)
assortativity(g, V(g)$friend)

inDegree <- degree(g, normalized = TRUE, mode = 'in')
outDegree <- degree(g, normalized = TRUE, mode = 'out')

set.seed(6)
png(filename = 'indegree.png', width = 1024, height = 1024)
plot(g, 
     vertex.size = inDegree*20, 
     vertex.color = delinqColor, 
     main = 'In Degree')
dev.off()


set.seed(6)
png(filename = 'outdegree.png', width = 1024, height = 1024)
plot(g, 
     vertex.size = outDegree*20, 
     vertex.color = delinqColor, 
     main = 'Out Degree')
dev.off()

print('Best 3 InDegree: ')
sort(inDegree, decreasing = TRUE, index.return = TRUE)$x[1:3]

print('Best 3 OutDegree: ')
sort(outDegree, decreasing = TRUE, index.return = TRUE)$x[1:3]

comps <- components(g, mode = 'strong')
gc <- subgraph(g, V(g)[comps$membership == which.max(comps$csize)])

inCloseness <- closeness(gc, normalized = TRUE, mode = 'in')
outCloseness <- closeness(gc, normalized = TRUE, mode = 'out')

set.seed(6)
png(filename = 'incloseness.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = inCloseness*20, 
     vertex.color = delinqColor, 
     main = 'Closeness In')
dev.off()

set.seed(6)
png(filename = 'outcloseness.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = outCloseness*20, 
     vertex.color = delinqColor, 
     main = 'Closeness Out')
dev.off()

print('Best 3 InCloseness: ')
sort(inCloseness, decreasing = TRUE, index.return = TRUE)$x[1:3]

print('Best 3 OutCloseness: ')
sort(outCloseness, decreasing = TRUE, index.return = TRUE)$x[1:3]


betw <- betweenness(gc, normalized = TRUE)


set.seed(6)
png(filename = 'betweenness.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = betw*20, 
     vertex.color = delinqColor, 
     main = 'Betweenness')
dev.off()

print('Best 3 Betweenness: ')
sort(betw, decreasing = TRUE, index.return = TRUE)$x[1:3]

geigen <- eigen_centrality(gc, scale = TRUE)$vector

set.seed(6)
png(filename = 'eigen.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = geigen*20, 
     vertex.color = delinqColor, 
     main = 'Eigenvector')
dev.off()

print('Best 3 Eigenvector: ')
sort(geigen, decreasing = TRUE, index.return = TRUE)$x[1:3]

centr_degree(gc, loops = FALSE, mode = 'in')$centralization
centr_degree(gc, loops = FALSE, mode = 'out')$centralization

centr_clo(gc, mode = 'in')$centralization
centr_clo(gc, mode = 'out')$centralization

centr_betw(gc)$centralization

centr_eigen(gc)$centralization




# Homogeneous simple random graph model in ERGM flavor
m0 <- ergm(net ~ 
              edges, 
              control = control.ergm(seed = 0))

# Non-homogeneous simple random graph model in ERGM flavor
m1 <- ergm(net ~ 
              edges 
              + sender 
              + receiver, 
              control = control.ergm(seed = 0))

# dyad independent model
m2Fail <- ergm(net ~ 
                   edges 
                 + sender 
                 + receiver
                 + mutual,
                 control = control.ergm(seed = 0))

# dyad independent model
m2NoAicBic <- ergm(net ~ 
               edges 
               + receiver
               + mutual,
               control = control.ergm(seed = 0))

# dyad independent model
m2 <- ergm(net ~ 
              edges 
              + mutual,
              control = control.ergm(seed = 0))

# dyad independent model with nodal attributes (too many)
m3TooManyAttrs <- ergm(net ~ 
                        edges 
                        + mutual 
                        + nodecov("delinq") 
                        + nodecov("friend") 
                        + nodefactor("gender") 
                        + absdiff("delinq")
                        + absdiff("friend")
                        + nodematch("gender"), 
                        control = control.ergm(seed = 0)) 

# dyad independent model with nodal attributes
m3 <- ergm(net ~ 
              edges 
            + mutual 
            + nodecov("delinq") 
            + absdiff("delinq")
            + nodematch("gender"), 
            control = control.ergm(seed = 0))

# markov model
m4Fail <- ergm(net ~ 
                edges
                + istar(2)
                + ostar(2)
                + triangle,
                control = control.ergm(seed = 0))


# markov model without triangle
m4StillFail <- ergm(net ~ 
                      edges
                      + istar(2)
                      + ostar(2),
                      control = control.ergm(seed = 0))

# markov model without triangle alternating k-star terms
m4TooSimple <- ergm(net ~ 
                     edges
                   + gwidegree(decay = 1, fixed = TRUE)
                   + gwodegree(decay = 1, fixed = TRUE),
                   control = control.ergm(seed = 0))

# markov model without triangle alternating k-star terms with nodal attributes
m4TooManyParams <- ergm(net ~ 
                     edges
                   + mutual
                   + nodecov("delinq") 
                   + absdiff("delinq")
                   + nodematch("gender")
                   + gwidegree(decay = 1, fixed = TRUE)
                   + gwodegree(decay = 1, fixed = TRUE),
                   control = control.ergm(seed = 0))

# markov model without triangle alternating k-star terms
m4 <- ergm(net ~ 
               edges
             + mutual
             + nodecov("delinq") 
             + absdiff("delinq")
             + nodematch("gender")
             + gwodegree(decay = 1, fixed = TRUE),
             control = control.ergm(seed = 0))

# Social circuit model with alternating k-triangles and k-2 paths
m5Fail <- ergm(net ~ 
                edges
              + gwesp(decay = 1, fixed = TRUE) 
              + gwdsp(decay = 1, fixed = TRUE),
              control = control.ergm(seed = 0))

# Social circuit model with alternating k-triangles
m5StillFail <- ergm(net ~ 
                 edges
               + gwesp(decay = 1, fixed = TRUE),
               control = control.ergm(seed = 0))

# Social circuit model with alternating k-2 paths
m5 <- ergm(net ~ 
                edges
              + mutual
              + nodecov("delinq") 
              + absdiff("delinq")
              + nodematch("gender")
              + gwodegree(decay = 1, fixed = TRUE)
              + gwdsp(decay = 1, fixed = TRUE),
              control = control.ergm(seed = 0))

# Social circuit model with alternating k-2 paths remove gender
m6 <- ergm(net ~ 
             edges
           + mutual
           + nodecov("delinq") 
           + absdiff("delinq")
           + gwodegree(decay = 1, fixed = TRUE)
           + gwdsp(decay = 1, fixed = TRUE),
           control = control.ergm(seed = 0))

mcmc.diagnostics(m6)

AIC(m0, m1, m2, m3, m4, m5, m6)
BIC(m0, m1, m2, m3, m4, m5, m6)

sdSim = c()
meanSim = c()
tranSim = c()

for (b in 1:100) {
  ig <- asIgraph(simulate(m0, burnin = 1000, nsim = 1, verbose = TRUE))
  sdSim[b] <- sd(degree(ig))
  meanSim[b] <- mean(degree(ig))
  tranSim[b] <- transitivity(ig)
}

par(mfrow = c(1,3))
hist(sdSim, main = 'sd Degree')
abline(v = sd(degree(g)), col = 'red', lty = 2, lwd = 2)
hist(meanSim, main = 'mean Degree')
abline(v = mean(degree(g)), col = 'red', lty = 2, lwd = 2)
hist(tranSim, main = 'transitivity', xlim = c(0, 0.5))
abline(v = transitivity(g), col = 'red', lty = 2, lwd = 2)

