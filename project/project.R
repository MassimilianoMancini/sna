# Project by Massimiliano Mancini

rm(list = ls())

library(igraph)
library(ergm)
#library(network)
library(intergraph)

# be sure to have the correct working directory
Y <- as.matrix(read.table('HN34_10.DAT'))
attrs <- read.table('CBE10.DAT')

g <- graph_from_adjacency_matrix(Y, mode = 'directed')

V(g)$gender <- attrs$V1
V(g)$delinq <- attrs$V2
V(g)$friend <- attrs$V3

net <- asNetwork(g)

# Y is the matrix
# g is the graph object
# net is the network object

# Density
rho <- graph.density(g)
oddsRho <- rho/(1-rho)

cat('Density: ', rho)
cat('Density odds: ', oddsRho)

# Type and number of dyads
dyad.census(g)

# Recioprocity
rec <- reciprocity(g)
oddsRec <- rec/(1-rec)

cat('Reciprocity: ', rec)
cat ('Reciprocity odds: ', oddsRec)

# Type and number of triads
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

# Transitivity
tran <- transitivity(g, type = 'global')
oddsTran <- tran/(1-tran)

standardizeTrans <- log(oddsTran) - log(oddsRho) 

cat('Transitivity: ', tran)
cat('Transitivity odds: ', oddsTran)
cat('Standardized transitivity: ', standardizeTrans)

# Plot graph and save it to file
fine = 500

# Gender colors
genderColor <- ifelse(V(g)$gender == 1, 'pink', 'lightblue')

# Delinq color
pal = colorRampPalette(c('white','red'))
delinqColor = pal(fine)[as.numeric(cut(V(g)$delinq,breaks = fine))]

# Friend color 
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

# Assortativity
assortativity(g, V(g)$gender)
assortativity(g, V(g)$delinq)
assortativity(g, V(g)$friend)

# Centrality
# In and out degree centrality
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

# For remaining centralities we need the biggest connected components
comps <- components(g, mode = 'strong')
gc <- subgraph(g, V(g)[comps$membership == which.max(comps$csize)])
gcnet <- asNetwork(gc)

# In and out closeness centrality
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

# Betweeness centrality
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


# Eigenvector centrality 
eigen <- eigen_centrality(gc, scale = TRUE)$vector

set.seed(6)
png(filename = 'eigen.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = eigen*20, 
     vertex.color = delinqColor, 
     main = 'Eigenvector')
dev.off()

print('Best 3 Eigenvector: ')
sort(eigen, decreasing = TRUE, index.return = TRUE)$x[1:3]

# Network centralization
centr_degree(gc, loops = FALSE, mode = 'in')$centralization
centr_degree(gc, loops = FALSE, mode = 'out')$centralization

centr_clo(gc, mode = 'in')$centralization
centr_clo(gc, mode = 'out')$centralization

centr_betw(gc)$centralization

centr_eigen(gc)$centralization

# Modeling
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

# dyad independent model, some nodes have -Inf value
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

png(filename = 'diagnostic-m6.png', width = 1024, height = 1024)
mcmc.diagnostics(m6)
dev.off()

aic <- AIC(m0, m1, m2, m3, m4, m5, m6)
bic <- BIC(m0, m1, m2, m3, m4, m5, m6)
data.frame(aic, bic$BIC )



meanDegreeSim <- c()
sdInDegreeSim <- c()
sdOutDegreeSim <- c()
recipSim <- c()
tranSim <- c()

for (b in 1:1000) {
  simg <- asIgraph(simulate(m6, burnin = 1000, nsim = 1, verbose = TRUE))
  meanDegreeSim[b] <- mean(degree(simg))
  sdInDegreeSim[b] <- sd(degree(simg, mode = 'in'))
  sdOutDegreeSim[b] <- sd(degree(simg, mode = 'out'))
  recipSim[b] <- reciprocity(simg)
  tranSim[b] <- transitivity(simg)
}


meanReal <- mean(degree(g))
sdIn <- sd(degree(g, mode = 'in'))
sdOut <- sd(degree(g, mode = 'out'))

png(filename = 'meanDegreeSim.png', width = 1024, height = 1024)
hist(meanDegreeSim, main = 'mean Degree')
abline(v = meanReal, col = 'red', lty = 2, lwd = 2)
dev.off()

png(filename = 'sdInDegreeSim.png', width = 1024, height = 1024)
hist(sdInDegreeSim, main = 'sd InDegree')
abline(v = sdIn, col = 'red', lty = 2, lwd = 2)
dev.off()

png(filename = 'sdOutDegreeSim.png', width = 1024, height = 1024)
hist(sdOutDegreeSim, main = 'sd OutDegree')
abline(v = sdOut, col = 'red', lty = 2, lwd = 2)
dev.off()

png(filename = 'recipSim.png', width = 1024, height = 1024)
hist(recipSim, main = 'Reciprocity')
abline(v = rec, col = 'red', lty = 2, lwd = 2)
dev.off()

png(filename = 'tranSim.png', width = 1024, height = 1024)
hist(tranSim, main = 'Transitivity', xlim = c(0,0.5))
abline(v = tran, col = 'red', lty = 2, lwd = 2)
dev.off()

t.test(meanDegreeSim, mu = meanReal)
t.test(sdInDegreeSim, mu = sdIn)
t.test(sdOutDegreeSim, mu = sdOut)
t.test(recipSim, mu = rec)
t.test(tranSim, mu = tran)



