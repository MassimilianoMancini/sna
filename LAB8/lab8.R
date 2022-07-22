rm(list = ls())

library(ergm)
library(igraph)
library(intergraph)

data(florentine)

marriages <- as.undirected(asIgraph(flomarriage))
V(marriages)$priorates <- network::get.vertex.attribute(flomarriage, 'priorates')
V(marriages)$totalties <- network::get.vertex.attribute(flomarriage, 'totalties')
V(marriages)$wealth <- network::get.vertex.attribute(flomarriage, 'wealth')
V(marriages)$name <- network::get.vertex.attribute(flomarriage, 'vertex.names')

par(mar = c(0,0,0,0))
plot(marriages, vertex.size = V(marriages)$wealth/5)
plot(marriages, vertex.size = V(marriages)$totalties/5)
plot(marriages, vertex.size = V(marriages)$priorates/5)

vcount(marriages)
ecount(marriages)
is.directed(marriages)

# network density
graph.density(marriages)

# dyads census
dyad.census(marriages)

# reciprocity
reciprocity(marriages)

# transitivity
transitivity(marriages)

# homophily
assortativity(marriages, V(marriages)$wealth)

# degree centrality
degree(marriages, normalized = TRUE)

# closeness centrality
closeness(marriages, normalized = TRUE)

# betweeness centrality
betweenness(marriages, normalized = TRUE)

# eigenvector centrality
eigen_centrality(marriages, scale = TRUE)$vector

brg <- ergm(flomarriage ~ edges)
summary(brg)

plogis(brg$coefficients)

markov <- ergm(flomarriage ~ edges + triangles)
summary(markov)

markov2 <- ergm(flomarriage ~ edges + triangles + gwdegree(1, fixed = TRUE))
summary(markov2)

partial <- ergm(flomarriage ~ edges + triangles + gwdegree(1, fixed = TRUE) + gwesp(1, fixed = TRUE) + gwdsp(1, fixed = TRUE))
summary(partial)

# usign vertex attributes
brgCov <- ergm(flomarriage ~ edges + nodecov('wealth'))
summary(brgCov)

brgCov2 <- ergm(flomarriage ~ edges + absdiff('wealth'))
summary(brgCov2)

brgCov3 <- ergm(flomarriage ~ edges + nodecov('wealth') + absdiff('wealth'))
summary(brgCov3)

sdSim = c()
meanSim = c()
tranSim = c()

for (b in 1:1000) {
  ig <- asIgraph(simulate(brgCov2, burnin = 1000, nsim = 1, verbose = TRUE))
  sdSim[b] <- sd(degree(ig))
  meanSim[b] <- mean(degree(ig))
  tranSim[b] <- transitivity(ig)
}

par(mfrow = c(1,3))
hist(sdSim, main = 'sd Degree')
abline(v = sd(degree(marriages)), col = 'red', lty = 2, lwd = 2)
hist(meanSim, main = 'mean Degree')
abline(v = mean(degree(marriages)), col = 'red', lty = 2, lwd = 2)
hist(tranSim, main = 'transitivity')
abline(v = transitivity(marriages), col = 'red', lty = 2, lwd = 2)

