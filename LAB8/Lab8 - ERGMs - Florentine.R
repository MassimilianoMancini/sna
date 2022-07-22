library(ergm)
library(igraph)

data(florentine) # loads flomarriage and flobusiness data
?florentine


# data on marriages between Renaissance Florentine families
# nodal attributes
# 1) wealth -- in thousands of lira
# 2) priorates -- number of seats in the civil council
# 3) totalities -- tolal number of marriage ties in the original dataset made by 116 families
# 4) names -- family name

flomarriage

# let us create the igraph object
adj = as.matrix.network(flomarriage)
priorates = get.vertex.attribute(flomarriage, "priorates")
totalties = get.vertex.attribute(flomarriage, "totalties")
wealth = get.vertex.attribute(flomarriage, "wealth")
names = get.vertex.attribute(flomarriage, "vertex.names")

library(igraph)
floM = graph_from_adjacency_matrix(adj, mode = "undirected")
V(floM)$priorates = priorates
V(floM)$totalities = totalties
V(floM)$wealth = wealth
V(floM)$names = names


# how does the network look like? 
par(mar = c(0,0,0,0))
plot(floM, vertex.size = 0.8, vertex.label = V(floM)$names)

# let us also represent the wealth of the family
plot(floM, vertex.size = V(floM)$wealth/5)

# let us also represent the total number of marriages per family
plot(floM, vertex.size = V(floM)$totalities/5)



# let us also represent the total number of priorates per family
plot(floM, vertex.size = V(floM)$priorates/5)


# network properties
# -------------------
vcount(floM)
ecount(floM)
is.directed(floM)


# 1) descriptive analysis of the network
# -----------------------------------------

# network density
graph.density(floM)
# the prob of observing a tie is quite low

# dyads (does it make sense?)
dyad.census(floM)
# what can we say?

# reciprocity (does it make sense?)
reciprocity(floM)
# R = 1 --  perfect reciprocity
# R = 0 --  no reciprocity
# what can we say?

# transitivity
transitivity(floM)
# C = 1 --  perfect transitivity
# C = 0 --  no transitivity
# what can we say?

# homophily and assortative coefficient
assortativity(floM, V(floM)$wealth)
# r = 1 -- perfect homophily
# r = -1 -- perfect eterophily
# what can we say?

# what about the other nodal attributes?


# 2) centrality analysis
# -----------------------

# 1) degree centrality
# a node is central if it is connected to many others
deg = degree(floM, normalized = T)

# which are the most central nodes? 
ord = order(deg, decreasing = T)
data.frame(Family = V(floM)$names[ord], Degree = deg[ord])[1:3,]


# 2) closeness centrality
# a node is central if it is close to many others
clo = closeness(floM, normalized = T)
# there is one isolated node (Pucci's family) -- geodesic distance is infinity

# which are the most central nodes? 
ord = order(clo, decreasing = T)
data.frame(Family = V(floM)$names[ord], Closeness = clo[ord])[1:3,]


# 3) betweeness centrality
# a node is central if it lies in between many others
bet = betweenness(floM, normalized = T)

# which are the most central nodes? 
ord = order(bet, decreasing = T)
data.frame(Family = V(floM)$names[ord], Betweenness = bet[ord])[1:3,]


# 4) eighenvector centrality
# a node is central if connected to other central nodes
eig = eigen_centrality(floM, scale = T)$vector

# which are the most central nodes? 
ord = order(eig, decreasing = T)
data.frame(Family = V(floM)$names[ord], Eigen = eig[ord])[1:3,]


# let us represent graphically centrality
pdf("file.pdf")
# par(mfrow = c(2,2), mar = c(0,0,3,0))
plot(floM, vertex.size = deg*50, main = "Degree", vertex.label.cex = 1.2, asp = 0)
plot(floM, vertex.size = clo*50, main = "Closeness", vertex.label.cex = 1.2, asp = 0)
plot(floM, vertex.size = bet*50, main = "Betweenees", vertex.label.cex = 1.2, asp = 0)
plot(floM, vertex.size = eig*50, main = "Eigenvector", vertex.label.cex = 1.2 , asp = 0)
dev.off()


# what can we say?
# Mecidi is the most central (important) family in the network 
# besides the centrality measure we consider


# ERG modelling 
# ------------------
# let us start from the BRG model -- use the original network (object of class "network")
brg = ergm(flomarriage ~ edges) 
summary(brg) 

# We obtain ML estimates
# edge term is the log-odds of observing a tie under the independence assumption
# the negative sign indicates the P(Y_ij = 1) <0.5

plogis(brg$coef)
# which was the density of the network?


# Let us consider a Markov dependence model
# ------------------------------------------

# Include a triangle term 
# network statistics
summary(flomarriage ~ edges + triangle) 

markov = ergm(flomarriage ~ edges + triangle, control = control.ergm(seed = 0))
summary(markov)
# the triangle parameter does not seem to be significant

# let's add the gwdegree statistic
markov2 = ergm(flomarriage ~ edges + triangle + 
            gwdegree(1, fixed = T), control = control.ergm(seed = 0))
summary(markov2)
# clustering doesn' seem to have an effect on this network


# let us try with the partial conditional dependence
# ---------------------------------------------------
partial = ergm(flomarriage ~ edges + triangle + 
                 gwdegree(1, fixed = T) + gwesp(1, fixed = T) + gwdsp(1, fixed = T),
                 control = control.ergm(seed = 0))
summary(partial)
# none of these effects seem to be significantly different from zero

# something else is affecting the observed relations


# maybe covariates?
# --------------------

brg.cov = ergm(flomarriage ~ edges + nodecov('wealth'))
summary(brg.cov)
exp(brg.cov$coef)
# positive and significant sign for the paramerter wealth: when family's wealth increase 
# the risk of observing a tie increases by 1% when a the family wealth increases by 1000 lire

# let us consider instead an homophily effect 
brg.cov2 = ergm(flomarriage ~ edges +  absdiff("wealth"))
summary(brg.cov2)
# homophily seems to play a role in determining the ties: 
# marriages between members of the same family are more likely than those among members of different families
exp(brg.cov2$coef)
# the risk or observing a tie between members of the same family is 1% higher that corresponding to members of
# different families


# what about both effects? 
brg.cov3 = ergm(flomarriage ~ edges +  absdiff("wealth") + nodecov("wealth"))
summary(brg.cov3)
# some overappling exists -- both effects are now non-significant

# which model should we prefer?
BIC(brg.cov, brg.cov2)
AIC(brg.cov, brg.cov2)

# check model adequacy 
library(intergraph)
sdDeg.sim = meanDeg.sim = tr.sim = c()
for(b in 1:100){
  sim = simulate(brg.cov2, burnin = 1000, nsim = 1, verbose = TRUE)
  ig = asIgraph(sim)
  sdDeg.sim[b] = sd(degree(ig))
  meanDeg.sim[b] = mean(degree(ig))
  tr.sim[b] = transitivity(ig)
}

# graphical comparison
par(mfrow = c(1,3))
hist(sdDeg.sim, main = "sd Degree")
abline(v = sd(degree(floM)), col = "red", lty = 2, lwd = 2)

hist(meanDeg.sim, main = "mean Degree")
abline(v = mean(degree(floM)), col = "red", lty = 2, lwd = 2)

hist(tr.sim, main = "Transitivity")
abline(v = transitivity(floM), col = "red", lty = 2, lwd = 2)

# The model seems appropriate to describe all these three features of the network.
