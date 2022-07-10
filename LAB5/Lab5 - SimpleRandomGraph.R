#--------------------------------------------#
#--------------------------------------------#
#                  Lab 5                     #
#           Simple random graph model        #
#--------------------------------------------#
#--------------------------------------------#

rm(list = ls())

# ----------
# 1. Set-up 
# ----------
library(igraph)




# ---- Generating a DIRECTED graph from Binomial random graph model G(n, p0) -----
par(mfrow = c(1,3))
set.seed(5)

# BRG 1
n = 50
p0 = 0.05
Y = matrix(0, n, n)
tmp = rbinom(n*(n-1), 1, p0)
Y[row(Y) != col(Y)] = tmp
diag(Y) = NA
graph.directed005 = graph_from_adjacency_matrix(Y, mode = "directed")

plot(graph.directed005, vertex.size = 10,  edge.arrow.size = 0.05)

# BRG 2
p0 = 0.15
Y = matrix(0, n, n)
tmp = rbinom(n*(n-1), 1, p0)
Y[row(Y) != col(Y)] = tmp
diag(Y) = NA
graph.directed015 = graph_from_adjacency_matrix(Y, mode = "directed")

plot(graph.directed015, vertex.size = 10,  edge.arrow.size = 0.05)

# BRG 3
p0 = 0.25
Y = matrix(0, n, n)
tmp = rbinom(n*(n-1), 1, p0)
Y[row(Y) != col(Y)] = tmp
diag(Y) = NA
graph.directed025 = graph_from_adjacency_matrix(Y, mode = "directed")

plot(graph.directed025, vertex.size = 10,  edge.arrow.size = 0.05)

# Clearly, the higher p the higher the number of observed links


# ---- Generating a DIRECTED graph from Uniform random graph model G(n, m) -----
set.seed(2)

n = 50
# let us fix m to the observed number of ties in the last directed network
m = sum(Y, na.rm = T)
m

N = n*(n-1)
ones = c(rep(1, m))
zeros = c(rep(0, N - m))
all = c(zeros, ones)

# let us exploit a simple random sampling without replacement
?sample
tmp = sample(all, N, replace = F)
# ATT! replace = F is the default option
sum(tmp)
# let us build the adjacency matrix
Y = matrix(NA, n, n)
Y[row(Y) != col(Y)] = tmp


# Maximum likelihood estimate of p
# ---------------------------------
p.MLE = mean(Y, na.rm = T)
p.MLE
p0
# the BRG and the URG models are equivalent  -- SIMPLE RANDOM GRAPH


# ---- Generating an UNDIRECTED graph from Binomial random graph model G(n, p0) -----
par(mfrow = c(1,3))
set.seed(2)

# BRG 1
n = 50
p0 = 0.05
Y = matrix(0, n, n)
tmp = rbinom(n*(n-1)/2, 1, p0)
Y[lower.tri(Y)] = tmp 
Y = Y + t(Y)
diag(Y) = NA

graph.undirected005 = graph_from_adjacency_matrix(Y, mode = "undirected")

plot(graph.undirected005, vertex.size = 10)


# BRG 2
p0 = 0.15
Y = matrix(0, n, n)
tmp = rbinom(n*(n-1)/2, 1, p0)
Y[lower.tri(Y)] = tmp 
Y = Y + t(Y)
diag(Y) = NA
graph.undirected015 = graph_from_adjacency_matrix(Y, mode = "undirected")

plot(graph.undirected015, vertex.size = 10)


# BRG 3
p0 = 0.25
Y = matrix(0, n, n)
tmp = rbinom(n*(n-1)/2, 1, p0)
Y[lower.tri(Y)] = tmp 
Y = Y + t(Y)
diag(Y) = NA
graph.undirected025 = graph_from_adjacency_matrix(Y, mode = "undirected")

plot(graph.undirected025, vertex.size = 10)

# As before, the higher p, the higher the number of observed links



# ---- Generating an UNDIRECTED graph from Uniform random graph model G(n, m) -----
set.seed(2)

n = 50
# let us fix m to the observed number of ties in the last network
m = sum(Y[upper.tri(Y)], na.rm = T)
m

N = n*(n-1)/2
ones = c(rep(1, m))
zeros = c(rep(0, N - m))
all = c(zeros, ones)

# let us exploit a simple random sampling without replacement
tmp = sample(all, replace = F)
# ATT! replace = F is the default option
Y = matrix(0, n, n)
Y[lower.tri(Y)] = tmp 
Y = Y + t(Y)
diag(Y) = NA


# Maximum likelihood estimate of p
# ---------------------------------
p.MLE = mean(Y[upper.tri(Y)], na.rm = T)
p.MLE
p0


# -----------------------------------
# ----- Assessing significance ------
# -----------------------------------


# let us consider advice relations from the high-tech worker network

# --------------------
# Load the data
# -------------------- 
load("HighTechNetworks.Rdata")


# let us extract the adjacency matrix
Y = get.adjacency(advice.net, sparse = F)
diag(Y) = NA


# number of nodes
n = nrow(Y)


# 1. Is the observed network coherent with a Binomial random graph model G(n, 0.42)?
# ATT: fixed value for p0 
# ----------------------------------------------------------------------------------
p0 = 0.42

# let us compare the observed network with the null model via the density statistic
rho.obs = mean(Y, na.rm = T)
rho.obs

B = 1000
rho.sim = c()
for(b in 1:B){
  tmp = rbinom(n^2,1,p0)  
  Y.sim = matrix(tmp, n,n)
  diag(Y.sim) = NA
  rho.sim[b] = mean(Y.sim, na.rm = TRUE)
}

# Graphical comparison
par(mfrow = c(1,1))
hist(rho.sim, col = "lightgray", main = "Null distribution")
abline(v = rho.obs, col = "red", lwd=2)
# what can we say? 

# compute an approximate p-value
mean(rho.sim >= rho.obs)
# is the observed network coherent with a BRG model G(n, 0.42)?



# let us consider the in-degree centralization and the clustering coefficient statistic
CId.obs = centr_degree(advice.net, mode = "in", loops = F)$centralization
CId.obs  
C.obs = transitivity(advice.net)
C.obs

B = 1000
CId.sim = C.sim = c()
for(b in 1:B){
  tmp = rbinom(n^2,1,p0)  
  Y.sim = matrix(tmp, n,n)
  diag(Y.sim) = NA
  g.sim = graph_from_adjacency_matrix(Y.sim)
  CId.sim[b] = centr_degree(g.sim, mode = "in", loops = F)$centralization
  C.sim[b] = transitivity(g.sim)
}

# Graphical comparison
par(mfrow = c(1,2))
hist(CId.sim, col = "lightgray", main = "Null distribution")
abline(v = CId.obs, col = "red", lwd=2)

hist(C.sim, col = "lightgray", main = "Null distribution")
abline(v = C.obs, col = "red", lwd=2)


# compute an approximate p-value
mean(CId.sim >= CId.obs)
mean(C.sim >= C.obs)
# is the network coherent with a BRG model G(n, 0.42) in terms of the above statistics?




# 2. Is the observed network coherent with the family of Binomial random graph models G(n, p)?
# a naive approach based on the best case scenario
# ------------------------------------------------------------------------------------

# maximum likelihood estimate of p
p.MLE = mean(Y, na.rm = T)

B = 1000
rho.sim = CId.sim = C.sim = c()
for(b in 1:B){
  tmp = rbinom(n^2,1,p.MLE)  
  Y.sim = matrix(tmp, n,n)
  diag(Y.sim) = NA
  rho.sim[b] = mean(Y.sim, na.rm = TRUE)
  g.sim = graph_from_adjacency_matrix(Y.sim)
  CId.sim[b] = centr_degree(g.sim, mode = "in", loops = F)$centralization
  C.sim[b] = transitivity(g.sim)
}


# Graphical comparison
par(mfrow = c(1,3))
low = pmin(min(rho.sim), rho.obs) - 0.05
up = pmax(max(rho.sim), rho.obs) + 0.05
hist(rho.sim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = rho.obs, col = "red", lwd=2)


low = pmin(min(CId.sim), CId.obs) - 0.05
up = pmax(max(CId.sim), CId.obs) + 0.05
hist(CId.sim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = CId.obs, col = "red", lwd=2)


low = pmin(min(C.sim), C.obs) - 0.05
up = pmax(max(C.sim), C.obs) + 0.05
hist(C.sim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = C.obs, col = "red", lwd=2)


# compute an approximate p-value
mean(rho.sim >= rho.obs)
mean(CId.sim >= CId.obs)
mean(C.sim >= C.obs)
# is the network coherent with the family of BRG models [G(n, p)] in terms of the above statistics?




# 3. Is the observed network coherent with a Binomial random graph model G(n, p)?
# a more formal approach based on the conditional uniform distribution
# ------------------------------------------------------------------------------------

# let us consider the degree centralization index and the transitivity coefficient
# NB. the density statistic makes no sense -- by fixing m, then rho = p is constant 
B = 5000
m =  sum(Y, na.rm = TRUE)
rho.sim = CId.sim = C.sim = c()
for(b in 1:B){
  Y.sim = matrix(, n, n) 
  ones = rep(1, m)
  zeros = rep(0, n*(n-1) - m)
  all = c(ones, zeros)
  Y.sim[col(Y.sim) != row(Y.sim)] = sample(all, n*(n-1))
  g.sim = graph_from_adjacency_matrix(Y.sim)
  CId.sim[b] = centr_degree(g.sim, mode = "in", loops = F)$centralization
  C.sim[b] = transitivity(g.sim)
}


# Graphical comparison
par(mfrow = c(1,2))
hist(CId.sim, col = "lightgray", main = "Null distribution")
abline(v = CId.obs, col = "red", lwd=2)


hist(C.sim, col = "lightgray", main = "Null distribution")
abline(v = C.obs, col = "red", lwd=2)


# compute an approximate p-value
mean(CId.sim >= CId.obs)
mean(C.sim >= C.obs)
# is the network coherent with a BRG model G(n, p) in terms of the above statistics?

# REPEAT THE ABOVE ANALYSIS FOR THE FRIENDSHIP NETWORK