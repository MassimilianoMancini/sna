rm(list = ls())

library(igraph)

# Generating a directed graph from binomial random graph model G(n, p0)
par(mfrow = c(1, 3))

# Binomial Random Graphs
for (p0 in c(0.05, 0.15, 0.25)) {
  n <- 50
  Y <- matrix(NA, n, n)
  Y[row(Y) != col(Y)] <- rbinom(n*(n-1), 1, p0)
  brg <- graph_from_adjacency_matrix(Y, mode = 'directed')
  plot (brg, vertex.size = 15, edge.arrow.size = 0.005, main = paste('p = ', p0))
  box()
}

# Uniform Random Graph
m <- sum(Y, na.rm = TRUE)
N <- n*(n-1)
Y <- matrix(NA, n, n)
ones <- c(rep(1, m))
zeros <- c(rep(0, N - m))
Y[row(Y) != col(Y)] <- sample(c(zeros, ones), N)
brg <- graph_from_adjacency_matrix(Y, mode = 'directed')
plot (brg, vertex.size = 15, edge.arrow.size = 0.005, main = paste('p = ', p0))
box()

# maximum likelihood estimate of p
pmle <- mean(Y, na.rm = TRUE)
pmle
p0


# Undirected Binomial Random Graphs
par(mfrow = c(1, 3))

for (p0 in c(0.05, 0.15, 0.25)) {
  n <- 50
  Y <- matrix(0, n, n)
  Y[lower.tri(Y)] <- rbinom(n*(n-1)/2, 1, p0)
  Y <- Y + t(Y)
  diag(Y) <- NA
  brg <- graph_from_adjacency_matrix(Y, mode = 'undirected')
  plot (brg, vertex.size = 15, edge.arrow.size = 0.005, main = paste('p = ', p0))
  box()
}

# Undirected Uniform Random Graph
m <- sum(Y[lower.tri(Y)], na.rm = TRUE)
n <- 50
Y <- matrix(0, n, n)
N <- n*(n-1)/2
ones <- c(rep(1, m))
zeros <- c(rep(0, N - m))
Y[lower.tri(Y)] <- sample(c(zeros, ones))
Y <- Y + t(Y)
diag(Y) <- NA
brg <- graph_from_adjacency_matrix(Y, mode = 'undirected')
plot (brg, vertex.size = 15, edge.arrow.size = 0.005, main = paste('p = ', p0))
box()

pmle <- mean(Y[lower.tri(Y)], na.rm = TRUE)
pmle
p0


# Experiment: is the observed network coherent with a BRG model G(n, 0.42)?
load('LAB 5/lab5.Rdata')
Y <- get.adjacency(advice, sparse = FALSE)
diag(Y) <- NA
n <- nrow(Y)
p0 <- 0.42

rhoObs <- mean(Y, na.rm = TRUE)

B <- 1000
rhoSim <- c()
for (b in 1:B) {
  Ysim <- matrix(rbinom(n^2, 1, p0), n, n)
  diag(Ysim) <- NA
  rhoSim[b] <- mean(Ysim, na.rm = TRUE)
}

# graphical comparison
par(mfrow = c(1,1))
hist(rhoSim, col = "lightgray", main = "Null distribution")
abline(v = rhoObs, col = "red", lwd = 2)

# p-value
mean(rhoSim >= rhoObs)

# Experiment: in-degree centralitazion and transitivity
zidObs <- centr_degree(advice, mode = 'in', loops = FALSE)$centralization
traObs <- transitivity(advice)

zidSim <- c()
traSim <- c()
B <- 1000
for (b in 1:B) {
  Ysim <- matrix(rbinom(n^2, 1, p0), n, n)
  diag(Ysim) <- NA
  brg <- graph_from_adjacency_matrix(Ysim)
  zidSim[b] <- centr_degree(brg, mode = 'in', loops = FALSE)$centralization
  traSim[b] <- transitivity(brg)
}

par(mfrow = c(1, 2))
hist(zidSim, col = 'lightgray', main = 'Null distribution')
abline(v = zidObs, col = 'red', lwd = 2)
hist(traSim, col = 'lightgray', main = 'Null distribution')
abline(v = traObs, col = 'red', lwd = 2)

# p-values
mean(zidSim >= zidObs)
mean(traSim >= traObs)

# Using a maximum likelihood estimate of p
pmle <- mean(Y, na.rm = TRUE)
B <- 1000
rhoSim <- c()
zidSim <- c()
traSim <- c()
for (b in 1:B) {
  Ysim <- matrix(rbinom(n^2, 1, pmle), n, n)
  diag(Ysim) <- NA
  brg <- graph_from_adjacency_matrix(Ysim)
  rhoSim[b] <- mean(Ysim, na.rm = TRUE)
  zidSim[b] <- centr_degree(brg, mode = 'in', loops = FALSE)$centralization
  traSim[b] <- transitivity(brg)
}

# graphical comparison
par(mfrow = c(1,3))
low <- pmin(min(rhoSim), rhoObs) - 0.05
up <- pmax(max(rhoSim), rhoObs) + 0.05
hist(rhoSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = rhoObs, col = "red", lwd = 2)

low <- pmin(min(zidSim), zidObs) - 0.05
up <- pmax(max(zidSim), zidObs) + 0.05
hist(zidSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = zidObs, col = "red", lwd = 2)

low = pmin(min(traSim), traObs) - 0.05
up = pmax(max(traSim), traObs) + 0.05
hist(traSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = traObs, col = "red", lwd = 2)

# p-values
mean(rhoSim >= rhoObs)
mean(zidSim >= zidObs)
mean(traSim >= traObs)

# using conditional uniform distribution
B <- 5000
m <- sum(Y, na.rm = TRUE)
zidSim <- c()
traSim <- c()
for (b in 1:B) {
  Ysim <- matrix(NA, n, n)
  ones <- rep(1, m)
  zeros <- rep(0, n*(n-1) - m)
  Ysim[col(Ysim) != row(Ysim)] <- sample(c(zeros, ones), n*(n-1))
  brg <- graph_from_adjacency_matrix(Ysim)
  zidSim[b] <- centr_degree(brg, mode = 'in', loops = FALSE)$centralization
  traSim[b] <- transitivity(brg)
}

low <- pmin(min(zidSim), zidObs) - 0.05
up <- pmax(max(zidSim), zidObs) + 0.05
hist(zidSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = zidObs, col = "red", lwd = 2)

low = pmin(min(traSim), traObs) - 0.05
up = pmax(max(traSim), traObs) + 0.05
hist(traSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = traObs, col = "red", lwd = 2)

# p-values
mean(zidSim >= zidObs)
mean(traSim >= traObs)

# friend network
Y <- get.adjacency(friend, sparse = FALSE)
diag(Y) <- NA
n <- nrow(Y)
p0 <- 0.42

rhoObs <- mean(Y, na.rm = TRUE)

B <- 1000
rhoSim <- c()
for (b in 1:B) {
  Ysim <- matrix(rbinom(n^2, 1, p0), n, n)
  diag(Ysim) <- NA
  rhoSim[b] <- mean(Ysim, na.rm = TRUE)
}

# graphical comparison
par(mfrow = c(1,1))
hist(rhoSim, col = "lightgray", main = "Null distribution")
abline(v = rhoObs, col = "red", lwd = 2)

# p-value
mean(rhoSim >= rhoObs)

# Experiment: in-degree centralitazion and transitivity
zidObs <- centr_degree(advice, mode = 'in', loops = FALSE)$centralization
traObs <- transitivity(advice)

zidSim <- c()
traSim <- c()
B <- 1000
for (b in 1:B) {
  Ysim <- matrix(rbinom(n^2, 1, p0), n, n)
  diag(Ysim) <- NA
  brg <- graph_from_adjacency_matrix(Ysim)
  zidSim[b] <- centr_degree(brg, mode = 'in', loops = FALSE)$centralization
  traSim[b] <- transitivity(brg)
}

par(mfrow = c(1, 2))
hist(zidSim, col = 'lightgray', main = 'Null distribution')
abline(v = zidObs, col = 'red', lwd = 2)
hist(traSim, col = 'lightgray', main = 'Null distribution')
abline(v = traObs, col = 'red', lwd = 2)

# p-values
mean(zidSim >= zidObs)
mean(traSim >= traObs)

# Using a maximum likelihood estimate of p
pmle <- mean(Y, na.rm = TRUE)
B <- 1000
rhoSim <- c()
zidSim <- c()
traSim <- c()
for (b in 1:B) {
  Ysim <- matrix(rbinom(n^2, 1, pmle), n, n)
  diag(Ysim) <- NA
  brg <- graph_from_adjacency_matrix(Ysim)
  rhoSim[b] <- mean(Ysim, na.rm = TRUE)
  zidSim[b] <- centr_degree(brg, mode = 'in', loops = FALSE)$centralization
  traSim[b] <- transitivity(brg)
}

# graphical comparison
par(mfrow = c(1,3))
low <- pmin(min(rhoSim), rhoObs) - 0.05
up <- pmax(max(rhoSim), rhoObs) + 0.05
hist(rhoSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = rhoObs, col = "red", lwd = 2)

low <- pmin(min(zidSim), zidObs) - 0.05
up <- pmax(max(zidSim), zidObs) + 0.05
hist(zidSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = zidObs, col = "red", lwd = 2)

low = pmin(min(traSim), traObs) - 0.05
up = pmax(max(traSim), traObs) + 0.05
hist(traSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = traObs, col = "red", lwd = 2)

# p-values
mean(rhoSim >= rhoObs)
mean(zidSim >= zidObs)
mean(traSim >= traObs)

# using conditional uniform distribution
B <- 5000
m <- sum(Y, na.rm = TRUE)
zidSim <- c()
traSim <- c()
for (b in 1:B) {
  Ysim <- matrix(NA, n, n)
  ones <- rep(1, m)
  zeros <- rep(0, n*(n-1) - m)
  Ysim[col(Ysim) != row(Ysim)] <- sample(c(zeros, ones), n*(n-1))
  brg <- graph_from_adjacency_matrix(Ysim)
  zidSim[b] <- centr_degree(brg, mode = 'in', loops = FALSE)$centralization
  traSim[b] <- transitivity(brg)
}

low <- pmin(min(zidSim), zidObs) - 0.05
up <- pmax(max(zidSim), zidObs) + 0.05
hist(zidSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = zidObs, col = "red", lwd = 2)

low = pmin(min(traSim), traObs) - 0.05
up = pmax(max(traSim), traObs) + 0.05
hist(traSim, col = "lightgray", main = "Null distribution", xlim = c(low, up))
abline(v = traObs, col = "red", lwd = 2)

# p-values
mean(zidSim >= zidObs)
mean(traSim >= traObs)


  
  

