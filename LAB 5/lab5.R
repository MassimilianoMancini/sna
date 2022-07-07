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

# Experiment: in-degree centralitazion

