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
}


