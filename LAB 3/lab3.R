# clean up
rm(list = ls())

# load igraph library
library(igraph)

# load data
load('LAB 3/lab3.Rdata')

# create the adjacency matrices
adviceY <- get.adjacency(advice, sparse = FALSE)
friendY <- get.adjacency(friend, sparse = FALSE)
reportY <- get.adjacency(report, sparse = FALSE)

diag(adviceY) <- NA
diag(friendY) <- NA
diag(reportY) <- NA

# standardized degree centrality
adviceZ <- rowSums(adviceY, na.rm = TRUE)/(vcount(advice) - 1)
friendZ <- rowSums(friendY, na.rm = TRUE)/(vcount(friend) - 1)
reportZ <- rowSums(reportY, na.rm = TRUE)/(vcount(report) - 1)

zid <- cbind(
  advice = adviceZ[order(as.numeric(names(adviceZ)))], 
  friend = friendZ[order(as.numeric(names(friendZ)))],
  report = reportZ[order(as.numeric(names(reportZ)))])

summary(zid)





