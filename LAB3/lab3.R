# clean up
rm(list = ls())

# load igraph library
library(igraph)

# load data
load('LAB 3/lab3.Rdata')

# transform graphs in undirected
adviceUn <- as.undirected(advice, mode = 'collapse')
friendUn <- as.undirected(friend, mode = 'collapse')
reportUn <- as.undirected(report, mode = 'collapse')

# create the adjacency matrices
adviceY <- get.adjacency(adviceUn, sparse = FALSE)
friendY <- get.adjacency(friendUn, sparse = FALSE)
reportY <- get.adjacency(reportUn, sparse = FALSE)

diag(adviceY) <- NA
diag(friendY) <- NA
diag(reportY) <- NA

# degree centrality
adviceD <- rowSums(adviceY, na.rm = TRUE)
friendD <- rowSums(friendY, na.rm = TRUE)
reportD <- rowSums(reportY, na.rm = TRUE)

# standardized degree centrality
adviceND <- rowSums(adviceY, na.rm = TRUE)/(vcount(adviceUn) - 1)
friendND <- rowSums(friendY, na.rm = TRUE)/(vcount(friendUn) - 1)
reportND <- rowSums(reportY, na.rm = TRUE)/(vcount(reportUn) - 1)

zid <- cbind(
  advice = adviceD[order(as.numeric(names(adviceD)))], 
  friend = friendD[order(as.numeric(names(friendD)))],
  report = reportD[order(as.numeric(names(reportD)))])

# more easily
degree(adviceUn)
degree(friendUn)
degree(reportUn)

# standardized
zind <- cbind(
  advice = adviceND[order(as.numeric(names(adviceND)))], 
  friend = friendND[order(as.numeric(names(friendND)))],
  report = reportND[order(as.numeric(names(reportND)))])

# more easily
degree(adviceUn, normalized = TRUE)
degree(friendUn, normalized = TRUE)
degree(reportUn, normalized = TRUE)

# closeness centrality
# standardized closeness centrality
adviceC <- 1/rowSums(distances(adviceUn))
friendC <- 1/rowSums(distances(friendUn))
reportC <- 1/rowSums(distances(reportUn))

zic <- cbind(
  advice = adviceC[order(as.numeric(names(adviceC)))],
  friend = friendC[order(as.numeric(names(friendC)))],
  report = reportC[order(as.numeric(names(reportC)))]
)

# more easily
closeness(adviceUn)
closeness(friendUn)
closeness(reportUn)

# standardized closeness centrality
adviceNC <- (vcount(adviceUn) - 1)/rowSums(distances(adviceUn))
friendNC <- (vcount(friendUn) - 1)/rowSums(distances(friendUn))
reportNC <- (vcount(reportUn) - 1)/rowSums(distances(reportUn))

zinc <- cbind(
  advice = adviceNC[order(as.numeric(names(adviceNC)))],
  friend = friendNC[order(as.numeric(names(friendNC)))],
  report = reportNC[order(as.numeric(names(reportNC)))]
)

# more easily
closeness(adviceUn, normalized = TRUE)
closeness(friendUn, normalized = TRUE)
closeness(reportUn, normalized = TRUE)

# betweenness centrality
adviceB <- betweenness(adviceUn, directed = FALSE)
friendB <- betweenness(friendUn, directed = FALSE)
reportB <- betweenness(reportUn, directed = FALSE)

zib <- cbind(
  advice = adviceB[order(as.numeric(names(adviceB)))],
  friend = friendB[order(as.numeric(names(friendB)))],
  report = reportB[order(as.numeric(names(reportB)))]
)

# standardized betweenness centrality
adviceNB <- betweenness(adviceUn, directed = FALSE, normalized = TRUE)
friendNB <- betweenness(friendUn, directed = FALSE, normalized = TRUE)
reportNB <- betweenness(reportUn, directed = FALSE, normalized = TRUE)

zinb <- cbind(
  advice = adviceNB[order(as.numeric(names(adviceNB)))],
  friend = friendNB[order(as.numeric(names(friendNB)))],
  report = reportNB[order(as.numeric(names(reportNB)))]
)


# eigenvector centrality
eigen_centrality(adviceUn)$vector
eigen_centrality(friendUn)$vector
eigen_centrality(reportUn)$vector

# network centralization

# degree centralization
dfun <- function(v) {
  n <- length(v)
  m <- max(v)
  num <- sum(m - v)
  den <- (n - 1)*(n - 2)
  return (num/den)
}

# closeness centralization
cfun <- function(v) {
  n <- length(v)
  m <- max(v)
  num <- sum(m - v)
  den <- (n - 2)/(2 * n - 3)
  return (num/den)
}

# betweenness centralization
bfun <- function(v) {
  n <- length(v)
  m <- max(v)
  num <- sum(m - v)
  den <- ((n - 1)^2)*(n - 2)/2
  return (num/den)
}

networkCentrDegree <- apply(zid, 2, dfun)
networkCentrClosne <- apply(zic, 2, cfun)
networkCentrBetwee <- apply(zib, 2, bfun)

# more easily
centr_degree(adviceUn, loops = FALSE)$centralization
centr_degree(friendUn, loops = FALSE)$centralization
centr_degree(reportUn, loops = FALSE)$centralization

centr_clo(adviceUn)$centralization
centr_clo(friendUn)$centralization
centr_clo(reportUn)$centralization

centr_betw(adviceUn)$centralization
centr_betw(friendUn)$centralization
centr_betw(reportUn)$centralization

# in and out degree 
degree(advice, mode = "in")
degree(advice, mode = "out")
degree(friend, mode = "in")
degree(friend, mode = "out")
degree(report, mode = "in")
degree(report, mode = "out")

# in and out closeness centrality
closeness(advice, mode = "in")
closeness(advice, mode = "out")
closeness(friend, mode = "in")
closeness(friend, mode = "out")
closeness(report, mode = "in")
closeness(report, mode = "out")


# betweenness centrality
betweenness(advice, directed = TRUE)
betweenness(friend, directed = TRUE)
betweenness(report, directed = TRUE)

# eigenvector centrality
eigen_centrality(advice, directed = TRUE)$vector
eigen_centrality(friend, directed = TRUE)$vector
eigen_centrality(report, directed = TRUE)$vector





