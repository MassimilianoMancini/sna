rm(list = ls())

library(igraph)

load('LAB6/lab6.Rdata')

# advice network
Y <- get.adjacency(advice, sparse = FALSE)
diag(Y) <- NA

n <-nrow(Y)

degreeIn <- 
degreeOut <- degree(advice, mode = 'out')
rec <- reciprocity(advice)

# best case scenario
p <- graph.density(advice)
simulatedInSD <- c()
simulatedOutSD <- c()
simulatedRec <- c()

for (b in 1:1000) {
  simulatedY <- matrix(rbinom(n^2, 1, p), n, n)
  diag(simulatedY) <- NA
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedInSD[b] <- sd(degree(simulatedGraph, mode = 'in'))
  simulatedOutSD[b] <- sd(degree(simulatedGraph, mode = 'out'))
  simulatedRec[b] <- reciprocity(simulatedGraph)
}

# let's take a look
par(mfrow = c(1, 3))
hist(simulatedInSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(in-degree)', xlim = c(2, 5))
abline(v = sd(degree(advice, mode = 'in')), col = 'red', lwd = 2)
hist(simulatedOutSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(out-degree)', xlim = c(2, 6))
abline(v = sd(degree(advice, mode = 'out')), col = 'red', lwd = 2)
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'reciprocity')
abline(v = reciprocity(advice), col = 'red', lwd = 2)

# formal approach
m <- ecount(advice)
simulatedInSD <- c()
simulatedOutSD <- c()
simulatedRec <- c()
for (b in 1:1000) {
  simulatedY <- matrix(NA, n, n)
  ones <- rep(1, m)
  zeros <- rep(0, n*(n-1) - m)
  simulatedY[col(simulatedY) != row(simulatedY)] <- sample(c(zeros, ones), n*(n-1))
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedInSD[b] <- sd(degree(simulatedGraph, mode = 'in'))
  simulatedOutSD[b] <- sd(degree(simulatedGraph, mode = 'out'))
  simulatedRec[b] <- reciprocity(simulatedGraph)
}

# let's take a look
par(mfrow = c(1, 3))
hist(simulatedInSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(in-degree)', xlim = c(2, 5))
abline(v = sd(degree(advice, mode = 'in')), col = 'red', lwd = 2)
hist(simulatedOutSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(out-degree)', xlim = c(2, 6))
abline(v = sd(degree(advice, mode = 'out')), col = 'red', lwd = 2)
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'reciprocity')
abline(v = reciprocity(advice), col = 'red', lwd = 2)

# non homogeneous simple random graph model
Y <- get.adjacency(advice, sparse = FALSE)
diag(Y) <- NA
y <- c(Y)

ri <- c(row(Y))
ci <- c(col(Y))

model1 <- glm (y ~ factor(ri) + factor(ci), family = 'binomial')

# manually compute tie probability
mu <- model1$coefficients[1]
ai <- c(0, model1$coefficients[2:n])
bi <- c(0, model1$coefficients[(n+1):(2*n-1)])

muij <- mu + rep(ai, n) + rep(bi, each = n)
pij <- exp(muij)/(1+exp(muij))

# easier
pij <- model1$fitted.values

# null distribution
simulatedInSD <- c()
simulatedOutSD <- c()
simulatedRec <- c()

for (b in 1:1000) {
  simulatedY <- matrix(rbinom(n^2, 1, pij), n, n)
  diag(simulatedY) <- NA
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedInSD[b] <- sd(degree(simulatedGraph, mode = 'in'))
  simulatedOutSD[b] <- sd(degree(simulatedGraph, mode = 'out'))
  simulatedRec[b] <- reciprocity(simulatedGraph)
}

# let's take a look
par(mfrow = c(1, 3))
hist(simulatedInSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(in-degree)')
abline(v = sd(degree(advice, mode = 'in')), col = 'red', lwd = 2)
hist(simulatedOutSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(out-degree)', xlim = c(2, 6))
abline(v = sd(degree(advice, mode = 'out')), col = 'red', lwd = 2)
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'reciprocity')
abline(v = reciprocity(advice), col = 'red', lwd = 2)

# p-values
mean(simulatedInSD < sd(degree(advice, mode = 'in')))
mean(simulatedOutSD < sd(degree(advice, mode = 'out')))
mean(simulatedRec > reciprocity(advice))

# formal approach

alternateRectangle <- function(Y, k) {
  Y1 <- matrix(c(0,1,1,0), 2, 2)
  Y2 <- matrix(c(1,0,0,1), 2, 2)
  
  n <- nrow(Y)
  
  for (s in 1:k) {
    ij <- sample(1:n, 4)
    rows <- ij[1:2]
    cols <- ij[3:4]
    
    if (all(Y[rows, cols] == Y1)) {
      Y[rows, cols] <- Y2
    } else if (all(Y[rows, cols] == Y2)) {
      Y[rows, cols] <- Y1
    }
  }
  return (Y)
}

simulatedRec <- c()

for (b in 1:1000) {
  simulatedY <- alternateRectangle(Y, 100)
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedRec[b] <- reciprocity(simulatedGraph)
}

# let's take a look
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'reciprocity')
abline(v = reciprocity(advice), col = 'red', lwd = 2)

# transitivity

simulatedTran <- c()

for (b in 1:1000) {
  simulatedY <- alternateRectangle(Y, 100)
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedTran[b] <- transitivity(simulatedGraph)
}

# let's take a look
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'transitivity')
abline(v = reciprocity(advice), col = 'red', lwd = 2)


# friend network
Y <- get.adjacency(friend, sparse = FALSE)
diag(Y) <- NA

n <-nrow(Y)

degreeIn <- 
  degreeOut <- degree(friend, mode = 'out')
rec <- reciprocity(friend)

# best case scenario
p <- graph.density(friend)
simulatedInSD <- c()
simulatedOutSD <- c()
simulatedRec <- c()

for (b in 1:1000) {
  simulatedY <- matrix(rbinom(n^2, 1, p), n, n)
  diag(simulatedY) <- NA
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedInSD[b] <- sd(degree(simulatedGraph, mode = 'in'))
  simulatedOutSD[b] <- sd(degree(simulatedGraph, mode = 'out'))
  simulatedRec[b] <- reciprocity(simulatedGraph)
}

# let's take a look
par(mfrow = c(1, 3))
hist(simulatedInSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(in-degree)', xlim = c(2, 5))
abline(v = sd(degree(friend, mode = 'in')), col = 'red', lwd = 2)
hist(simulatedOutSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(out-degree)', xlim = c(2, 6))
abline(v = sd(degree(friend, mode = 'out')), col = 'red', lwd = 2)
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'reciprocity')
abline(v = reciprocity(friend), col = 'red', lwd = 2)

# formal approach
m <- ecount(friend)
simulatedInSD <- c()
simulatedOutSD <- c()
simulatedRec <- c()
for (b in 1:1000) {
  simulatedY <- matrix(NA, n, n)
  ones <- rep(1, m)
  zeros <- rep(0, n*(n-1) - m)
  simulatedY[col(simulatedY) != row(simulatedY)] <- sample(c(zeros, ones), n*(n-1))
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedInSD[b] <- sd(degree(simulatedGraph, mode = 'in'))
  simulatedOutSD[b] <- sd(degree(simulatedGraph, mode = 'out'))
  simulatedRec[b] <- reciprocity(simulatedGraph)
}

# let's take a look
par(mfrow = c(1, 3))
hist(simulatedInSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(in-degree)', xlim = c(2, 5))
abline(v = sd(degree(friend, mode = 'in')), col = 'red', lwd = 2)
hist(simulatedOutSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(out-degree)', xlim = c(2, 6))
abline(v = sd(degree(friend, mode = 'out')), col = 'red', lwd = 2)
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'reciprocity')
abline(v = reciprocity(friend), col = 'red', lwd = 2)

# non homogeneous simple random graph model
Y <- get.adjacency(friend, sparse = FALSE)
diag(Y) <- NA
y <- c(Y)

ri <- c(row(Y))
ci <- c(col(Y))

model1 <- glm (y ~ factor(ri) + factor(ci), family = 'binomial')

# manually compute tie probability
mu <- model1$coefficients[1]
ai <- c(0, model1$coefficients[2:n])
bi <- c(0, model1$coefficients[(n+1):(2*n-1)])

muij <- mu + rep(ai, n) + rep(bi, each = n)
pij <- exp(muij)/(1+exp(muij))

# easier
pij <- model1$fitted.values

# null distribution
simulatedInSD <- c()
simulatedOutSD <- c()
simulatedRec <- c()

for (b in 1:1000) {
  simulatedY <- matrix(rbinom(n^2, 1, pij), n, n)
  diag(simulatedY) <- NA
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedInSD[b] <- sd(degree(simulatedGraph, mode = 'in'))
  simulatedOutSD[b] <- sd(degree(simulatedGraph, mode = 'out'))
  simulatedRec[b] <- reciprocity(simulatedGraph)
}

# let's take a look
par(mfrow = c(1, 3))
hist(simulatedInSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(in-degree)')
abline(v = sd(degree(friend, mode = 'in')), col = 'red', lwd = 2)
hist(simulatedOutSD, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'sd(out-degree)', xlim = c(2, 6))
abline(v = sd(degree(friend, mode = 'out')), col = 'red', lwd = 2)
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'reciprocity')
abline(v = reciprocity(friend), col = 'red', lwd = 2)

# p-values
mean(simulatedInSD < sd(degree(friend, mode = 'in')))
mean(simulatedOutSD < sd(degree(friend, mode = 'out')))
mean(simulatedRec > reciprocity(friend))

# formal approach

alternateRectangle <- function(Y, k) {
  Y1 <- matrix(c(0,1,1,0), 2, 2)
  Y2 <- matrix(c(1,0,0,1), 2, 2)
  
  n <- nrow(Y)
  
  for (s in 1:k) {
    ij <- sample(1:n, 4)
    rows <- ij[1:2]
    cols <- ij[3:4]
    
    if (all(Y[rows, cols] == Y1)) {
      Y[rows, cols] <- Y2
    } else if (all(Y[rows, cols] == Y2)) {
      Y[rows, cols] <- Y1
    }
  }
  return (Y)
}

simulatedRec <- c()

for (b in 1:1000) {
  simulatedY <- alternateRectangle(Y, 100)
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedRec[b] <- reciprocity(simulatedGraph)
}

# let's take a look
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'reciprocity')
abline(v = reciprocity(friend), col = 'red', lwd = 2)

# transitivity

simulatedTran <- c()

for (b in 1:1000) {
  simulatedY <- alternateRectangle(Y, 100)
  simulatedGraph <- graph_from_adjacency_matrix(simulatedY)
  simulatedTran[b] <- transitivity(simulatedGraph)
}

# let's take a look
hist(simulatedRec, breaks = 20, col = 'lightgray', main = '', prob = TRUE, xlab = 'transitivity')
abline(v = reciprocity(friend), col = 'red', lwd = 2)




