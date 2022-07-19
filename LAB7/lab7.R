rm(list = ls())

library(igraph)
library(ergm)
library(intergraph)


load('LAB7/lab7.Rdata')
friend
attr <- read.table('LAB7/attributes.txt', sep = '', head = TRUE)
head(attr)

Y <- get.adjacency(friend, sparse = FALSE)

diag(Y) = NA

# homogeneous binomial random graph model
pmle = mean(Y, na.rm = TRUE)

net <- network(Y, directed = TRUE)

net %v% "Dept" <- attr$Dept
net %v% "Age" <- attr$Age

# mod0 null model SRG
mod0 <- ergm(net ~ edges)

# Same stuff
exp(coef(mod0))
graph.density(friend)
pmle
exp(coef(mod0))/(1+exp(coef(mod0)))

# Non homogeneous SRG
mod1 <- ergm(net ~ edges + sender + receiver)
BIC(mod0, mod1)
AIC(mod0, mod1)

# dyad independence model
mod2 <- ergm(net ~ edges + sender + receiver + mutual, control = control.ergm(seed = 1))

mod3 <- ergm(net ~ edges + receiver + mutual, control = control.ergm(seed = 1))

mod4 <- ergm(net ~ edges + mutual, control = control.ergm(seed = 1))

# model diagnostic
mcmc.diagnostics(mod4)

# include nodal attribute
mod5 <- ergm(net ~ edges + mutual + nodecov('Age') + nodefactor('Dept') + absdiff('Age') + nodematch('Dept'), control = control.ergm(seed = 1))

mod6 <- ergm(net ~ edges + mutual + nodefactor('Dept') + nodematch('Dept'), control = control.ergm(seed = 1))

mod7 <- ergm(net ~ edges + mutual + nodefactor('Dept'), control = control.ergm(seed = 1))

deptNew <- rep(0, 21)
deptNew[attr$Dept == 2] <- 1

net %v% 'deptNew' <- deptNew

mod8 <- ergm(net ~ edges + mutual + nodefactor('deptNew'), control = control.ergm(seed = 1))

mcmc.diagnostics(mod8)

# markov graph model

mod9 <- ergm(net ~ edges + mutual + nodefactor('Dept') + istar(2) + ostar(2) + triangle, control = control.ergm(seed = 1))


mod10 <- ergm(net ~ edges + mutual + nodefactor('Dept') + istar(2) + ostar(2), control = control.ergm(seed = 1))

mod11 <- ergm(net ~ edges + mutual + triangle + nodefactor('Dept') + gwodegree(decay = 1, fixed = TRUE), control = control.ergm(seed = 1))

mod12 <- ergm(net ~ edges + mutual + triangle + nodefactor('Dept') + gwidegree(decay = 1, fixed = TRUE), control = control.ergm(seed = 1))

# social circuit
mod13 <- ergm(net ~ edges + mutual + nodefactor('Dept') + gwesp(decay = 1, fixed = TRUE) + gwdsp(decay = 1, fixed = TRUE), control = control.ergm(seed = 1))


mod14 <- ergm(net ~ edges + mutual + nodefactor('Dept') + gwdsp(decay = 1, fixed = TRUE), control = control.ergm(seed = 1))

sim <- simulate(mod14,  burnin = 1000, nsim = 100, verbose = TRUE, seed = 1)

getStats <- function(g) {
  ig <- asIgraph(g)
  tr <- transitivity(ig)
  ideg <- sd(degree(ig, mode = 'in'))
  odeg <- sd(degree(ig, mode = 'out'))
}

nullDistr <- matrix(NA, 100, 3)

for (b in 1:100) {
  nullDistr[b,] <- getStats(sim[[b]])
}

dev.new()
par(mfrow = c(3,1))


hist(unlist(nullDistr[,1]), xlab = 'transitivity')
abline(v = transitivity(friend), col = 'red')

hist(unlist(nullDistr[,2]), xlab = 'in-degree')
abline(v = sd(degree(friend, mode = 'in')), col = 'red')

hist(unlist(nullDistr[,3]), xlab = 'out-degree')
abline(v = sd(degree(friend, mode = 'out')), col = 'red')

