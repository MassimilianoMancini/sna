#--------------------------------------------#
#--------------------------------------------#
#                  Lab 6                     #
#           Simple random graph model:       #
#             modelling heterogeneity        #
#--------------------------------------------#
#--------------------------------------------#

rm(list = ls())

# -----------------------#
# ------ 1. Set-up ----- #
# -----------------------#
library(igraph)

# ----- SRG model for the advice network  ------

#--- Load the data ----

load("HighTechNetworks.Rdata")
advice.net

# let us extract the adjacency matrix
Y = get.adjacency(advice.net, sparse = F)
diag(Y) = NA


# number of nodes
n = nrow(Y)
n

#---- assessing significance in terms of in-degree, out-degree, and reciprocity ---   

# observed statistics
plot(advice.net, edge.arrow.size = 0.1, vertex.size = 20)

din = degree(advice.net, mode = "in")
din
dout = degree(advice.net, mode = "out")
dout

par(mfrow = c(1,2))
hist(din, breaks = 10, col = "lightgray", main = "", prob = TRUE, xlab = "in-degree")
hist(dout, breaks = 10, col = "lightgray", main = "", prob = TRUE, xlab = "out-degree")

reciprocity(advice.net)

# ----- best-case scenario -----

p = graph.density(advice.net)
sdIn.sim = sdOut.sim = recip.sim =  c()
for(b in 1:1000){
  tmp = rbinom(n^2,1,p)  
  Y.sim = matrix(tmp, n,n); diag(Y.sim) = NA
  g.sim = graph_from_adjacency_matrix(Y.sim)
  sdIn.sim[b] = sd(degree(g.sim, mode = "in"))
  sdOut.sim[b] = sd(degree(g.sim, mode = "out"))
  recip.sim[b] = reciprocity(g.sim)
}

# graphical comparison
par(mfrow = c(1,3))
hist(sdIn.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(in-degree)", xlim = c(2,6))
abline(v = sd(din), col = "red", lwd = 2)


hist(sdOut.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(out-degree)", xlim = c(2,9))
abline(v = sd(dout), col = "red", lwd = 2)

hist(recip.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "reciprocity")
abline(v = reciprocity(advice.net), col = "red", lwd = 2)


# is the G(n,p) model coherent with the observed network in terms of in- and out-degree variability?
# And in terms of reciprocity?



# ---- formal approach -- conditional uniform distribution ----
m = ecount(advice.net)
B = 1000
sdIn.sim = sdOut.sim = recip.sim = c()
for(b in 1:B){
  ones = rep(1, m)
  zeros = rep(0, n*(n-1) - m)
  all = c(ones, zeros)
  tmp = sample(all, n*(n-1))
  Y.sim = matrix(tmp, n,n)
  diag(Y.sim) = NA
  g.sim = graph_from_adjacency_matrix(Y.sim)
  sdIn.sim[b] = sd(degree(g.sim, mode = "in"))
  sdOut.sim[b] = sd(degree(g.sim, mode = "out"))
  recip.sim[b] = reciprocity(g.sim)
  
}

# graphical comparison
par(mfrow = c(1,3))
hist(sdIn.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(in-degree)", xlim = c(2,6))
abline(v = sd(din), col = "red", lwd = 2)


hist(sdOut.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(out-degree)", xlim = c(2,9))
abline(v = sd(dout), col = "red", lwd = 2)


hist(recip.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(out-degree)")
abline(v = reciprocity(advice.net), col = "red", lwd = 2)

# is the G(n,p) model coherent with the observed network in terms of in- and out-degree variability?
# And in terms of reciprocity?


# ------ non-homogeneous SRG model ----

#-------- model parameter estimation -------
?glm
# how to use the function? 
# formula, family, data
# fitted values



# response variable 
Y = as.matrix(get.adjacency(advice.net), sparse = F)
diag(Y) = NA
y = c(Y)

# covariates -- sender and receiver effects
rowIdx = row(Y)
colIdx = col(Y)
rowIdx[1:4, 1:4]
colIdx[1:4, 1:4]

rowidx = c(rowIdx)
colidx = c(colIdx)

# estimate the parameters
mod = glm(y ~ factor(rowidx) + factor(colidx), family = "binomial")
mod

# more..
summary(mod)


# compute tie probability
mu = mod$coefficients[1]
ai = mod$coefficients[2:n]
ai = c(0, ai)
bj = mod$coefficients[(n+1):(2*n-1)]
bj = c(0, bj)

ai.all = rep(ai, n)
bj.all = rep(bj, each = n)

muij = mu + ai.all + bj.all
pij = exp(muij)/(1+exp(muij))

# or more simply 
pij = mod$fitted.values
# ATT: when using this latter approach, only off diagonal elements are given!


# ---------- assessing significance for the non-homogeneous SRG model --------

# null distribution
sdIn.sim = sdOut.sim = recip.sim = c()
for(b in 1:1000){
  tmp = rbinom(n*(n-1),1,pij)
  Y.sim = matrix(, n,n); Y.sim[row(Y.sim) != col(Y.sim)] = tmp
  g.sim = graph_from_adjacency_matrix(Y.sim)
  sdIn.sim[b] = sd(degree(g.sim, mode = "in"))
  sdOut.sim[b] = sd(degree(g.sim, mode = "out"))
  recip.sim[b] = reciprocity(g.sim)
}

# graphical comparison
par(mfrow = c(1,3))
hist(sdIn.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(in-degree)")
abline(v = sd(din), col = "red", lwd = 2)

hist(sdOut.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(out-degree)")
abline(v = sd(dout), col = "red", lwd = 2)

# graphical comparison
hist(recip.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "reciprocity")
abline(v = reciprocity(advice.net), col = "red", lwd = 2)



# p-value
mean(sdIn.sim < sd(din))
mean(sdOut.sim < sd(dout))
mean(recip.sim > reciprocity(advice.net))
# is the G(n,p_ij) model coherent with the observed network in terms of in- and out-degree variability?
# is the G(n,p_ij) model coherent with the observed network in terms of reciprocity?


# --------- formal approach -- conditional uniform distribution ---------
# we need to condition to the sufficient statistics
# 1) m = y..
# 2) y_i., i = 1, ..., n
# 3) y_.j, j = 1, ..., n

# simulate from the conditional uniform distribution
# simulating function -- alternating rectangles
aRect.fnc = function(Y, k){
  
  # Y = adjacency matrix
  # k = n. of steps in the alternating rectangles algorithm 
  
  Y1 = matrix(c(0,1,1,0), 2, 2)
  Y2 = 1 - Y1
  
  n = nrow(Y)
  
  for(s in 1:k){
    # draw 4 distinct indexes
    # two rows and two columns
    ij = sample(1:n,4,replace = F)
    
    # select the corresponding sub-matrix
    rows = ij[1:2]
    cols = ij[3:4]
    Yij = Y[rows, cols]
    
    # perturbation
    if(all(Yij == Y1)) Yij = Y2 else if(all(Yij == Y2))  Yij = Y1
      
    Y[rows, cols] = Yij
  }
    
  return(Y)
}


sdIn.sim = sdOut.sim = recip.sim = c()
for(b in 1:1000){
  Y.sim = aRect.fnc(Y, 100)
  # print number of perturbed elements in Y.sim 
  cat(sum(Y != Y.sim, na.rm = T), "*", sep="")
  g.sim = graph_from_adjacency_matrix(Y.sim)
  sdIn.sim[b] = sd(degree(g.sim, mode = "in"))
  sdOut.sim[b] = sd(degree(g.sim, mode = "out"))
  recip.sim[b] = reciprocity(g.sim)

}



# graphical comparison
par(mfrow = c(1,3))
hist(sdIn.sim)
abline(v = sd(din), col = "red", lty = 2)
sdIn.sim

hist(sdOut.sim)
abline(v = sd(dout), col = "red", lty = 2)
sdOut.sim

hist(recip.sim)
abline(v = reciprocity(advice.net), col = "red", lty = 2)

# p-value
mean(sdIn.sim > sd(din))
mean(sdOut.sim > sd(dout))
mean(recip.sim > reciprocity(advice.net))

# let us also look at the transitivity of the network
cl.sim = c()
for(b in 1:1000){
  Y.sim = aRect.fnc(Y, 100)
  # print number of perturbed elements in Y.sim 
  cat(sum(Y != Y.sim, na.rm = T), "*", sep="")
  g.sim = graph_from_adjacency_matrix(Y.sim)
  cl.sim[b] = transitivity(g.sim)
}


# graphical comparison
par(mfrow = c(1,1))
hist(cl.sim)
abline(v = transitivity(advice.net), col = "red", lty = 2)

# p-value
mean(cl.sim < transitivity(advice.net))
# is the G(n,p_ij) model coherent with the observed network in terms of transitivity?


# REPEAT THE ABOVE ANALYSIS FOR THE FRIENDSHIP NETWORK