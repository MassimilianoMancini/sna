#--------------------------------------------#
#--------------------------------------------#
#                  Lab 7                     #
#     Exponential Random Graph Models        #
#--------------------------------------------#
#--------------------------------------------#

rm(list = ls())

# -----------------------#
# ------ 1. Set-up ----- #
# -----------------------#
library(igraph)

#--- Load the data ----

# network
load("friend.Rdata")
friend.net

# nodal attributes
attr = read.table("attributes.txt", sep = "", head = TRUE)
head(attr)

# how does the network look like?
plot(friend.net, vertex.size = 20, edge.arrow.size = 0.5)


# AIM: build an adequate ERGM
# ---------------------------------
Y = get.adjacency(friend.net, sparse = F)
diag(Y)  = NA


# -------- mod1: homogeneous binomial random graph model --------
Y = get.adjacency(friend.net, sparse = F)
diag(Y) = NA
p.MLE = mean(Y, na.rm = T)
p.MLE


# how can we obtain an equivalent result by considering the  ERGM specification?
# ergm package
library(ergm)
?ergm
# what should we note? 
# 1) the function allows to estimate an ergm 
# 2) estimation methods which are implemented are: 
#    - maximum-pseudo likelihood (MPLE)
#    - approximate Monte Carlo Markov Chain ML (MLE)
#    - CD -- another approximate estimating procedure (experimental)
# 3) In the formula, the response is a "network" object 


# to start, let us transform the adjacency matrix in a 
# network object
?network
net = network(Y, directed = T)
class(net)
net

# add network and edge attributes
net %v% "Dept" = attr$Dept
net %v% "Age" = attr$Age

# let us give a look at the properties
net

# let us assume, we also have some edge attributes --  
# let us create a "fake" attribute
attr.test = rnorm(21*20, 0, 1)
net %e% "Test" = attr.test

# let us give a look at the properties
net

#************************************
# Let's estimate the NULL model: BRG
#************************************
# let us estimate parameters of a SRG model via the ergm function
# which is the sufficient statistic for theta? 
# y.. -> the number of observed edges
mod0 = ergm(net ~ edges)
mod0
# look in more depth
summary(mod0)

# which information do we have? 
# 1) formula: recall the estimated model
# 2) how many iterations for convergence? 
# 3) estimates -- in this case, as independence holds, 
#    a standar ML approach is used (besides the title!)
# 4) as this model simply corresponds to a logistic regression
#    we can interpret the results as usual

# edge parameter
# 1) is it significant? Yes -> there is a significant difference 
# between Pr(Y_ij = 1) and Pr(Y_ij = 0)  
# That is, these quantities are significantly different from 0.5
# 2) sign? Negative ->  Pr(Y_ij = 1)< 0.5

# indeed...
exp(coef(mod0))
# the odds of observing a relation between two randomly 
# selected nodes is about 70% lower than that of not observing it

# indeed... 
graph.density(friend.net)
# or
p.MLE
# or, via the expit function
exp(coef(mod0))/(1+exp(coef(mod0)))


# -------------------------------------------
# let us move towards the non-homogeneous BRG
# -------------------------------------------
# which sufficient statistics for theta?
# n. of edges
# in-degrees -> receiver effects
# out-degrees -> sender effects
mod1 = ergm(net ~ edges + sender + receiver)
summary(mod1)

# how do we interpret the results? 
# 1) significant parameters? 
# 2) sign of the parameters? 

# which model is more appropriate? 
BIC(mod0, mod1)
AIC(mod0, mod1)
# as expected, bic is more conservative

# -----------------------------------------------
# let us now consider the dyad independence model
# -----------------------------------------------
# which sufficient statistics for theta?
# n. of edges
# in-degrees -> receiver effects
# out-degrees -> sender effects
# n. of mutual relations 

# ! ATT -- MCMCMLE used for estimation 
# let us fix a seed to replicate results
mod2 = ergm(net ~ edges + sender + receiver + mutual, 
            control = control.ergm(seed = 1))

summary(mod2)
# cannot use BIC and AIC for model selection because of some infinite sender coefficients
#let us remove this term from the model 
mod3 = ergm(net ~ edges + receiver + mutual, control = control.ergm(seed = 1))
summary(mod3)
# many parameters not significantly different from zero
# let us also exclude the receiver parameter
mod4 = ergm(net ~ edges + mutual)
summary(mod4)
# which model should we prefer? let us use the BIC
BIC(mod0, mod3, mod4)

# mod4 seems the best
summary(mod4)
# interpret the results
# 1) edge parameter? 
#    negative and significant - less ties than expected by chance
#    -> non-ties are more frequent than ties
# 2) mutuality parameter? 
#    positive and significant --- mutuality plays a role ->
#    there is a positive tendency to return ties -> 
#    more mutual ties than expected by chance

# ----------------
# model diagnostic
# ----------------
# has the model converged?
mcmc.diagnostics(mod4)
# Look at standard results you consider for any MCMC estimation: 
# well-mixed, stationary chains
# results look pretty ok
# 1) the chains explore the parameter space
# 2) posterior distributions are almost bell-shaped

# If the MCMC procedure does not converge
# A) the model is ill-specified
# B) try to improve things by increasing the length of the MCMC routine

# ---------------------------------------------
# let us include in the model nodal attributes
# ---------------------------------------------
# main effects 
# quantitative attributes -- nodecov(attr)
# qualitative attributes -- nodefactor(attr)

# homophily effects
# quantitative attributes -- absdiff(attr) - sum_ij[abs(attr_i - attr_j)y_ij]
# qualititave attributes -- nodematch(attr) 
mod5 = ergm(net ~ edges + mutual + nodecov("Age") + nodefactor("Dept") + 
              absdiff("Age") + nodematch("Dept"), control = control.ergm(seed = 1)) 

# let us look at the results
summary(mod5)
# age doesn't play a role, let us exclude it from the model

mod6 = ergm(net ~ edges + mutual + nodefactor("Dept") + nodematch("Dept"), control = control.ergm(seed = 1)) 
summary(mod6)
# only the main effect plays a role, let us exclude the homophily effect
mod7 = ergm(net ~ edges + mutual + nodefactor("Dept"), control = control.ergm(seed = 1)) 
summary(mod7)
# how do we interpret results?
# negative and significant edges
# positive and significant mutuality
# positive and significant effect for the main effect of Dep2 

# let us build a new attribute
Dept.new = rep(0, 21)
Dept.new[attr$Dept == 2] = 1

# add the attribute to the network
net %v% "Dept.new" = Dept.new

# estimate again the model
mod8 = ergm(net ~ edges + mutual + nodefactor("Dept.new"), control = control.ergm(seed = 1)) 
summary(mod8)
# how do we interpret results?
# negative and significant edges
# positive and significant mutuality
# positive and significant effect for the main effect of Dep2 

# let us verify convergence
mcmc.diagnostics(mod8)
# again, convergence seems to be acheaved


# -------------------------------------------
# let us move towards the Markov graph model
# -------------------------------------------
# Let us add to the model the triangle term, the in- and the out-stars of order 2
# (indicating the tendency to form clusters in the network)

mod9 = ergm(net ~ edges + mutual + nodefactor("Dept.new") + istar(2) + ostar(2) + triangle, 
             control = control.ergm(seed = 1))
# the model cannot be estimated due to degeneracy issues

# let us try to remove triangles
mod10 = ergm(net ~ edges + mutual + nodefactor("Dept.new") + istar(2) + ostar(2), 
             control = control.ergm(seed = 1))
# still, the model suffer from some estimation issues

# let us try to solve the issue by considering the alternating k-star term
# --------------------------------------------------------------------------
# ergm implements some modified forms of such statistics

# for undirected networks 
# gwdegree(decay, fixed = FALSE) -- decay = log(lambda/(lambda-1))
# gwdegree -- geometrically weighted degree distribution

# for directed networks 
# gwidegree(decay, fixed = FALSE) -- decay = log(lambda/(lambda-1))
# gwodegree(decay, fixed = FALSE) -- decay = log(lambda/(lambda-1))
# gwidegree/gwodegree -- geometrically weighted in-/out- degree distribution

# positive estimates -- centralized network --  few nodes with high degree
# negative estimates -- non centralized network

# a standard choice for decay is 1, but model selection can be used!
mod10 = ergm(net ~ edges + triangle + nodefactor("Dept.new") + gwidegree(decay = 1, fixed = TRUE), 
             control = control.ergm(seed = 1))
mod11 = ergm(net ~ edges + triangle + nodefactor("Dept.new") + gwodegree(decay = 1, fixed = TRUE), 
             control = control.ergm(seed = 1))
mod12 = ergm(net ~ edges + triangle + nodefactor("Dept.new") + gwidegree(decay = 1, fixed = TRUE)
             + gwodegree(decay = 1, fixed = TRUE), 
             control = control.ergm(seed = 1))
# again, the model cannot be estimated

# ----------------------------------------
# let us consider the social circuit model 
# ----------------------------------------
# alternating k-triangles --> gwesp(decay = 0, fixed = FALSE) 
# geometrically weighted edge-wise shared partners
# the corresponding parameter expresses the tendency for tied nodes 
# to have multiple shared partners

# alternating k-2-paths --> gwdsp(decay = 0, fixed = FALSE)
# geometrically weighted dyad-wise shared partners
# the corresponding parameter expresses the tendency for dyads 
# (whether tied or not) to have multiple shared partners

mod13 = ergm(net ~ edges + mutual + nodefactor("Dept.new") + 
               gwesp(decay = 1, fixed = T) + 
                gwdsp(decay = 1, fixed = T), control = control.ergm(seed=1))
summary(mod13)

# let us remove the gwesp term
mod14 = ergm(net ~ edges + mutual + nodefactor("Dept.new") + 
               gwdsp(decay = 1,fixed = T), control = control.ergm(seed=1))
summary(mod14)


# how can we interpret the parameters?
# 1) the edge parameter is significantly different from zero and is negative 
#    --> less ties those expected by chance 
# 2) the mutuality parameter is significantly different from zero and is positive
#   --> positive tendency to return ties in the network
# 3) the parameter associated to the nodal attribute is significant and positive
#   --> nodes in department 2 are more active than others
# 3) the gwdsp parameter is significant and negative
#    -> lower tendency to form clusters than that expected by chance

# check convergence
mcmc.diagnostics(mod14)


# still look at model selection via BIC
BIC(mod14, mod8)
# mod14 seems to be the optimal choice

# let us evaluate the goodness of fit
# ----------------------------------
# simulate from the model
sim = simulate(mod14, burnin = 1000, nsim = 100, verbose = TRUE, seed = 1)

# let us assume we want to verify whether the model is appropriate to represent the degree and the 
# transitivity in the network

# install.packages("intergraph")
library(intergraph)
?asIgraph

fnc = function(xx){
  ig = asIgraph(xx)
  tr = transitivity(ig)
  ideg = sd(degree(ig, mode = "in"))
  odeg = sd(degree(ig, mode = "out"))
  return(c(tr, ideg, odeg))
}

null.distr = matrix(,100,3)
for(b in 1:100){
  null.distr[b,]  = fnc(sim[[b]])
}
dev.new()
par(mfrow = c(3,1))
hist(unlist(prova[,1]), xlab = "transitivity"); abline(v = transitivity(friend.net), col = "red")
hist(unlist(prova[,2]), xlab = "in-degree"); abline(v = sd(degree(friend.net, mode = "in")), col = "red")
hist(unlist(prova[,3]), xlab = "out-degree"); abline(v = sd(degree(friend.net, mode = "out")), col = "red")
