library('igraph')
load('lab2_networks.Rdata')

advice
friend
report

V(advice)
V(friend)
V(report)

E(advice)
E(friend)
E(report)

vcount(advice)
vcount(friend)
vcount(report)

ecount(advice)
ecount(friend)
ecount(report)

is.directed(advice)
is.directed(friend)
is.directed(report)

# density
ecount(advice)/(vcount(advice)*(vcount(advice)-1))
graph.density(advice)                

dyad.census(advice)
unlist(dyad.census(advice))
