#-------------------------------------------------#
#-------------------------------------------------#
#                      Lab 1                      #
#  Create, manipulate, and visualize network data #
#-------------------------------------------------#
#-------------------------------------------------#


# Let us install the igraph package
# install.packages("igraph")
library(igraph)
?igraph


# -----------------------------------
# Creating a network object 
# --------------------------

# "By hand"
# ----------

# UNDIRECTED NETWORK
# ------------------ 
?graph_from_literal
g = graph_from_literal(1-2, 1-3, 2-3, 2-4, 3-5, 4-5, 4-6, 4-7, 5-6, 6-7)
g
# which information? 
# first line 
# a) object type -- IGRAPH object 
# b) type of network --  UN: undirected, named, network (name of vertices specified) 
# c) number of vertices and edges -- 7 and 10, respectively
# second line
# a) attributes -- name - vertex, character, attribute
# third line 
# a) edges 

# or more simply
# g = graph_from_literal(1-2:3, 2-3:4, 3-5, 4-5:6:7, 5-6, 6-7)
# g

# graphical representation
plot(g)

# which are the nodes/vertices?
V(g)
# which are the edges?
E(g)


# DIRECTED NETWORK
# -----------------
# - denotes the sending node
# + denotes the receiving node

gd = graph_from_literal(1-+2, 1-+3, 2++3, 4++1, 4+-3, 5-+1)
gd
# which information? 
# first line 
# a) object type -- IGRAPH object 
# b) type of network --  DN: directed, named, network (name of vertices is specified)
# c) number of vertices and edges -- 5 and 8, respectively
# second line
# a) attributes -- name - vertex, character, attribute
# third line 
# a) edges 

# graphical representation
plot(gd)
set.seed(1)
plot(gd)


# how to give names to vertices? 
gd.lab = graph_from_literal(A-+B, A-+C, B++C, D++A, D+-C, E-+A)
plot(gd.lab)
# or defining labels "a posteriori"
gd.lab = graph_from_literal(1-+2, 1-+3, 2++3, 4++1, 4+-3, 5-+1)
set.seed(1)
plot(gd.lab)
V(gd.lab)$name = c("A", "B", "C", "D", "E")
set.seed(1)
plot(gd.lab)


# networks with isoleted nodes
# ----------------------------
gd.iso = graph_from_literal(1-+2, 1-+3, 2++3, 4++1, 4+-3, 5-+1, 6)
gd.iso
# we are in the presence of a directed network with 6 vertices and 8 edges 
set.seed(1)
plot(gd.iso)


# network with 2 components
# -------------------------
# a) 1->2, 1->3, 1->4, 
# b) 2<->3 
# c) 3->4 
# d) 5->1
# e) 6<->7, 6<->8

gd.2comp = graph_from_literal(1-+2, 1-+3, 1-+4, 2++3, 3-+4, 5-+1, 6++7, 6++8)
gd.2comp
# we are in the presence of a directed network with 8 vertices and 11 edges 

# or more simply
# gd.2comp = graph_from_literal(1-+2:3:4, 2++3, 3-+4, 5-+1, 6++7:8)
# gd.2comp

# graphical representation
set.seed(1)
plot(gd.2comp)

# how to verify?
components(gd.2comp)


# how many nodes? 
# ----------------
vcount(gd.2comp)

# how many edges? 
# -------------------
ecount(gd.2comp)

# is the network directed?
is.directed(gd.2comp)

# moving from a directed to an undirected network
?as.undirected
par(mfrow = c(2,2), mar = c(1,1,1,1))
plot(gd.2comp)

# an undirected relation is created if it exists at least one directed 
# relation between nodes
plot(as.undirected(gd.2comp, mode = "collapse"))
# an undirected relation is created if it exists a bidirectional 
# relation between nodes
plot(as.undirected(gd.2comp, mode = "mutual"))
# an undirected relation is created for each directed relation 
# multiple edges can be created
plot(as.undirected(gd.2comp, mode = "each"))


# moving from an undirected to a directed network
?as.directed
par(mfrow = c(1,3), mar = c(0,0,0,0))
plot(g)
# bidirectional relations are created for each undirected one
plot(as.directed(g, mode = "mutual"))
# a directed relation with arbitrary direction is created for each undirected one
plot(as.directed(g, mode = "arbitrary"))



# *********************************************************************************************

# MATRIX REPRESENTATION
# ----------------------

# from an igraph object
# ---------------------
Y = get.adjacency(g)
Y
# the function retun an object of class dgCMatrix
# This is useful to deal with sparse numeric matrices to save memory

# to obtain a standard adjacency matrix
Y = get.adjacency(g, sparse = FALSE)
Y
# diagonal elements are set to 0 -- transform them in NAs
diag(Y) = NA
Y

# "by hand"
# ---------
Y = matrix(0, 7, 7)
# for UNDIRECTED networks, we may consider the upper/lower triangle of Y only
tmp = c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1)
Y[lower.tri(Y)] = tmp
Y
# what about the upper/lower triangle?
Y = Y + t(Y)
diag(Y) = NA
Y
# is it symmetric? Check!
sum(Y != t(Y), na.rm=TRUE)

# for DIRECTED networks, we need to consider both the upper/lower triangle of Y
tmp = c(0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0)
Yd = matrix(tmp, 5, 5)
Yd
diag(Yd) = NA
# is it symmetric?
sum(Yd != t(Yd), na.rm = TRUE)
# != 0 -- not symmetric


# from the matrix to the igraph object
# ------------------------------------
?graph_from_adjacency_matrix
g2 = graph_from_adjacency_matrix(Y, mode = "undirected")
gd2 = graph_from_adjacency_matrix(Yd, mode = "directed")

par(mfrow = c(2,2), mar = c(1,1,1,1))
# let us set a seed to fix position of the nodes on the graph
set.seed(1)
plot(g)
set.seed(1)
plot(g2)

set.seed(1)
plot(gd)
set.seed(1)
plot(gd2)

par(mfrow = c(1,1))


# create a weighted graph from the adjacency matrix
tmp = c(0, 0, 3, 1,0, 1, 0, 2, 0, 0, 1, 3, 0, 0, 0, 4, 0, 2, 0, 0, 5, 0, 0, 0, 1)
Yd = matrix(tmp, 5, 5)
Yd
diag(Yd) = NA

wg = graph_from_adjacency_matrix(Yd, mode = "directed", weighted = T, diag = F)
wg
# which information?
# first line 
# a) object type -- IGRAPH object 
# b) type of network --  D-W: directed weighted
# c) number of vertices and edges -- 5 and 9, respectively
# second line
# a) attributes -- weight - edge (e), numeric (n), attribute
# third line 
# a) edges 

# let us give a label to the vertices of wg
V(wg)$name = c("Pippo", "Pluto", "Paperino", "Cip", "Ciop")
wg
# when printing graph information, two attributes are now present
# one for the vertices (name), one for edges

# let us give a name to the graph
wg$name = "Disney"
wg
# when printing graph information, three attributes are now present
# one for the graph (name), one for vertices (name), and one for edges

plot(wg)


# ******************************************************************************

# ------------------------
# improve network drawing
# ------------------------

# let us consider a simple graph
gd = graph_from_literal(1-+2:3, 2++3)
V(gd)$name = c("Francesca", "Paola", "Giovanni")
plot(gd)


# modify the nodes 
# -----------------
# increase the size
plot(gd, vertex.size = 35)

# change the color 
plot(gd, vertex.size = 35, vertex.color = "cyan")

# change the shape 
plot(gd, vertex.size = 35, vertex.color = "cyan", vertex.shape = "square")

# change the label color
plot(gd, vertex.size = 35, vertex.color = "cyan", vertex.shape = "square",
     vertex.label.color = "red")

# increase the font
plot(gd, vertex.label.cex = 1.3, 
     vertex.size = 50, vertex.color = "cyan", vertex.shape = "square",
     vertex.label.color = "red")



# modify the ties
# ----------------
# change the color
plot(gd, vertex.size = 35, edge.color = "blue")

# change the arrow size
plot(gd, vertex.size = 35, edge.color = "blue", edge.arrow.size = 0.5)

# change the edge width
plot(gd, vertex.size = 35, edge.color = "blue", edge.arrow.size = 0.5, 
     edge.width = 2)


# for more options
?igraph.plotting


# ------------------------------
# Decorating a graph: attributes 
# ------------------------------
# Let us consider a simple graph
gd = graph_from_literal(1-+2:3, 2++3)
plot(gd)

# nodal attributes
# ----------------
# give a color to the nodes
V(gd)$color = c("red", "green", "cyan")
gd
plot(gd)


# give a color to the nose based on some other attributes 
V(gd)$name = c("Francesca", "Paola", "Giovanni")
V(gd)$gender = c("F", "F", "M")
gd

# based on the gender
V(gd)[gender == "F"]$color = "cyan"
V(gd)[gender == "M"]$color = "green"
V(gd)$color

# graphical representation
plot(gd, vertex.size = 45)


# tie attributes 
# ----------------
# give a width to edges
E(gd)$width = c(1,5,1,2)
gd
# give a different width to ties based on the width  
plot(gd, vertex.size = 45)


# give a width based on another attribute
E(gd)$other = 1:4
plot(gd, vertex.size = 45, edge.width = E(gd)$other)


# weighted adjacency matrix from gd
get.adjacency(gd, attr = "other")
# or 
get.adjacency(gd, attr = "other", sparse = FALSE)


# network attributes 
# we may also define attributes for the whole graph
gd$name = "My first directed graph"
plot(gd, vertex.size = 45, main = gd$name, edge.width = E(gd)$other)


#---------------------------------#
#         Real-data examples       #
#---------------------------------#
library(igraphdata)
data(package='igraphdata')

# Zachary's karate club network
data("karate")

?karate
# structure of the dataset
karate
# how to read the output?
# -----------------------
# The first line always starts with IGRAPH, showing you that the object is an igraph object
# The first letter distinguishes between directed  networks ('DN') and undirected networks ('UN'). 
# The second letter is 'N' for named graphs, i.e. graphs having NAME among vertex attribute. 
# The third letter is 'W' for weighted graphs, i.e. graphs having WEIGHT among the edge attributes

# Then, the name of the graph 
# From the second line, the attributes of the graph are listed, separated by a comma. 
# After the attribute names, the kind of the attribute - graph ('g'), vertex ('v') or edge ('e') - is reported
# and the type of the attribute as well, character ('c'), numeric ('n'), logical ('l'), or other ('x').



# how to list attributes?
# --------------------------------
# graph attributes
# -----------------
graph.attributes(karate)
# or only the names
names(graph.attributes(karate))

# vertex attributes
# -----------------
vertex.attributes(karate)
names(vertex.attributes(karate))

# extract members' names
V(karate)$name
# extract members' faction
V(karate)$Faction

# edge attributes
# -----------------
edge.attributes(karate)
names(edge.attributes(karate))

# extract edges' weights
E(karate)$weight


# visualize the graph
par(mar = c(1,1,1,1))
plot(karate)
# two separate colors, based on the first vertex attribute (by default)
plot(karate, vertex.size = 20)
# use a different shape for two main vertices: John A and Mr Hi
V(karate)$shape = "circle"
V(karate)[c("Mr Hi", "John A")]$shape = "square"
plot(karate, vertex.size = 20)


# represent tie weights into the graph
E(karate)$weight
plot(karate, vertex.size = 20, edge.width = E(karate)$weight)

# use a different color for the two leaders
V(karate)[c("Mr Hi", "John A")]$color = "cyan"
plot(karate, vertex.size = 20, edge.width = E(karate)$weight)

