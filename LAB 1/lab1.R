# Load igraph library
library('igraph')

# Create a simple undirected graph
g <-
  graph_from_literal(1 - 2, 1 - 3, 2 - 3, 2 - 4, 3 - 5, 4 - 5, 4 - 6, 4 -
                       7, 5 - 6, 6 - 7)

# Plot the graph
plot(g)

# Display vertices
V(g)

# Display edges
E(g)

# Create simple directed graph
gd <- graph_from_literal(1-+2, 1-+3, 2++3, 4++1, 4+-3, 5-+1)

# Plot the directed graph with reasonable arrow dimension
plot(gd, edge.arrow.size = 0.1)

# another directed graph example
gd_lab <- graph_from_literal(A-+B, A-+C, B++C, D++A, D+-C, E-+A)
plot(gd_lab, edge.arrow.size = 0.1)

# yet another example using labels
gd_lab <- graph_from_literal(1-+2, 1-+3, 2++3, 4++1, 4+-3, 5-+1)
plot(gd_lab, edge.arrow.size = 0.1)

# Set vertices labels
V(gd_lab)$name <- c("A", "B", "C", "D", "E")
plot(gd_lab, edge.arrow.size = 0.12)

# a graph with isolated vertices
gd_iso <- graph_from_literal(1-+2, 1-+3, 2++3, 4++1, 4+-3, 5-+1, 6)
plot(gd_iso, edge.arrow.size = 0.12)

# a graph with 2 different components
gd_2comp <- graph_from_literal(1-+2:4, 2++3, 3-+4, 5-+1, 6++7:8)
plot(gd_2comp, edge.arrow.size = 0.12)

# Enumerate the components
components(gd_2comp)

# Count number of vertices
vcount(gd_2comp)

# count number of edges
ecount(gd_2comp)

# Check if graph is directed
is.directed(gd_2comp)

# Set a 2x2 layout for aligned plots
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))

# Plot graph as is, i.e. directed
plot(gd_2comp, edge.arrow.size = 0.12, main = "directed")
box()

# Plot directed graph as undirected by collapse
plot(as.undirected(gd_2comp, mode = "collapse"), main = "collapse")
box()

# Plot directed graph as undirected by mutual
plot(as.undirected(gd_2comp, mode = "mutual"), main = "mutual")
box()

# Plot directed graph as undirected by each
plot(as.undirected(gd_2comp, mode = "each"), main = "each")
box()

# Set a 1x3 layout for aligned plots
par(mfrow = c(1, 3), mar = c(1, 1, 1, 1))

# Plot undirected graph as is
plot(g, main = 'undirected')
box()

# Plot undirected graph as directed by mutual
plot(as.directed(g, mode = 'mutual'),
     edge.arrow.size = 0.13,
     main = 'mutual')
box()

# Plot undirected graph as directed by arbitrary
plot(as.directed(g, mode = 'arbitrary'),
     edge.arrow.size = 0.13,
     main = 'arbitrary')
box()

# Matrix representation
Y <- get.adjacency(g)

# Set diagonal as NA
diag(Y) <- NA

# Create a matrix from scratch
Y <- matrix(0, 7, 7)
tmp = c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1)

# For undirected matrix only a lower (or upper) triangle is needed
Y[lower.tri(Y)] <- tmp

# Then replicate on the upper triangle
Y <- Y + t(Y)
diag(Y) <- NA

# Is the matrix simmetric?
sum(Y != t(Y), na.rm = TRUE) == 0

# Now let's buil a directed graph by a asimmetric matrix
tmp <-
  c(0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0)
Yd <- matrix(tmp, 5, 5)
diag(Yd) = NA

# Is the matrix simmetric?
sum(Yd != t(Yd), na.rm = TRUE) == 0

# From matrix to graph
g2 <- graph_from_adjacency_matrix(Y, mode = "undirected")
gd2 <- graph_from_adjacency_matrix(Yd, mode = "directed")

# Plot the resulting graphs
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))

plot(g, main = 'undirected 1')
box()

plot(g2, main = 'undirected 2')
box()

plot(gd, edge.arrow.size = 0.13, main = 'directed 1')
box()

plot(gd2, edge.arrow.size = 0.13, main = 'directed 2')
box()

# Create a weighted graph from the adjacency matrix
tmp <-
  c(0, 0, 3, 1, 0, 1, 0, 2, 0, 0, 1, 3, 0, 0, 0, 4, 0, 2, 0, 0, 5, 0, 0, 0, 1)
Yd <- matrix(tmp, 5, 5)
diag(Yd) <- NA
wg <-
  graph_from_adjacency_matrix(Yd,
                              mode = "directed",
                              weighted = T,
                              diag = F)

# let us give a label to the vertices of wg
V(wg)$name <- c("Pippo", "Pluto", "Paperino", "Cip", "Ciop")

# let us give a name to the graph
wg$name <- "Disney"

# plot
par(mfrow = c(1, 1))
plot(wg, edge.arrow.size = 0.13)

# Let's play with network drawing
gd <- graph_from_literal(1-+2:3, 2++3)
V(gd)$name <- c("Francesca", "Paola", "Giovanni")
plot(gd, edge.arrow.size = 0.3)

# increase the size
plot(gd, vertex.size = 35, edge.arrow.size = 0.5)

# change the color 
plot(gd, vertex.size = 35, vertex.color = "cyan", edge.arrow.size = 0.5)

# change the shape 
plot(gd, vertex.size = 35, vertex.color = "cyan", vertex.shape = "square", edge.arrow.size = 0.5)

# change the label color
plot(gd, vertex.size = 35, vertex.color = "cyan", vertex.shape = "square",
     vertex.label.color = "red", edge.arrow.size = 0.5)

# increase the font
plot(gd, vertex.label.cex = 1.3, 
     vertex.size = 50, vertex.color = "cyan", vertex.shape = "square",
     vertex.label.color = "red", edge.arrow.size = 0.5)

# change the color of edge
plot(gd, vertex.size = 35, edge.color = "blue", edge.arrow.size = 0.5)

# change the edge width
plot(gd, vertex.size = 35, edge.color = "blue", edge.arrow.size = 0.5, 
     edge.width = 2)

# Let us consider a simple graph
gd = graph_from_literal(1-+2:3, 2++3)
plot(gd, edge.arrow.size = 0.5)

# give a color to the nodes
V(gd)$color <- c("red", "green", "cyan")
plot(gd, edge.arrow.size = 0.5)

# give a color to the nose based on some other attributes 
V(gd)$name <- c("Francesca", "Paola", "Giovanni")
V(gd)$gender <- c("F", "F", "M")

# based on the gender
V(gd)[gender == "F"]$color <- "cyan"
V(gd)[gender == "M"]$color <- "green"
V(gd)$color

# graphical representation
plot(gd, vertex.size = 45, edge.arrow.size = 0.5)

# give a width to edges
E(gd)$width <- c(2,3,1,1)

# give a different width to ties based on the width  
plot(gd, vertex.size = 45, edge.arrow.size = 1)

# give a width based on another attribute
E(gd)$other <- 1:4
plot(gd, vertex.size = 45, edge.width = E(gd)$other)

# weighted adjacency matrix from gd
get.adjacency(gd, attr = "other")
# or 
get.adjacency(gd, attr = "other", sparse = FALSE)

# we may also define attributes for the whole graph
gd$name <- "My first directed graph"
plot(gd, vertex.size = 45, main = gd$name, edge.width = E(gd)$other)


# Real-data examples
library(igraphdata)
data(package = 'igraphdata')

# Zachary's karate club network
data("karate")

# structure of the dataset
karate

# graph attributes
graph.attributes(karate)

# or only the names
names(graph.attributes(karate))

# vertex attributes
vertex.attributes(karate)
names(vertex.attributes(karate))

# extract members' names
V(karate)$name
# extract members' faction
V(karate)$Faction

# edge attributes
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
