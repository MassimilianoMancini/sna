library(igraph)
for (i in c(1, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23)) {
  for (j in c('N', 'HN')) {
    filename <- paste(j, '34_', i, '.DAT', sep = '')
    Y <- as.matrix(read.table(filename))
    g <- graph_from_adjacency_matrix(Y, mode = 'directed')
    plot(g, vertex.size = 10, edge.arrow.size = 0.1, vertex.label.cex = 0.3, main = filename)
    readline(prompt="Press [enter] to continue")
  }
}


par(mfrow = c(1,3), mar = c(0,0,2,0))

Y7 <- as.matrix(read.table('N34_7.DAT'))
g7 <- graph_from_adjacency_matrix(Y7, mode = 'directed')
plot(g7, vertex.size = 10, edge.arrow.size = 0.1, vertex.label.cex = 0.3, main = '7')
box()

Y10 <- as.matrix(read.table('HN34_10.DAT'))
g10 <- graph_from_adjacency_matrix(Y10, mode = 'directed')
plot(g10, vertex.size = 10, edge.arrow.size = 0.1, vertex.label.cex = 0.3, main = '10')
box()

Y19 <- as.matrix(read.table('HN34_19.DAT'))
g19 <- graph_from_adjacency_matrix(Y19, mode = 'directed')
plot(g19, vertex.size = 10, edge.arrow.size = 0.1, vertex.label.cex = 0.3, main = '19')
box()
graph.density(g7)
graph.density(g10)
graph.density(g19)
dyad.census(g7)
dyad.census(g10)
dyad.census(g19)
attrs <- read.table('CBE10.DAT')
V(g10)$gender <- attrs$V1
V(g10)$delinq <- attrs$V2
V(g10)$friend <- attrs$V3
genderColor <- ifelse(V(g10)$gender == 1, 'pink', 'lightblue')
plot(g10, vertex.size = 10, edge.arrow.size = 0.1, vertex.label.cex = 0.3, main = '10', vertex.color = genderColor)

fine = 500 # this will adjust the resolving power.
pal = colorRampPalette(c('white','red'))

#this gives you the colors you want for every point
delinqCol = pal(fine)[as.numeric(cut(V(g10)$delinq,breaks = fine))]

# now you just need to plot it with those colors
plot(g10, vertex.size = 10, edge.arrow.size = 0.1, vertex.label.cex = 0.3, main = '10', vertex.color = delinqCol)
legend('bottomright', legend = c('0-1', '1-2', '2-3'), fill = c("#FFFFFF", "#FF7777", "#FF5858"))

fine = 500 # this will adjust the resolving power.
pal = colorRampPalette(c('white','green'))

#this gives you the colors you want for every point
friendCol = pal(fine)[as.numeric(cut(V(g10)$friend,breaks = fine))]

# now you just need to plot it with those colors
plot(g10, vertex.size = 10, edge.arrow.size = 0.1, vertex.label.cex = 0.3, main = '10', vertex.color = friendCol)


assortativity(g10, V(g10)$gender)
assortativity(g10, V(g10)$delinq)
assortativity(g10, V(g10)$friend)

inDegree <- degree(g10, normalized = TRUE, mode = 'in')
outDegree <- degree(g10, normalized = TRUE, mode = 'out')

plot(g10, vertex.size = inDegree*70, main = '10', vertex.color = genderColor)
plot(g10, vertex.size = outDegree*70, edge.arrow.size = 0.3, main = '10', vertex.color = genderColor)



inCloseness <- closeness(g10, normalized = TRUE, mode = 'in')
outCloseness <- closeness(g10, normalized = TRUE, mode = 'out')

g11 <- g10[is.finite(inCloseness)]

plot(g11, vertex.size = inCloseness*100, main = '10', vertex.color = genderColor)
plot(g11, vertex.size = outCloseness, edge.arrow.size = 0.3, main = '10', vertex.color = genderColor)


inCloseness[is.finite(inCloseness)]
V(g10)$delinq
delinqCol

V(g10)$gender
genderColor

count_components(g10)
components(g10)

summary(attrs)
hist(attrs$V3)
?plot

plot(V(g10)$delinq[order(V(g10)$delinq)])
order(attrs$V2)
     