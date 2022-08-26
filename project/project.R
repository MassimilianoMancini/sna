# Project by Massimiliano Mancini

rm(list = ls())

library(igraph)

setwd("~/unifi/sna/project")

Y <- as.matrix(read.table('HN34_10.DAT'))
attrs <- read.table('CBE10.DAT')

g <- graph_from_adjacency_matrix(Y, mode = 'directed')

V(g)$gender <- attrs$V1
V(g)$delinq <- attrs$V2
V(g)$friend <- attrs$V3

rho <- graph.density(g)
oddsRho <- rho/(1-rho)

cat('Density: ', rho)
cat('Density odds: ', oddsRho)

dyad.census(g)

rec <- reciprocity(g)
oddsRec <- rec/(1-rec)

cat('Reciprocity: ', rec)
cat ('Reciprocity odds: ', oddsRec)


censusLabels = c('empty',
                  'A->B, C',
                  'A<->B, C',
                  'A<-B->C',
                  'A->B<-C',
                  'A->B->C',
                  'A<->B<-C',
                  'A<->B->C',
                  'A->B<-C, A->C',
                  'A<-B<-C, A->C',
                  'A<->B<->C',
                  'A<-B->C, A<->C',
                  'A->B<-C, A<->C',
                  'A->B->C, A<->C',
                  'A->B<->C, A<->C',
                  'A<->B<->C, A<->C')

data.frame(censusLabels, triad.census(g))

tran <- transitivity(g, type = 'global')
oddsTran <- tran/(1-tran)

standardizeTrans <- log(oddsTran) - log(oddsRho) 

cat('Transitivity: ', tran)
cat('Transitivity odds: ', oddsTran)
cat('Standardized transitivity: ', standardizeTrans)

fine = 500

genderColor <- ifelse(V(g)$gender == 1, 'pink', 'lightblue')

pal = colorRampPalette(c('white','red'))
delinqColor = pal(fine)[as.numeric(cut(V(g)$delinq,breaks = fine))]

pal = colorRampPalette(c('blue','white'))
friendColor = pal(fine)[as.numeric(cut(V(g)$friend,breaks = fine))]

set.seed(6)
png(filename = 'gender.png', width = 1024, height = 1024)

plot(g, 
     vertex.size = 10, 
     vertex.color = genderColor, 
     main = 'Gender')
dev.off()

png(filename = 'dt.png', width = 1024, height = 1024)

set.seed(6)
plot(g, 
     vertex.size = 10, 
     vertex.color = delinqColor, 
     main = 'Delinquent tendency')
dev.off()

set.seed(6)
png(filename = 'fi.png', width = 1024, height = 1024)
plot(g, 
     vertex.size = 10, 
     vertex.color = friendColor, 
     main = 'Friendship importance')
dev.off()


assortativity(g, V(g)$gender)
assortativity(g, V(g)$delinq)
assortativity(g, V(g)$friend)

inDegree <- degree(g, normalized = TRUE, mode = 'in')
outDegree <- degree(g, normalized = TRUE, mode = 'out')

set.seed(6)
png(filename = 'indegree.png', width = 1024, height = 1024)
plot(g, 
     vertex.size = inDegree*70, 
     vertex.color = delinqColor, 
     main = 'In Degree')
dev.off()


set.seed(6)
png(filename = 'outdegree.png', width = 1024, height = 1024)
plot(g, 
     vertex.size = outDegree*70, 
     vertex.color = delinqColor, 
     main = 'Out Degree')
dev.off()

print('Best 3 InDegree: ')
sort(inDegree, decreasing = TRUE, index.return = TRUE)$x[1:3]

print('Best 3 OutDegree: ')
sort(outDegree, decreasing = TRUE, index.return = TRUE)$x[1:3]

comps <- components(g, mode = 'strong')
gc <- subgraph(g, V(g)[comps$membership == which.max(comps$csize)])

inCloseness <- closeness(gc, normalized = TRUE, mode = 'in')
outCloseness <- closeness(gc, normalized = TRUE, mode = 'out')

set.seed(6)
png(filename = 'incloseness.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = inCloseness*40, 
     vertex.color = delinqColor, 
     main = 'Closeness In')
dev.off()

set.seed(6)
png(filename = 'outcloseness.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = outCloseness*40, 
     vertex.color = delinqColor, 
     main = 'Closeness Out')
dev.off()

print('Best 3 InCloseness: ')
sort(inCloseness, decreasing = TRUE, index.return = TRUE)$x[1:3]

print('Best 3 OutCloseness: ')
sort(outCloseness, decreasing = TRUE, index.return = TRUE)$x[1:3]


betw <- betweenness(gc, normalized = TRUE)


set.seed(6)
png(filename = 'betweenness.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = betw*70, 
     vertex.color = delinqColor, 
     main = 'Betweenness')
dev.off()

print('Best 3 Betweenness: ')
sort(betw, decreasing = TRUE, index.return = TRUE)$x[1:3]

geigen <- eigen_centrality(gc, scale = TRUE)$vector

set.seed(6)
png(filename = 'eigen.png', width = 1024, height = 1024)
plot(gc, 
     vertex.size = geigen*20, 
     vertex.color = delinqColor, 
     main = 'Eigenvector')
dev.off()

print('Best 3 Eigenvector: ')
sort(geigen, decreasing = TRUE, index.return = TRUE)$x[1:3]

centr_degree(gc, loops = FALSE, mode = 'in')$centralization
centr_degree(gc, loops = FALSE, mode = 'out')$centralization

centr_clo(gc, mode = 'in')$centralization
centr_clo(gc, mode = 'out')$centralization

centr_betw(gc)$centralization

centr_eigen(gc)$centralization

