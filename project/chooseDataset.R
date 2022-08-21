for (i in c(1, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23)) {
  for (j in c('N', 'HN')) {
    filename <- paste(j, '34_', i, '.DAT', sep = '')
    Y <- as.matrix(read.table(filename))
    g <- graph_from_adjacency_matrix(Y, mode = 'directed')
    plot(g, vertex.size = 10, edge.arrow.size = 0.1, vertex.label.cex = 0.3, main = filename)
    readline(prompt="Press [enter] to continue")
  }
}



