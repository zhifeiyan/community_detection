#' Spectral clustering to the adjacency matrix
#' 
#' For details, see paper "CONSISTENCY OF SPECTRAL CLUSTERING IN STOCHASTIC BLOCK MODELS" by Jing Lei and Alessandro Rinaldo
#' 
#' @export
adjsc <- function(adj, nclust, kmeans_iter = 1e3, kmeans_rep = 1e2) {
  eig <- eigen(adj)
  # Extract leading nclust eigenvectors (ordered by absolute eigenvalues)
  eigenmat <- eig$vectors[, order(abs(eig$values), 
                                  decreasing = TRUE)[1:nclust]]
  kmout <- kmeans(eigenmat, nclust, 
                  iter.max = kmeans_iter, nstart = kmeans_rep)
  list(eigenmat = eigenmat, label = kmout$cluster)
}