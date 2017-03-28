#' Regularized spectral clustering
#'
#' Regularized spectral clustering where entries of the
#' adjacency matrix are inflated by a constant small number. 
#' See paper "PSEUDO-LIKELIHOOD METHODS FOR COMMUNITY DETECTION IN LARGE SPARSE NETWORKS" by ARASH A. AMINI, AIYOU CHEN, PETER J. BICKEL AND ELIZAVETA LEVINA for details"
#'
#' @export
rsc <- function(adj, n, nclust, const = 0.25,
                kmeans_iter = 1e3, kmeans_rep = 1e2){
  # Adjust adjacency matrix
  adj <- adj + const * sum(adj) / n^2
  d <- rowSums(adj)^(-0.5)
  lap <- sweep(adj, 1, d, '*')
  lap <- sweep(lap, 2, d, '*')
  eig <- eigen(lap)
  eigenmat <- eig$vectors[, order(abs(eig$values), 
                                  decreasing = TRUE)[1:nclust]]
  kmout <- kmeans(eigenmat, nclust, 
                  iter.max = kmeans_iter, nstart = kmeans_rep)

  list(eigenmat = eigenmat, label = kmout$cluster)
}
