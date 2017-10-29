# Return a sequence of estimated probability matrices from Zhang's neighborhood smoothing

#' @export
nbs_realnet <- function(adj, h_seq) {
  probmat_est <- vector('list', length = length(h_seq))

  n <- nrow(adj)
  d <- compute_dist_nbs(adj)

  # For the ease of latter comparison, set diag to be infinity
  diag(d) <- Inf

  for (j in seq_along(h_seq)) {
    h <- h_seq[j]

    membermat <- matrix(0, n, n)
    for (i in 1:n) {
      q_h <- quantile(d[i, -i], h)
      nb_id <- which(d[i, ] <= q_h)
      membermat[i, nb_id] <- 1 / length(nb_id)
    }
    probmat_est[[j]] <- 0.5 * (membermat %*% adj + adj %*% t(membermat))
  }
  
  probmat_est
}