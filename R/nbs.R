# Created by Zhifei Yan
# Last update: 2017-5-12

# MSE of neighborhood smoothing method proposed by Yuan Zhang
# This function accepts a sequence of tuning parameter h - hth sample quantile

#' @export
nbs <- function(adj, probmat_true, h_seq) {
  n <- nrow(adj)
  d <- compute_dist_nbs(adj)
  # For the ease of latter comparison, set diag to be infinity
  diag(d) <- Inf

  h_best <- NA
  mse_best <- Inf
  probmat_est_best <- NA

  for (h in h_seq) {
    membermat <- matrix(0, n, n)
    for (i in 1:n) {
      q_h <- quantile(d[i, -i], h)
      nb_id <- which(d[i, ] <= q_h)
      membermat[i, nb_id] <- 1 / length(nb_id)
    }
    probmat_est <- 0.5 * (membermat %*% adj + adj %*% t(membermat))
    mse_cur <- norm(probmat_est - probmat_true, 'f')^2 / n^2

    if (mse_cur < mse_best) {
      mse_best <- mse_cur
      h_best <- h
      probmat_est_best <- probmat_est
    }

  }
  
  return(list(mse = mse_best, h = h_best, probmat_est = probmat_est_best))

}