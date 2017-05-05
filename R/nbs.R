# Created by Zhifei Yan
# Last update: 2017-4-23

# MSE of neighborhood smoothing method proposed by Yuan Zhang
# This function accepts a sequence of tuning parameter h - hth sample quantile

#' @export
nbs <- function(A, expect_A, h_seq) {
  n <- nrow(A)
  d <- compute_dist_nbs(A)
  # For the ease of latter comparison, set diag to be infinity
  diag(d) <- Inf
  mse_seq <- rep(NA, length(h_seq))

  for (k in seq_along(h_seq)) {
    h <- h_seq[k]
    membermat <- matrix(0, n, n)
    for (i in 1:n) {
      q_h <- quantile(d[i, -i], h)
      nb_id <- which(d[i, ] <= q_h)
      membermat[i, nb_id] <- 1 / length(nb_id)
    }
    est_prob <- 0.5 * (membermat %*% A + A %*% t(membermat))
    mse_seq[k] <- norm(est_prob - expect_A, 'f')^2 / n^2
  }
  
  mse_seq
}