
################################################################################
# Generate input matrices for the proximal gradient descent
################################################################################

#' @title Get H matrix giving Graph Penalization
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @return The matrix H defined in the reference.
#' @references Smoothing Proximal Gradient Method for General Structured Sparse Regressoin
#' @export
get_H <- function(R) {
  K <- nrow(R)
  R[lower.tri(R, diag = TRUE)] <- 0
  R_obs <- which(R != 0, arr.ind = T)
  H <- matrix(0, K, max(nrow(R_obs), 1)) # protect against no-edges case
  dimnames(H) <- list(1:K, apply(R_obs, 1, paste0, collapse = "-"))

  for(i in seq_len(nrow(R_obs))) {
    cur_row <- R_obs[i, 1]
    cur_col <- R_obs[i, 2]
    rij <- R[cur_row, cur_col]
    tau_rij <- abs(rij)
    cur_h_col <- paste(R_obs[i, ], collapse = "-")
    H[cur_row, cur_h_col] <- tau_rij
    H[cur_col, cur_h_col] <- -sign(rij) * tau_rij
  }
  H
}

#' @title Get automatic step-sizes
#' @param X The data matrix.
#' @param H  The matrix H defined in the reference.
#' @param lambda The l1 regularization parameter.
#' @param gamma The graph regularization parameter.
#' @param mu The smoothing parameter.
#' @return Lu The automatically chosen step sizes defined in the reference.
#' @references Smoothing Proximal Gradient Method for General Structured Sparse Regressoin
#' @export
get_Lu <- function(X, C, lambda, gamma, mu) {
  eigen(t(X) %*% X)$values[1] +
    (1 / mu) * (lambda ^ 2 + 2 * gamma ^ 2 * max(rowSums(C)))
}
