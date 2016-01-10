
################################################################################
# Generate input matrices for the proximal gradient descent
################################################################################

#' @title Get H matrix giving Graph Penalization
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @return The matrix H defined on page 6 of the reference.
#' @references Graph-Structured Multi-task Regression adn Efficient Optimization for General Fused Lasso
#' @export
get_H <- function(R) {
  K <- nrow(R)
  R[lower.tri(R)] <- 0
  R_obs <- which(R != 0, arr.ind = T)
  H <- matrix(0, K, nrow(R_obs))
  dimnames(H) <- list(1:K, apply(R_obs, 1, paste0, collapse = "-"))

  for(i in seq_len(nrow(R_obs))) {
    cur_row <- R_obs[i, 1]
    cur_col <- R_obs[i, 2]
    tau_rij <- abs(R[cur_row, cur_col])
    cur_h_col <- paste(R_obs[i, ], collapse = "-")
    H[cur_row, cur_h_col] <- tau_rij
    H[cur_col, cur_h_col] <- -sign(rij) * tau_rij
  }
  H
}

#' @title Get C matrix giving L1 penalization
#' @param H The H matrix output by \code{get_H()}.
#' @param lambda The l1 regularization parameter.
#' @param gamma The graph regularization parameter.
#' @return THe C matrix defined on page 6 of the reference.
#' @references Graph-Structured Multi-task Regression adn Efficient Optimization for General Fused Lasso
#' @export
get_C <- function(H, lambda, gamma) {
  cbind(lambda * diag(nrow(H)), gamma * H)
}

#' @title Get optimal A matrix from the dual norm representation
#' @param B The J x K regression coefficient estimates.
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param mu The smoothing parameter.
#' @return The optimal A matrix defined by Lemma 1.
#' @references Graph-Structured Multi-task Regression adn Efficient Optimization for General Fused Lasso
#' @export
get_A_opt <- function(B, C, mu) {
  A <- (1 / mu) * B %*% C
  A[A > 1] <- 1
  A[A < -1] <- -1
  A
}

#' @title Get automatic step-sizes
#' @param X The data matrix.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @param lambda The l1 regularization parameter.
#' @param gamma The graph regularization parameter.
#' @param mu The smoothing parameter.
#' @return Lu The automatically chosen step sizes defined by equation (12)
#' @references Graph-Structured Multi-task Regression adn Efficient Optimization for General Fused Lasso
#' @export
get_Lu <- function(X, R, lambda, gamma, mu) {
  tau_R <- abs(R)
  eigen(t(X) %*% X)$values[1] +
    (1 / mu) * (lambda ^ 2 + 2 * gamma ^ 2 * max(colSums(tau_R ^ 2)))
}
