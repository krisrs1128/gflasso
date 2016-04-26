
################################################################################
# Nesterov's accelerated gradient descent, for graph fused lasso
################################################################################

# utils ------------------------------------------------------------------------

#' @title Standard soft-thresholding operator
#' @param v The vector to soft threshold
#' @param lambda The amount to soft-threshold by
#' @export
soft_threshold <- function(v, lambda) {
  res <- v
  res[v > lambda] <- v[v > lambda] - lambda
  res[v < -lambda] <- v[v < -lambda] + lambda
  res[v > -lambda & v < lambda] <- 0
  res
}

#' @title Get the objective of the current estimates
#' @param X The data matrix.
#' @param B The J x K regression coefficient estimates.
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @return The error of the reconstruction according to the original (nonsmooth)
#' objective.
#' @export
objective <- function(X, B, Y, C, lambda) {
  (1 / 2) * sum( (Y - X %*% B) ^ 2) + lambda * sum(abs(B)) + sum(abs(B %*% t(C)))
}

# gradient-descent-funs --------------------------------------------------------

#' @title Get optimal alpha matrix from the dual norm representation
#' @param B The J x K regression coefficient estimates.
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param mu The smoothing parameter.
#' @return The optimal A matrix defined by Lemma 1.
#' @references Smoothing Proximal Gradient Method for General Structured Sparse Regressoin
#' @export
get_alpha_opt <- function(C, W, mu) {
  alpha <- (1 / mu) * W %*% t(C)
  alpha[alpha > 1] <- 1
  alpha[alpha < -1] <- -1
  alpha
}

#' @title Calculate the gradient for one step
#' @param Y The matrix of regression responses.
#' @param X The data matrix.
#' @param B The J x K regression coefficient estimates.
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param mu The smoothing parameter.
#' @return The gradient descent direction for B, given the current estimates.
#' See equation (11) in the reference.
#' @references Smoothing Proximal Gradient Method for General Structured Sparse Regressoin
#' @export
get_grad_f <- function(X, Y, W, C, mu) {
  alpha <- get_alpha_opt(C, W, mu)
  t(X) %*% (X %*% W - Y) + alpha %*% C
}

#' @title Get next value of B in gradient descent
#' @param W The W matrix in Algorithm 1 of the reference
#' @param grad_f The gradient computed in Algorithm 1
#' @param L The lipshitz constant, which determines the step size
#' @param lambda The l1 regularization parameter.
#' @export
get_B_next <- function(W, grad_f, L, lambda) {
  B_next <- W - (1 / L) * grad_f
  soft_threshold(B_next, lambda / L)
}

#' @title Nesterov's Accelerated Gradient Descent
#' @param X The data matrix.
#' @param Y The matrix of regression responses.
#' @param H  The matrix H defined in the reference.
#' @param opts A list of gradient descent tuning parameters, \cr
#'   $iter_max: The maximum number of iterations.
#'   $mu: The smoothing parameter
#'   $lambda: The l1 regularization parameter
#'   $delta_conv: The convergence criterion for delta
#' @return A list containing the following elements \cr
#'   $B The final estimates of the multitask regression coefficients
#'   $obj The objective function across gradient descent iterations
#' @export
accgrad <- function(X, Y, C, opts) {
  # initialize results
  W <- matrix(0, ncol(X), ncol(Y))
  B <- matrix(0, ncol(X), ncol(Y))
  obj <- vector(length = opts$iter_max)

  theta <- 1
  if (opts$verbose) cat("\titer\t|\tobj\t|\t|B(t + 1) - B(t)| \n")
  for (iter in seq_len(opts$iter_max)) {
    # make a step
    grad_f <- get_grad_f(X, Y, W, C, opts$mu)
    B_next <- get_B_next(W, grad_f, opts$L, opts$lambda)

    theta_next <- 2 / (iter + 1)
    W <- B_next + ((1 - theta) / theta) * (theta_next) * (B_next - B)

    # check convergence, and update counter
    delta <- sum(abs(B_next - B))
    B <- B_next
    obj[iter]  <- objective(X, B, Y, C, opts$lambda)
    if (iter %% 10 == 0 & opts$verbose) {
      cat(sprintf("%d \t | %f \t | %f \n", iter, obj[iter], delta))
    }
    if (delta < opts$delta_conv | iter > opts$iter_max) break
    theta <- theta_next
  }
  list(B = B, obj = obj[seq_len(iter)])
}

