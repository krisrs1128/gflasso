
################################################################################
# proximal gradient descent step
################################################################################

# descent-utils ----------------------------------------------------------------
#' @title Calculate the gradient for one step
#' @param Y The matrix of regression responses.
#' @param X The data matrix.
#' @param B The J x K regression coefficient estimates.
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param mu The smoothing parameter.
#' @return The gradient descent direction for B, given the current estimates.
#' See equation (11) in the reference.
#' @references Smoothing proximal gradient method for general structured sparse regression
get_grad_f <- function(X, Y, B, C, mu) {
  A_opt <- get_A_opt(B, C, mu)
  t(X) %*% (X %*% B - Y) + A_opt %*% t(C)
}

#' @title Get the objective of the current estimates
#' @param X The data matrix.
#' @param B The J x K regression coefficient estimates.
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @return The error of the reconstruction according to the original (nonsmooth)
#' objective.
#' @export
objective <- function(X, B, R, C) {
  (1 / 2) * sum( (Y - X %*% B) ^ 2) + sum(abs(B %*% C))
}

#' @title Get the Z update in one step of the gradient descent
#' @param Lu The automatically chosen step sizes defined by equation (12)
#' @param grad_f A list of gradient descent directions for B at each iteration.
#' @param t The current iteration of the gradient descent.
update_Z <- function(Lu, grad_f, t) {
    weights <- (1 / (2 * Lu)) * seq_len(t)
    Reduce("+", lapply(seq_len(t), function(i) {
      weights[i] * grad_f[[i]]
    }))
}

# descent ----------------------------------------------------------------------
#' @title Merge default proximal gradient descent options
#' @param opts A potentially partially specified list of regularization and
#' gradient descent parameters. The currently supported options are,
#' $delta_conv How small does the change in B need to be to declare convergence?
#' $eps A tolerance used in calculating the degree of smoothing in the proximal
#' gradient objective.
#' $gamma The graph regularization parameter.
#' $iter_max What is the maximum number of iterations we should run?
#' $lambda The l1 regularization parameter.
#' @return A modified version of opts with defaults filled in.
#' @export
merge_proxgrad_opts <- function(opts) {
  default_opts <- list()
  default_opts$delta_conv <- 1e-2
  default_opts$eps <- 0.0001
  default_opts$gamma <- 0.1
  default_opts$iter_max <- 1e3
  default_opts$lambda <- 0.1
  modifyList(default_opts, opts)
}

#' @title Proximal Gradient Descent for Graph Fused Lasso
#' @param Y The matrix of regression responses.
#' @param X The data matrix.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @param opts A potentially partially specified list of regularization and
#' gradient descent parameters. See merge_proxgrad_opts().
#' @return A list containing the estimated parameters B, as well as various
#' intermediate gradient descent quantities described in Algorithm 1 of the
#' reference.
#' @references Graph-Structured Multi-task Regression adn Efficient Optimization for General Fused Lasso
#' @export
proxgrad <- function(Y, X, R, opts = list()) {
  # get opts
  J <- ncol(X)
  K <- ncol(Y)
  opts <- merge_proxgrad_opts(opts)

  # get L1 penalty matrices
  H <- get_H(R)
  C <- get_C(H, opts$lambda, opts$gamma)

  # calculate automatic step size
  D <- (1 / 2) * J * (K + ncol(H) / 2)
  mu <- eps / (2 * D)
  Lu <- get_Lu(X, R, lambda, gamma, mu)

  # initialize results
  W <- list(matrix(0, J, K))
  B <- matrix(0, J, K)
  grad_f <- list()
  Z <- list()
  obj <- list()

  t <- 1
  while(TRUE) {
    # make a step
    grad_f[[t]] <- get_grad_f(X, Y, B, C, mu)
    B_new <- W[[t]] - (1 / Lu) * grad_f[[t]]
    Z[[t]] <- update_Z(Lu, grad_f, t)
    W[[t + 1]] <- (t + 1) / (t + 3) * B_new + (2 / (t + 3)) * Z[[t]]

    # check convergence, and update counter
    delta <- sum(abs(B_new - B))
    obj[[t]]  <- objective(X, B, R, C)
    if(delta < opts$delta_conv | t > opts$iter_max) break
    B <- B_new
    t <- t + 1
  }
  list(B = B, obj = obj, Z = Z, W = W, grad_f = grad_f, Lu = Lu)
}

