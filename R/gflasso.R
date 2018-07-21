
################################################################################
# Graph fused lasso
################################################################################

#' @title Merge default proximal gradient descent options
#' @param opts A potentially partially specified list of regularization and
#' gradient descent parameters. The currently supported options are,
#' $delta_conv How small does the change in B need to be to declare convergence?
#' $eps A tolerance used in calculating the degree of smoothing in the proximal
#' gradient objective.
#' $gamma The graph regularization parameter.
#' $iter_max What is the maximum number of iterations we should run?
#' $lambda The l1 regularization parameter.
#' $verbose Should the gradient descent print its progress?
#' @return A modified version of opts with defaults filled in.
#' @export
merge_proxgrad_opts <- function(opts, J, K) {
  default_opts <- list()
  default_opts$delta_conv <- 1e-2
  default_opts$eps <- 0.005
  default_opts$gamma <- 1
  default_opts$iter_max <- 1e3
  default_opts$lambda <- 1
  default_opts$verbose <- FALSE
  default_opts$B0 <- matrix(0, J, K)
  modifyList(default_opts, opts)
}

#' @title Graph Fused Lasso via Smoothed Proximal Gradient Descent
#' @param Y The matrix of regression responses, scaled and centered as necessary.
#' @param X The data matrix, scaled and centered as necessary.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @param opts A potentially partially specified list of regularization and
#' gradient descent parameters. See merge_proxgrad_opts().
#' @return A list containing the following quantities: \cr
#'   $B The estimated beta coefficient matrix. Its rows correspond columns of
#'    X, while its columns correspond to columns of Y. \cr
#'   $obj The graph fused lasso objective function, over iterations.
#'   $Lu The automatically calculated step size. \cr
#' reference.
#' @references Smoothing Proximal Gradient Method for General Structured Sparse Regression
#' @export
gflasso <- function(Y, X, R, opts = list()) {
  # get opts
  opts <- merge_proxgrad_opts(opts, ncol(X), ncol(Y))

  # get L1 penalty matrix
  C <- opts$gamma * t(get_H(R))

    # calculate automatic step size
  D <- (1 / 2) * ncol(X) * (ncol(Y) + ncol(C) / 2)
  mu <- opts$eps / (2 * D)
  Lu <- get_Lu(X, C, opts$lambda, opts$gamma, mu)

  accgrad_opts <- list(lambda = opts$lambda, L = Lu, mu = mu,
                       iter_max = opts$iter_max, delta_conv = opts$delta_conv,
                       verbose = opts$verbose, B0 = opts$B0)
  optim_result <- accgrad(X, Y, C, accgrad_opts)
  list(B = optim_result$B, obj = optim_result$obj, Lu = Lu)
}


#' @title GFLASSO prediction
#'
#' @param model The model object from `gflasso`.
#' @param new.data The data from which to predict.
#'
#' @return A n x k matrix carrying all k predicted responses across all n samples
#' @export
predict_gflasso <- function(model, new.data){
      # Simple matrix multiplication
      return(new.data %*% model$B)
}
