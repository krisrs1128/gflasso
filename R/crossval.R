#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Cross-validation functions for Graph-Fused lasso. This can't just be done
## directly by caret since that package doesn't support multitask regression.
##
## author: Francisco Lima (https://github.com/monogenea), with revisions by
## sankaran.kris@gmail.com
## date: 11/26/2017

#' Cross Validation for GFLasso
#'
#' Cross validation for GFLasso, supports crossval function.
#'
#' @param Y The matrix of regression responses.
#' @param X The data matrix.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @param opts A potentially partially specified list of regularization and
#' gradient descent parameters. See merge_proxgrad_opts().
#' @param cvIndex [list] each element of which carries the indices of the rows
#'   for a particular fold
#' @return rmse The average RMSE per element in the lambda*gamma grid
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' u <- matrix(rnorm(10), 10, 1)
#' B <- u %*% t(u) + matrix(rnorm(10 * 10, 0, 0.1), 10, 10)
#' Y <- X %*% B + matrix(rnorm(100 * 10), 100, 10)
#' cvIndex <- caret::createFolds(1:nrow(Y))
#' cv_gflasso(Y, X, cor(Y), list(), cvIndex)
#' @export
cv_gflasso <- function(Y, X, R, opts, cvIndex) {
  rmse <- rep(NA, length(cvIndex))
  for(i in seq_along(cvIndex)) {
    mod <- gflasso(Y = Y[-cvIndex[[i]],], X = X[-cvIndex[[i]],], R = R, opts = opts)
    pred <- X[cvIndex[[i]],] %*% mod$B
    rmse[i] <- sqrt(mean((pred - Y[cvIndex[[i]],])**2))
  }

  rmse
}

#' Cross Validate on a prespecified grid
#'
#' Wrapper function for cv_gflasso, applying it along a grid of lambda and gamma
#' values.
#'
#' @param Y The matrix of regression responses.
#' @param X The data matrix.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @param params The grid of lambda and gamma values to try
#' @param cvIndex [list] each element of which carries the indices of the rows
#'   for a particular fold
#' @return cvMatrix A matrix of errors across a grid of lambda (row) and gamma
#'   (column) values.
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' u <- matrix(rnorm(10), 10, 1)
#' B <- u %*% t(u) + matrix(rnorm(10 * 10, 0, 0.1), 10, 10)
#' Y <- X %*% B + matrix(rnorm(100 * 10), 100, 10)
#' cvIndex <- caret::createFolds(1:nrow(Y))
#' crossval(X, Y, cor(Y), seq(0, 1, by = 0.1), cvIndex)
#' @export
crossval <- function(X, Y, R, additionalOpts = NULL, k = 5, times = 1,
                     params = seq(0,1,by=0.1)) {
      cvIndex <- caret::createMultiFolds(1:nrow(Y), k = k, times = times)    
      cvArray <- array(NA, dim = c(rep(length(params), 2), k * times))
      dimnames(cvArray) <- list(params, params, names(cvIndex))
      grid <- expand.grid(lambda = params, gamma = params)
      for(i in 1:nrow(grid)){
            cv <- cv_gflasso(X = X, Y = Y, R = R, opts = list(lambda = grid[i,1], gamma = grid[i,2], additionalOpts),
                             cvIndex = cvIndex)
            cvArray[as.character(grid[i,1]),as.character(grid[i,2]),] <- cv
            message(paste(round((i/(nrow(grid)))*100, 2), "% completion", collapse = " "))
            }
      cvMean <- apply(cvArray, 1:2, mean)
      cvSE <- apply(cvArray, 1:2, function(x) sd(x) / sqrt(k * times))
      optimal <- grid[which.min(cvMean),]
      return(list("mean" = cvMean, "SE" = cvSE, "optimal" = optimal))
}
