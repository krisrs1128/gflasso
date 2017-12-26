#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Cross-validation functions for Graph-Fused lasso. This can't just be done
## directly by caret since that package doesn't support multitask regression.
##
## author: Francisco Lima (https://github.com/monogenea), with revisions by
## sankaran.kris@gmail.com
## date: 12/26/2017

#' Cross Validation for GFLasso
#'
#' @param Y The matrix of regression responses.
#' @param X The data matrix.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @param additionalOpts Additional options to pass alongside lambda and gamma (cf. gflasso)
#' @param k Number of folds
#' @param times Number of repetitions (Note: Total number of RMSE estimates = k x times)
#' @param params The grid of lambda and gamma values to try
#' @param nCores The number of CPU cores to be used, >1 represents parallelized executions
#' @return cvMatrix A matrix of errors across a grid of lambda (row) and gamma
#'   (column) values.
#' @importFrom parallel mclapply
#' @importFrom caret createMultiFolds
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' u <- matrix(rnorm(10), 10, 1)
#' B <- u %*% t(u) + matrix(rnorm(10 * 10, 0, 0.1), 10, 10)
#' Y <- X %*% B + matrix(rnorm(100 * 10), 100, 10)
#' R <- ifelse(cor(Y) > .8, 1, 0)
#' system.time(testCV <- cv.gflasso(X, Y, R, nCores = 1))
#' system.time(testCV <- cv.gflasso(X, Y, R, nCores = 2))
#' cvPlot.gflasso(testCV)
#' @export
cv.gflasso <- function(X, Y, R, additionalOpts = NULL, k = 5, times = 1,
                     params = seq(0,1,by=0.1), nCores = NULL) {

  if (is.null(nCores)) {
    nCores <- detectCores() - 1
  }

      cvFUN <- function(Y, X, R, opts, cvIndex) {
            rmse <- lapply(as.list(seq_along(cvIndex)), function(i){
                  mod <- gflasso(Y = Y[-cvIndex[[i]],], X = X[-cvIndex[[i]],], R = R, opts = opts)
                  pred <- X[cvIndex[[i]],] %*% mod$B
                  return(sqrt(mean((pred - Y[cvIndex[[i]],])**2)))
            })
            return(unlist(rmse))
      }
      cvIndex <- caret::createMultiFolds(1:nrow(Y), k = k, times = times)
      cvArray <- array(NA, dim = c(rep(length(params), 2), k * times))
      dimnames(cvArray) <- list(params, params, names(cvIndex))
      grid <- expand.grid(lambda = params, gamma = params)
      allCV <- mclapply(as.list(1:nrow(grid)), function(x){
            cv <- cvFUN(X = X, Y = Y, R = R, opts = list(lambda = grid[x,1], gamma = grid[x,2], additionalOpts),
                             cvIndex = cvIndex)
      }, mc.cores = nCores)
      for(i in 1:nrow(grid)){
            cvArray[as.character(grid[i,1]),as.character(grid[i,2]),] <- allCV[[i]]
      }
      cvMean <- apply(cvArray, 1:2, mean)
      cvSE <- apply(cvArray, 1:2, function(x) sd(x) / sqrt(k * times))
      optimal <- grid[which.min(cvMean),]

      return(list("mean" = cvMean, "SE" = cvSE, "optimal" = optimal))
}
