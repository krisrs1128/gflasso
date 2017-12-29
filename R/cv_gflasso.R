#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Cross-validation functions for Graph-Fused lasso. This can't just be done
## directly by caret since that package doesn't support multitask regression.
##
## author: Francisco Lima (https://github.com/monogenea), with revisions by
## sankaran.kris@gmail.com
## date: 12/26/2017

#' Compute the root mean squared error
rmse <- function(pred, y) {
  sqrt(mean( (pred - y) ^ 2 ))
}

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
#' @param err_fun A function that computes the error between vectors of
#'   predicted and true responses. Defaults to rmse(pred, y) = sqrt(mean( (pred - y) ^ 2)).
#' @return cvMatrix A matrix of errors across a grid of lambda (row) and gamma
#'   (column) values.
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @importFrom caret createMultiFolds
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' u <- matrix(rnorm(10), 10, 1)
#' B <- u %*% t(u) + matrix(rnorm(10 * 10, 0, 0.1), 10, 10)
#' Y <- X %*% B + matrix(rnorm(100 * 10), 100, 10)
#' R <- ifelse(cor(Y) > .8, 1, 0)
#' system.time(testCV <- cv_gflasso(X, Y, R, nCores = 1))
#' system.time(testCV <- cv_gflasso(X, Y, R, nCores = 2))
#' cv_plot_gflasso(testCV)
#' @export
cv_gflasso <- function(X, Y, R, additionalOpts = list(), k = 5, times = 1,
                       params = seq(0, 1, by = 0.1), nCores = NULL,
                       err_fun = rmse) {

  additionalOpts <- merge_proxgrad_opts(additionalOpts, ncol(X), ncol(Y))
  if (is.null(nCores)) {
    nCores <- detectCores() - 1
  }

  cvFUN <- function(Y, X, R, opts, cvIndex) {
    sapply(seq_along(cvIndex), function(i){
      mod <- gflasso(Y = Y[-cvIndex[[i]],], X = X[-cvIndex[[i]],], R = R, opts = opts)
      pred <- X[cvIndex[[i]],] %*% mod$B
      err_fun(pred, Y[cvIndex[[i]], ])
    })
  }

  cvIndex <- caret::createMultiFolds(1:nrow(Y), k = k, times = times)
  cvArray <- array(NA, dim = c(rep(length(params), 2), k * times))
  dimnames(cvArray) <- list(params, params, names(cvIndex))

  grid <- expand.grid(lambda = params, gamma = params)
  allCV <- mclapply(
    as.list(1:nrow(grid)),
    function(x) {
      if (additionalOpts$verbose && x %% 10 == 0) {
        cat(sprintf("CV grid %s/%s \n", x, nrow(grid)))
      }

      cvFUN(
        X = X, Y = Y, R = R,
        opts = list(lambda = grid[x, 1], gamma = grid[x, 2], additionalOpts),
        cvIndex = cvIndex
      )
    },
    mc.cores = nCores
  )

  print(allCV[[1]])
  for(i in 1:nrow(grid)){
    cvArray[as.character(grid[i,1]),as.character(grid[i,2]),] <- allCV[[i]]
  }

  cvMean <- apply(cvArray, 1:2, mean)
  list(
    "mean" = cvMean,
    "SE" = apply(cvArray, 1:2, sd) / sqrt(k * times),
    "optimal" = grid[which.min(cvMean), ]
  )
}

#' Plot Results from Cross Validation
#'
#' @importFrom pheatmap pheatmap
#' @export
cv_plot_gflasso <- function(cv.gflasso){
  pheatmap(cv.gflasso$mean, cluster_rows = F, cluster_cols = F,
           main = paste("CV mean RMSE\nOptimal pars:", "lambda =", cv.gflasso$optimal$lambda,
                        ",", "gamma =", cv.gflasso$optimal$gamma))
}
