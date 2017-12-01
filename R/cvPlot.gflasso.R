

cvPlot.gflasso <- function(cv.gflasso){
      pheatmap(cv.gflasso$mean, cluster_rows = F, cluster_cols = F,
               main = paste("CV mean RMSE\nOptimal pars:", "lambda =", cv.gflasso$optimal$lambda,
                             ",", "gamma =", cv.gflasso$optimal$gamma))
}