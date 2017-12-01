

cvPlot.gflasso <- function(cv.gflasso){
      if(require("pheatmap")){
      }else{
            install.packages("pheatmap")
            library(pheatmap)
      }
      pheatmap(cv.gflasso$mean, cluster_rows = F, cluster_cols = F,
               main = paste("CV mean RMSE\nOptimal pars:", "lambda =", cv.gflasso$optimal$lambda,
                             ",", "gamma =", cv.gflasso$optimal$gamma))
}