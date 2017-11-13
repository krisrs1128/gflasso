

# CV GFLASSO (RMSE-based, can be easily changed to R2)
# cvIndex must be a list, each element of which carrying the indices of the rows for a particular fold (I used caret)
# Kris, can you make the code more efficient?
# Note that cvGFLASSO is supporting crossval
# For now, crossval provides the average RMSE per element in the lambda*gamma grid

cvGFLASSO <- function(Y, X, R, opts, cvIndex){
      rmse <- rep(NA, length(cvIndex))
      for(i in 1:length(cvIndex)){
            mod <- gflasso(Y = Y[-cvIndex[[i]],], X = X[-cvIndex[[i]],], R = R, opts = opts)
            pred <- X[cvIndex[[i]],] %*% mod$B
            error <- sqrt(mean((pred - Y[cvIndex[[i]],])**2))
            rmse[i] <- error
      }
      return(rmse)
}

crossval <- function(X, Y, R, params = seq(0,1,by=0.1), cvIndex){
      cvMatrix <- matrix(NA, length(params),length(params))
      dimnames(cvMatrix) <- list(params, params)
      grid <- expand.grid(lambda = params, gamma = params)
      
      for(i in 1:nrow(grid)){
            cv <- cvGFLASSO(X = X, Y = Y, R = R, opts = list(lambda = grid[i,1], gamma = grid[i,2]),
                            cvIndex = myfolds)
            cvMatrix[as.character(grid[i,1]),as.character(grid[i,2])] <- mean(cv)
            message(paste(round((i/(nrow(grid)))*100, 2), "% completion", collapse = " "))
      }
      return(cvMatrix)
}