
################################################################################
# R script to accompany simulated_data.Rmd
################################################################################

## ---- sparse-sim-packages ----
# List of packages for session
.packages <- c(
  "gflasso",
  "plyr",
  "dplyr",
  "reshape2",
  "glmnet",
  "Rcpp",
  "ggplot2",
  "data.table"
)

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.packages[!.inst], repos='http://cran.rstudio.com/')
}

# Load packages into session
lapply(.packages, require, character.only=TRUE)
set.seed(04032016)

cat("\014")  # Clear console

rm(list=ls()) # Delete all existing variables
graphics.off() # Close all open plots
theme_set(theme_bw())
sim_methods_dir <- "~/Documents/programming/simFuns/simMethods/"


## ---- functions ----
generate_data <- function(X, mu, eta, K, J, J0, n) {
  beta <- rnorm(J, -mu, 1)
  beta[sample(J)[1:(J / 2)]] <- rnorm(J / 2, -mu, 1)
  B <- beta %*% t(rep(1, K)) + matrix(rnorm(J * K, 0, sqrt(mu)), J, K)
  B[-sample(J, J0), 1 : (K / 2)] <- 0
  B[-sample(J, J0), (K / 2 + 1) : K] <- 0
  Y <- X %*% B + matrix(rnorm(n * K, 0, eta * mu), n, K)
  list(B = B, Y = Y)
}

fit_models <- function(X, Y, R, lambda, gamma, eps, n_iter) {
  gflasso_res <- gflasso(Y, X, R,
                         list(delta_conv = 1, lambda = lambda, gamma = gamma,
                              iter_max = n_iter, verbose = T, eps = eps))
  J <- ncol(X)
  K <- ncol(Y)

  b_lasso <- matrix(0, J, K)
  for (k in seq_len(K)) {
    glmnet_fit <- cv.glmnet(x = X, y = Y[, k], intercept = F)
    b_lasso[, k] <- coef(glmnet_fit)[-1]
  }

  list(gflasso = gflasso_res$B, lasso = b_lasso)
}

save_output <- function(b, b_hat, params, output_dir) {
  timestamp <- gsub("[^0-9]", "", Sys.time())
  feather::write_feather(
    data.frame(type = "lasso", k = seq_len(nrow(b)), params, b_hat$lasso),
    file.path(output_dir, paste0(timestamp, "1.feather"))
  )
  feather::write_feather(
    data.frame(type = "gflasso", k = seq_len(nrow(b)), params, b_hat$gflasso),
    file.path(output_dir, paste0(timestamp, "2.feather"))
  )
  feather::write_feather(
    data.frame(type = "truth", k = seq_len(nrow(b)), params, b),
    file.path(output_dir, paste0(timestamp, "3.feather"))
  )
}

## ---- simdata-params ----
K <- 200 # number of columns of Y
J <- 30 # number of columns of X
J0 <- 5 # number of columns of X actually contribute
n <- 40 # number of samples

R <- matrix(1, K, K)
R[(K / 2 + 1) : K, 1:(K / 2)] <- 0
R[1:(K / 2), (K / 2 + 1) : K] <- 0
n_iter <- 200

lambda <- .15
gamma <- .05
eps <- .1
X <- matrix(rnorm(n * J), n, J)

## ---- run_sims ----
output_dir <- "/Users/krissankaran/Documents/programming/gflasso/sims/"
dir.create(output_dir)

params <- expand.grid(
  eta = 10 ^ (seq(-3, 3, length.out = 10)),
  mu = 10 ^ (seq(-3, 3, length.out = 10))
)

for (i in seq_len(nrow(params))) {
  cat(sprintf("starting params %d \n", i))
  data <- generate_data(X, params$mu[i], params$eta[i], K, J, J0, n)
  b_hat <- fit_models(X, data$Y, R, lambda, gamma, eps, n_iter)
  save_output(data$B, b_hat, params[i,, drop = F], output_dir)
}

## ---- postprocess_results ----
files <- list.files(output_dir, full.names = TRUE)
results <- list()
for (i in seq_along(files)) {
  results[[i]] <- data.table(feather::read_feather(files[i]))
}

results <- rbindlist(results)
head(results)

mresults <- results %>%
  melt(id.vars = c("type", "eta", "k", "mu")) %>%
  dcast(eta + mu + k + variable ~ type, value.var = "value")

stats <- mresults %>%
  group_by(eta, mu) %>%
  summarise(
    mae_lasso = mean(abs(lasso - truth)),
    mae_gflasso = mean(abs(gflasso - truth))
  )
stats
