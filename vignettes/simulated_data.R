
################################################################################
# R script to accompany simulated_data.Rmd
################################################################################

## ---- simulate-data ----
K <- 50 # number of columns of Y
J <- 100 # number of columns of X
J0 <- 10 # number of columns of X actually contribute
n <- 200 # number of samples

B <- matrix(rnorm(J * K), J, K)
B[-sample(J, J0), ] <- 0
X <- matrix(rnorm(n * J), n, J)

Y <- X %*% B + matrix(rnorm(n * K), n, K)

## ---- gflasso ----
R <- matrix(1, K, K) # tie all the responses together
res <- gflasso(Y, X, R, list(delta_conv = 1e-10, lambda = 5, gamma = 5, iter_max = 1000, eps = .5))

## ---- vis-objective ----
plot(res$obj)

## ---- coef-hat-pred ----
plot(res$B, B, asp = 1)
abline(a = 0, b = 1, col = "red")

## ---- vis-coefs ----
image(B == 0)
image(res$B == 0)
table(B == 0, res$B == 0)
