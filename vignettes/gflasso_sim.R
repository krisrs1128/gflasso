
################################################################################
# R script to accompany simulated_data.Rmd
################################################################################

## ---- sparse-sim-packages ----
# List of packages for session
.packages = c("gflasso",
              "reshape2",
              "ggplot2")

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

## ---- simdata-params ----
K <- 50 # number of columns of Y
J <- 100 # number of columns of X
J0 <- 20 # number of columns of X actually contribute
n <- 200 # number of samples

## ---- gflasso-simulate-data ----
beta <- rnorm(J)
B <- beta %*% t(rep(1, K)) + matrix(rnorm(J * K, 0, .5), J, K)
B[-sample(J, J0), 1 : (K / 2)] <- 0
B[-sample(J, J0), (K / 2 + 1) : K] <- 0
X <- matrix(rnorm(n * J), n, J)
Y <- X %*% B + matrix(rnorm(n * K), n, K)

## ---- vis-B ----
m_b <- melt(B, varnames = c("feature", "task"), value.name = "beta")
ggplot(m_b) +
  geom_tile(aes(x = task, y = feature, fill = beta)) +
  scale_fill_gradient2(midpoint = 0, high = "#90ee90", low = "#000080") +
  ggtitle("True Regression Coefficients") +
  theme(panel.grid = element_blank())

## ---- vis-reg-data ----
m_y <- melt(Y, varnames = c("sample", "task"), value.name = "y")
m_x <- melt(X, varnames = c("sample", "feature"), value.name = "x")
reg_data <- m_y %>%
  left_join(m_x) %>%
  left_join(m_b)
cur_data <- reg_data %>%
  filter(task %in% 20:30,
         feature %in% 10:18)

## ---- vis-reg-lines ----
ggplot(cur_data) +
  geom_point(aes(x = x, y = y, col = as.factor(beta == 0)), 
             size = .3, alpha = 0.3) +
  geom_abline(data = cur_data %>% filter(beta != 0),
              aes(slope = beta, intercept = 0, col = as.factor(beta == 0))) +
  scale_color_manual(values = c("#8068ab", "#d9bad8")) +
  facet_grid(feature ~ task) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

## ---- gflasso ----
R <- matrix(1, K, K)
R[(K / 2 + 1) : K, 1:(K / 2)] <- 0
R[1:(K / 2), (K / 2 + 1) : K] <- 0
res <- gflasso(Y, X, R, list(delta_conv = 1e-10, lambda = 10, gamma = 10, iter_max = 400, eps = .5))

## ---- vis-objective ----
plot(res$obj)

## ---- coef-hat-pred ----
mB_compare <- melt(list(truth = B, fit = res$B)) %>%
  dcast(Var1 + Var2 ~ L1)
mB_compare$group <- ifelse(mB_compare$Var2 <= K / 2, "group_1", "group_2")
ggplot(mB_compare) +
  geom_point(aes(x = truth, y = fit, col = group)) +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  ggtitle("True vs. Fitted Coefficients")

## ---- vis-coefs ----
mres_B <- melt(res$B)
mres_B$small <- mres_B$value < 5e-4

ggplot(mres_B) +
  geom_tile(aes(x = Var2, y = Var1, fill = small)) +
  scale_fill_brewer(palette = "RdBu")

table(B == 0, abs(res$B) < 1e-5)
