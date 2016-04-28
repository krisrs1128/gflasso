
################################################################################
# R script to accompany simulated_data.Rmd
################################################################################

## ---- sparse-sim-packages ----
# List of packages for session
.packages = c("gflasso",
              "plyr",
              "dplyr",
              "reshape2",
              "glmnet",
              "Rcpp",
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
theme_set(theme_bw())
sim_methods_dir <- "~/Documents/programming/simFuns/simMethods/"

## ---- simdata-params ----
K <- 200 # number of columns of Y
J <- 30 # number of columns of X
J0 <- 5 # number of columns of X actually contribute
n <- 40 # number of samples

## ---- gflasso-simulate-data ----
beta <- rnorm(J, .5, 1)
beta[sample(J)[1:(J / 2)]] <- rnorm(J / 2, -.5, 1)
B <- beta %*% t(rep(1, K)) + matrix(rnorm(J * K, 0, .1), J, K)
B[-sample(J, J0), 1 : (K / 2)] <- 0
B[-sample(J, J0), (K / 2 + 1) : K] <- 0
X <- matrix(rnorm(n * J), n, J)
Y <- X %*% B + matrix(rnorm(n * K, 0, .5), n, K)

## ---- vis-B ----
m_b <- melt(B, varnames = c("feature", "task"), value.name = "beta")
ggplot(m_b) +
  geom_tile(aes(x = task, y = feature, fill = beta)) +
  scale_fill_gradient2(midpoint = 0, high = "#90ee90", low = "#000080") +
  theme(panel.grid = element_blank())

## ---- vis-reg-data ----
m_y <- melt(Y, varnames = c("sample", "task"), value.name = "y")
m_x <- melt(X, varnames = c("sample", "feature"), value.name = "x")
reg_data <- m_y %>%
  left_join(m_x) %>%
  left_join(m_b)
cur_data <- reg_data %>%
  filter(task %in% 95:105,
         feature %in% 1:10)

## ---- vis-reg-lines ----
ggplot(cur_data) +
  geom_point(aes(x = x, y = y, col = beta == 0),
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
n_iter <- 200

lambda <- .15
gamma <- .05
eps <- .1
gflasso_res <- gflasso(Y, X, R,
                       list(delta_conv = 1, lambda = lambda, gamma = gamma,
                            iter_max = n_iter, verbose = T, eps = eps))

## ---- vis-objective ----
obj_data <- data.frame(iteration = seq_along(gflasso_res$obj),
                       objective = gflasso_res$obj)
ggplot(obj_data) +
  geom_line(aes(x = iteration, y = objective)) +
  theme(panel.grid = element_blank())

## ---- coef-hat-pred ----
b_compare <- melt(list(truth = B, fit = gflasso_res$B)) %>%
  dcast(Var1 + Var2 ~ L1)
b_compare$group <- ifelse(b_compare$Var2 <= K / 2, "task_group_1", "task_group_2")
ggplot(b_compare) +
  geom_abline(slope = 1, intercept = 0, col = '#696969', alpha = 0.8, size = .3) +
  geom_vline(xintercept = 0, col = '#696969', alpha = 0.5, size = .05) +
  geom_hline(yintercept = 0, col = '#696969', alpha = 0.5, size = .05) +
  geom_point(aes(x = truth, y = fit, col = group), size = .5, alpha = 0.6) +
  scale_color_manual(values = c("#008080", "#cd5c5c")) +
  coord_fixed() +
  theme(panel.grid = element_blank())

## ---- vis-Bhat-gflasso ----
gflasso_mb <- melt(gflasso_res$B, varnames = c("feature", "task"),
               value.name = "beta_hat")

ggplot(gflasso_mb) +
  geom_tile(aes(x = task, y = feature, fill = beta_hat)) +
  scale_fill_gradient2(midpoint = 0, high = "#90ee90", low = "#000080") +
  theme(panel.grid = element_blank())

## ---- gflasso-reg-bhat ----
reg_fit_data <- m_y %>%
  left_join(m_x) %>%
  left_join(m_b) %>%
  left_join(gflasso_mb)
reg_fit_data$zero <- reg_fit_data$beta == 0
cur_fit_data <- reg_fit_data %>%
  filter(task %in% 95:105,
         feature %in% 1:10) %>%
  melt(measure.vars = c("beta", "beta_hat"), variable.name = "coef")

## ---- vis-reg-fit ----
ggplot(cur_fit_data) +
  geom_point(aes(x = x, y = y, col = zero), 
             size = .3, alpha = 0.05) +
  geom_abline(data = cur_fit_data,
              aes(slope = value, intercept = 0, col = zero, linetype = coef),
              alpha = 0.9) +
  scale_color_manual(values = c("#8068ab", "#d9bad8")) +
  facet_grid(feature ~ task) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

## ---- class-accuracy ----
table(B == 0, abs(gflasso_res$B) < 5e-2)

## ---- lasso-baseline ----
b_lasso <- matrix(0, J, K)
for (k in seq_len(K)) {
  glmnet_fit <- cv.glmnet(x = X, y = Y[, k], intercept = F)
  b_lasso[, k] <- coef(glmnet_fit)[-1]
}

## ---- lasso-bhat-b ----
b_compare <- melt(list(truth = B, fit = b_lasso)) %>%
  dcast(Var1 + Var2 ~ L1)
b_compare$group <- ifelse(b_compare$Var2 <= K / 2, "task_group_1", "task_group_2")
ggplot(b_compare) +
  geom_abline(slope = 1, intercept = 0, col = '#696969', alpha = 0.8, size = .3) +
  geom_vline(xintercept = 0, col = '#696969', alpha = 0.5, size = .05) +
  geom_hline(yintercept = 0, col = '#696969', alpha = 0.5, size = .05) +
  geom_point(aes(x = truth, y = fit, col = group), size = .5, alpha = 0.6) +
  scale_color_manual(values = c("#008080", "#cd5c5c")) +
  coord_fixed() +
  theme(panel.grid = element_blank())

## ---- vis-Bhat-lasso ----
lasso_mb <- melt(b_lasso, varnames = c("feature", "task"),
                 value.name = "beta_hat_lasso")

ggplot(lasso_mb) +
  geom_tile(aes(x = task, y = feature, fill = beta_hat_lasso)) +
  scale_fill_gradient2(midpoint = 0, high = "#90ee90", low = "#000080") +
  theme(panel.grid = element_blank())

## ---- lasso-reg-bhat ----
reg_fit_data <- m_y %>%
  left_join(m_x) %>%
  left_join(m_b) %>%
  left_join(lasso_mb)

reg_fit_data$zero <- reg_fit_data$beta == 0
cur_fit_data <- reg_fit_data %>%
  filter(task %in% 95:105,
         feature %in% 1:10) %>%
  melt(measure.vars = c("beta", "beta_hat_lasso"), variable.name = "coef")

## ---- vis-lasso-fit ----
ggplot(cur_fit_data) +
  geom_point(aes(x = x, y = y, col = zero), 
             size = .3, alpha = 0.05) +
  geom_abline(data = cur_fit_data,
              aes(slope = value, intercept = 0, col = zero, linetype = coef),
              alpha = 0.9) +
  scale_color_manual(values = c("#8068ab", "#d9bad8")) +
  facet_grid(feature ~ task) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
