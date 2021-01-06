rm(list = ls())

library(grf)
library(ggplot2)
library(data.table)
library(ggpubr)
library(gridExtra)
library(dplyr)


sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}
relu <- function(x) {
  return(ifelse(x > 0, x, 0))
}

n1 <- 5000
# calibration
n2 <- 1000
num_trees <-  4000
# validation
p_continuous <- 4
p_discrete <- 1
p_demog <- 1
n <- n1 + n2


set.seed(123456789)

X_cont <- matrix(rnorm(n * p_continuous), n, p_continuous)
X_disc <- matrix(rbinom(n * p_discrete, 1, 0.3), n, p_discrete)
Z <- sapply(1:1, function(x) as.matrix(rbinom(n, 1, sigmoid(X_cont[, 2]*10))))
Z <- cbind(Z, rbinom(n, 1, 0.2), rbinom(n, 1, 0.7), rnorm(n, 0, 1))

X <- cbind(X_cont, X_disc)
p <- p_continuous + p_discrete

W <- rbinom(n, 1, 0.5)

# Simulate 'tau' -- TAU depends on X[2] but no on Z
tau_function <- function(x) {
  temp <- (-1 + pmax(x[, 1], 0) + x[, 2] + abs(x[, 3]) + x[, 5])
  return(temp)
}
tau <- tau_function(X)

# Simulate 'Y'
y_function <- function(x, z) {
  temp <- x[, 1] - 2 * x[, 2] + x[, 4] + z[, 1]
  return(temp)
}
noise_y <- runif(n)
Y <- y_function(X, Z) + tau * W + noise_y


# For now, let's use only train data
train_data <- data.table(Y = Y[c(1:n1)], Z = Z[c(1:n1), ], W = W[c(1:n1)], tau = tau[c(1:n1)], X = X[c(1:n1), ])
x_cols <- grep("X", names(train_data), value = TRUE)
X_train <- train_data[, .SD, .SDcols = x_cols]
Y_train <- train_data$Y
W_train <- train_data$W
Zcols <- grep("Z", names(train_data), value = TRUE)
Z_train <- train_data[, .SD, .SDcols = Zcols]


## Some validity plots to make sure that tau discriminates in Z
if (FALSE) {
  # Tau condition on  W and Tau
  data_to_plot <- data.table(Y = Y, W = W, Z = Z[, 1], tau)
  data_to_plot[, W := as.character(W)]
  data_to_plot[, Z := as.character(Z)]
  p_ate <- ggdensity(data = data_to_plot, x = "Y", color = "W", fill = "W", alpha = 0.2, add = "mean")
  p_bias <- ggdensity(data = data_to_plot, x = "tau", color = "Z", fill = "Z", alpha = 0.2, add = "mean",
                      title='Z_1')
  grid.arrange(p_ate, p_bias, nrow = 1)
  ggsave("Y_tau_density.png", width=12, height=6)
  # tau ~ x
  colnames(Z) = paste0("Z_", 1:dim(Z)[2])
  tau_z = data.table(tau, Z)
  tau_z[, tau_bin := as.numeric(cut(tau, 256))]

  tau_z = melt(tau_z, id.vars = c('tau','tau_bin'))

  tau_z = tau_z[, .(avg_z = mean(value)), by=.(tau_bin , variable)]
  ggscatter(data=tau_z, x='tau_bin', y='avg_z', facet.by='variable',
            add = 'reg.line', add.params = list(color='red')) +
    labs(x='Bins (Evenly Distributed 256 Bins for Tau)', y='Average Z',
         title='Average Z conditioned on Tau',
         subtitle = 'Red line is simple regression; bin is the same range for all')
  ggsave("average_z_bin_tau.png", width=12, height=8)
  # Y ~ Z
  Y_Z = data.table(Z, Y)
  Y_Z[, y_bin := as.numeric(cut(Y, 256))]

  Y_Z = melt(Y_Z, id.vars = c('Y','y_bin'))

  Y_Z = Y_Z[, .(avg_y = mean(value)), by=.(y_bin , variable)]
  ggscatter(data=Y_Z, x='y_bin', y='avg_y', facet.by='variable',
            add = 'reg.line', add.params = list(color='red')) +
    labs(x='Bins (Evenly Distributed 256 Bins for Y)', y='Average Z',
         title='Average Z conditioned on Y',
         subtitle = 'Red line is simple regression; bin is the same range for all')
  ggsave("average_z_bin_y.png", width=12, height=8)

}

# fit_grf_v2 <- causal_forest(X_train, Y_train, W_train,
#                             Y.hat = 0, W.hat=0,
#                             target.weights = as.matrix(Z_train),
#                             target.weight.penalty = 10,
#                             target.weight.bins.breaks = 256,
#                             #target.weights.hat=0,
#                             honesty = TRUE,
#                             num.trees = num_trees, #10, #num_trees,
#                             # tune.parameters = "all",
#                             num.threads = 1,
#                             seed=1
# )
## ---------------------------------------------------
##   Estimate 'tau' using GRF
## ---------------------------------------------------


# X_hat = regression_forest(X = as.matrix(Z_train), Y = as.matrix(X_train[, 1]), num.trees = num_trees, seed=1,
                          # tune.parameters = 'all')
# X_org = X_train - predict(X_hat)$prediction

time1 = Sys.time()
fit_grf <- causal_forest(X_train, Y_train, W_train,
                         honesty = TRUE,
                         #Y.hat = 0, W.hat=0,
                         #target.weights = as.matrix(Z_train) ,
                         target.weight.penalty=0,
                         num.trees = num_trees,
                         # num.threads=1,
                         tune.parameters='all',
                         seed=1
)
print(Sys.time() - time1)

time1 = Sys.time()
fit_grf_v2 <- causal_forest(X_train, Y_train, W_train,
                            target.weights = as.matrix(Z_train),
                            target.weight.penalty = 0,
                            target.weight.standardize = TRUE,
                            # target.weights.hat = 0,
                            target.weight.bins.breaks = 256,
                            honesty = TRUE,
                            num.trees = num_trees, #10, #num_trees,
                            tune.parameters = "all",
                            # num.threads = 1,
                            seed=1
)
print(Sys.time() - time1)



tau_train.grf <- predict(fit_grf, estimate.variance = TRUE) # newdata=X_validate,
tau.pred <- tau_train.grf$predictions

tau_train.grf_v2 <- predict(fit_grf_v2, estimate.variance = TRUE) # newdata=X_validate,
tau.pred_v2 <- tau_train.grf_v2$predictions

data_to_plot <- data.table(
  Z = Z_train[, 1],
  tau_grf_no_weight = tau_train.grf$predictions,
  tau_grf_fair_tree = tau_train.grf_v2$predictions
) #
names(data_to_plot)[1] <- "Z"
data_to_plot[, Z := as.character(Z)]

data_to_plot2 <- melt(data_to_plot, id.vars = "Z")

ggdensity(
  data = data_to_plot2, x = "value", color = "Z", fill = "Z", alpha = 0.2, add = "mean",
  facet.by = "variable", ncol = 1,
  scales = "free_y",
  title='Tau Density',
  subtitle='Weight 10'
)
ggsave('tau_density_z_weight_10.png', width = 8, height=8)


## -------------------------------------------
## Generate Marginals (moderating effects)
## -------------------------------------------

generate_prediction_per_dim= function(x_dim) {
  if(x_dim<5){
    x <- seq(-2, 2, length.out = 101)
    X.marginals <- matrix(0, 101, p)
  }else if (x_dim==5){
    x <- c(0, 1)
    X.marginals <- matrix(0, 2, p)
  }
  X.marginals[, x_dim] <- x
  X.marginals <- as.data.table(X.marginals)
  # X.marginals.hat = predict(X_hat, newdata=X.marginals)$prediction
  tau.hat <- predict(fit_grf, X.marginals, estimate.variance = TRUE)
  sigma.hat <- sqrt(tau.hat$variance.estimates)

  tau.weight.hat <- predict(fit_grf_v2, X.marginals, estimate.variance = TRUE)
  sigma.weight.hat <- sqrt(tau.weight.hat$variance.estimates)

  # out = data.table(x=x, x_dim=x_dim)
  # out[, tau.true :=tau_function(X.marginals)]

  out.no.weight = data.table(
    x = x, x_dim = x_dim,tau.true = tau_function(as.matrix(X.marginals)),
    tau.pred = tau.hat$predictions,
    tau.pred.upper=  tau.hat$predictions + 1.96 * sigma.hat,
    tau.pred.lower =  tau.hat$predictions - 1.96 * sigma.hat,
    tau.type = 'no weight'
  )
  out.with.weight = data.table(
    x = x, x_dim = x_dim, tau.true = tau_function(as.matrix(X.marginals)),
    tau.pred = tau.weight.hat$predictions,
    tau.pred.upper=  tau.weight.hat$predictions + 1.96 * sigma.weight.hat,
    tau.pred.lower =  tau.weight.hat$predictions - 1.96 * sigma.weight.hat,
    tau.type = 'fair tree'
  )
  out = rbindlist(list(out.no.weight,out.with.weight))
  return(out)
}

dat_pred <- rbindlist(lapply(1:p, generate_prediction_per_dim))

ggline(dat_pred, x='x',
       y='tau.pred',
       short.panel.labs = FALSE,
       facet.by = 'x_dim',
       scale='free_y',
       title='Marginal Tau',
       subtitle='Weight: 10',
       color='tau.type',numeric.x.axis = TRUE) +
  geom_line(mapping = aes(x = x, y = tau.true), color = "purple") +
  geom_line(mapping = aes(x = x, y = tau.pred.upper, color = tau.type), linetype='dashed') +
  geom_line(mapping = aes(x = x, y = tau.pred.lower, color = tau.type), linetype='dashed')
ggsave("tau_x_weight_10.png", width=14, height=8)

