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
Z <- sapply(1:1, function(x) as.matrix(sigmoid(X_cont[, 2]) ))
Z <- cbind( rbinom(n, 1, sigmoid(X_cont[, 2] * 10)) , rbinom(n, 1, 0.2), rpois(n, sigmoid(X_cont[,4])))
#

X <- cbind(X_cont, X_disc)
p <- p_continuous + p_discrete


# Simulate 'Y'
y_function <- function(x) {
  temp <- x[, 1] - 2 * x[, 2] + x[, 4]
  return(temp)
}
noise_y <- runif(n)
Y <- y_function(X)  + noise_y


# For now, let's use only train data
train_data <- data.table(Y = Y[c(1:n1)], Z = Z[c(1:n1), ], X = X[c(1:n1), ])
x_cols <- grep("X", names(train_data), value = TRUE)
X_train <- train_data[, .SD, .SDcols = x_cols]
Y_train <- train_data$Y
Zcols <- grep("Z", names(train_data), value = TRUE)
Z_train <- train_data[, .SD, .SDcols = Zcols]


## Some validity plots to make sure that tau discriminates in Z
if (FALSE) {
  # Tau condition on  W and Tau
  data_to_plot <- data.table(Y = Y, Z = Z[, 2])
  data_to_plot[, Z := as.character(Z)]
  ggdensity(data = data_to_plot, x = "Y", color = "Z", fill = "Z", alpha = 0.2, add = "mean")

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
fit_grf <- regression_forest(X_train, Y_train,
                         honesty = TRUE,
                         num.trees = num_trees,
                         # num.threads=1,
                         #tune.parameters='all',
                         seed=1
)
print(Sys.time() - time1)

time1 = Sys.time()
fit_grf_v2 <- balanced_regression_forest(X_train, Y_train,
                            target.weights = as.matrix(Z_train),
                            target.weight.penalty = 3,
                            target.weight.standardize = TRUE,
                            # target.weights.hat = 0,
                            target.weight.bins.breaks = 256,
                            honesty = TRUE,
                            num.trees = num_trees, #10, #num_trees,
                            #tune.parameters = "all",
                            # num.threads = 1,
                            seed=1
)
print(Sys.time() - time1)

fit_grf_v2$tuning.output

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
  title='Regression Y Density',
  subtitle='Weight 3'
)
ggsave('regression_y_density_z_weight_3.png', width = 8, height=8)


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
  y_hat <- predict(fit_grf, X.marginals, estimate.variance = TRUE)
  sigma.hat <- sqrt(y_hat$variance.estimates)

  y.weight.hat <- predict(fit_grf_v2, X.marginals, estimate.variance = TRUE)
  sigma.weight.hat <- sqrt(y.weight.hat$variance.estimates)

  # out = data.table(x=x, x_dim=x_dim)
  # out[, tau.true :=tau_function(X.marginals)]

  out.no.weight = data.table(
    x = x, x_dim = x_dim,
    y.true = y_function(X.marginals)$V1,
    y.pred = y_hat$predictions,
    y.pred.upper=  y_hat$predictions + 1.96 * sigma.hat,
    y.pred.lower =  y_hat$predictions - 1.96 * sigma.hat,
    y.type = 'no weight'
  )
  out.with.weight = data.table(
    x = x, x_dim = x_dim,
    y.true = y_function(X.marginals )$V1,
    y.pred = y.weight.hat$predictions,
    y.pred.upper=  y.weight.hat$predictions + 1.96 * sigma.weight.hat,
    y.pred.lower =  y.weight.hat$predictions - 1.96 * sigma.weight.hat,
    y.type = 'fair tree'
  )
  out = rbindlist(list(out.no.weight,out.with.weight))
  return(out)
}

dat_pred <- rbindlist(lapply(1:p, generate_prediction_per_dim))

ggline(dat_pred, x='x',
       y='y.pred',
       short.panel.labs = FALSE,
       facet.by = 'x_dim',
       scale='free_y',
       title='Regression Y',
       xlab='Marginal X',
       subtitle='Weight: 3',
       color='y.type',numeric.x.axis = TRUE) +
  geom_line(mapping = aes(x = x, y = y.true), color = "purple") +
  geom_line(mapping = aes(x = x, y = y.pred.upper, color = y.type), linetype='dashed') +
  geom_line(mapping = aes(x = x, y = y.pred.lower, color = y.type), linetype='dashed')

ggsave("reg_y_x_weight_10.png", width=14, height=8)

