
# install.packages("grf")
library(grf)



library(ggplot2)
library(data.table)
library(ggpubr)
library(gridExtra)
library(dplyr)

# Empty workspace
rm(list = ls())

sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}
relu <- function(x) {
  return(ifelse(x > 0, x, 0))
}

n1 <- 5000
# calibration
n2 <- 1000
# validation
p_continuous <- 4
p_discrete <- 1
p_demog <- 1
n <- n1 + n2
num_trees <- 1000


set.seed(123456789)

# Covariates = X_cont + Binary + Demog (binary)
X_cont <- matrix(rnorm(n * p_continuous), n, p_continuous)
X_disc <- matrix(rbinom(n * p_discrete, 1, 0.3), n, p_discrete)

# Demographic variable correlated with X[,2]
# Z = as.matrix(rbinom(n, 1, 0.2)) # This is to test whether grf_v2 gives same fit when Z is independent from X.
# Z = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,2]))))

## Higher dimension for demographics -- 1 of them related to X_cont[,2], two other binary and a continuous, all independent
Z <- sapply(1:1, function(x) as.matrix(rbinom(n, 1, 1 / (1 + exp(-X_cont[, 2])))))
# Z <- cbind(Z, rbinom(n, 1, 0.2), rbinom(n, 1, 0.7), rnorm(n, 0, 1)) #

# print(c("validity_check for X,Z dependency",mean(X_cont[Z[,1]==0,2]),mean(X_cont[Z[,1]==1,2])))


# For now, we assume that Z is not a covariate in the model
X <- cbind(X_cont, X_disc)
p <- p_continuous + p_discrete

# Random assigment
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

# calculate_bin_avg = function(condition_var, target_var, num_breaks=256){
#   df = data.table(x = c(condition_var), y= c(target_var))
#   df[, bin := cut(x, breaks = num_breaks)]
#   df[, bin_avg  := mean(y), by=bin]
#   # return(unique(df[, .(bin, bin_avg)]))
#   return(df$bin_avg)
# }
# a = calculate_bin_avg(X[,1], Z[, 1], num_breaks = 256)
# uniqueN(a)
# ggbarplot(data=a, x='bin', y='bin_avg')
#

# validate_data = data.frame(Y=Y[c(n1:n2)],Z=Z[c(n1:n2)],W=W[c(n1:n2)],tau = tau[c(n1:n2)],X=X[c(n1:n2),])
# X_validate = as.matrix(validate_data[,c(5:dim(validate_data)[2])])
# Y_validate = validate_data$Y
# W_validate = validate_data$W
# Z_validate = validate_data$Z


## Some validity plots to make sure that tau discriminates in Z
if (FALSE) {
  data_to_plot <- data.table(Y = Y, W = W, Z = Z[, 1], tau)
  data_to_plot[, W := as.character(W)]
  data_to_plot[, Z := as.character(Z)]
  p_ate <- ggdensity(data = data_to_plot, x = "Y", color = "W", fill = "W", alpha = 0.2, add = "mean")
  p_bias <- ggdensity(data = data_to_plot, x = "tau", color = "Z", fill = "Z", alpha = 0.2, add = "mean")
  grid.arrange(p_ate, p_bias, nrow = 1)
}

calculate_bin_avg = function(condition_var, target_var, num_breaks=256){
  df = data.table(x = c(condition_var), y= c(target_var))
  df[, bin := cut(x, breaks = num_breaks)]
  df[, bin_avg  := mean(y), by=bin]
  # return(unique(df[, .(bin, bin_avg)]))
  return(df$bin_avg)
}
target.weight.X = sapply(X_train, calculate_bin_avg, target_var = Z_train$Z)

fit_grf_v2 <- causal_forest(X_train, Y_train, W_train,
                            target.weights = as.matrix(Z_train),
                            target.weight.penalty = 50,
                            target.weight.X = target.weight.X,
                            honesty = TRUE,
                            #W.hat = W.hat, Y.hat = Y.hat,
                            num.trees = num_trees,
                            # tune.parameters ='target.weight.penalty',#c('target.weight.penalty',my_tuning_param), # "target.weight.penalty",# my_tuning_param, # "none",
                            # sample.fraction = fit_grf$tuning.output$params$sample.fraction,
                            # mtry = fit_grf$tuning.output$params$mtry,
                            # min.node.size = fit_grf$tuning.output$params$min.node.size,
                            # honesty.fraction = fit_grf$tuning.output$params$honesty.fraction,
                            # honesty.prune.leaves = fit_grf$tuning.output$params$honesty.prune.leaves,
                            # alpha = fit_grf$tuning.output$params$alpha,
                            # imbalance.penalty = fit_grf$tuning.output$params$imbalance.penalty,
                            # num.threads = 1,
                            seed=1
)

tau_train.grf_v2 <- predict(fit_grf_v2, estimate.variance = TRUE) # newdata=X_validate,
plot(tau[1:n1],tau_train.grf_v2$predictions )
data_to_plot2 <-data.table(Z=as.character(c(Z_train$Z)), tau = tau_train.grf_v2$predictions)

ggdensity(
  data = data_to_plot2, x = "tau", color = "Z", fill = "Z", alpha = 0.2, add = "mean")



## ---------------------------------------------------
##   Estimate 'tau' using GRF
## ---------------------------------------------------
all_tuning_param <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty")


# propensity scores
forest.W <- regression_forest(X_train, W_train)
W.hat <- predict(forest.W)$predictions


# Y_hat, for GRF
forest.Y <- regression_forest(X_train, Y_train)
Y.hat <- predict(forest.Y)$predictions

# Set spect for regular forest
# my_tuning_param = c("none")


# Estimate REGULAR Causal forest (using 'grf' package)
start.time <- Sys.time()

fit_grf <- causal_forest(X_train, Y_train, W_train,
                         honesty = TRUE,
                         W.hat = W.hat, Y.hat = Y.hat,
                         num.trees = num_trees,
                         # tune.parameters = my_tuning_param,
                         seed=1
)
end.time <- Sys.time()
time_cf_estim <- end.time - start.time
print(time_cf_estim)
tau_train.grf <- predict(fit_grf, estimate.variance = TRUE) # newdata=X_validate,


# Estimate Causal forest with TARGET WEIGHTS, using the tuneable parameters from the regular forest
start.time <- Sys.time()


calculate_bin_avg = function(condition_var, target_var, num_breaks=256){
  df = data.table(x = c(condition_var), y= c(target_var))
  df[, bin := cut(x, breaks = num_breaks)]
  df[, bin_avg  := mean(y), by=bin]
  # return(unique(df[, .(bin, bin_avg)]))
  return(df$bin_avg)
}
target.weight.X = sapply(X_train, calculate_bin_avg, target_var = Z_train$Z)

fit_grf_v2 <- causal_forest(X_train, Y_train, W_train,
                            target.weights = as.matrix(Z_train),
                            target.weight.penalty = 100,
                            target.weight.X = target.weight.X,
                            honesty = TRUE,
                            #W.hat = W.hat, Y.hat = Y.hat,
                            num.trees = num_trees,
                            # tune.parameters ='target.weight.penalty',#c('target.weight.penalty',my_tuning_param), # "target.weight.penalty",# my_tuning_param, # "none",
                            # sample.fraction = fit_grf$tuning.output$params$sample.fraction,
                            # mtry = fit_grf$tuning.output$params$mtry,
                            # min.node.size = fit_grf$tuning.output$params$min.node.size,
                            # honesty.fraction = fit_grf$tuning.output$params$honesty.fraction,
                            # honesty.prune.leaves = fit_grf$tuning.output$params$honesty.prune.leaves,
                            # alpha = fit_grf$tuning.output$params$alpha,
                            # imbalance.penalty = fit_grf$tuning.output$params$imbalance.penalty,
                            # num.threads = 1,
                            seed=1
)

tau_train.grf_v2 <- predict(fit_grf_v2, estimate.variance = TRUE) # newdata=X_validate,

data_to_plot <- data.table(
  Z = Z_train[, 1],
  tau_grf_no_weight = tau_train.grf$predictions,
  tau_grf_with_weight = tau_train.grf_v2$predictions
) #
names(data_to_plot)[1] <- "Z"
data_to_plot[, Z := as.character(Z)]

data_to_plot2 <- melt(data_to_plot, id.vars = "Z")

ggdensity(
  data = data_to_plot2, x = "value", color = "Z", fill = "Z", alpha = 0.2, add = "mean",
  facet.by = "variable", ncol = 1,
  scales = "free_y"
)



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

  tau.hat <- predict(fit_grf_v1, X.marginals, estimate.variance = TRUE)
  sigma.hat <- sqrt(tau.hat$variance.estimates)

  tau.weight.hat <- predict(fit_grf_v2, X.marginals, estimate.variance = TRUE)
  sigma.weight.hat <- sqrt(tau.weight.hat$variance.estimates)

  out = data.table(x=x, x_dim=x_dim)
  out[, tau.true :=tau_function(X.marginals)]

  out.no.weight = data.table(
    x = x, x_dim = x_dim,
    tau.pred = tau.hat$predictions,
    tau.pred.upper=  tau.hat$predictions + 1.96 * sigma.hat,
    tau.pred.lower =  tau.hat$predictions - 1.96 * sigma.hat,
    tau.type = 'no weight'
  )
  out.with.weight = data.table(
    x = x, x_dim = x_dim,
    tau.pred = tau.weight.hat$predictions,
    tau.pred.upper=  tau.weight.hat$predictions + 1.96 * sigma.weight.hat,
    tau.pred.lower =  tau.weight.hat$predictions - 1.96 * sigma.weight.hat,
    tau.type = 'with weight'
  )
  out = rbindlist(list(out.no.weight,out.with.weight))
  out[, tau.true :=tau_function(X.marginals)]
  return(out)
}

dat_pred <- rbindlist(lapply(1:p, generate_prediction_per_dim))
dat_pred = melt(dat_pred, id.vars = c('x_dim', 'x', grep("upper|lower", names(dat_pred), value=TRUE)))


ggline(dat_pred, x='x',
       y='value',
       short.panel.labs = FALSE,
       facet.by = 'x_dim',
       scale='free_y',
       color='variable',numeric.x.axis = TRUE)


generate_prediction_per_dim_cont <- function(x_dim) {
  x <- seq(-2, 2, length.out = 101)
  X.marginals <- matrix(0, 101, p)
  X.marginals[, x_dim] <- x
  X.marginals <- as.data.frame(X.marginals)

  tau.hat <- predict(fit_grf_v1, X.marginals, estimate.variance = TRUE)
  mean.hat <- tau.hat$predictions
  sigma.hat <- sqrt(tau.hat$variance.estimates)

  tau.true <- tau_function(X.marginals)
  upper <- mean.hat + 1.96 * sigma.hat
  lower <- mean.hat - 1.96 * sigma.hat
  out <- data.table(x = x, tau_pred = mean.hat, tau_true = tau.true, tau_upper = upper, tau_lower = lower, x_dim = x_dim)
  return(out)
}


dat_pred_cont <- rbindlist(lapply(1:p_continuous, generate_prediction_per_dim_cont))

ggline(
  data = dat_pred_cont, x = "x", y = "tau_pred", facet.by = "x_dim",
  scale='free_y',
  plot_type = "l", xlab = "x", ylab = "tau", title = "Marginal Effects - Continuous var.",
  numeric.x.axis = TRUE, short.panel.labs = FALSE
) +
  geom_line(mapping = aes(x = x, y = tau_true), color = "red") +
  geom_line(mapping = aes(x = x, y = tau_upper), linetype = "dashed") +
  geom_line(mapping = aes(x = x, y = tau_lower), linetype = "dashed")

