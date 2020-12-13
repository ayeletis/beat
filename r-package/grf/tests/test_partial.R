
# install.packages("grf")
library(grf)


library(ggplot2)
library(data.table)
library(ggpubr)
library(gridExtra)
library(dplyr)

# Empty workspace
rm(list = ls())
dev.off()
n1 <- 5000
# calibration
n2 <- 1000
num_trees <- 4000
# validation
p_continuous <- 4
p_discrete <- 1
p_demog <- 1
n <- n1 + n2


set.seed(123456789)

# Covariates = X_cont + Binary + Demog (binary)
X_cont <- matrix(rnorm(n * p_continuous), n, p_continuous)
X_disc <- matrix(rbinom(n * p_discrete, 1, 0.3), n, p_discrete)

## Higher dimension for demographics -- 1 of them related to X_cont[,2], two other binary and a continuous, all independent
Z <- sapply(1:1, function(x) as.matrix(rbinom(n, 1, 1 / (1 + exp(-X_cont[, 2])))))
X <- cbind(X_cont, X_disc)
p <- p_continuous + p_discrete

# Random assigment
W <- rbinom(n, 1, 0.5)

# Simulate 'tau' -- TAU depends on X[2] but no on Z
tau_function <- function(x) {
  temp <- (-1 + pmax(x[, 1], 0) + 5* x[, 2] + abs(x[, 3]) + x[, 5])
  # temp <- 5*x[, 2]
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
if (TRUE) {
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



forest.W <- regression_forest(X_train, W_train  )
W.hat <- predict(forest.W)$predictions

# Y_hat, for GRF
forest.Y <- regression_forest(X_train, Y_train )
Y.hat <- predict(forest.Y)$predictions

fit_grf_v1 <- causal_forest(X_train, Y_train, W_train,
                            target.weights = as.matrix(Z_train) * 0,
                            target.weight.penalty = 0,
                            target.weight.X = target.weight.X* 0,
                            honesty = TRUE,
                            W.hat = W.hat, Y.hat = Y.hat,
                            num.trees = num_trees,
                            # num.threads = 1,
                            seed=1
)
tau_train.grf_v1 <- predict(fit_grf_v1, estimate.variance = TRUE) # newdata=X_validate,

fit_grf_v2 <- causal_forest(X_train, Y_train, W_train,
                            target.weights = as.matrix(Z_train),
                            target.weight.penalty =100,
                            target.weight.X = target.weight.X,
                            honesty = TRUE,
                            W.hat = W.hat, Y.hat = Y.hat,
                            num.trees = num_trees,
                            # num.threads = 1,
                            seed=1
)

tau_train.grf_v2 <- predict(fit_grf_v2, estimate.variance = TRUE) # newdata=X_validate,


data_to_plot <- data.table(
  Z = Z_train[, 1],
  tau_grf_no_weight = tau_train.grf_v1$predictions,
  tau_grf_with_weight = tau_train.grf_v2$predictions
) #
names(data_to_plot)[1] <- "Z"
data_to_plot[, Z := as.character(Z)]

data_to_plot2 <- melt(data_to_plot, id.vars = "Z")

ggdensity(
  data = data_to_plot2, x = "value", color = "Z",
  facet.by ='variable',ncol=1,
  fill = "Z", alpha = 0.2, add = "mean")

dat_plot = copy(X_train)
dat_plot[, true_tau := tau[1:n1]]
dat_plot[, tau_pred_1 := tau_train.grf_v1$predictions]
dat_plot[, tau_pred_2 := tau_train.grf_v2$predictions]
dat_plot[, Z:= as.character(Z_train$Z)]

get_one_plot = function(x){
  df = copy(dat_plot)
  df = df[, bins := cut(get(x),breaks = 256)]
  df = df[, .(tau_mean_true = mean(true_tau),
                   tau_hat_with_weight = mean(tau_pred_2),
                   tau_hat_no_weight=  mean(tau_pred_1)), by=.(bins, Z)]
  df[,bins:=as.numeric(bins)]
  df = melt(df, id.vars = c('Z','bins','tau_mean_true'))
    p = ggline(data=df, x='tau_mean_true',
             facet.by = 'Z', ncol=1,short.panel.labs = FALSE,
             y='value', color='variable', fill='variable', alpha=0.2,
             numeric.x.axis = TRUE) +
      geom_abline(slope = 1, intercept = 0, color='purple', linetype='dashed') +
    labs(title=x)
  return(p)
}
# a = lapply(names(X_train), get_one_plot)
# ggarrange(plotlist = a)



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
    tau.type = 'with weight'
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
       color='tau.type',numeric.x.axis = TRUE) +
  geom_line(mapping = aes(x = x, y = tau.true), color = "purple") +
  geom_line(mapping = aes(x = x, y = tau.pred.upper, color = tau.type), linetype='dashed') +
  geom_line(mapping = aes(x = x, y = tau.pred.lower, color = tau.type), linetype='dashed')


# plot(tau[1:n1], tau_train.grf_v1$predictions)
# abline(0, 1)


