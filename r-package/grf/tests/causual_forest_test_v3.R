rm(list = ls())
# install.packages("grf")
library(grf)
#
construct_target_weight_mean = function(x, z, num_breaks = 256) {
  calculate_avg = function(dat, x_col, z_col, num_breaks) {
    df = dat[, .SD, .SDcols = c(x_col, z_col)]
    df[, cut_bins := cut(get(x_col), breaks=num_breaks)]

    df[, mean_value := mean(get(z_col)), by = cut_bins]
    # out = df[, mean_value] # if z is orthogonal (z hat)
    out = df[, mean_value - mean(mean_value)] # if they are independent, then error will center around 0
    return(out)
  }

  stopifnot(is.matrix(z))
  stopifnot(dim(x)[1] == dim(z)[1])
  dat = as.data.table(x)
  x_cols = paste0('x_', 1:dim(x)[2])
  names(dat) = x_cols
  dat_z = as.data.table(z)
  z_cols = paste0("z_", 1:dim(z)[2])
  names(dat_z) = z_cols
  dat = cbind(dat, dat_z)
  # output matrix: 3D [var, n, z]
  out = vector('list', length = length(x_cols))
  for (i in 1:length(x_cols)) {
    out[[i]] = sapply(z_cols, calculate_avg, x_col = x_cols[i], num_breaks = num_breaks, dat=dat)
  }

  return(out)
}

library(ggplot2)
library(data.table)
library(ggpubr)
library(gridExtra)
library(dplyr)

# Empty workspace


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

# Covariates = X_cont + Binary + Demog (binary)
X_cont <- matrix(rnorm(n * p_continuous), n, p_continuous)
X_disc <- matrix(rbinom(n * p_discrete, 1, 0.3), n, p_discrete)

# Demographic variable correlated with X[,2]
# Z = as.matrix(rbinom(n, 1, 0.2)) # This is to test whether grf_v2 gives same fit when Z is independent from X.
# Z = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,2]))))

## Higher dimension for demographics -- 1 of them related to X_cont[,2], two other binary and a continuous, all independent
Z <- sapply(1:1, function(x) as.matrix(rbinom(n, 1, 1 / (1 + exp(-X_cont[, 2]*10)))))
Z <- cbind(Z, rbinom(n, 1, 0.2), rbinom(n, 1, 0.7), rnorm(n, 0, 1))

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


## Some validity plots to make sure that tau discriminates in Z
if (FALSE) {
  data_to_plot <- data.table(Y = Y, W = W, Z = Z[, 1], tau)
  data_to_plot[, W := as.character(W)]
  data_to_plot[, Z := as.character(Z)]
  p_ate <- ggdensity(data = data_to_plot, x = "Y", color = "W", fill = "W", alpha = 0.2, add = "mean")
  p_bias <- ggdensity(data = data_to_plot, x = "tau", color = "Z", fill = "Z", alpha = 0.2, add = "mean")
  grid.arrange(p_ate, p_bias, nrow = 1)
}


## ---------------------------------------------------
##   Estimate 'tau' using GRF
## ---------------------------------------------------

# target.weights = construct_target_weight_mean(X_train, as.matrix(Z_train), num_breaks = 256)
fit_grf <- causal_forest(X_train, Y_train, W_train,
                         honesty = TRUE,
                         #Y.hat = 0, W.hat=0,
                         target.weights = as.matrix(Z_train) ,
                         target.weight.penalty=0,
                         num.trees = num_trees,
                         # num.threads=1,
                         # tune.parameters='all'
                         seed=1
)


fit_grf_v2 <- causal_forest(X_train, Y_train, W_train,
                            target.weights = as.matrix(Z_train),
                            target.weight.penalty = 10,
                            target.weight.bins.breaks = 128,
                            honesty = TRUE,
                            num.trees = num_trees, #10, #num_trees,
                            # tune.parameters = "all",
                            # num.threads = 1,
                            seed=1
)



tau_train.grf <- predict(fit_grf, estimate.variance = TRUE) # newdata=X_validate,
tau.pred <- tau_train.grf$predictions

tau_train.grf_v2 <- predict(fit_grf_v2, estimate.variance = TRUE) # newdata=X_validate,
tau.pred_v2 <- tau_train.grf_v2$predictions

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


## ---------------------------------------------------------
##   Explore the possible bias and whether new GRF got rid of it
## ---------------------------------------------------------

get_ecdf_values <- function(x) {
  v <- seq(0, 1, 0.01)
  return(list(v, ecdf(x)(v)))
}

# ggsave("result_with_z_50_correlation.png", width = 8, height=6)

## Hengyu, in the plot above, I want to have the x-axis fixed (to compare values) and the y-axis free, how can i do that?

W.hat <- predict(forest.W)$predictions

# Y_hat, for GRF
forest.Y <- regression_forest(X_train, Y_train, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions

# Set spect for regular forest

all_tuning_param <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty")
my_tuning_param <- all_tuning_param[1:7]
# my_tuning_param = c("none")


# Estimate REGULAR Causal forest (using 'grf' package)
start.time <- Sys.time()


end.time <- Sys.time()
time_cf_estim <- end.time - start.time
print(time_cf_estim)


# Estimate Causal forest with TARGET WEIGHTS, using the tuneable parameters from the regular forest
start.time <- Sys.time()

end.time <- Sys.time()
time_cf_estim_v2 <- end.time - start.time
print(time_cf_estim_v2)


## ---------------------------------------------------
##   Explore GRF ability at recovering 'tau'
## ---------------------------------------------------

# Regular GRF predictions
plot(tau[c(1:n1)], tau.pred, main = "in-sample - not weigthed")
abline(0, 1)
# (Good fit)


# Target.Weighted predictions
plot(tau[c(1:n1)], tau.pred_v2, main = "in-sample - weigthed")
abline(0, 1)
# (Directionally OK)


# Weighted vs. not
plot(tau.pred, tau.pred_v2, main = "in-sample - weigthed vs. not.")
abline(0, 1)


## -------------------------------------------
## Generate Marginals (moderating effects)
## -------------------------------------------
generate_prediction_per_dim_cont <- function(x_dim) {
  x <- seq(-2, 2, length.out = 101)
  X.marginals <- matrix(0, 101, p)
  X.marginals[, x_dim] <- x
  X.marginals <- as.data.frame(X.marginals)

  tau.hat <- predict(model_fit, X.marginals, estimate.variance = TRUE)
  mean.hat <- tau.hat$predictions
  sigma.hat <- sqrt(tau.hat$variance.estimates)

  tau.true <- tau_function(X.marginals)
  upper <- mean.hat + 1.96 * sigma.hat
  lower <- mean.hat - 1.96 * sigma.hat
  out <- data.table(x = x, tau_pred = mean.hat, tau_true = tau.true, tau_upper = upper, tau_lower = lower, x_dim = x_dim)
  return(out)
}

generate_prediction_per_dim_discrete <- function(x_dim) {
  x <- c(0, 1)
  X.marginals <- matrix(0, 2, p)
  X.marginals[, x_dim] <- x

  tau.hat <- predict(model_fit, X.marginals, estimate.variance = TRUE)
  mean.hat <- tau.hat$predictions
  sigma.hat <- sqrt(tau.hat$variance.estimates)

  tau.true <- tau_function(X.marginals)
  upper <- mean.hat + 1.96 * sigma.hat
  lower <- mean.hat - 1.96 * sigma.hat
  out <- data.table(x = x, tau_pred = mean.hat, tau_true = tau.true, tau_upper = upper, tau_lower = lower, x_dim = x_dim)
  return(out)
}

### Hengyu, I am getting an error in these functions and i don't understand why.  Could you please check.


# generate tau marginal per *continuos* dimension using Regular GRF
model_fit <- fit_grf # we use regular GRF for now

dat_pred_cont <- rbindlist(lapply(1:p_continuous, generate_prediction_per_dim_cont))

p_marginal_cont <- ggline(
  data = dat_pred_cont, x = "x", y = "tau_pred", facet.by = "x_dim",
  plot_type = "l", xlab = "x", ylab = "tau", title = "Marginal Effects - Continuous var.",
  numeric.x.axis = TRUE, short.panel.labs = FALSE
) +
  geom_line(mapping = aes(x = x, y = tau_true), color = "red") +
  geom_line(mapping = aes(x = x, y = tau_upper), linetype = "dashed") +
  geom_line(mapping = aes(x = x, y = tau_lower), linetype = "dashed")


# generate tau marginal per *discrete* dimension using Regular GRF
dat_pred_disc <- rbindlist(lapply((p_continuous + 1):p, generate_prediction_per_dim_discrete))
dat_pred_disc[, x_dim := as.character(x_dim)]
dat_pred_disc[, x := as.character(x)]
dat_pred_disc

p_marginal_disc <- ggbarplot(
  data = dat_pred_disc, x = "x", y = "tau_pred", fill = "x_dim",
  facet.by = "x_dim",
  palette = "jco", short.panel.labs = FALSE,
  legend = "none"
) +
  geom_errorbar(mapping = aes(ymin = tau_lower, ymax = tau_upper), width = .2) +
  geom_point(mapping = aes(x = x, y = tau_true), color = "red") +
  labs(x = "X", y = "Prediction") +
  geom_text(
    data = data.table(x = 1, tau_pred = 3, x_dim = factor("5", levels = c("5", "6"))),
    label = "Red point: True value", color = "red"
  )


## ---------------------------------------------------------
##   Explore the possible bias and whether new GRF got rid of it
## ---------------------------------------------------------

get_ecdf_values <- function(x) {
  v <- seq(0, 1, 0.01)
  return(list(v, ecdf(x)(v)))
}

# ggsave("result_with_z_50_correlation.png", width = 8, height=6)

## Hengyu, in the plot above, I want to have the x-axis fixed (to compare values) and the y-axis free, how can i do that?
