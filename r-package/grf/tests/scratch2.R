
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

n1 <- 1000
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
X <- cbind(X_cont, X_disc)

# Demographic variable correlated with X[,2]
# Z = as.matrix(rbinom(n, 1, 0.2)) # This is to test whether grf_v2 gives same fit when Z is independent from X.
# Z = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,2]))))

## Higher dimension for demographics -- 1 of them related to X_cont[,2], two other binary and a continuous, all independent
Z <- sapply(1:1, function(x) as.matrix(rbinom(n, 1, 1 / (1 + exp( -X_cont[, 2])))))
Z <- cbind(Z, rbinom(n, 1, 0.2), rbinom(n, 1, 0.7),sigmoid(X[,2]*10) )

# print(c("validity_check for X,Z dependency",mean(X_cont[Z[,1]==0,2]),mean(X_cont[Z[,1]==1,2])))


# For now, we assume that Z is not a covariate in the model
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
train_data <- data.table(Y = Y[c(1:n1)], Z = Z[c(1:n1), ], W = W[c(1:n1)], tau = tau[c(1:n1)], X = X[c(1:n1), ])

x = X
z= Z
I = 1
z_hat = predict(multi_regression_forest(X, z))$prediction
z_org  = z-z_hat


dat1 = data.table(x)
# dat1[, z1 := z_org[, I]]
dat1[, z1 := z[, I]]
dat1[, bins := as.numeric(cut(V2, 256))]
A = as.matrix(dat1[, mean(z1), by=bins])
plot(A[,1], A[,2]-mean(A[,2]))
sum(( A[,2]-mean(A[,2]))^2)

construct_target_weight_mean = function(x, z, num_breaks = 256) {

  calculate_avg = function(x_col, z_col, num_breaks) {
    df = dat[, .SD, .SDcols = c(x_col, z_col)]
    df[, bins := cut(get(x_col), breaks = num_breaks)]
    df[, mean_value := mean(get(z_col)), by = bins]
    out = df[, mean_value]
    return(out)
  }

  stopifnot(is.matrix(x))
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
    out[[i]] = sapply(z_cols, calculate_avg, x_col = x_cols[i], num_breaks = num_breaks)
  }
  return(out)
}
a = construct_target_weight_mean(x, z_org)
plo
