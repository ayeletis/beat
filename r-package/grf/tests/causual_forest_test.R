
# install.packages("grf")
library(grf)

library(ggplot2)
library(data.table)
library(ggpubr)
library(gridExtra)
library(dplyr)

# Empty workspace
rm(list=ls())

sigmoid = function(x){
  return( 1/(1+exp(-x))   )
}
relu = function(x){
  return(ifelse(x>0, x, 0))
}


num_trees = 4000

n1 = 5000; #calibration
n2 = 1000; #validation
p_continuous = 4
p_discrete = 1
p_demog = 1
n = n1 + n2


set.seed(123456789)
# Covariates = X_cont + Binary + Demog (binary)
X_cont = matrix(rnorm(n*p_continuous), n, p_continuous)
X_disc = matrix(rbinom(n*p_discrete, 1, 0.3),n,p_discrete)
# Demographic variable correlated with X[,2]
# Z = rbinom(n, 1, 0.2) # This is to test whether grf_v2 gives same fit when Z is independent from X.
Z = as.matrix(rbinom(n, 1, sigmoid(X_cont[,2])))

## high dimensions

# Z =matrix(nrow=n, ncol=2)
# for(i in 1:2){
#   Z[, i] = rbinom(n, 1, sigmoid(sample(c( exp, relu, sin, cos, tanh), 1)[[1]](X_cont[,2])))
# }

#  another correlation
# X_Z = rnorm(n, mean=4, sd=10)
# Z = sigmoid(X_cont[,2] + exp(X_Z) + X_Z * X_cont[,2])
# Z = ifelse(Z>0.5, 1, 0)


print(c("validity_check",mean(X_cont[Z[,1]==0,2]),mean(X_cont[Z[,1]==1,2])))

# For now, we assume that Z is not a covariate in the model
X = cbind(X_cont,X_disc)
p = p_continuous + p_discrete


# Random assigment
W = rbinom(n, 1, 0.5)

# Simulate 'tau' -- TAU depends on X[2] but no on Z
tau_function <- function(x) {
  temp <- (-1 + pmax(x[,1], 0) + x[,2] + abs(x[,3]) + x[,5])
  return(temp)
}
tau = tau_function(X)

# Simulate 'Y'
y_function <- function(x,z) {
  temp <-  x[,1] - 2*x[,2] +  x[,4] + z #+ rowMeans( z)
  return(temp)
}
noise_y = runif(n)
Y =  y_function(X,Z) + tau*W + noise_y
Y_n = X[,1] - 2*X[,2] +  X[,4] + tau*W + noise_y


# For now, let's use only train data
train_data = data.table(Y=Y[c(1:n1)],Z=Z[c(1:n1), ],W=W[c(1:n1)],tau = tau[c(1:n1)],X=X[c(1:n1),])
x_cols = grep("X", names(train_data), value=TRUE)
X_train = train_data[, .SD, .SDcols=x_cols]
Y_train = train_data$Y
W_train = train_data$W
Zcols =grep('Z', names(train_data), value=TRUE)
Z_train = train_data[, .SD, .SDcols=Zcols]


# validate_data = data.frame(Y=Y[c(n1:n2)],Z=Z[c(n1:n2)],W=W[c(n1:n2)],tau = tau[c(n1:n2)],X=X[c(n1:n2),])
# X_validate = as.matrix(validate_data[,c(5:dim(validate_data)[2])])
# Y_validate = validate_data$Y
# W_validate = validate_data$W
# Z_validate = validate_data$Z


## Some validity plots
if(FALSE) {
data_to_plot = data.table(Y=Y,W=W,Z = Z[,1], tau)
data_to_plot[, W:=as.character(W)]
data_to_plot[, Z:=as.character(Z)]
p_ate = ggdensity(data=data_to_plot, x='Y', color='W', fill='W', alpha=0.2, add = "mean")
p_bias = ggdensity(data=data_to_plot, x='tau', color='Z', fill='Z', alpha=0.2, add = "mean")
grid.arrange(p_ate, p_bias, nrow = 1)

# ggsave("density_diff.png",p, width=12, height=6)
}


## ---------------------------------------------------
##   Estimate 'tau' using GRF
## ---------------------------------------------------

# propensity scores
# forest.W = regression_forest(X_train, W_train, tune.parameters = "all")
# W.hat = predict(forest.W)$predictions

# Y_hat, for GRF
# forest.Y = regression_forest(X_train, Y_train, tune.parameters = "all")
# Y.hat = predict(forest.Y)$predictions

# Estimate Causal forest (using 'grf' package)
start.time <- Sys.time()

fit_grf = causal_forest(X_train, Y_train, W_train,
                               honesty = TRUE,
                               # W.hat = W.hat, Y.hat = Y.hat,
                              # tune.parameters = "all",
                               num.trees = num_trees,
                        seed=1)
end.time <- Sys.time()
time_cf_estim <- end.time - start.time
print(time_cf_estim)

### Hengyu, here you would test your causal_forest_v2 function, which will also depend on Z_train
if (TRUE) {

# forest.W = regression_forest(X_train, W_train, target.weights=as.matrix(Z_train),
#                              tune.parameters = "all")
# W.hat = predict(forest.W)$predictions

# Y_hat, for GRF
# forest.Y = regression_forest(X_train, Y_train,target.weights=as.matrix(Z_train),
#                              tune.parameters = "all")
# Y.hat = predict(forest.Y)$predictions
start.time <- Sys.time()
fit_grf_v2 = causal_forest(X_train, Y_train, W_train,
                        target.weights = as.matrix(Z_train),
                        honesty = TRUE,
                        # W.hat = W.hat, Y.hat = Y.hat,
                       #tune.parameters = "all",
                        num.trees = num_trees,
                       seed = 1)
end.time <- Sys.time()
time_cf_estim <- end.time - start.time
print(time_cf_estim)

}


## ---------------------------------------------------
##   Explore GRF ability at recovering 'tau'
## ---------------------------------------------------

# GRF predictions
tau_train.grf = predict(fit_grf, estimate.variance = TRUE) # newdata=X_validate,
tau.pred = tau_train.grf$predictions
# plot(tau[c(1:n1)], tau.pred,main="in-sample"); abline(0,1)
#(Good fit)



# V2 predictions
if(TRUE) {
  tau_train.grf_v2 = predict(fit_grf_v2, estimate.variance = TRUE)  # newdata=X_validate,
  tau.pred_v2 = tau_train.grf_v2$predictions
  # plot(tau[c(1:n1)], tau.pred_v2,main="in-sample"); abline(0,1)
}
#(Hengyu, here the fit can be poor, especially when Z and X are heavily correlated.
# In fact, one way to test teh v2 is to make Z and X *not* correlated (like in line 50 above "Z = rbinom(n, 1, 0.5)"). In that case, your tau predictions should be as good as the regular GRF



## ---------------------------------------------------------
##   Explore the possible bias and whether V@ got rid of it
## ---------------------------------------------------------
rsq <- function (x, y) cor(x, y) ^ 2

get_ecdf_values = function(x){
  v = seq(0, 1, 0.01)
  return(list(v, ecdf(x)(v)))
}

data_to_plot =  data.table(Z= Z_train[, 1],
                          tau_grf_no_weight = tau_train.grf$predictions,
                          tau_grf_with_weight = tau_train.grf_v2$predictions) #
names(data_to_plot)[1] = 'Z'
data_to_plot[, Z:=as.character(Z)]

data_to_plot2 = melt(data_to_plot, id.vars = 'Z')

ggdensity(data=data_to_plot2, x='value',
          color='Z', fill='Z',
          alpha=0.2, add = "mean",
          facet.by = 'variable', ncol=1
          ,scales='free'
          )
# ggsave("result_with_z_50_correlation.png", width = 8, height=6)
dat = data.table(
  tau_true = tau,
  tau_pred_no_weights = tau_train.grf$predictions,
  tau_pred_weights= tau_train.grf_v2$predictions
)
dat[, idx:=1:.N]
dat = melt(dat, id.vars = 'idx')
ggline(data=dat, x='idx', y='value', color='variable')

ggecdf(data=dat, x='value', color='variable')
ggsave("prediction_ecdf.png", width = 8, height=6)
