
# install.packages("grf")
library(grf)

library(ggplot2)
library(data.table)
library(ggpubr)
library(gridExtra)
library(dplyr)

# Empty workspace
rm(list=ls())


# Working directory
dir <- "/Applications/Dropbox/research/heterog_treat_effects/"
setwd(dir)
dir_code <- paste(getwd(),"/code/",sep="")
dir_data <- paste(getwd(),"/data/",sep="")
dir_output <- paste(getwd(),"/output/",sep="")
dir_plots <- paste(getwd(),"/plots",sep="")

setwd(dir_code)

num_trees = 4000
# install.packages('grf')
## -------------------------------------------
## Simulate data.
#       X has continuous and discrete variables
#       Z is a demographic (dummy) correlated with X
#       tau depends on  X, but not on Z
#       W is random assigment
#       Y depends on X, Z, and tau*W
## -------------------------------------------

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
#Z = rbinom(n, 1, 0.5) # This is to test whether grf_v2 gives same fit when Z is independent from X.
Z = rbinom(n, 1, 1/(1+exp(-X_cont[,2])))
print(c("validity_check",mean(X_cont[Z==0,2]),mean(X_cont[Z==1,2])))

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
  temp <- x[,1] - 2*x[,2] + x[,4] + z
  return(temp)
}
noise_y = runif(n)
Y =  y_function(X,Z) + tau*W + noise_y


# For now, let's use only train data
train_data = data.frame(Y=Y[c(1:n1)],Z=Z[c(1:n1)],W=W[c(1:n1)],tau = tau[c(1:n1)],X=X[c(1:n1),])
X_train = as.matrix(train_data[,c(5:dim(train_data)[2])])
Y_train = train_data$Y
W_train = train_data$W
Z_train = train_data$Z #Hengyu, I started with Z being dim=1, but we can start with higher dimensionality

## Some validity plots
if(FALSE) {
  data_to_plot = data.table(Y,W,Z, tau)
data_to_plot[, W:=as.character(W)]
data_to_plot[, Z:=as.character(Z)]

p_ate = ggdensity(data=data_to_plot, x='Y', color='W', fill='W', alpha=0.2, add = "mean")
p_bias = ggdensity(data=data_to_plot, x='tau', color='Z', fill='Z', alpha=0.2, add = "mean")
p_test = ggdensity(data=data_to_plot, x='Y', color='W', fill='W', alpha=0.2, facet.by = 'Z',
          short.panel.labs = FALSE, add = "mean")

grid.arrange(p_ate, p_bias, nrow = 1)
}


## ---------------------------------------------------
##   Estimate 'tau' using GRF
## ---------------------------------------------------

# propensity scores
forest.W = regression_forest(X_train, W_train, tune.parameters = "all")
W.hat = predict(forest.W)$predictions

# Y_hat, for GRF
forest.Y = regression_forest(X_train, Y_train, tune.parameters = "all")
Y.hat = predict(forest.Y)$predictions

# Estimate Causal forest (using 'grf' package)
start.time <- Sys.time()
fit_grf = causal_forest(X_train, Y_train, W_train,
                               honesty = TRUE,
                               W.hat = W.hat, Y.hat = Y.hat,
                               tune.parameters = "all",
                               num.trees = num_trees)
end.time <- Sys.time()
time_cf_estim <- end.time - start.time
print(time_cf_estim)


### Hengyu, here you would test your causal_forest_v2 function, which will also depend on Z_train
if (TRUE) {

forest.W = regression_forest(X_train, W_train, target.weights=as.matrix(Z_train),
                             tune.parameters = "all")
W.hat = predict(forest.W)$predictions

# Y_hat, for GRF
forest.Y = regression_forest(X_train, Y_train,target.weights=as.matrix(Z_train),
                             tune.parameters = "all")
Y.hat = predict(forest.Y)$predictions

start.time <- Sys.time()
fit_grf_v2 = causal_forest(X_train, Y_train, W_train,
                                         target.weights = as.matrix(Z_train),
                        honesty = TRUE,
                        W.hat = W.hat, Y.hat = Y.hat,
                        tune.parameters = "all",
                        num.trees = num_trees)
end.time <- Sys.time()
time_cf_estim <- end.time - start.time
print(time_cf_estim)
}


## ---------------------------------------------------
##   Explore GRF ability at recovering 'tau'
## ---------------------------------------------------

# GRF predictions
tau_train.grf = predict(fit_grf, estimate.variance = TRUE)
tau.pred = tau_train.grf$predictions
# plot(tau[c(1:n1)], tau.pred,main="in-sample"); abline(0,1)
#(Good fit)



# V2 predictions
if(TRUE) {
  tau_train.grf_v2 = predict(fit_grf_v2, estimate.variance = TRUE)
  tau.pred_v2 = tau_train.grf_v2$predictions
  # plot(tau[c(1:n1)], tau.pred_v2,main="in-sample"); abline(0,1)
}
#(Hengyu, here the fit can be poor, especially when Z and X are heavily correlated.
# In fact, one way to test teh v2 is to make Z and X *not* correlated (like in line 50 above "Z = rbinom(n, 1, 0.5)"). In that case, your tau predictions should be as good as the regular GRF



## ---------------------------------------------------------
##   Explore the possible bias and whether V@ got rid of it
## ---------------------------------------------------------


data_to_plot = data.table(Z = Z_train, tau_grf = tau_train.grf$predictions, tau_grf_v2 = tau_train.grf_v2$predictions) #
data_to_plot[, Z:=as.character(Z)]

data_to_plot = melt(data_to_plot, id.vars = 'Z')

ggdensity(data=data_to_plot, x='value', color='Z', fill='Z', alpha=0.2, add = "mean",
          facet.by = 'variable', ncol=1)

p_grf = ggdensity(data=data_to_plot, x='tau_grf', color='Z', fill='Z', alpha=0.2, add = "mean")
p_grf_v2 = ggdensity(data=data_to_plot, x='tau_grf_v2', color='Z', fill='Z', alpha=0.2, add = "mean")


grid.arrange(p_grf,p_grf_v2, nrow = 1)
## HENGYU, here the second plot should not have two separated densities. They should be almost total overlap.
