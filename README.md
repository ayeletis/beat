# BEAT (Bias-Eliminating Adapted Trees)
Forked from [https://github.com/grf-labs/grf](https://github.com/grf-labs/grf)

## Installation

This package is not available on CRAN and must be built from source. On Windows you will
need RTools from https://cran.r-project.org/bin/windows/Rtools/ and on OS X you will 
need the developer tools documented in https://mac.r-project.org/tools/.

First install the `devtools` package, then use it to install `beat`:

```R
install.packages(c("devtools", "ggplot2")) ## ggplot2 only needed for example
devtools::install_version("RcppEigen", "0.3.3.7.0") ## beat does not work with newer RcppEigen
devtools::install_github("ayeletis/beat")  ## do not update the RcppEigem package if prompted
```
## Changes relative to GRF
The package offers three new functions: balanced_causal_forest, balanced_regression_forest, 
balanced_probability_forest. 

All arguments are the same as the original package, but there are two new inputs: 
target.weight.penalty indicates the penalty assigned to the protected attributes.
target.weights is a matrix that includes the protected characteristics. X should not 
inlcude the protected characteristics.

See full details about the BEAT method in: XXX link.

## Sample Usage

```R
library(beat)
library(ggplot2)


## ----------------------------------------------
##   Simulate some data
## ----------------------------------------------

n1 = 1000; #calibration 
n2 = 1000; #validation
p_continuous = 4  # number of continuous features (unprotected)
p_discrete = 3  # number of discrete features (unprotected)
p_demog = 1 # number of protected attributes
n = n1 + n2

                                        # Features (unprotected)
X_cont = matrix(rnorm(n*p_continuous), n, p_continuous)
X_disc = matrix(rbinom(n*p_discrete, 1, 0.3),n,p_discrete)
X = cbind(X_cont,X_disc)

                                        # Protected attributes, discrete and continuous, where the first one is correlated with X[,2]
Z = rbinom(n, 1, 1/(1+exp(-X_cont[,2])))

                                        # Tau -- in this example in depends on X[2] but no on Z
tau <- (-1 + pmax(X[,1], 0) + X[,2] + abs(X[,3]) + X[,5]) 

                                        # Random assignment
W = rbinom(n, 1, 0.5)

                                        # Output for causal forest
Y =  X[,1] - 2*X[,2] + X[,4] + 3*Z + tau*W + runif(n)   # Y is function of X, Z(demo), tau*W

train_data = data.frame(Y=Y[c(1:n1)], Z=Z[c(1:n1)], W=W[c(1:n1)], X=X[c(1:n1),], tau = tau[c(1:n1)])
test_data = data.frame(Y=Y[c((n1+1):(n1+n2))], Z=Z[c((n1+1):(n1+n2))], W=W[c((n1+1):(n1+n2))], X=X[c((n1+1):(n1+n2)),], tau = tau[c((n1+1):(n1+n2))])

Xcols = grep("X", names(train_data), value=TRUE)
Zcols =grep('Z', names(train_data), value=TRUE)


## train
X_train = train_data[,c(4:10)]
Y_train = train_data$Y
W_train = train_data$W
Z_train = train_data[,2]

## test
X_test = test_data[,c(4:10)]
Z_test = test_data$Z

## model specs
num_trees = 2000
my_penalty = 10 # When penalty = 0 it corresponds to GRF

## ----------------------------------------------
##   Estimate Balanced Causal Forest 
## ----------------------------------------------
fit_causal_beat <- balanced_causal_forest(X_train, Y_train, W_train,
                                          target.weights = as.matrix(Z_train),
                                          target.weight.penalty = my_penalty,
                                          num.trees = num_trees)


## Predict CBT causal scores
cbt_causal_train = predict(fit_causal_beat)$predictions
cbt_causal_test = predict(fit_causal_beat, X_test)$predictions


## ----------------------------------------------
##   Estimate Balanced Regression Forest 
## ----------------------------------------------
fit_regression_beat <- balanced_regression_forest(X_train, Y_train,
                                                  target.weights = as.matrix(Z_train),
                                                  target.weight.penalty = my_penalty,
                                                  num.trees = num_trees)

## Predict CBT regression scores
cbt_regression_train = predict(fit_regression_beat)$predictions
cbt_regression_test = predict(fit_regression_beat, X_test)$predictions


## ----------------------------------------------
##   Check balance in test scores
## ----------------------------------------------
dat.plot = data.frame(cbt_causal = cbt_causal_test,
                      cbt_regr = cbt_regression_test,
                      Z = as.factor(Z_test))

ggdensity(data=dat.plot,
          x='cbt_causal', color='Z', fill='Z', alpha=0.2, add = "mean")


ggdensity(data=dat.plot,
          x='cbt_regr', color='Z', fill='Z', alpha=0.2, add = "mean")



 
