library(beat)
rm(list=ls())
library(data.table)
library(ggpubr)

n = 4000
p = 10
num.trees <- 2000

sigma.tau = 0.1
sigma.noise = 1
sigma.m = 1

sigmoid = function(x){1/(1+exp(-x))}

set.seed(1)

# modify from generate_causal_data(dgp == "aw2")
X <- matrix(runif(n * p), n, p)
zeta1 <- 1 + 1 / (1 + exp(-20 * (X[, 1] - (1 / 3))))
zeta2 <- 1 + 1 / (1 + exp(-20 * (X[, 2] - (1 / 3))))
weights <- rbinom(n, 1, 1 / (1 + exp(-20 * (X[, 3] - (1 / 3))))) # new

tau <- zeta1 * zeta2 * weights # new
e <- rep(0.5, n)  # Randomized trial (no confounding)
W <- rbinom(n = n, size = 1, prob = e)
m <- e * tau
V <- 1
if (!is.na(sd(m)) & !(sd(m) == 0)) {
  m <- m / sd(m) * sigma.m
}
if (!is.na(sd(tau)) & !(sd(tau) == 0)) {
  tau <- tau / sd(tau) * sigma.tau
}

V <- V / mean(V) * sigma.noise^2
Y <- m + (W - e) * tau + sqrt(V) * rnorm(n)


# test
data = data.table(Y, Z = as.character(weights), W=as.character(W), tau)
p_ate = ggdensity(data=data, x='Y', color='W', fill='W', alpha=0.2, add = "mean")
p_bias = ggdensity(data=data, x='tau', color='Z', fill='Z', alpha=0.2, add = "mean")
ggarrange(p_ate, p_bias, nrow = 1)


forest.tau.noweight <- causal_forest(X, Y, W,
                            num.trees = num.trees, seed=1)

forest.tau <- causal_forest(X, Y, W,target.weights = as.matrix(weights),
                            num.trees = num.trees, seed=1)


mean((predict(forest.tau.noweight)$predictions - tau)^2)
mean((predict(forest.tau)$predictions - tau)^2)


data_to_plot =  data.table(Z= weights,
                           tau_grf_no_weight = predict(forest.tau.noweight)$predictions,
                           tau_grf_with_weight = predict(forest.tau)$predictions) #

data_to_plot[, Z:=as.character(Z)]

data_to_plot2 = melt(data_to_plot, id.vars = 'Z')

ggdensity(data=data_to_plot2, x='value',
          color='Z', fill='Z',
          alpha=0.2, add = "mean",
          facet.by = 'variable', ncol=1
          ,scales='free'
)
