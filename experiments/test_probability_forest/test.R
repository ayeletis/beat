library(grf)
library(ggpubr)
library(data.table)
p <- 5
n <- 2000
X1 = as.matrix(rbinom(n=n, size=1, prob=0.3))
X2 = matrix(rnorm(n*p), n, p)
X = cbind(X1, X2)

W1 = rbinom(n=n, size=1, prob=0.3)
W2 = runif(n=n, min=-10, max=10)
W = cbind(W1, W2)
prob <- 1 / (1 + exp(-X[, 1]*W1*3 + X[,2] - X[, 3] ))
Y <- as.factor(rbinom(n, 5, prob))

dat_plot = data.table(x=X[,2], Y=Y, w=W1)
dat_plot[, w:= as.character(w)]
dat_plot = dat_plot[,.N, by=.(Y, w)]
dat_plot[, rate:= N/sum(N), by=w]
ggbarplot(data=dat_plot,
          x="Y",y='rate',
          position = position_dodge(),
          fill='w')

# standard
p.forest <- probability_forest(X, Y, seed = 1)
p.hat <- predict(p.forest, X, estimate.variance = TRUE)
p.pred = apply(p.hat$predictions, 1, which.max)-1

dat_plot2 = data.table(y=p.pred, x=X[,1], w=W1)
dat_plot2 = dat_plot2[,.N, by=.(Y, w)]

dat_plot2[, rate:= N/sum(N), by=w]
dat_plot2[, w:= as.character(w)]
ggbarplot(data=dat_plot2,
          x="Y",y='rate',
          position = position_dodge(),
          fill='w')


p.balanced = balanced_probability_forest(X, Y,
                                         target.weights = W,
                                         target.weight.penalty = 100,
                                         target.weight.bins.breaks = 256,
                                         target.weight.standardize = TRUE,
                                         target.weight.penalty.metric = "split_l2_norm_rate"
)

p.hat.balanced <- predict(p.balanced, X, estimate.variance = TRUE)
p.pred.balanced = apply(p.hat.balanced$predictions, 1, which.max)-1

dat_plot3 = data.table(y=p.pred.balanced, x=X[,1], w=W1)
dat_plot3 = dat_plot3[,.N, by=.(Y, w)]

dat_plot3[, rate:= N/sum(N), by=w]
dat_plot3[, w:= as.character(w)]
ggbarplot(data=dat_plot3,
          x="Y",y='rate',
          position = position_dodge(),
          fill='w')


a = p.hat.balanced$predictions
b = p.hat$predictions
p.pred.balanced[1:3]
p.pred[1:3]
variable_importance(p.balanced)
variable_importance(p.forest)
