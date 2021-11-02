rm(list=ls())
library(beat)
library(data.table)
library(ggpubr)
r_2 = function(y_true, y_hat){
  mean_y_true = mean(y_true)
  top = (y_hat - mean_y_true)^2
  bottom = (y_true - mean_y_true)^2
  return(sum(top)/sum(bottom))
}
set.seed(1)
n <- 50
p <- 3
X <- matrix(rnorm(n * p), n, p)
Z = matrix(c(sample(x=c(0,1), size=n, replace=TRUE, prob=c(0.3, 0.7)),
             sample(x=c(0,1), size=n, replace=TRUE)),
           n, 2)
W = rbinom(n, 1, 0.5)

# Z = matrix(c(sample(x=seq(-100,100), size=n, replace=TRUE ),
#              sample(x=seq(-1000,1000), size=n, replace=TRUE )),
#            n, 2)
Y <- X[, 1] * rnorm(n) +  Z[, 1] * rnorm(n)  +  Z[, 2] * rnorm(n)

r.forest <- causal_forest(X, Y,W=W, num.trees = 2, target.weights=NULL) #


# regression_forest

y_hat = predict(r.forest)

out = data.table(y_true = Y, y_hat = y_hat$predictions, Z=Z)
out[, idx := 1:.N]
out = melt(out, measure.vars = c("y_true", 'y_hat'))
out[, Z.V1:= as.character(Z.V1)]
out[, Z.V2:= as.character(Z.V2)]

p1 = ggpubr::ggdensity(data=out, x='value', color='Z.V1',facet.by = 'variable',
                       title='With Weights',
                       subtitle = sprintf("R2: %.4f", r_2(Y,y_hat$predictions)))

p2 = ggpubr::ggdensity(data=out, x='value', color='Z.V2',facet.by = 'variable'  )
p3 = ggpubr::ggdensity(data=out, x='value', color='variable'  )
ggarrange(p1,p2,p3)
ggsave("test_weights.png", width = 14, height=10)



r.forest <- regression_forest(X, Y, num.trees = 1000, target.weights=NULL) #

y_hat = predict(r.forest)
out = data.table(y_true = Y, y_hat = y_hat$predictions, Z=Z)
out[, idx := 1:.N]
out = melt(out, measure.vars = c("y_true", 'y_hat'))
out[, Z.V1:= as.character(Z.V1)]
out[, Z.V2:= as.character(Z.V2)]

p1 = ggpubr::ggdensity(data=out, x='value', color='Z.V1',facet.by = 'variable',
                       title='Without Weights',
                       subtitle = sprintf("R2: %.4f", r_2(Y,y_hat$predictions)))
p2 = ggpubr::ggdensity(data=out, x='value', color='Z.V2',facet.by = 'variable'  )
p3 = ggpubr::ggdensity(data=out, x='value', color='variable'  )
ggarrange(p1,p2,p3 )
ggsave("test_no_weights.png", width = 14, height=10)
