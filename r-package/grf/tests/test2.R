library(grf)
library(data.table)
library(ggpubr)

r_2 = function(y_true, y_hat){
  mean_y_true = mean(y_true)
  top = (y_hat - mean_y_true)^2
  bottom = (y_true - mean_y_true)^2
  return(sum(top)/sum(bottom))
}

set.seed(1)
n <- 100
p <- 4
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5)

Z = matrix(c(sample(x=c(-20, 20), size=n, replace=TRUE, prob=c(0.3, 0.7)),
             sample(x=c(-10,10), size=n, replace=TRUE),
             sample(x=c(-100,100), size=n, replace=TRUE)),
           n, 3)
print(colMeans(Z, na.rm=TRUE))

Y = pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n) +  Z[, 1] * rnorm(n)  +  Z[, 2] * rnorm(n)


c.forest <- causal_forest(X, Y, W, num.trees =100, target.weights=NULL)
y_hat = predict(c.forest)


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



r.forest <- causal_forest(X, Y, W, num.trees = 100, target.weights=NULL) #

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
# ggsave("test_no_weights.png", width = 14, height=10)
