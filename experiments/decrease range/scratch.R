library(ggpubr)
library(data.table)
library(scales)

n = 10000
set.seed(1)
w = sample(c(0,1), size=n, replace = TRUE)
w_mean = mean(w)
y = runif(n, min=2, max=10) # rbeta(n, 1,4) #
y = rescale(y, to=c(0,1))

y_mean = mean(y)

beta = sum((w - w_mean)*(y -y_mean)) / sum((w - w_mean)^2)
response = (w  - w_mean)*((y - y_mean) - beta*(w - w_mean))
response = response^2


cal_decrease = function(i){
  out = sum(response[1:i])^2/i + sum(response[-c(1:i)])^2/(n-i)
  return(out)
}

split_point = 2:(n-1)

decrease = sapply(split_point, cal_decrease)
dat = data.table(split = split_point,
                 dec = decrease)

ggline(dat, x="split", y="dec",
       plot_type = 'l',
       numeric.x.axis = TRUE,
       xlab="split point",
       ylab="sum decreaes",
       title='Y is randomlyi sample in [2,10], set to 0 when w=0, and rescaled to [0,1]'
       )

ggsave('random_sample.png', width=8,height=6)
