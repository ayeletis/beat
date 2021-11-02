library(ggpubr)
library(data.table)


construct_target_weight_mean = function(x, z, num_breaks = 256, demean=TRUE) {
  calculate_avg = function(dat, x_col, z_col, num_breaks) {
    df = dat[, .SD, .SDcols = c(x_col, z_col)]
    df[, cut_bins := cut(get(x_col), breaks=num_breaks)]
    if(isTRUE(demean)){
      df[, (z_col) := get(z_col) - mean(get(z_col))]
    }
    df[, mean_value := mean(get(z_col)), by = cut_bins]
    out = df[, mean_value]
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


sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}
relu <- function(x) {
  return(ifelse(x > 0, x, 0))
}

x_dim = 4
n = 1000
set.seed(1)
X = matrix(rnorm(n * x_dim), n, x_dim)
colnames(X) = paste0("X", 1:x_dim)
Z1 = rbinom(n, 1, sigmoid(X[, 1]*10))
Z2 = sigmoid(relu(X[, 2]))
Z3 = sigmoid(X[, 3] + rgamma(n, 3) )
Z4 = rbeta(n , 2, 4)
dat = data.table(X, Z1,Z2, Z3, Z4)

setwd('./plots/')

for(i in 1:4){
  x_var = paste0("X", i)
  dat_x = melt(dat[, lapply(.SD, mean), by=c(x_var),
                   .SDcols=paste0("Z", 1:4)], id.vars = x_var)
  dat_x[, corre := cor(get(x_var), value), by=variable]
  dat_x[, variable := sprintf("%s, Corr(%s, %s): %.2f", variable, variable, x_var,  corre )]
  p = ggscatter(data=dat_x, x_var, y='value',
            plot_type = 'p',add='loess',
            add.params = list(color='red', linetype='dashed'),
            facet.by = 'variable',
            numeric.x.axis = TRUE,
            title=paste0('Correlation of Z and X', i) )
  ggsave(sprintf("correlation_z_x_%s.png", i),p, width=8, height=6)
}

Z = data.table(Z1, Z2, Z3,Z4)
Z = as.matrix(Z)
Z_new = construct_target_weight_mean(X, Z, 256, demean = TRUE )

for(i in 1:4){
  x_var = paste0("X", i)
  z_i = Z_new[[i]]
  colnames(z_i) = paste0('Z', 1:4)
  dat_new = data.table(x = dat[, get(x_var)], z_i)
  dat_new = melt(dat_new, id.vars='x')
  dat_new[, corre := cor(x, value), by=variable]
  dat_new[, variable := sprintf("%s, Corr(%s, %s): %.2f", variable, variable, x_var,  corre )]
  p = ggscatter(data=dat_new, x='x', y='value',
                plot_type = 'p',add='loess',
                add.params = list(color='red', linetype='dashed'),
                facet.by = 'variable',
                numeric.x.axis = TRUE,
                title=paste0('Correlation of Z and X', i),
                subtitle='Transformed' )
  ggsave(sprintf("correlation_z_new_x_%s.png", i),p, width=8, height=6)
}
