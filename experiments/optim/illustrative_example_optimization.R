library(data.table)

rm(list=ls())

# Working directory
dir <- "/Applications/Dropbox (Harvard University)/research/heterog_treat_effects/SHARED_BEAT_unintended_bias"
setwd(dir)
dir_code <- paste(getwd(),"/code/",sep="")
dir_data <- paste(getwd(),"/data/",sep="")
dir_output <- paste(getwd(),"/output/",sep="")
dir_plots <- paste(getwd(),"/plots",sep="")
setwd(dir_code)


# ================================================================ 
#            Get data
# ----------------------------------------------------------------

setwd(dir_output)
# save(list = c('data'), file = "data_illustrative_example.RData") 
load("data_illustrative_example.RData")
train_data = data$train_data
test_data = data$test_data

Xcols = grep("X", names(train_data), value=TRUE)
Zcols =grep('Z', names(train_data), value=TRUE)

## train
X_train = train_data[, .SD, .SDcols=Xcols]
Y_train = train_data$Y
W_train = train_data$W
Z_train = train_data[, .SD, .SDcols=Zcols]
tau_train = train_data$tau
colMeans(Z_train)
cor(X_train,Z_train)
#mean(tau_train)
#var(tau_train)

## test
X_test = test_data[, .SD, .SDcols=Xcols]
Y_test = test_data$Y
W_test = test_data$W
Z_test = test_data[, .SD, .SDcols=Zcols]
tau_test = test_data$tau

## with demographic variables
X_train_demog = cbind(X_train,Z_train)
X_test_demog = cbind(X_test,Z_test)


# ================================================================ 
#            Get taus
# ----------------------------------------------------------------

load("estimated_models_illustrative_example.RData")
## Predict TAUs
set.seed(my_seed) # we fix seed for the random prediction
{
  ## Get score predictions
  tau.pred = list(
    'grf_demog' = list(
      'train' = predict(fit_grf_demog)$predictions,
      'test' = predict(fit_grf_demog, X_test_demog)$predictions
    ),
    'grf' = list(
      'train' = predict(fit_grf)$predictions,
      'test' = predict(fit_grf, X_test)$predictions
    ),
    'grf_fair' =  list(
      'train' = predict(fit_fair)$predictions,
      'test' = predict(fit_fair, X_test)$predictions
    ),
    'random' =  list(
      'train' = runif(length(Y_train),min(predict(fit_grf)$predictions),max(predict(fit_grf)$predictions)),
      'test' = runif(length(Y_test),min(predict(fit_grf)$predictions),max(predict(fit_grf)$predictions))
    ),
    'grf_residual_rf' = list(
      'train' = predict(fit_rfr, X_rfr_train)$predictions,
      'test' = predict(fit_rfr, X_rfr_test)$predictions
    )
  )
} # end of predict taus


# ================================================================ 
#            Set allocation
# ----------------------------------------------------------------

#setwd(dir_output)
#save(list = c('tau.pred','Z_test'), file = "cache_for_optimization.RData") 
load("cache_for_optimization.RData")

setwd(dir_code)
source("optimization_select_targets.R")

tau_vector = tau.pred$grf_demog$test
demographic_data = Z_test
min_tau = 0


# optimization -- all data 
selected_targets = optimization_search(tau_vector=tau.pred$grf_demog$test,
                              demographic_data=Z_test ,
                              min_tau=min_tau,
                              max_target_rate=0.2)

#

is_targeted = selected_targets


## Summarize results

output = data.table(id = c(1:length(tau_vector)),
                    tau = tau_vector, 
                    is_targeted = is_targeted,
                    Z = demographic_data)

colnames(output)[4:7] = c('z1','z2','z3','z4')
output[, lapply(.SD, mean), by=is_targeted][order(is_targeted)]

output[, sum(is_targeted)]

output[, tau:= tau.pred$grf_demog$test]

# origin 
dat_origin = copy(Z_test)
dat_origin[, tau := tau.pred$grf_demog$test]
dat_origin[, is_targeted :=as.integer(tau >=min_tau)]
dat_origin[, lapply(.SD, mean), by=is_targeted][order(is_targeted)]
dat_origin[, sum(is_targeted)]

z_col_names = grep("Z", names(dat_origin), value=TRUE)
z_origin = dat_origin[is_targeted==1,  lapply(.SD, mean), .SDcols=z_col_names] - dat_origin[selected==0,  lapply(.SD, mean), .SDcols=z_col_names]
dat_origin[is_targeted==1, sum(tau)] - sqrt(sum(z_origin^2))

z_col_names = grep("z", names(output), value=TRUE)
z_out = output[is_targeted ==1,  lapply(.SD, mean), .SDcols=z_col_names] - output[is_targeted ==0,  lapply(.SD, mean), .SDcols=z_col_names]
output[is_targeted==1,sum(tau)] - sqrt(sum(z_out^2))
