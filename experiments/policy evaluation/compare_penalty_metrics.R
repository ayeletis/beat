  library(grf)
  library(data.table)
  library(ggpubr)
  library(ggnewscale)

  library(fst)
  library(directlabels)
  library(ggforce)


  source("functions_off_policy_evaluation.R")
  ## Functions for data generation for the MOTIVATION example
  {
    # TAU, which depends on X and on Z, normalized to ~ N(0,1) when which_tau>10
    tau_function <- function(x,z) {
      aux_baseline <- (pmax(x[,1], 0) - 0.5*x[,3])
      aux_bias <- z[,1]
      # aux_bias <- z[,1] + x[,2]
      aux = aux_baseline + aux_bias
      temp <- (aux - mean(aux))/sd(aux)
      return(temp)
    }

    # Z = Demographics, which depend on the correlation between Z and X
    z_function <- function(n,X_cont) {
      Z = as.matrix(rbinom(n, 1, 1/(1+exp(-6*X_cont[,2]))))
      Z = cbind(Z,rbinom(n, 1, 0.2),rbinom(n, 1, 0.7),rnorm(n,0,1))
      return(Z)
    }

    # DATA generation process -- uses the dimensions of the data, the the "version" which is the dependence between X and Z
    generate_data = function(n1, n2,
                             p_continuous, p_discrete,
                             sigma_y,seed=my_seed){
      set.seed(seed)
      n = n1 + n2

      X_cont = matrix(rnorm(n*p_continuous), n, p_continuous)
      X_disc = matrix(rbinom(n*p_discrete, 1, 0.3),n,p_discrete)
      X = cbind(X_cont,X_disc)
      p = p_continuous + p_discrete

      # Protected vars
      Z = z_function(n,X_cont)

      # Random assignment
      W = rbinom(n, 1, 0.5)

      # Tau = Treatment effect
      tau = tau_function(X,Z)

      # Y = Outcome (continuous)
      noise_y = sigma_y*runif(n)
      Y =  tau*W + noise_y

      train_data = data.table(Y=Y[c(1:n1)],
                              Z=Z[c(1:n1),],
                              W=W[c(1:n1)],
                              tau = tau[c(1:n1)],
                              X=X[c(1:n1),])

      test_data = data.table(Y=Y[c((n1+1):(n1+n2))],
                             Z=Z[c((n1+1):(n1+n2)), ],
                             W=W[c((n1+1):(n1+n2))],
                             tau = tau[c((n1+1):(n1+n2))],
                             X=X[c((n1+1):(n1+n2)),])

      return(list(train_data=train_data,
                  test_data=test_data))
    }

  } # end of functions

  ## Model specs
  num_trees = 2000
  my_penalty = 4
  my_seed = 1
  my_tunable_params = c("sample.fraction","mtry","min.node.size",
                        "honesty.fraction","honesty.prune.leaves",
                        "alpha","imbalance.penalty")
  # THis list has the tunable params from regular GRF.
  # We could tune "target.weight.penalty" too, using [tune.parameters='all'], but chose not too because to play with it

  ## Data specs
  n1 = 10000  # training
  n2 = 5000; # testing
  p_continuous = 4
  p_discrete = 1
  sigma_y = 1

# ================================================================
#            Simulate DATA
# ----------------------------------------------------------------
data = generate_data(n1=n1, n2=n2, p_continuous = p_continuous,
                     p_discrete = p_discrete, sigma_y=sigma_y,
                     seed=my_seed)

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
#            Estimate TAUs
# ----------------------------------------------------------------

a = Sys.time()
out <- balanced_causal_forest(X_train, Y_train, W_train,
                              target.weights = as.matrix(Z_train ) ,
                              target.weight.penalty = 1,
                              target.weight.standardize = TRUE,
                              target.weight.bins.breaks = 256,
                              #target.weight.penalty.metric = 'cosine_similarity_rate',
                              honesty = TRUE,
                              # num.threads = 1,
                              num.trees = num_trees,
                              tune.parameters=my_tunable_params,
                              seed=my_seed)
print(Sys.time()-a)
# out = balanced_causal_forest(X_train, Y_train, W_train,
#                                 target.weights = as.matrix(Z_train ) ,
#                                 target.weight.penalty = 2,
#                                 target.weight.standardize = TRUE,
#                                 target.weight.bins.breaks = 256,
#                                 target.weight.penalty.metric = 'cosine_similarity_rate',
#                                 honesty = TRUE,
#                                 # num.threads = 1,
#                                 num.trees = num_trees,
#                                 tune.parameters=my_tunable_params,
#                                 seed=my_seed)
# out_pred = predict(out)



# library(microbenchmark)
# microbenchmark(
#   out <- balanced_causal_forest(X_train, Y_train, W_train,
#                                 target.weights = as.matrix(Z_train ) ,
#                                 target.weight.penalty = 1,
#                                 target.weight.standardize = TRUE,
#                                 target.weight.bins.breaks = 256,
#                                 target.weight.penalty.metric = 'cosine_similarity_rate',
#                                 honesty = TRUE,
#                                 # num.threads = 1,
#                                 num.trees = num_trees,
#                                 tune.parameters=my_tunable_params,
#                                 seed=my_seed),
#   times = 10
# )


available_metrics = c("split_l2_norm_rate", # left, right: l2 norm(colmean target weight)* penalty rate * node decrease
                      "euclidean_distance_rate", # (left+right decrease) *  Euclidean distance (column mean target weight left, right ) * penalty rate
                      "cosine_similarity_rate", # (left+right decrease) *  (1-cos_sim(column mean target weight left, right )) * penalty rate

                      "split_l2_norm", #  sum(left,right l2 norm(colmean target weight))* penalty rate
                      "euclidean_distance", #  Euclidean distance (column mean target weight left, right ) * penalty rate
                      "cosine_similarity" #  (1-cos_sim(column mean target weight left, right )) * penalty rate
)

# pre build the array to save time
# target.avg.weights = construct_target_weight_mean(x = X_train,
#                                                   z = as.matrix(Z_train),
#                                                   num_breaks = 256)


# Function that selects targets based on tau_pred
target_tau_proportion  = function(tau_pred, target_rate){
  return(as.integer(tau_pred > quantile(tau_pred, c(1-target_rate))))
}

# Function that runs BEAT for different levels of penalty -- returns TAUs
estimate_taus_per_penalty = function(penalty_value, metric){
# this following copied from codes above and recycle the same parameters
  print(sprintf("Penalty rate: %s", penalty_value))
  fit_fair <- balanced_causal_forest(X_train, Y_train, W_train,
                                     target.weights = as.matrix(Z_train),
                                     target.weight.penalty = penalty_value,
                                     target.weight.standardize = TRUE,
                                     target.weight.bins.breaks = 256,
                                     target.weight.penalty.metric = metric,
                                     #target.avg.weights =  target.avg.weights,
                                     honesty = TRUE,
                                     num.trees = num_trees,
                                     tune.parameters=my_tunable_params,
                                     seed=my_seed)

  tau_train = predict(fit_fair)$predictions
  tau_test = predict(fit_fair, X_test)$predictions

  out_train = data.table(row_id = 1: dim(X_train)[1],
                         tau_predict = tau_train,
                         penalty_rate =penalty_value,
                         data_type = "train")
  out_test = data.table(row_id = (dim(X_train)[1]+1): (dim(X_train)[1]+dim(X_test)[1]),
                         tau_predict = tau_test,
                         penalty_rate =penalty_value,
                         data_type = "test")
  out = rbind(out_train, out_test)
  out[, penalty_metric := metric]
  return(out)

}

## Set range for penalty
penalty_vect =c( seq(0.1, 2,0.2), 2.5, 3)

available_metrics = c("split_l2_norm_rate", # left, right: l2 norm(colmean target weight)* penalty rate * node decrease
                      "euclidean_distance_rate", # (left+right decrease) *  Euclidean distance (column mean target weight left, right ) * penalty rate
                      "cosine_similarity_rate", # (left+right decrease) *  (1-cos_sim(column mean target weight left, right )) * penalty rate
                      "split_l2_norm", #  sum(left,right l2 norm(colmean target weight))* penalty rate
                      "euclidean_distance", #  Euclidean distance (column mean target weight left, right ) * penalty rate
                      "cosine_similarity" #  (1-cos_sim(column mean target weight left, right )) * penalty rate
)

taus_per_penalty = list()
for(i in available_metrics){
  a = Sys.time()
  print(i)
  out = rbindlist(lapply(penalty_vect, estimate_taus_per_penalty, metric=i))
  taus_per_penalty[[i]] = out
  print(Sys.time()-a)
}


taus_per_penalty = rbindlist(taus_per_penalty)

write_fst(taus_per_penalty, "cache_tau_per_penalty_metric.fst")
taus_per_penalty = read_fst('cache_tau_per_penalty.fst', as.data.table = TRUE)

# ================================================================
#            Benchmark
# ----------------------------------------------------------------

# demographic tau
fit_grf_demog <- causal_forest(X = X_train_demog, Y = Y_train, W = W_train,
                               honesty = TRUE,
                               num.trees = num_trees,
                               tune.parameters='all',
                               seed=my_seed)


# regular tau -- no demographics

fit_grf <- causal_forest(X = X_train, Y = Y_train, W = W_train,
                         honesty = TRUE,
                         num.trees = num_trees,
                         tune.parameters='all',
                         seed=my_seed)
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
    'random' =  list(
      'train' = runif(length(Y_train),min(predict(fit_grf)$predictions),max(predict(fit_grf)$predictions)),
      'test' = runif(length(Y_test),min(predict(fit_grf)$predictions),max(predict(fit_grf)$predictions))
    )
  )
} # end of predict taus


# ================================================================
#   Collect policy results per target_rate levels -- both for the benchmarks and for BEAT penalty
# ----------------------------------------------------------------
W.hat_train = mean(W_train) # We use prop(treatment) because the data comes from an experiment
W.hat_test = mean(W_test)
# Function that collects policy results per target vector using BEAT by penalty
collect_penalty_results = function(target_vector, data_type="train"){
  if(data_type=="train"){
    stopifnot(length(target_vector)==dim(X_train)[1])
    out = off_policy_evalution_v3(X = X_train,
                                  W = W_train,
                                  Y = Y_train,
                                  Z = Z_train,
                                  policy_targets = target_vector,
                                  policy_name = "balanced_grf" ,
                                  data_type = 'train',
                                  W.hat = W.hat_train)
    out = rbind(unlist(out))
    out = as.data.table(out)
  }else if(data_type=="test"){
    stopifnot(length(target_vector)==dim(X_test)[1])
    out = off_policy_evalution_v3(X = X_test,
                                  W = W_test,
                                  Y = Y_test,
                                  Z = Z_test,
                                  policy_targets = target_vector,
                                  policy_name = "balanced_grf" ,
                                  data_type = 'test',
                                  W.hat = W.hat_test)
    out = rbind(unlist(out))
    out = as.data.table(out)
  }else{
    stop("data_type is either train or test")
  }
  return(out)
}

# Function that collects policy results per target vector for BENCHMARKS
collect_benchmark_results = function(target_rate){

  evaluate_policy_results = function(policy_list, policy_name){
    policy = policy_list[[policy_name]]
    out_train = off_policy_evalution_v3(X_train, W_train, Y_train, Z_train,policy$train, policy_name, 'train', W.hat_train)
    out_test = off_policy_evalution_v3(X_test, W_test, Y_test, Z_test, policy$test, policy_name, 'test', W.hat_test)
    return(data.table(rbind(unlist(out_train),unlist(out_test) )))
  }

  my_targeting_policies = list(
    'grf_demog' = list(
      'train' = as.integer(tau.pred$grf_demog$train > quantile(tau.pred$grf_demog$train, c(1-target_rate))),
      'test' = as.integer(tau.pred$grf_demog$test > quantile(tau.pred$grf_demog$test, c(1-target_rate)))
    ),
    'grf' = list(
      'train' = as.integer(tau.pred$grf$train > quantile(tau.pred$grf$train, c(1-target_rate))),
      'test' = as.integer(tau.pred$grf$test > quantile(tau.pred$grf$test, c(1-target_rate)))
    ),
    'random' = list(
      'train' = as.integer(tau.pred$random$train > quantile(tau.pred$random$train, c(1-target_rate))),
      'test' = as.integer(tau.pred$random$test > quantile(tau.pred$random$test, c(1-target_rate)))
    )
  )

  results = rbindlist(lapply(names(my_targeting_policies),
                             evaluate_policy_results,
                             policy_list= my_targeting_policies))

}


dat_result_list = list()
count = 1
for(target_rate in c(0.2,0.4,0.5,0.6,0.8)){

  ## 1.- Collect results of "per_penalty", by different target rates
  if ('is_targeted' %in%  names(taus_per_penalty)){
    taus_per_penalty[, is_targeted:= NULL]
  }

  # Add "is_targeted" column per penalty rate
  taus_per_penalty[, is_targeted := target_tau_proportion(tau_predict,
                                                     target_rate=target_rate),
              by=.(penalty_rate, data_type)]

  # Collect results per data type
  results_penalty_train = taus_per_penalty[data_type=="train",
                                  collect_penalty_results(is_targeted,"train"),
                            by=.(penalty_rate, penalty_metric )]
  results_penalty_test = taus_per_penalty[data_type=="test",
                                  collect_penalty_results(is_targeted,"test"),
                                  by=.(penalty_rate, penalty_metric )]

  results_penalty = rbind(results_penalty_train, results_penalty_test)

  results_penalty_m = melt(results_penalty, id.vars = c("penalty_rate",'data_type','policy_name', "penalty_metric"),
                       measure.vars = c("ipw",
                                        "ate_targeted",
                                        "ate_non_targeted",
                                        'target_weight_imbalance',
                                        "target_weight_imbalance_per_dim.Z.V1"
                       ))
  results_penalty_m[, value := as.numeric(value)]


  ## 2.- Collect results of benchmarks, by different target rates
  results_benchmarks = collect_benchmark_results(target_rate)
  results_benchmarks[,penalty_rate:=10]

  results_benchmarks_m = melt(results_benchmarks, id.vars = c("penalty_rate",'data_type','policy_name'),
                               measure.vars = c("ipw",
                                                "ate_targeted",
                                                "ate_non_targeted",
                                                'target_weight_imbalance',
                                                "target_weight_imbalance_per_dim.Z.V1"
                               ))

  ## 3.- Combine and prepare results
  results_all = rbindlist(list(results_penalty,results_benchmarks), use.names = TRUE, fill = TRUE)
  results_all_m = melt(results_all, id.vars = c("penalty_rate",'data_type','policy_name','penalty_metric'),
                              measure.vars = c("ipw",
                                               "ate_targeted",
                                               "ate_non_targeted",
                                               'target_weight_imbalance',
                                               "target_weight_imbalance_per_dim.Z.V1"
                              ))

  results_all_m[, value := as.numeric(value)]
  results_all[, target_rate := target_rate]
  dat_result_list[[count]] = results_all
  count = count+1
  # cache output data
  # out_target_rate = results_all[,.(target_rate = target_rate,
  #                          data_type, penalty_rate,policy_name, ipw,ate_targeted,
  #
  #                          imbalance_Z1 = target_weight_imbalance_per_dim.Z.V1)]

  # foo = rbind(out_target_rate,foo)
}

dat_reulst = rbindlist(dat_result_list, use.names = TRUE)
write_fst(dat_reulst, "cache_policy_results_metric.fst")



# ================================================================
#            Plots
# ----------------------------------------------------------------
# 1. results vs penalty (color metric)

# when penalty rate > 1.5 for cosine similarity rate , predicted tau in [0.01, 0.015]
# so some points are NA
target_rate_plot = 0.2
dat_reulst = read_fst("cache_policy_results_metric.fst", as.data.table = TRUE)
dat_reulst_m = melt(dat_reulst, id.vars = c("penalty_rate",'data_type','policy_name','penalty_metric','target_rate'),
                    measure.vars = c("ipw",
                                     "ate_targeted",
                                     "ate_non_targeted",
                                     'target_weight_imbalance',
                                     "target_weight_imbalance_per_dim.Z.V1"
                    ))
dat_reulst_m[, value := as.numeric(value) ]

ggline(data=dat_reulst_m[!is.na(penalty_metric) & target_rate ==target_rate_plot & data_type=='train' &
                           penalty_metric=='split_l2_norm_rate'&variable=="ipw"],
       x='penalty_rate',
       y="value",
       numeric.x.axis = TRUE,
       #add='loess',
       # plot_type = 'p',
       xlab="Penalty Rate",
       color="penalty_metric",
       title='Comparison of Distance Metrics Per Penalty Rate On Test Data',
       subtitle = sprintf('Target rate: %s', target_rate_plot),
       facet.by = c("variable" ),
       scales='free_y')


ggsave("benchmark_metrics_penalty_rate_test_data.png", width=14, height=10)

# 2. ipw vs penalty rate (color metric and target rate)

dat_reulst[, imbalance_z1:= as.numeric(target_weight_imbalance_per_dim.Z.V1)]
dat_reulst[, ate_targeted:= as.numeric(ate_targeted)]
setorder(dat_reulst, imbalance_z1)

ggline(data=dat_reulst[!is.na(penalty_metric) & target_rate == target_rate_plot],
          x = "imbalance_z1",
          y = "ate_targeted",
       numeric.x.axis = TRUE,
          color="data_type",
          #add = 'loess',
          title='Comparison of Distance Metrics On All Penalty Rate',
          subtitle = sprintf('Target rate: %s; points are smoothed by loess', target_rate_plot),
          xlab = "Imbalance Z1",
          ylab = 'Ate Targeted',
       facet.by = 'penalty_metric',
       scales='free'
       )
ggsave("benchmark_metrics_imbalance_ate.png", width=14, height=10)
