  library(grf)
  library(data.table)
  library(ggpubr)
  library(ggnewscale)



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

# demographic tau
# fit_grf_demog <- causal_forest(X = X_train_demog, Y = Y_train, W = W_train,
#                                honesty = TRUE,
#                                num.trees = num_trees,
#                                tune.parameters='all',
#                                seed=my_seed)
#
#
# # regular tau -- no demographics
# a = Sys.time()
# fit_grf <- causal_forest(X = X_train, Y = Y_train, W = W_train,
#                          honesty = TRUE,
#                          num.trees = num_trees,
#                          tune.parameters='all',
#                          seed=my_seed)
#
# print(Sys.time()-a)

a = Sys.time()

# regular tau -- no demographics
fit_grf = balanced_causal_forest(X_train, Y_train, W_train,
                       target.weights = as.matrix(Z_train ) ,
                       target.weight.penalty = 1,
                       target.weight.standardize = TRUE,
                       target.weight.bins.breaks = 256,
                       # target.weight.penalty.metric = 'cosine_similarity_rate',
                       honesty = TRUE,
                       num.trees = num_trees,
                       tune.parameters=my_tunable_params,
                       seed=my_seed)
print(Sys.time()-a)

# ================================================================
#    Estimate BEAT for different PENALTY and different TARGET_RATE
# ----------------------------------------------------------------
library(fst)
library(ggnewscale)
library(directlabels)
library(ggforce)


available_metrics = c("split_l2_norm_rate", # left, right: l2 norm(colmean target weight)* penalty rate * node decrease
                      "euclidean_distance_rate", # (left+right decrease) *  Euclidean distance (column mean target weight left, right ) * penalty rate
                      "cosine_similarity_rate", # (left+right decrease) *  (1-cos_sim(column mean target weight left, right )) * penalty rate

                      "split_l2_norm", #  sum(left,right l2 norm(colmean target weight))* penalty rate
                      "euclidean_distance", #  Euclidean distance (column mean target weight left, right ) * penalty rate
                      "cosine_similarity" #  (1-cos_sim(column mean target weight left, right )) * penalty rate
)


# Function that selects targets based on tau_pred
target_tau_proportion  = function(tau_pred, target_rate){
  return(as.integer(tau_pred > quantile(tau_pred, c(1-target_rate))))
}

# Function that runs BEAT for different levels of penalty -- returns TAUs
estimate_taus_per_penalty = function(penalty_value, metric){
# this following copied from codes above and recycle the same parameters
  print(penalty_value)
  fit_fair <- balanced_causal_forest(X_train, Y_train, W_train,
                                     target.weights = as.matrix(Z_train),
                                     target.weight.penalty = penalty_value,
                                     target.weight.standardize = TRUE,
                                     target.weight.bins.breaks = 256,
                                     target.weight.penalty.metric = metric,
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
penalty_vect = c(0.1, 0.5, 1, 1.5, 2, 3)   #c(seq(0.1,0.7,0.1),seq(0.8, 1.4, 0.02),seq(1.5,3,0.5))

taus_per_penalty = list()

for(i in available_metrics){
  a = Sys.time()
    print(i)
  taus_per_penalty[[i]] = rbindlist(lapply(penalty_vect, estimate_taus_per_penalty, metric=i))
  print(Sys.time()-a)
}



write_fst(taus_per_penalty, "cache_tau_per_penalty.fst")
taus_per_penalty = read_fst('cache_tau_per_penalty.fst', as.data.table = TRUE)


# ================================================================
#   Collect policy results per target_rate levels -- both for the benchmarks and for BEAT penalty
# ----------------------------------------------------------------

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


foo = NULL
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
                            by=.(penalty_rate )]
  results_penalty_test = taus_per_penalty[data_type=="test",
                                  collect_penalty_results(is_targeted,"test"),
                                  by=.(penalty_rate )]

  results_penalty = rbind(results_penalty_train, results_penalty_test)

  results_penalty_m = melt(results_penalty, id.vars = c("penalty_rate",'data_type','policy_name'),
                       measure.vars = c("ipw",
                                        "ate_targeted",
                                        "ate_non_targeted",
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
                                                "target_weight_imbalance_per_dim.Z.V1"
                               ))

  ## 3.- Combine and prepare results
  results_all = rbind(results_penalty,results_benchmarks)
  results_all_m = melt(results_all, id.vars = c("penalty_rate",'data_type','policy_name'),
                              measure.vars = c("ipw",
                                               "ate_targeted",
                                               "ate_non_targeted",
                                               "target_weight_imbalance_per_dim.Z.V1"
                              ))

  results_all_m[, value := as.numeric(value)]


  # cache output data
  out_target_rate = results_all[,.(target_rate = target_rate,
                           data_type, penalty_rate,policy_name, ipw,ate_targeted,
                           imbalance_Z1 = target_weight_imbalance_per_dim.Z.V1)]

  foo = rbind(out_target_rate,foo)
}
policy_results = foo

setwd(dir_output)
write_fst(policty_results, "cache_policy_results.fst")
policy_results = read_fst('cache_policy_results.fst', as.data.table = TRUE)




# ================================================================
#    PLOTS
# ----------------------------------------------------------------
# converte to numeric
policy_results[, ipw:= as.numeric(ipw)]
policy_results[, ate_targeted:= as.numeric(ate_targeted)]
policy_results[, imbalance_Z1:= as.numeric(imbalance_Z1)]


policy_results_m = melt(policy_results,
                        id.vars = c("target_rate", 'penalty_rate',"data_type", "policy_name") )
policy_results_m[, value:= as.numeric(value)]
## by target rate
my_target_rate = 0.2
data_to_plot = policy_results_m[target_rate==my_target_rate]
results_benchmarks_m = data_to_plot[policy_name!="balanced_grf"]
results_benchmarks_m[, line_color := paste0(policy_name,'-',data_type)]

### HENGYU, Could you please recreate your plots, from here on, using policy_results??
# line plots for penalty rate

   p = ggline(data_to_plot[policy_name=="balanced_grf"],
         x = "penalty_rate",
         y = "value",
         facet.by=c( "variable"),
         color="data_type",
         shape = "data_type",
         numeric.x.axis = TRUE,
         scale = "free_y",
         legend='bottom',
         title="Targeting Policy Evaluation on Penalty Rate for Balanced GRF",
         subtitle=sprintf("Z.V1 is correlated with X; Target Rate: %s", my_target_rate),
         ncol=2
  ) + grids() +
    coord_cartesian(xlim=c(0.8, 1.5)) + scale_color_brewer(palette="Set1") +
    scale_x_continuous(breaks=seq(0.8, 2, 0.05))  +
    new_scale_color()  + scale_color_brewer(palette="Set1") +
    geom_hline(data=results_benchmarks_m,
               mapping=aes(yintercept=value,color=data_type, linetype=policy_name      )
               )
    # we include too much text and hard to read
    # geom_text(data=results_benchmarks_m,
    #           mapping=aes(x=1, y=value, label=policy_name),
    #           position = position_jitter(seed=1, width=0.1, height = 0))
    #

  ggsave(sprintf("target_policy_penalty_rate_balanced_grf_target_rate_%s.png", my_target_rate), p,
         width = 14, height=10)


  # plot per target rate
  # last point per group

  point_annot_dat = policy_results[target_rate==my_target_rate,
                                   lapply(.SD, tail, n=1),
                                   by=.(policy_name, data_type)]
  main_points = policy_results[ target_rate==my_target_rate & policy_name=="balanced_grf"]
  benchmark_points = policy_results[ target_rate==my_target_rate &  policy_name!="balanced_grf"]
  #
  p_cont =  ggline(data=main_points,
         x = "imbalance_Z1",
         y = "ate_targeted",
         color = "data_type",
         palette = "Set1",
         size=1,
         numeric.x.axis = TRUE,
         legend="top",
         xlab = "Imbalance Z1",
         ylab = 'Ate Targeted',
         title="Ate Targeted VS Imbalance Z1",
         subtitle=sprintf("Target rate: %s; Z1 is correlated to X", my_target_rate),
         shape='policy_name') +
    geom_point(data=point_annot_dat,
               mapping=aes(x=imbalance_Z1, y=ate_targeted,
                           color=data_type,
                           shape=policy_name ),
               size=3)   +
   geom_text(data=point_annot_dat[data_type=='train'],
             mapping=aes(x=imbalance_Z1, y=ate_targeted, label=policy_name),
             position = position_nudge())

  ggsave(sprintf("line_plot_ate_vs_imbalance_z1_target_rate_%s.png", my_target_rate), p_cont,
         width = 12, height=10)



