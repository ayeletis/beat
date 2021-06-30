  library(grf)
  library(data.table)
  library(ggpubr)
  library(ggnewscale)
  rm(list=ls())
  is_source=FALSE

  # Working directory
  dir <- "/Applications/Dropbox (Harvard University)/research/heterog_treat_effects/SHARED_BEAT_unintended_bias"
  setwd(dir)
  dir_code <- paste(getwd(),"/code/",sep="")
  dir_data <- paste(getwd(),"/data/",sep="")
  dir_output <- paste(getwd(),"/output/",sep="")
  dir_plots <- paste(getwd(),"/plots",sep="")
  setwd(dir_code)

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

# Explore TAU vs. X from simulation
if(FALSE) {
aux_data = data.table(tau_true = tau_train, Z = 1*(Z_train[,1]==1),
                      X1 = X_train[,1], X2 = X_train[,2], X3 = X_train[,3], X4 = X_train[,4])
names(aux_data)[2] = "Z"
aux_data$Z = as.factor((aux_data$Z))
names(aux_data)[3] = "X1"
names(aux_data)[4] = "X2"
names(aux_data)[5] = "X3"
names(aux_data)[6] = "X4"

# create plots
{
  ## Targeting plots = TAU by Z
  #true
  p0a = ggdensity(data=aux_data,
                  x='tau_true', color='Z', fill='Z', alpha=0.2, add = "mean")
  ## Plot TAUs vs. X
  # true
  p0b = ggplot(aux_data, aes(x=X1, y=tau_true)) +
    geom_point()+
    geom_smooth()
  p0c = ggplot(aux_data, aes(x=X2, y=tau_true)) +
    geom_point()+
    geom_smooth()
  p0d = ggplot(aux_data, aes(x=X3, y=tau_true)) +
    geom_point()+
    geom_smooth()
  p0e = ggplot(aux_data, aes(x=X4, y=tau_true)) +
    geom_point()+
    geom_smooth()
} # end of creating plots

# plot config.
{
  base_font_size = 14
  label_font_size = 14
  if(Sys.info()['sysname']=="Windows"){windowsFonts(Times=windowsFont("Times New Roman"))}
  plot_theme =  theme(
    text=element_text(family="Times", face="plain", size=base_font_size),
    axis.title=element_text(family="Times", face="plain",size=label_font_size),
    panel.background = element_rect(fill = "white", colour = "gray",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80",size = 0.1),
    panel.grid.minor.y = element_line(color = "grey80",size = 0.1))
  theme_set(plot_theme)

} # end of plot config.

# combining plots
p_tau_x = ggarrange(p0b+ylab("")+xlim(c(-2.5,2.5)),
                      p0c+ylab("")+xlim(c(-2.5,2.5)),
                      p0d+ylab("")+xlim(c(-2.5,2.5)),
                      p0e+ylab("")+xlim(c(-2.5,2.5)),
                      ncol = 2,nrow = 2)
} # end of explore TAU vs. X

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
{
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


  # balanced tau
  fit_fair <- balanced_causal_forest(X_train, Y_train, W_train,
                                     target.weights = as.matrix(Z_train),
                                     target.weight.penalty = my_penalty,
                                     target.weight.standardize = TRUE,
                                     target.weight.bins.breaks = 256,
                                     target.weight.penalty.metric = 'split_l2_norm_rate',
                                     honesty = TRUE,
                                     num.trees = num_trees,
                                     tune.parameters=my_tunable_params,
                                     seed=my_seed)

} # end of estimation methods


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
    )
  )
} # end of predict taus

## -------------------------------------------
##   Select target units (both in-sample and out-of-sample)
## -------------------------------------------
target_rate = 0.2 # target ## % of the population
{
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
      ),
      'grf_fair' =  list(
        'train' = as.integer(tau.pred$grf_fair$train > quantile(tau.pred$grf_fair$train, c(1-target_rate))),
        'test' = as.integer(tau.pred$grf_fair$test > quantile(tau.pred$grf_fair$test, c(1-target_rate)))
      )
    )
} # end of select targets

# -------------------------------------------------
# OFF-POLICY EVALUATION
# -------------------------------------------------
{
  W.hat_train = mean(W_train) # We use prop(treatment) because the data comes from an experiment
  W.hat_test = mean(W_test)
  # Otherwise, estimate propensities from the data.
  if (FALSE) {
    forest.W = regression_forest(X_train_demog, W_train, tune.parameters = "all")
    W.hat_train = predict(forest.W)$predictions
    W.hat_test = predict(forest.W,X_test_demog)$predictions
  }

  # Function that combines the results (train and test) for each targeting policy
  collect_results = function(policy_list, policy_name){
    policy = policy_list[[policy_name]]
    out_train = off_policy_evalution_v3(X_train, W_train, Y_train, Z_train,policy$train, policy_name, 'train', W.hat_train)
    out_test = off_policy_evalution_v3(X_test, W_test, Y_test, Z_test, policy$test, policy_name, 'test', W.hat_test)
    return(data.table(rbind(unlist(out_train),unlist(out_test) )))
  }

  # Get results for all the policies in "targeting
  results = rbindlist(lapply(names(my_targeting_policies),
                             collect_results,
                             policy_list= my_targeting_policies))


} # end of policy evaluation
# a = results[policy_name == "grf_fair"]
# b=readRDS("a.rds")
# all.equal(a,b)

## -----------------------------------------------------------------------------------------
##   Visuals for motivation and results
## -----------------------------------------------------------------------------------------
{
aux_data = data.table(tau_true = tau_train, tau_grf_d = tau.pred$grf_demog$train, tau_grf = tau.pred$grf$train,
                          tau_fair = tau.pred$grf_fair$train, Z = 1*(Z_train[,1]==1),
                          X1 = X_train[,1], X2 = X_train[,2], X3 = X_train[,3], X4 = X_train[,4])
names(aux_data)[5] = "Z"
aux_data$Z = as.factor((aux_data$Z))
names(aux_data)[6] = "X1"
names(aux_data)[7] = "X2"
names(aux_data)[8] = "X3"
names(aux_data)[9] = "X4"

aux_data[, Protected := factor(Z,
                                       labels = c("0", "1"),
                                       levels=c('0','1'))]

# create plots
{
    ## SCORES
    #true
    p0 = ggdensity(data=aux_data,
                    x='tau_true', fill = "lightgray", alpha=0.2, add = "mean")
    # grf + demog
    p1 = ggdensity(data=aux_data,
                    x='tau_grf_d',fill = "lightgray", alpha=0.2, add = "mean")
    # grf
    p2 = ggdensity(data=aux_data,
                    x='tau_grf', fill = "lightgray", alpha=0.2, add = "mean")
    # beat
    p3 = ggdensity(data=aux_data,
                    x='tau_fair', fill = "lightgray", alpha=0.2, add = "mean")

    ## SCORES by Z
    #true
    p0a = ggdensity(data=aux_data,
                      x='tau_true', color='Protected', fill='Protected', alpha=0.2, add = "mean")
    # grf + demog
    p1a = ggdensity(data=aux_data,
                      x='tau_grf_d', color='Protected', fill='Protected', alpha=0.2, add = "mean")
    # grf
    p2a = ggdensity(data=aux_data,
                      x='tau_grf', color='Protected', fill='Protected', alpha=0.2, add = "mean")
    # beat
    p3a = ggdensity(data=aux_data,
                      x='tau_fair', color='Protected', fill='Protected', alpha=0.2, add = "mean")

    ## Plot TAUs vs. X
    # true
    p0b = ggplot(aux_data, aes(x=X1, y=tau_true)) +
      geom_point()+
      geom_smooth()
    p0c = ggplot(aux_data, aes(x=X2, y=tau_true)) +
      geom_point()+
      geom_smooth()
    p0d = ggplot(aux_data, aes(x=X3, y=tau_true)) +
      geom_point()+
      geom_smooth()
    p0e = ggplot(aux_data, aes(x=X4, y=tau_true)) +
      geom_point()+
      geom_smooth()

    # grf + demog
    p1b = ggplot(aux_data, aes(x=X1, y=tau_grf_d)) +
      geom_point()+
      geom_smooth()
    p1c = ggplot(aux_data, aes(x=X2, y=tau_grf_d)) +
      geom_point()+
      geom_smooth()
    p1d = ggplot(aux_data, aes(x=X3, y=tau_grf_d)) +
      geom_point()+
      geom_smooth()
    p1e = ggplot(aux_data, aes(x=X4, y=tau_grf_d)) +
      geom_point()+
      geom_smooth()

    # grf
    p2b = ggplot(aux_data, aes(x=X1, y=tau_grf)) +
      geom_point()+
      geom_smooth()
    geom_smooth()
    p2c = ggplot(aux_data, aes(x=X2, y=tau_grf)) +
      geom_point()+
      geom_smooth()
    p2d = ggplot(aux_data, aes(x=X3, y=tau_grf)) +
      geom_point()+
      geom_smooth()
    p2e = ggplot(aux_data, aes(x=X4, y=tau_grf)) +
      geom_point()+
      geom_smooth()

    # fair
    p3b = ggplot(aux_data, aes(x=X1, y=tau_fair)) +
      geom_point()+
      geom_smooth()
    p3c = ggplot(aux_data, aes(x=X2, y=tau_fair)) +
      geom_point()+
      geom_smooth()
    p3d = ggplot(aux_data, aes(x=X3, y=tau_fair)) +
      geom_point()+
      geom_smooth()
    p3e = ggplot(aux_data, aes(x=X4, y=tau_fair)) +
      geom_point()+
      geom_smooth()

  } # end of creating plots

# plot config.
{
    base_font_size = 14
    label_font_size = 14
    if(Sys.info()['sysname']=="Windows"){windowsFonts(Times=windowsFont("Times New Roman"))}
    plot_theme =  theme(
      text=element_text(family="Times", face="plain", size=base_font_size),
      axis.title=element_text(family="Times", face="plain",size=label_font_size),
      panel.background = element_rect(fill = "white", colour = "gray",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major.y = element_line(color = "grey80",size = 0.1),
      panel.grid.minor.y = element_line(color = "grey80",size = 0.1))
    theme_set(plot_theme)

  } # end of plot config.

# combining plots
{
    my_xlim = c(-2.5,2.5)
    my_ylim = c(-2.5,2.5)

    p_example_grf_demog = ggarrange(p1+ylab("")+xlim(1.2*my_xlim)+xlab("(a) Predicted CATE"),
                            p1b+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(b) X1"),
                        p1c+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(c) X2 (corr. Z1)"),
                        p1d+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(d) X3"),
                        p1e+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(e) X4"),
                        p1a+ylab("")+xlim(1.2*my_xlim)+xlab("(f) Target score"),
                        ncol = 6,nrow = 1)

    p_example_grf = ggarrange(p2+ylab("")+xlim(1.2*my_xlim)+xlab("(a) Predicted CATE"),
                              p2b+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(b) X1"),
                              p2c+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(c) X2 (corr. Z1)"),
                              p2d+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(d) X3"),
                              p2e+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(e) X4"),
                              p2a+ylab("")+xlim(1.2*my_xlim)+xlab("(f) Target score"),
                              ncol = 6,nrow = 1)

    p_example_beat = ggarrange(p3+ylab("")+xlim(1.2*my_xlim)+xlab("(a) Predicted CBT"),
                            p3b+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(b) X1"),
                            p3c+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(c) X2 (corr. Z1)"),
                            p3d+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(d) X3"),
                            p3e+ylab("")+xlim(my_xlim)+ylim(my_ylim)+xlab("(e) X4"),
                            p3a+ylab("")+xlim(1.2*my_xlim)+xlab("(f) Target score"),
                            ncol = 6,nrow = 1)

  } # end of combining plots

p_example_grf_demog
p_example_grf
p_example_beat


setwd(dir_plots)
#ggsave("p_example_grf_demog.png",plot = p_example_grf_demog, width = 14, height=3)
#ggsave("p_example_grf.png",plot = p_example_grf, width = 14, height=3)
#ggsave("p_example_beat.png",plot = p_example_beat, width = 14, height=3)
}

## -----------------------------------------------------------------------------------------
##   Degree of over-targeting and output
## -----------------------------------------------------------------------------------------
{
ratio_grf_demo = c(mean((Z_train[,1]==1)*my_targeting_policies$grf_demog$train),
                   mean((Z_train[,1]==0)*my_targeting_policies$grf_demog$train))

ratio_grf = c(mean((Z_train[,1]==1)*my_targeting_policies$grf$train),
              mean((Z_train[,1]==0)*my_targeting_policies$grf$train))

ratio_grf_fair = c(mean((Z_train[,1]==1)*my_targeting_policies$grf_fair$train),
                   mean((Z_train[,1]==0)*my_targeting_policies$grf_fair$train))

prop_targeted = data.table(rbind(ratio_grf_demo,ratio_grf,ratio_grf_fair))
names(prop_targeted) = c("z_1","z_0")
prop_targeted[,ratio := z_1/z_0]
prop_targeted

results[data_type=='train',.(ate_targeted,policy_name)]
}

# ================================================================
#    Estimate BEAT for different PENALTY and different TARGET_RATE
# ----------------------------------------------------------------
library(fst)
library(ggnewscale)
library(directlabels)
library(ggforce)

# Function that selects targets based on tau_pred
target_tau_proportion  = function(tau_pred, target_rate){
  return(as.integer(tau_pred > quantile(tau_pred, c(1-target_rate))))
}

# Function that runs BEAT for different levels of penalty -- returns TAUs
estimate_taus_per_penalty = function(penalty_value){
# this following copied from codes above and recycle the same parameters
  print(penalty_value)
  fit_fair <- balanced_causal_forest(X_train, Y_train, W_train,
                                     target.weights = as.matrix(Z_train),
                                     target.weight.penalty = penalty_value,
                                     target.weight.standardize = TRUE,
                                     target.weight.bins.breaks = 256,
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

  return(out)

}

## Set range for penalty
penalty_vect =   c(seq(0.1,0.7,0.1),seq(0.8, 1.4, 0.02),seq(1.5,3,0.5))
#penalty_vect =   c(1,5)
penalty_vect = seq(0.1, 2.1, 0.3)
## Estimate BEAT per penalty and save TAUs
taus_per_penalty = rbindlist(lapply(penalty_vect, estimate_taus_per_penalty))
setwd(dir_output)
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



