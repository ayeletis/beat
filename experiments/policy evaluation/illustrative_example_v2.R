  library(grf)
  library(data.table)
  library(ggpubr)
  
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
  my_penalty = 10
  my_seed = 1
  
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

# Explore TAU vs. X
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
  
  
  # for BEAT, lets' use the same tunable parameters as in regular GRF
  my_min_node_size  = fit_grf$tunable.params$min.node.size 
  my_imbalance_penalty  = fit_grf$tunable.params$imbalance.penalty  
  my_alpha  = fit_grf$tunable.params$alpha 
  my_sample_fraction = fit_grf$tunable.params$sample.fraction 
  my_mtry = fit_grf$tunable.params$mtry 
  
  # balanced tau
  fit_fair <- balanced_causal_forest(X_train, Y_train, W_train,
                                     target.weights = as.matrix(Z_train),
                                     target.weight.penalty = my_penalty,
                                     target.weight.standardize = TRUE,
                                     target.weight.bins.breaks = 256,
                                     honesty = TRUE,
                                     num.trees = num_trees,
                                     min.node.size = my_min_node_size,
                                     imbalance.penalty = my_imbalance_penalty,
                                     alpha = my_alpha,
                                     sample.fraction = my_sample_fraction,
                                     mtry = my_mtry,
                                     seed=my_seed)
} # end of estim scores

## -------------------------------------------
##   Select target units (both in-sample and out-of-sample)
## -------------------------------------------
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
  
  ## Target based on score
  target_proportion = 0.2 # target half the population
  set.seed(my_seed) # we fix seed for the random targeting

    my_targeting_policies = list(
      'grf_demog' = list(
        'train' = as.integer(tau.pred$grf_demog$train > quantile(tau.pred$grf_demog$train, c(1-target_proportion))),
        'test' = as.integer(tau.pred$grf_demog$test > quantile(tau.pred$grf_demog$test, c(1-target_proportion)))
      ),
      'grf' = list(
        'train' = as.integer(tau.pred$grf$train > quantile(tau.pred$grf$train, c(1-target_proportion))),
        'test' = as.integer(tau.pred$grf$test > quantile(tau.pred$grf$test, c(1-target_proportion)))
      ),
      'random' = list(
        'train' = as.integer(tau.pred$random$train > quantile(tau.pred$random$train, c(1-target_proportion))),
        'test' = as.integer(tau.pred$random$test > quantile(tau.pred$random$test, c(1-target_proportion)))
      ),
      'grf_fair' =  list(
        'train' = as.integer(tau.pred$grf_fair$train > quantile(tau.pred$grf_fair$train, c(1-target_proportion))),
        'test' = as.integer(tau.pred$grf_fair$test > quantile(tau.pred$grf_fair$test, c(1-target_proportion)))
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

## -----------------------------------------------------------------------------------------
##   Visuals for motivation and results
## -----------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------
##   Degree of over-targeting and output
## -----------------------------------------------------------------------------------------

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

# ================================================================ 
#            Check Penalty Rate vs ATE vs Z
# ----------------------------------------------------------------
# Z1 is correlated with X, defined in  z_function
#
library(fst)
library(ggnewscale)

library(directlabels)
library(ggforce)

target_tau_proportion  = function(tau_pred, target_rate){
  return(as.integer(tau_pred > quantile(tau_pred, c(1-target_rate))))
}

# for balanced grf
gather_tau_results = function(penalty_value){
# this following copied from codes above and recycle the same parameters  
  print(penalty_value) 
  fit_fair <- balanced_causal_forest(X_train, Y_train, W_train,
                                     target.weights = as.matrix(Z_train),
                                     target.weight.penalty = penalty_value,
                                     target.weight.standardize = TRUE,
                                     target.weight.bins.breaks = 256,
                                     honesty = TRUE,
                                     num.trees = num_trees,
                                     min.node.size = my_min_node_size,
                                     imbalance.penalty = my_imbalance_penalty,
                                     alpha = my_alpha,
                                     sample.fraction = my_sample_fraction,
                                     mtry = my_mtry,
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

# get results per penalty  1.4
penalty_vect =   c(seq(0.8, 1.4, 0.02))


tau_penalty = rbindlist(lapply(penalty_vect, gather_tau_results))

# write_fst(tau_penalty, "cache_tau_per_penalty.fst")
# write_fst(tau_penalty, "cache_tau_per_penalty_combine.fst")
# tau_penalty = read_fst('cache_tau_per_penalty_combine.fst', as.data.table = TRUE)

# organize results into one row per penalty 
W.hat_train = mean(W_train) 

collect_balance_grf_results = function(target_vector, data_type="train"){
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
                                  W.hat = W.hat_train)
    out = rbind(unlist(out))
    out = as.data.table(out)
  }else{
    stop("data_type is either train or test")
  }
  
  return(out)
  
}

gather_benchmark_others = function(target_rate){
  # copy the same code for above 
  
  collect_results_other = function(policy_list, policy_name){
    policy = policy_list[[policy_name]]
    out_train = off_policy_evalution_v3(X_train, W_train, Y_train, Z_train,policy$train, policy_name, 'train', W.hat_train)
    return(data.table(rbind(unlist(out_train))))
  }
  
  my_targeting_policies = list(
    'grf_demog' = list(
      'train' = as.integer(tau.pred$grf_demog$train > quantile(tau.pred$grf_demog$train, c(1-target_rate)))
    ),
    'grf' = list(
      'train' = as.integer(tau.pred$grf$train > quantile(tau.pred$grf$train, c(1-target_rate)))
    ),
    'random' = list(
      'train' = as.integer(tau.pred$random$train > quantile(tau.pred$random$train, c(1-target_rate)))
    ) 
  )
  
  results = rbindlist(lapply(names(my_targeting_policies),
                             collect_results_other,
                             policy_list= my_targeting_policies))
}
 

target_data_list = list()
counter=1
for(target_rate in c(0.2,0.4,0.5, 0.6, 0.8)){
  
  if ('is_targeted' %in%  names(tau_penalty)){
    tau_penalty[, is_targeted:= NULL]
  }
  
  tau_penalty[, is_targeted := target_tau_proportion(tau_predict, 
                                                     target_rate=target_rate), 
              by=.(penalty_rate, data_type)]
  

  tau_results_train = tau_penalty[data_type=="train", 
                                  collect_balance_grf_results(is_targeted,"train"),  
                            by=.(penalty_rate )]
  tau_results_test = tau_penalty[data_type=="test", 
                                  collect_balance_grf_results(is_targeted,"test"),  
                                  by=.(penalty_rate )]
  
  tau_results = rbind(tau_results_train, tau_results_test)
  
  tau_results_m = melt(tau_results, id.vars = c("penalty_rate",'data_type'),
                       measure.vars = c("ipw",
                                        "ate_targeted",
                                        "ate_non_targeted",
                                        "target_weight_imbalance_per_dim.Z.V1" 
                                        
                       ))
  tau_results_m[, value := as.numeric(value)]
  
  # use results from above
  benchmark_others = gather_benchmark_others(target_rate)
  
  tau_benchmark = benchmark_others[data_type=="train" & policy_name %in% c("grf_demog", "grf", "random"),
                          .(benchmark = policy_name, ipw, 
                            ate_targeted, ate_non_targeted,
                            target_weight_imbalance_per_dim.Z.V1)]
  
  
  tau_benchmark_m = melt(tau_benchmark, id.vars = 'benchmark')
  tau_benchmark_m[, value := as.numeric(value)]

  p = ggline(tau_results_m,
         x = "penalty_rate",
         y = "value",
         facet.by="variable",
         color="data_type",
         shape = "data_type",
         numeric.x.axis = TRUE,
         scale = "free_y",
         legend='bottom',
         title="Targeting Policy Evaluation on Penalty Rate for Balanced GRF",
         subtitle=sprintf("Z.V1 is correlated with X; Target Rate: %s", target_rate),
         ncol=2
  ) + grids() + 
    coord_cartesian(xlim=c(0.8, 1.4)) + 
    scale_x_continuous(breaks=seq(0.8, 1.4, 0.05)) + scale_color_brewer(palette="Set1") + 
    new_scale_color()  + 
    geom_hline(data=tau_benchmark_m,
               mapping=aes(yintercept=value,colour=benchmark),
               linetype="longdash") + 
    # add label to hlines
    geom_text(data=tau_benchmark_m, 
              mapping=aes(x=1.4, y=value, label=benchmark))
  
  
  ggsave(sprintf("target_policy_penalty_rate_balanced_grf_target_rate_%s.png", target_rate), p,
         width = 14, height=10)

  
  # cache output data
  out_tau = tau_results[,.(target_rate = target_rate, 
                           data_type, penalty_rate,policy_name, ipw,ate_targeted, 
                 imbalance_Z1 = target_weight_imbalance_per_dim.Z.V1)]

  out_bench = tau_benchmark[,.(target_rate = target_rate,
                               data_type = "benchmark",
                               policy_name = benchmark, 
                               ipw,ate_targeted, 
                   imbalance_Z1 = target_weight_imbalance_per_dim.Z.V1)]
  
  
  fixed_cols= expand.grid(penalty_rate = out_tau[, penalty_rate], policy_name = out_bench[, policy_name])
  fixed_cols = as.data.table(fixed_cols)
  out_bench = fixed_cols[out_bench, on=.(policy_name)] 
  target_data_list[[counter]] = rbind(out_tau, out_bench )
  counter = counter  + 1
}

target_data = rbindlist(target_data_list)
fwrite(target_data, "target_data_penalty_benchmark.csv")

# counter plot
# x: imbalance z1, y: ate target, z: target proportion 
get_data = function(target_rate){
  tau_penalty[, is_targeted := target_tau_proportion(tau_predict, 
                                                     target_rate=target_rate), 
              by=.(penalty_rate, data_type)]
  
  
  tau_results_train = tau_penalty[data_type=="train", collect_balance_grf_results(is_targeted,'train'), 
                                  by=penalty_rate]
  tau_results_test = tau_penalty[data_type=="test", collect_balance_grf_results(is_targeted,"test"),
                                 by=penalty_rate]
  tau_results = rbind(tau_results_train, tau_results_test)
  
  out = tau_results[,.(imbalance_z1 = as.numeric(target_weight_imbalance_per_dim.Z.V1), 
                       ate_targeted = as.numeric(ate_targeted),
                       ipw = as.numeric(ipw),
                       data_type )]
  out[, target_rate:= target_rate]
  return(out)
}

contour_data = rbindlist(lapply(seq(0.1,0.9,0.1), get_data))



contour_data[, add_text := paste0("T.R: ",target_rate)]
 
p_min = ggline(contour_data[target_rate>0.6 & data_type == "train"],
               x="imbalance_z1",
               y="ate_targeted", 
               color="target_rate",
               legend="none",
               xlab="",
               ylab="",
               numeric.x.axis = TRUE) + 
  geom_dl(mapping = aes(label = add_text), 
          method = list("last.points", 
                        hjust=-0.1, vjust=1.5) ) + 
  coord_cartesian(xlim=c(0,0.02))+
  scale_color_gradient(low="red",high="blue",limits=c(0.1,0.9) )

ggline(contour_data[  data_type == "train"],
       x="imbalance_z1",
       y="ate_targeted", 
       xlab="Z1 Imbalance",
       ylab="ATE",
       title="Comparison on ATE VS Imbalance In The Targeted Group",
       subtitle="Z1 is correlated to X",
       color="target_rate",
       legend="right",
       numeric.x.axis = TRUE) + 
  labs(color="Target Rate") + 
  scale_color_gradient(low="red",high="blue",limits=c(0.1,0.9) )+
  coord_cartesian(xlim=c(0,0.09)) +
  # scale_x_log10()
  scale_x_continuous(brea=seq(0,0.09, 0.01)) +
  scale_y_continuous(breaks=seq(0, 2, 0.2)) +
  geom_dl(mapping = aes(label = add_text), 
          method = list("last.points", 
                        hjust=-0.1, vjust=1.5) )  +
  annotation_custom(ggplotGrob(p_min), xmin=0.05, xmax=0.09, ymin=0.1, ymax=0.8)
ggsave("ate_vs_z1_imbalance.png", width=10, height=8)


ggline(contour_data[  data_type == "train"],
       x="imbalance_z1",
       y="ipw", 
       xlab="Z1 Imbalance",
       ylab="IPW",
       title="Comparison on IPW VS Imbalance In The Targeted Group",
       subtitle="Z1 is correlated to X",
       color="target_rate",
       legend="right",
       numeric.x.axis = TRUE) + 
  labs(color="Target Rate") + 
  scale_color_gradient(low="red",high="blue",limits=c(0.1,0.9) )+
  coord_cartesian(xlim=c(0,0.09), ylim=c(0.63, 0.9)) +
  scale_x_continuous(brea=seq(0,0.09, 0.01)) +
  scale_y_continuous(breaks=seq(0, 1, 0.05)) +
  geom_dl(mapping = aes(label = add_text), 
          method = list("last.points", 
                        hjust=-0.1, vjust=1.5) )   
ggsave("ipw_vs_z1_imbalance.png", width=10, height=8)
