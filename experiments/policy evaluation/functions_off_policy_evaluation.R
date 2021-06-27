# Inverse probability scores (function) -- it needs Y, W, the propensities, and the policy to evaluate
v_ips_function <- function(Y,W,ehat,my_policy) {
  ips = (W==my_policy)*(W==0)/(1-ehat) + (W==my_policy)*(W==1)/ehat
    temp = mean(Y*ips)
   return(temp)
}


# Model based (function) -- it needs Y.hat and the policy to evaluate
v_mb_function <- function(Y_0,Y_1,my_policy) {
  rewards = Y_0*(my_policy==0) + Y_1*(my_policy==1)
  temp = mean(rewards)
  return(temp)
}


# Doubly robust Model-based (function) -- It needs Y.hat, Y, W, the propensities, and the polity to evaluate
v_dr_function <- function(Y_0,Y_1,Y,W,ehat,my_policy) {
  rewards = Y_0*(my_policy==0) + Y_1*(my_policy==1)
  dr_rewards = rewards + (Y-rewards)*((W==my_policy)*(W==0)/(1-ehat) + (W==my_policy)*(W==1)/ehat)
  temp = mean(dr_rewards)
  return(temp)
}


# Function that evaluates the empirical ATE mean(Y|treated) - mean(Y|non_treated) for each policy
empirical_ate <- function(Y,W,my_policy,is_targeted=1) {
  dt = data.table(Y=Y,W=W,my_policy=my_policy)
  ate = dt[my_policy==is_targeted & W==1,mean(Y)] - dt[my_policy==is_targeted & W==0,mean(Y)]
  return(ate)
}


# Function that evaluates the predicted TAU : = mean(Y|treated) - mean(Y|non_treated) for each policy
empirical_tau <- function(tau,my_policy,is_targeted=1) {
  dt = data.table(tau=tau,my_policy=my_policy)
  tau = dt[my_policy==is_targeted,mean(tau)] 
  return(tau)
}

# Function that evaluates (1) Reward (IPW), (2) % users targeted, (3) imbalance, for each policy and dataset
off_policy_evalution = function(X, W, Y, Z,
                                policy_targets,
                                policy_name,
                                data_type,
                                W.hat){
  
  
  # Reward (IPW estimator)
  ipw  = v_ips_function(Y, W, W.hat, policy_targets)
  
  # Proportion of users targeted
  target_ratio = mean(policy_targets)
  
  # Amount of imbalance
  aux_avg = colMeans(Z)
  aux_policy = colMeans(Z[policy_targets==1])
#  target_weight_imbalance  = sum(abs(aux_policy-aux_avg))
  target_weight_imbalance  = sum((aux_policy-aux_avg)^2)
  
  return(list('ipw' = ipw,
              'target_ratio' = target_ratio,
              'target_weight_imbalance' = target_weight_imbalance,
              'policy_name'=policy_name,
              'data_type'=data_type))
}



# Function that evaluates (1) Reward (IPW), (2) % users targeted, (3) imbalance, (4) empirical ATE for treated and for non-treated, for each policy and dataset
off_policy_evalution_v3 = function(X, W, Y, Z,
                                policy_targets,
                                policy_name,
                                data_type,
                                W.hat){
  
  
  # Reward (IPW estimator)
  ipw  = v_ips_function(Y, W, W.hat, policy_targets)
  
  # Proportion of users targeted
  target_ratio = mean(policy_targets)
  
  # Amount of imbalance
  aux_avg = colMeans(Z)
  aux_policy = colMeans(Z[policy_targets==1])
  #  target_weight_imbalance  = sum(abs(aux_policy-aux_avg))
  target_weight_imbalance  = sum((aux_policy-aux_avg)^2)
  
  target_weight_imbalance_per_dim = (aux_policy-aux_avg)^2
  
  
  # Empirical ATE
  ate_targeted = empirical_ate(Y,W,policy_targets,1)
  ate_non_targeted = empirical_ate(Y,W,policy_targets,0)
  
  return(list('ipw' = ipw,
              'target_ratio' = target_ratio,
              'target_weight_imbalance' = target_weight_imbalance,
              'target_weight_imbalance_per_dim' = target_weight_imbalance_per_dim,
              'ate_targeted' = ate_targeted,
              'ate_non_targeted' = ate_non_targeted,
              'policy_name'=policy_name,
              'data_type'=data_type))
}


