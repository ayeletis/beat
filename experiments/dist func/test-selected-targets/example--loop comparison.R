library(arrangements)
library(foreach)
library(parallel)
library(doParallel)
library(data.table)
library(progress)
source("generate_data.R")

# ================================================================ 
#            Setup Data
# ----------------------------------------------------------------
n1 = 5000  # training
n2 = 5000; # testing
p_continuous = 4
p_discrete = 1
sigma_y = 1
which_version = 3
which_tau = 0
my_seed = 1

data = generate_data(n1=n1, n2=n2, p_continuous = p_continuous,
                     p_discrete = p_discrete, sigma_y=sigma_y,
                     which_version = which_version,
                     which_tau = which_tau ,
                     seed=my_seed)
data = data$train_data
 

# ================================================================ 
#            Greedy Search
# ----------------------------------------------------------------

library(scales)
greedy_search = function(selected){
  
  ##############################################
  # functions for optimization
  tau_sum = function(x){
    # x is a binary vector of 1 or 0, same order with tau
    # stopifnot(all(x %in% c(0,1)))
    return(sum(tau_target*x))
  }
  group_mean_l2 = function(x){
    # x is a binary vector of 1 or 0, same order with tau
    # stopifnot(all(x %in% c(0,1)))
    
    # target groups
    group_1 = group_target[x==1]
    group_selected = colMeans(group_1, na.rm=TRUE)
    #setnafill(group_selected, type="const", fill=0)
    # rest: not selected + not in threshold
    group_2 = group_target[x==0]
    group_2_n = dim(group_2)[1]
    
    group_2_sum = colSums(group_2, na.rm=TRUE)
    
    group_nonselected = (group_2_sum + group_rest_sum)/(group_rest_n + group_2_n)
    
    out = (group_selected-group_nonselected)^2
    out = sqrt(sum(out, na.rm=TRUE)) # l2 norm
    out = out/dim(group_target)[2]
    return(out)
  }
  
  out_tau = tau_sum(selected)
  out_l2 = group_mean_l2(selected)
  return(list(tau = out_tau, 
              l2 = out_l2,
              N= sum(selected)))
}

## ================================================================ 
#            ROI
# ----------------------------------------------------------------
library(ROI)
main_roi=function(min_tau){
  init_tau_select = function(){
    if(group_rest_n>0){
      l2_init = sqrt(rowSums((group_target - group_rest_sum/group_rest_n)^2))
    }else{
      l2_init = sqrt(rowSums((group_target)^2))
    }
    l2_init = l2_init/dim(group_target)[2]
    
    out = tau_target - l2_init  
    # out = 1/(1+exp(-out))
    # out = rescale(out, to=c(0,1))
    out = as.integer(out>=0)  
    return(out)
  }
  
  tau_sum = function(x){
    # x is a binary vector of 1 or 0, same order with tau
    # stopifnot(all(x %in% c(0,1)))
    return(sum(tau_target * x))
  }
  group_mean_l2 = function(x){
    # x is a binary vector of 1 or 0, same order with tau
    # stopifnot(all(x %in% c(0,1)))
    
    # target groups
    group_selected = colMeans(group_target[x==1], na.rm=TRUE)
    
    # rest: not selected + not in threshold
    
    
    group_2 = group_target[x==0]
    group_2_n = dim(group_2)[1]
    group_2_sum = colSums(group_2, na.rm=TRUE)
    group_nonselected = (group_2_sum + group_rest_sum)/(group_rest_n + group_2_n)
    group_nonselected[is.na(group_nonselected)] = 0
    out = (group_selected - group_nonselected)^2
    out = sqrt(sum(out, na.rm=TRUE)) # l2 norm
    out = out/dim(group_target)[2]
    return(out)
  }
  opt_fun = function(x){
    selected = as.integer(x>=0.5)
    return(tau_sum(selected) - group_mean_l2(selected))
  }
  num_vars = length(tau_target)
  Objective = ROI::F_objective(F=opt_fun, n = num_vars)
  
  opt = ROI::OP(objective = Objective,
                maximum=TRUE, 
                bounds = V_bound(lb = rep(0.0, num_vars),
                                 ub = rep(1.0, num_vars))
  )
  # ROI_applicable_solvers(opt)
  
  
  # init_vars =  rnorm(num_vars, mean=0.5, sd=0.1) #
  # init_vars = rep(0.5, num_vars) #
  # init_vars = runif(num_vars, min=0.4, max=0.6)
  init_vars = init_tau_select()
  solution = ROI_solve(opt,
                       solver="nlminb", #"optimx",#"alabama",
                       #max_iter = 1000,
                       start =init_vars,
                       #method = "Nelder-Mead", 
                       verbose=FALSE)
  selected_targets = solution$solution
  selected_targets = as.integer(selected_targets>=0.5)
  out = data.table(tau=tau_sum(selected_targets), 
                   l2 = group_mean_l2(selected_targets),
                   N = sum(selected_targets))
  out[, criteria:= tau - l2]
  return(out)
}
 

# prepare data 
# data = fread("sample_data.csv")
setorder(data, -tau)
# dat_group = as.data.table(data[,-'tau', with=FALSE])
dat_group =  data[, .SD, .SDcols=grep("Z", names(data), value=TRUE)]
# rescale group values
dat_group = dat_group[, lapply(.SD,scales::rescale, to=c(0,1))]
names(dat_group) = paste0("group_",1:dim(dat_group)[2])
tau_vector = data$tau

# setting 
# cl = makeCluster(detectCores())
cl = makeCluster(12)
doParallel::registerDoParallel(cl)
TOTALN = 15
min_tau_vect= c(-2, -1,-0.5, 0,  0.5, 1, 2)
seed_vect = seq(1,40)

total_iter = length(min_tau_vect) * length(seed_vect)
results_greddy = vector(mode='list',length = total_iter)
results_opt = vector(mode='list',length =total_iter)
results_min_tau = vector(length = total_iter)
results_seed = vector( length = total_iter)


pb= progress_bar$new(
  format = " Remaining [:bar] :percent in :elapsed",
  total = total_iter,
  clear = FALSE, 
  width= 50)


counter = 1

for(min_tau in min_tau_vect){
  for(seed in seed_vect){
    pb$tick()
    results_min_tau[[counter]] = min_tau
    results_seed[[counter]] = seed
 
 
    target_idx = tau_vector>= min_tau
    tau_target = copy(tau_vector[target_idx] )
    group_target = copy(dat_group[target_idx] )
    
    set.seed(seed)
    s_idx = sample(1:length(tau_target), size=TOTALN, replace=FALSE)
    
    tau_target = tau_target[s_idx]
    group_target = group_target[s_idx]
    
    
    group_rest = dat_group[!target_idx]
    group_rest_sum = colSums(group_rest, na.rm=TRUE)
    group_rest_n = dim(group_rest)[1]
    
    
    # greedy search 
    all_pairs_iter = arrangements::ipermutations(k=dim(group_target)[1], v=c(0,1), replace=TRUE)
    all_results = foreach(x = all_pairs_iter,   
                          .inorder = FALSE,
                          .packages = c("data.table")) %dopar% {
                            greedy_search(x)
                          }
    
    greedy_results=rbindlist(all_results)
    greedy_results[, criteria := tau - l2]
    setorder(greedy_results, -criteria)
    results_greddy[[counter]] = greedy_results[1]
    # opt
    roi_results=main_roi(min_tau)
    results_opt[[counter]] = roi_results
    
    # counter
    counter = counter + 1
  }
}

results_greddy = rbindlist(results_greddy)
results_greddy[, type := "greedy"]
results_greddy[, min_tau := results_min_tau]
results_greddy[, seed := results_seed]

results_opt = rbindlist(results_opt)
results_opt[, type:= "optimization"]
results_opt[, min_tau := results_min_tau]
results_opt[, seed := results_seed]

results = rbind(results_opt, results_greddy)
saveRDS(results, 'results.rds')

results_w = dcast(results, min_tau + seed~type, value.var = c("tau","l2", "N", "criteria"))
results_w[, gap_rate_tau := (tau_optimization-tau_greedy )/tau_greedy ]
results_w[, gap_rate_l2 := (l2_optimization -l2_greedy  )/l2_greedy  ]
results_w[, gap_rate_N := (N_optimization   -N_greedy   )/N_greedy   ]

library(ggpubr)
results_w_m = melt(results_w, id.vars = c("min_tau","seed"), 
                   measure.vars = grep("gap", names(results_w), value=TRUE))
ggbarplot(data=results_w_m, 
          x="min_tau", 
          y='value',
          #color="variable",
          fill="variable",
          legend="none",
          add="mean_sd",
          title="Comparison -- Gap Rate: (Optimization-Greedy)/Greedy ",
          subtitle="Sample size: 15; 40 Reps Per Min Tau; Mean and SE",
          scales="free_y",
          facet.by = "variable")
ggsave("gap_rate.png", width=12, height=6)


results_rw = melt(results, id.vars = c("min_tau", "seed", "type"))
ggbarplot(data=results_rw, 
          x="min_tau", 
          y='value',
          #color="variable",
          fill="type",
          # legend="none",
          add="mean_sd",
          title="Comparison -- Gap: Optimization-Greedy",
          subtitle="Sample size: 15; 40 Reps Per Min Tau; Mean and SE",
          scales="free_y",
          position = position_dodge(width=0.4),
          facet.by = "variable")
ggsave("gap.png", width=12, height = 6)
