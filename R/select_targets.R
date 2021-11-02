## function that targets by diffrent threshold depending on whether is make/female
# finds the threshold that targets the same amount of that group as tau>0 are in the other group, or the same ratio as in the population

# The algorithm works in the following procedure:
# 1/ Calculate group ratio and normalized by the smallest group
# 2/ Select the first group (A) that has ratio of 1
# 3/ Denote target pool as P; denoted each group count as C_*
# 4/ Generate line sequence 1 to C_A
# 5/ For each group (denoted as I ) not equaling to A, do the following:
#   a/ Calculate expected count C_I’ by its ratio to A and round down the nearest to integer
#   b/ Set values to NA if C_I’ > C_I  (larger than the target group)
# 6/ Then we have a data frame, each row representing a possible combination
# 7/ Take the row that has the largest sum of all groups, so we selected the max possible top N per group
# 8/ Output the selected list


balanced_sampling = function(tau_vector, group_vector,
                             min_tau=0, target_type ='ratio'){

  # assume both vectors have the same order
  # select the same amount of targets for each group and the highest tau values
  # tau_vector: (float)  from model prediction
  # group_vector: (discrete values)
  # target_type: either 'equal' or 'ratio' that matches the same group ratio
  # output: selected target vector (1 or 0)

  stopifnot(length(tau_vector) == length(group_vector))
  stopifnot(target_type %in% c("ratio", 'equal'))
  dat_tau = data.table(group = group_vector, tau=tau_vector)
  dat_tau[, idx := 1:.N] # keep track order
  setorder(dat_tau, group, -tau)

  # count number of potential targets (tau > min tau)
  candidate_count = dat_tau[tau >min_tau, .N, by=group]

  if (target_type=='equal'){
    # set max targets per group by the min
    max_targets = candidate_count[, min(N)]

    # select targets per group from the highest tau
    # order tau by group, highest tau first
    target_idx = dat_tau[tau >min_tau, .(idx=head(idx, max_targets)), by=group]
  }else if (target_type=='ratio'){

    # get ratio per group,  baseline is the smallest group
    group_count = dat_tau[,.N, by=group]
    group_count[, ratio := N/min(N)]
    stopifnot(group_count[,min(ratio)==1])

    target_pool_count = dat_tau[tau >min_tau, .N, by=group]
    # list all possible combinations that have the same ratio, value round down
    smallest_group = group_count[ratio==1, head(group, 1)]
    smallest_group_count = target_pool_count[group==smallest_group, N]

    dat_target_rates = data.table(x= seq(1, smallest_group_count))

    for(i in group_count[group!=smallest_group, group]){
      new_group_name = paste0('group_', i)
      group_rate = group_count[group==i, ratio]
      max_available_target = target_pool_count[group==i, N]
      dat_target_rates[, y := as.integer(floor(x *group_rate))]
      dat_target_rates[y>max_available_target, y:= NA]
      setnames(dat_target_rates, 'y', new_group_name)
    }

    setnames(dat_target_rates, 'x', paste0('group_', smallest_group ))
    dat_target_rates = na.omit(dat_target_rates)
    dat_target_rates[,group_totals:= rowSums(.SD)]
    # find the count set such that total is largest
    targets_count = dat_target_rates[which.max(group_totals)]
    targets_count[, idx := 1:.N]
    targets_count[, group_totals := NULL]
    targets_count = melt(targets_count, 'idx')
    # number of targets per group
    targets_count[, group := gsub("group_", '', variable) ]


    target_idx = rbindlist(mapply(function(group_value, first_n){
      dat_tau[group==group_value, .(group=group_value,
                                    idx = head(idx, first_n))]
    },
    group_value=targets_count$group,
    first_n=targets_count$value,
    SIMPLIFY = FALSE))


  }else{
    stop("target_type must be 'ratio' or 'equal'. ")
  }


  # init a vector of 0 and set the selected idx to 1
  dat_out = data.table(idx = 1:length(tau_vector), selected = 0)
  dat_out[target_idx, on=.(idx), selected := 1]
  return(dat_out$selected)
}

# select targets that max tau and min user profile matrix dist(selected, unselected)
library(doParallel)
library(parallel)
library(arrangements)
library(foreach)
greedy_search = function(tau_vector,
                         demographic_data,
                         min_tau,
                         rescale_demographic=TRUE,
                         num_threads = -1,
                         max_target_rows=15,
                         seed=1
                         ){
  # return dataframe of sum(tau|selected), l2 norm for all possible solutions

  # helper functions
  tau_sum = function(x){return(sum(tau_target*x, na.rm=TRUE))}
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
  get_one_result = function(selected){
    out_tau = tau_sum(selected)
    out_l2 = group_mean_l2(selected)
    return(list(tau = out_tau,
                l2 = out_l2,
                N= sum(selected)))
  }

  # config
  if(num_threads==-1){
    cl = makeCluster(detectCores())
  }else if(num_threads>0){
    cl = makeCluster(num_threads)
  }else{
    stop("num_threads shall be -1 or positive")
  }

  doParallel::registerDoParallel(cl)

  # pre-process data

  ## rescale
  dat_group = as.data.table(demographic_data)
  dat_group = dat_group[, lapply(.SD, scales::rescale, to=c(0,1))]
  names(dat_group) = paste0("group_",1:dim(dat_group)[2])

  ## split data
  ### target group
  target_idx = tau_vector>= min_tau
  tau_target = copy(tau_vector[target_idx] )
  group_target = copy(dat_group[target_idx] )

  set.seed(seed)
  s_idx = sample(1:length(tau_target), size=max_target_rows, replace=FALSE)

  tau_target = tau_target[s_idx]
  group_target = group_target[s_idx]

  ### not in target group
  group_rest = dat_group[!target_idx]
  group_rest_sum = colSums(group_rest, na.rm=TRUE)
  group_rest_n = dim(group_rest)[1]


  all_pairs_iter = arrangements::ipermutations(k=dim(group_target)[1], v=c(0,1), replace=TRUE)
  all_results = foreach(x = all_pairs_iter,
                        .inorder = FALSE,
                        .packages = c("data.table")) %dopar% {
                          get_one_result(x)
                        }
  doParallel::stopCluster(cl)
  greedy_results=rbindlist(all_results)
  greedy_results[, criteria := tau - l2]
  setorder(greedy_results, -criteria)

  return(greedy_results)
}

library(ROI)
optimization_search = function(tau_vector,
                               demographic_data,
                               min_tau,
                               rescale_demographic=TRUE,
                               seed=1){
  # return binary vector, same length and order as tau_vector, 1 is selected
  stopifnot(length(tau_vector)==dim(demographic_data)[1])
  ## helper functions
  init_tau_select = function(){
    if(group_rest_n>0){
      l2_init = sqrt(rowSums((group_target - group_rest_sum/group_rest_n)^2))
    }else{
      l2_init = sqrt(rowSums((group_target)^2))
    }
    l2_init = l2_init/dim(group_target)[2]

    out = tau_target - l2_init
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

  main_roi = function(){
    num_vars = length(tau_target)
    Objective = ROI::F_objective(F=opt_fun, n = num_vars)

    opt = ROI::OP(objective = Objective,
                  maximum=TRUE,
                  bounds = V_bound(lb = rep(0.0, num_vars),
                                   ub = rep(1.0, num_vars))
    )
    init_vars = init_tau_select()
    set.seed(seed)
    solution = ROI_solve(opt,
                         solver="nlminb", #"optimx",#"alabama",
                         #max_iter = 100,
                         start =init_vars,
                         # method = "BFGS",
                         verbose=FALSE)
    selected_targets = solution$solution
    selected_targets = as.integer(selected_targets>=0.5)
    #out = data.table(tau=tau_sum(selected_targets),
    #                 l2 = group_mean_l2(selected_targets),
    #                 N = sum(selected_targets))
    #out[, criteria:= tau - l2]
    return(selected_targets)

  }

  ## rescale
  dat_group = as.data.table(demographic_data)
  dat_group = dat_group[, lapply(.SD, scales::rescale, to=c(0,1))]
  names(dat_group) = paste0("group_",1:dim(dat_group)[2])

  ## split data
  ### target group
  target_idx = tau_vector>= min_tau
  tau_target = copy(tau_vector[target_idx] )
  group_target = copy(dat_group[target_idx] )

  ### not in target group
  group_rest = dat_group[!target_idx]
  group_rest_sum = colSums(group_rest, na.rm=TRUE)
  group_rest_n = dim(group_rest)[1]
  roi_result=main_roi()

  ## so fill back in the roi_results for tau > min_tau
  target_idx = as.integer(target_idx)
  target_idx[target_idx==1] = roi_result # some are 1 but some are 0

  return(roi_result)
}

