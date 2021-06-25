source("select_targets.R")
source("generate_data.R")
library(data.table)

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

tau_vector = data$tau
demographic_data = data[, .SD, .SDcols=grep("Z", names(data), value=TRUE)]
min_tau = 0.1

# ================================================================ 
#            Comparison
# ----------------------------------------------------------------

# greedy search-- only sample 15 rows (need 8 seconds). 17 rows need around 6 minutes 
dat_greedy = greedy_search(tau_vector=tau_vector,
                           demographic_data=demographic_data,
                           min_tau=min_tau,
                           rescale_demographic=TRUE,
                           num_threads = -1,
                           max_target_rows=15,
                           seed=1)
best_result_greedy = dat_greedy[1]

# optimization -- all data 
selected_targets = optimization_search(tau_vector=tau_vector,
                              demographic_data=demographic_data,
                              min_tau=min_tau,
                              rescale_demographic=TRUE)
selected_targets
