library(data.table)
library(grf)
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
## test
X_test = test_data[, .SD, .SDcols=Xcols]
Y_test = test_data$Y
W_test = test_data$W
Z_test = test_data[, .SD, .SDcols=Zcols]
tau_test = test_data$tau

## with demographic variables
X_train_demog = cbind(X_train,Z_train)
X_test_demog = cbind(X_test,Z_test)


AA  = c("split_l2_norm_rate", # left, right: l2 norm(colmean target weight)* penalty rate * node decrease
  "euclidean_distance_rate", # (left+right decrease) *  Euclidean distance (column mean target weight left, right ) * penalty rate
  "cosine_similarity_rate", # (left+right decrease) *  (1-cos_sim(column mean target weight left, right )) * penalty rate

  "split_l2_norm", #  sum(left,right l2 norm(colmean target weight))* penalty rate
  "euclidean_distance", #  Euclidean distance (column mean target weight left, right ) * penalty rate
  "cosine_similarity" #  (1-cos_sim(column mean target weight left, right )) * penalty rate
)

a = Sys.time()

fit_fair <- balanced_regression_forest(X_train, Y_train, W_train,
                                   target.weights = as.matrix(Z_train),
                                   target.weight.penalty = 1,
                                   target.weight.standardize = TRUE,
                                   target.weight.bins.breaks = 256,
                                   target.weight.penalty.metric = AA[3],
                                   honesty = TRUE,
                                   num.trees = num_trees,
                                   tune.parameters= "none",#my_tunable_params,
                                   seed=my_seed)
print(Sys.time()-a)
