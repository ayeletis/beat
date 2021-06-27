library(data.table)
sigmoid = function(x){return( 1/(1+exp(-x)))}
relu = function(x){return(ifelse(x>0, x, 0))}

tau_function <- function(x,z,which_tau) {
  if(which_tau==0) {
    temp <- (pmax(x[,1], 0) + 0.5*abs(x[,3]) - x[,4] - z[,1])
  } else if (which_tau==1) {
    temp <- (pmax(x[,1], 0) + 0.5*abs(x[,3]) - x[,4] - z[,1] + 2*x[,2])
  } else if (which_tau==2) {
    temp <- (pmax(x[,1], 0) + 0.5*abs(x[,3]) - x[,4] + z[,1] - 2*x[,2])
  } else {print("ERROR .- which y tau?")}
  return(temp)
}

y_function <- function(x,z) {
  temp <- x[,1] - 2*x[,2] + x[,4] + z[,1]
  return(temp)
}


# DATA generation process -- uses the dimensions of the data, the the "version" which is the dependence between X and Z
generate_data = function(n1, n2,
                         p_continuous, p_discrete,
                         sigma_y,  which_version, which_tau, seed=my_seed){
  set.seed(seed)
  n = n1 + n2
  
  X_cont = matrix(rnorm(n*p_continuous), n, p_continuous)
  X_disc = matrix(rbinom(n*p_discrete, 1, 0.3),n,p_discrete)
  X = cbind(X_cont,X_disc)
  p = p_continuous + p_discrete
  
  # Z = Demographics
  # v3 := 4 demographic variables -- 1 of them related to X_cont[,2], two other binary and a continuous, all independent
  # v10 := like v3, with stronger correlation
  # v34 := like v3, with z1 continuous and cor(z1,x2) = 0.4
  # v38 := like v3, with z1 continuous and cor(z1,x2) = 0.8
  # (other versions I tested in the past)
  # v0 := No correlation. This is to test whether "fair" gives same fit when Z is independent from X.
  # v1 := One demographic variable correlated with X[,2], which is linear on TAU
  # v2 := Two demographic variables, only one correlated with X[,2]
  # v4 := 10 demographic variables -- 2 of them related to X_cont[,2], two other binary and a continuous, all independent
  # v5 := One demographic variable correlated with X[,3], which is "abs()" on TAU
  # v6 := One demographic variable correlated with X[,1], which is "pmax" on TAU
  # v7 := One demographic variable correlated with pmax(X[,2],0), which is linear on TAU
  # v8 := One demographic variable correlated with abs(X[,2]), which is linear on TAU
  # v9 := One demographic variable correlated with abs(X[,2]) + X[,1]
  # v32 := like v3, with z1 continuous and cor(z1,x2) = 0.2
  
  if (which_version==0) {Z = as.matrix(rbinom(n, 1, 0.2))
  } else if (which_version==1) {Z = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,2]))))
  } else if (which_version==2) {Z = cbind(as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,2])))),as.matrix(rbinom(n, 1, 0.2)))
  } else if (which_version==3) {Z = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,2]))))
  Z = cbind(Z,rbinom(n, 1, 0.2),rbinom(n, 1, 0.7),rnorm(n,0,1))
  } else if (which_version==4) {Z1 = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,2]))))
  Z2 = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,4]))))
  Z = cbind(Z1,Z2,rbinom(n, 1, 0.2),rbinom(n, 1, 0.7),rnorm(n,0,1),rbinom(n, 1, 0.2),rbinom(n, 1, 0.7),rnorm(n,0,1),rbinom(n, 1, 0.2),rbinom(n, 1, 0.7))
  } else if (which_version==5) {Z = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,3]))))
  } else if (which_version==6) {Z = as.matrix(rbinom(n, 1, 1/(1+exp(-X_cont[,1]))))
  } else if (which_version==7) {Z = as.matrix(rbinom(n, 1, 1/(1+exp(-(pmax(X_cont[,2],0))))))
  } else if (which_version==8) {Z = as.matrix(rbinom(n, 1, 1/(1+exp(-(abs(X_cont[,2]))))))
  } else if (which_version==9) {Z = as.matrix(rbinom(n, 1, 1/(1+exp(-(abs(X_cont[,2])+X[,1])))))
  } else if (which_version==10) {Z = as.matrix(rbinom(n, 1, 1/(1+exp(-6*X_cont[,2]))))
  Z = cbind(Z,rbinom(n, 1, 0.2),rbinom(n, 1, 0.7),rnorm(n,0,1))
  } else if (which_version==32) {Z = X_cont[,2] + rnorm(n,0,5)
  Z = cbind(Z,rbinom(n, 1, 0.2),rbinom(n, 1, 0.7),rnorm(n,0,1))
  } else if (which_version==34) {Z = X_cont[,2] + rnorm(n,0,2.5)
  Z = cbind(Z,rbinom(n, 1, 0.2),rbinom(n, 1, 0.7),rnorm(n,0,1))
  } else if (which_version==38) {Z = X_cont[,2] + rnorm(n,0,0.75)
  Z = cbind(Z,rbinom(n, 1, 0.2),rbinom(n, 1, 0.7),rnorm(n,0,1))
  } else {print("ERROR .- which version?")}
  
  # Random assignment
  W = rbinom(n, 1, 0.5)
  
  # Tau = Treatment effect
  tau = tau_function(X,Z,which_tau)
  
  # Y = Outcome (continuous)
  noise_y = sigma_y*runif(n)
  Y =  y_function(X,Z) + tau*W + noise_y
  
  # For now, let's assume that Z is not part of the features
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
