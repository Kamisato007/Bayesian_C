



library(survival)
library(mvtnorm)

# Compute theta
THETA <- function(A,beta){
  n = dim(A)[1]
  Abeta <- A%*%beta
  Abeta_bar = A%*%beta- mean(Abeta)
  theta = Abeta_bar/as.numeric(t(Abeta_bar)%*%Abeta_bar)
  theta[is.nan(theta)] = 0
  return(theta)
}

# Compute C-index type loss function
HarrellC_Wmat <- function(Y, delta, tau) {
  n <- length(Y)
  Wmat <- matrix(0, n, n)
  for(i in 1:n) {
    if(delta[i] != 0) {
      Wmat[i,] <- delta[i]*(Y[i] < Y)*(Y[i] < tau)
    }
  }
  Wmat <- Wmat/sum(Wmat)
  return(Wmat)
}


C_index <- function(theta, Wmat) {
  theta <- c(theta)
  Thetamat <- (outer(theta, theta, FUN="-") > 0)
  # Case that theta are all 0
  if (sum(Wmat*Thetamat) == 0){
    return (0.0001)
  }
  
  return(sum(Wmat*Thetamat))
}

# Wight matrix for Uno's C-statistics

UnoC_Wmat <- function(Y, delta, tau) {
  
  # An index indicates whether the observation is censored
  censor = ifelse(delta==0,1,0)
  
  # Censoring Distribution Estimate using Kaplan-Meier Estimator
  KM_G = survfit(formula = Surv(Y ,censor) ~ 1)
  
  # Get G(Y) for each observation 
  # (Since G(Y) in KM_G is ordered，we want each G(Y) to match original Y)
  G_y = KM_G$surv[match(Y,KM_G$time)]
  
  
  n <- length(Y)
  Wmat <- matrix(0, n, n)
  for(i in 1:n) {
    if(delta[i] != 0) {
      Wmat[i,] <- delta[i]*(Y[i] < Y)*(Y[i] < tau)*G_y[i]^(-2)
    }
  }
  Wmat <- Wmat/sum(Wmat)
  return(Wmat)
}




matrix_K <- function(X,lambda){
  
  return(exp(-0.5*(outer(rowSums(t(t(X)^2/lambda)),
                         rowSums(t(t(X)^2/lambda)),'+')-
                     2*(t(t(X)/sqrt(lambda)))%*%(t(X)/sqrt(lambda)))))
}











library(invgamma)
MH_GP_Sampling <- function(tti,Y,Y.test,delta,delta.test,tau,
                           A,A.all,beta0,alpha0,v0,kappa,
                           m,B,eta,K.all,n,
                           Wmat_option=0){
  
  
  accept_beta = 0
  accept_lambda = 0
  beta = beta0
  lambda = lambda0
  
  # What we want to record
  BETA = matrix(0,m,dim(A)[1])
  BETA.test = matrix(0,m,(dim(A.all)[1]-dim(A)[1]))
  LAMBDA = matrix(0,m,dim(A)[2])
  C_stat = c()
  C_stat.test = c()
  
  # For safety m>B
  if (B>m){
    B = 0
  }
  
  
  # 0 means we use Harrell C statistics
  # 1 means we use Uno C statistics
  if (Wmat_option==0){
    Wmat <- HarrellC_Wmat(Y, delta, tau)
    Wmat.test <- HarrellC_Wmat(Y.test, delta.test, tau)
  }else if (Wmat_option==1){
    Wmat <- UnoC_Wmat(Y, delta, tau)
    Wmat.test <- UnoC_Wmat(Y.test, delta.test, tau)
  }else{  # Other Possible C index...
    Wmat <- HarrellC_Wmat(Y, delta, tau)
    Wmat.test <- HarrellC_Wmat(Y.test, delta.test, tau)
  }
  
  
  for (i in 1:m){
    # Get covariance matrix for training set
    K = K.all[1:tti,1:tti]
    K.test = K.all[(tti+1):n,(tti+1):n]
    KK = K.all[1:tti,(tti+1):n]
    
    # Sample beta from proposal distribution
    beta.p = t(rmvnorm(1,beta,kappa*K))
    
    
    # Compute theta from current and last iteration
    theta.p = beta.p
    theta = beta
    
    
    
    # Compute C-statistics from current and last iteration
    
    HC.p = C_index(theta.p, Wmat)
    HC = C_index(theta, Wmat)
    
    # Record C-statistics from last iteration
    C_stat = c(C_stat,HC)
    
    
    
    #############################################
    # Compute log of MH ratio
    
    lrMH = eta*log(HC.p) +
      dmvnorm(as.numeric(beta.p),beta0,K,log=T)-
      eta*log(HC) -
      dmvnorm(as.numeric(beta),beta0,K,log=T)
    
    if (log(runif(1))<lrMH){
      beta = beta.p
      accept_beta = accept_beta + 1
    }
    BETA[i,] = beta
    #############################################
    
    # Calculate the C_index for the testing data    
    interim = t(KK)%*%solve(K)
    
    
    
    mu = interim%*%beta
    sig = K.test - interim%*%KK
    beta.test = mu
    theta.test = beta.test
    
    
    BETA.test[i,] = beta.test
    
    HC.test = C_index(theta.test,Wmat.test)
    C_stat.test = c(C_stat.test,HC.test)
    
    ###############################################################
    # Compute log of MH_lambda ratio
    lambda.p = exp(t(rnorm(dim(A)[2],log(lambda),rep(1,dim(A)[2]))))
    lambda.p = as.vector(lambda.p)
    K.p = matrix_K(A,lambda.p)
    
    lrMH_lambda = dmvnorm(as.numeric(beta),beta0,K.p,log=T)+
      sum(dinvgamma(lambda.p,alpha0,v0,log = T))  -
      dmvnorm(as.numeric(beta),beta0,K,log=T) -
      sum(dinvgamma(lambda,alpha0,v0,log = T))
    
    if (log(runif(1))<lrMH_lambda){
      lambda = lambda.p
      K.all = matrix_K(A.all,lambda)
      accept_lambda = accept_lambda + 1
    }
    LAMBDA[i,] = lambda
  }
  #####################################################################
  
  
  
  
  
  if (B == 0){
    return(list(BETA=BETA,
                BETA_test=BETA.test,
                LAMBDA = LAMBDA,
                accept_beta=accept_beta/m,
                C_stat = C_stat,
                C_stat_test = C_stat.test))
  }else{
    return(list(BETA=BETA[-c(1:B),],
                BETA_test=BETA.test[-c(1:B),],
                LAMBDA = LAMBDA[-c(1:B),],
                accept_beta=accept_beta/m,
                C_stat = C_stat[-c(1:B)],
                C_stat_test = C_stat.test[-c(1:B)]))
  }
  
  
  
}



