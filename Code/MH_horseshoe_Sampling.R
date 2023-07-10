
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

C_index <- function(theta, Wmat) {
  theta <- c(theta)
  Thetamat <- (outer(theta, theta, FUN="-") > 0)
  # Case that theta are all 0
  if (sum(Wmat*Thetamat) == 0){
    return (0.0001)
  }
  
  return(sum(Wmat*Thetamat))
}


MH_horseshoe_Sampling <- function(Y,delta,tau,
                                   A,beta0,sigma0,var.prop,
                                   m,B,eta,
                                   Wmat_option=0){
  
  
  accept_beta = 0
  accept_lambda = 0
  beta = beta0
  lambda = lambda0
  sigma = sigma0
  
  # What we want to record
  BETA = matrix(0,m,dim(A)[2])
  LAMBDA = matrix(0,m,dim(A)[2])
  V = matrix(0,m,dim(A)[2])
  ThetaRecord <- matrix(0, m, length(Y))
  C_stat = c()
  
  
  # For safety m>B
  if (B>m){
    B = 0
  }
  
  
  # 0 means we use Harrell C statistics
  # 1 means we use Uno C statistics
  if (Wmat_option==0){
    Wmat <- HarrellC_Wmat(Y, delta, tau)
  }else if (Wmat_option==1){
    Wmat <- UnoC_Wmat(Y, delta, tau)
  }else{  # Other Possible C index...
    Wmat <- HarrellC_Wmat(Y, delta, tau)
  }
  
  
  for (i in 1:m){
    
    # Sample beta from proposal distribution
    beta.p = t(rmvnorm(1,beta,var.prop))
    
    
    # Compute theta from current and last iteration
    theta.p = THETA(A,beta.p)
    theta = THETA(A,beta)
    
    # Record theta from last iteration
    ThetaRecord[i,] <- theta
    
    # Compute C-statistics from current and last iteration
    HC.p = C_index(theta.p, Wmat)
    HC = C_index(theta, Wmat)
    
    # Record C-statistics from last iteration
    C_stat = c(C_stat,HC)
    
    
    
    lrMH = eta*log(HC.p) +
      sum(dnorm(beta.p,beta0,sigma,log=T))- 
      eta*log(HC) -
      sum(dnorm(beta,beta0,sigma,log=T)) 
    
    
    
    if (log(runif(1))<lrMH){
      beta = beta.p
      accept_beta = accept_beta + 1
    }
    BETA[i,] = beta
    
    
    ##############################################################################
    ########## For sampling lambda, I try Gibbs sampling this time ###############
    lambda = sqrt(1/rgamma(dim(A)[2],shape=1,rate=(1/v)+(beta^2/(2*beta_tau^2))))
    v = 1/rgamma(dim(A)[2],shape=1,rate=1+1/lambda^2)
    sigma = lambda*beta_tau
    
    
    LAMBDA[i,] = lambda
    V[i,] = v
  }
  ######################################################################
  
  
  
  if (B == 0){
    return(list(BETA=BETA,
                LAMBDA = LAMBDA,
                V = V,
                accept_beta=accept_beta/m,
                THETA = ThetaRecord,
                C_stat = C_stat))
  }else{
    return(list(BETA=BETA[-c(1:B),],
                LAMBDA = LAMBDA[-c(1:B),],
                V = V[-c(1:B)],
                accept_beta=accept_beta/m,
                THETA = ThetaRecord[-c(1:B),],
                C_stat = C_stat[-c(1:B)]))
  }
  
  
}
