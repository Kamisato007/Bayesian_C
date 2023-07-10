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

RiskC_Wmat1 <- function(Y,delta,tau){
  n <- length(Y)
  Wmat <- matrix(0, n, n)
  for(i in 1:n) {
    if(delta[i] != 0) {
      Wmat[i,] <- 
        (delta[i]==1)*
        ((Y[i] < Y)+(Y[i]>= Y & delta == 2))*
        (Y[i] < tau)
    }
  }
  Wmat <- Wmat/sum(Wmat)
  return(Wmat)
}

RiskC_Wmat2 <- function(Y,delta,tau,Cox_G=NULL){
  
  
  
  if (is.null(Cox_G)){
    # An index indicates whether the observation is censored
    censor = ifelse(delta==0,1,0)
    # Censoring Distribution Estimate using Kaplan-Meier Estimator
    KM_G = survfit(formula = Surv(Y ,censor) ~ 1)
    G_y = KM_G$surv[match(Y,KM_G$time)]
  }else{
    # Censoring Distribution Estimate using Kaplan-Meier Estimator
    G_y = Cox_G$surv[match(Y,Cox_G$time)]
  }
  
  # Get G(Y) for each observation 
  # (Since G(Y) in KM_G is ordered，we want each G(Y) to match original Y)
  
  
  n <- length(Y)
  Wmat <- matrix(0, n, n)
  for(i in 1:n) {
    if(delta[i] != 0) {
      Wmat[i,] <- 
        (delta[i]==1)*
        ((Y[i] < Y)*G_y[i]^(-2)+
           (Y[i]>= Y & delta == 2)*G_y[i]^(-1)*
           ifelse(G_y!=0,G_y^(-1),G_y[max(which(G_y != 0))]^(-1)))*
        (Y[i] < tau)
    }
  }
  Wmat <- Wmat/sum(Wmat)
  return(Wmat)
  
  
}

MH_CR_Sampling <- function(Y,delta,tau,
                           A,beta0,sigma0,var.prop,
                           m,B,eta,
                           Wmat_option,Cox_G=NULL){
  
  
  accept = 0
  beta = beta0
  
  # What we want to record
  BETA = matrix(0,m,dim(A)[2])
  ThetaRecord <- matrix(0, m, length(Y))
  C_stat = c()
  
  
  # For safety m>B
  if (B>m){
    B = 0
  }
  
  
  # 0 means we use Harrell C statistics
  # 1 means we use Uno C statistics
  if (Wmat_option==1){
    Wmat <- RiskC_Wmat1(Y, delta, tau)
  }
  if (Wmat_option==2){
    if (is.null(Cox_G)){
      Wmat <- RiskC_Wmat2(Y, delta, tau) 
    }else{
      Wmat <- RiskC_Wmat2(Y, delta, tau,Cox_G)
    }
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
    
    
    # Compute log of MH ratio
    
    lrMH = eta*log(HC.p) +
      sum(dnorm(beta.p,beta0,sigma0,log=T))- 
      eta*log(HC) -
      sum(dnorm(beta,beta0,sigma0,log=T)) 
    
    if (log(runif(1))<lrMH){
      beta = beta.p
      accept = accept + 1
    }
    BETA[i,] = beta
    
  }
  
  if (B == 0){
    return(list(BETA=BETA,
                accept_rate=accept/m,
                THETA = ThetaRecord,
                C_stat = C_stat))
  }else{
    return(list(BETA=BETA[-c(1:B),],
                accept_rate=accept/m,
                THETA = ThetaRecord[-c(1:B),],
                C_stat = C_stat[-c(1:B)]))
  }
  
  
}
