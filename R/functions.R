#================================================================================================================
# Functions in the algorithms
#
# author: Fatih Kızılaslan (fatih.kizilaslan@medisin.uio.no)
# date: 15-January-2024
#================================================================================================================

library(glmnet)
library(lbfgs)


# scaling of the data matrix which includes some dummy variables 
scale.dummy.matrix <- function(X) {
  Xscaled = matrix(NA, nrow = nrow(X), ncol = ncol(X))
  for (j in 1:ncol(X)) {
    if (sum(X[,j]==0)> nrow(X)/2 | sum(X[,j]==1)> nrow(X)/2) {
      Xscaled[,j] = X[,j]
    } else{
      Xscaled[,j] = scale(X[,j])
    }
  }
  return(Xscaled)
}


# It is similar to the "expsafe" function in SCinCRM paper by in Shi et al. (2020)
# Xingjie Shi, Shuangge Ma, and Yuan Huang. Promoting sign consistency in the cure model estimation and
# selection. Statistical Methods in Medical Research, 29(1):15–28, 2020.
expsafe<- function(x) {
  for (i in 1:length(x)) {
    if (x[i] < -150.0) x[i] = -150;
    if (x[i] >  50.0)  x[i] = 50;  
  }
  return(x)		 
}


# lambda path based on the Cox regression of glmnet R package
lambda.path.glmnet.cox <- function(Xu, Xp, Time, Delta, nlambda, alpha_enet){
  coxfit = glmnet(x = cbind(scale.dummy.matrix(Xu), scale.dummy.matrix(Xp)) , y = Surv(time = Time, event = Delta), family = "cox",  alpha=alpha_enet, 
                  penalty.factor = c( rep(0,ncol(Xu) ), rep(1,ncol(Xp)) ), nlambda = nlambda )
  lambdapath = coxfit$lambda
  return(lambdapath)
}

# Likelihood functions in the EM algorithm

# negative of the observed log-likelihood function based on the right censoring data
observed_LogLik_survival_data <- function( par, b0, b_u, b_p, beta_u, beta_p, alpha, gammaa, X_u, X_p, Z_u, Z_p, Time, Delta){
  theta = par
  n = nrow(X_u)
  bZ = b0 + Z_u %*% b_u + Z_p %*% b_p 
  betaX = X_u %*% beta_u + X_p %*% beta_p 
  uncure.rate = exp(bZ)/(1+exp(bZ))   
  
  term.fr = 1+ ( (alpha*(Time^gammaa)*exp(betaX)) /theta )
  f.uncured = alpha*gammaa*(Time^(gammaa-1))*exp(betaX)* (term.fr^(-theta-1))
  F.uncured = 1- ( term.fr^(-theta) )
  S.pop = 1-(uncure.rate*F.uncured)
  f.pop = uncure.rate*f.uncured
  loglik = Delta* log(pmax(f.pop,rep(1e-20,n))) + (1-Delta)*log(pmax(S.pop,rep(1e-20,n))) 
  return( (-1)*sum(loglik) )
}

# Weighted likelihood functions for the penalized parameters b_p and beta_p version of the cost function for the adaptive elastic-net penalty
# Note: weight=1/par_int where par_int = the elastic-net estimates in the first step 
lc1_weighted_elastic_penalized_min <- function(par, par_int, b0, b_u, Z_u, Z_p, pir, alpha_enet, lambda_enet) {
  n = nrow(Z_p)
  n_z = ncol(Z_p)
  bT_p = par # b_tilda
  bp_nonzero = which(bT_p!=0)
  Zp_new = matrix(NA, nrow=n, ncol=n_z)
  for (i in 1:n_z) { Zp_new[,i] = Z_p[,i]*par_int[i]  }
  bZ = b0 + Z_u %*% b_u + Zp_new[,bp_nonzero, drop=FALSE] %*% bT_p[bp_nonzero]
  ll = (-1/n)* sum( pir*bZ - log(1+exp(expsafe(bZ))) ) + ( 0.5*lambda_enet*(1-alpha_enet)*sum((bT_p*par_int)^2) )   # lambda*alpha*sum(abs(bT_p)) is adding in the "lbfgs" function using "orthantwise_c"
  return(ll)
}

grad_lc1_weighted_elastic_penalized_min  <- function(par, par_int, b0, b_u, Z_u, Z_p, pir, alpha_enet, lambda_enet) {
  n = nrow(Z_p)
  n_z = ncol(Z_p)
  bT_p = par
  bp_nonzero = which(bT_p!=0)
  Zp_new = matrix(NA, nrow=n, ncol=n_z)
  for (i in 1:n_z) { Zp_new[,i] = Z_p[,i]*par_int[i]  }
  bZ = b0 + Z_u %*% b_u + Zp_new[,bp_nonzero, drop=FALSE] %*% bT_p[bp_nonzero]
  piZ = 1/ (1+exp(expsafe(-bZ)) )
  grad = (-1/n)*matrix(pir-piZ,1) %*% Zp_new + (lambda_enet*(1-alpha_enet)*(bT_p*(par_int^2)) )
  return( grad) }



lc2_weighted_elastic_penalized_min <- function(par, par_int, alpha, gammaa, beta_u, X_u, X_p, c_i, alpha_enet, lambda_enet, Time, Delta) {
  n = nrow(X_p)
  n_x = ncol(X_p)
  betaT_p = par # beta_p tilda 
  betap_nonzero = which(betaT_p!=0)
  Xp_new = matrix(NA, nrow=n, ncol=n_x)
  for (i in 1:n_x) { Xp_new[,i] = X_p[,i]*par_int[i]  }
  betaX = X_u %*% beta_u + Xp_new[,betap_nonzero, drop=FALSE] %*% betaT_p[betap_nonzero]
  ll1 = Delta *( log(alpha) + log(gammaa) + (gammaa-1)*log( pmax(Time,1e-15) ) + betaX )
  ll2 = c_i*alpha*(Time^gammaa)*exp(betaX)
  ll = (-1/n)* sum( ll1 - ll2 ) + ( 0.5*lambda_enet*(1-alpha_enet)*sum((betaT_p*par_int)^2) )   # lambda*alpha*sum(abs(bT_p)) is adding in the "lbfgs" function using "orthantwise_c"
  return(ll)
}

grad_lc2_weighted_elastic_penalized_min <- function(par, par_int, alpha, gammaa, beta_u, X_u, X_p, c_i, alpha_enet, lambda_enet, Time, Delta) {
  n = nrow(X_p)
  n_x = ncol(X_p)
  betaT_p = par
  betap_nonzero = which(betaT_p!=0)
  Xp_new = matrix(NA, nrow=n, ncol=n_x)
  for (i in 1:n_x) { Xp_new[,i] = X_p[,i]*par_int[i]  }
  betaX = X_u %*% beta_u + Xp_new[,betap_nonzero, drop=FALSE] %*% betaT_p[betap_nonzero]
  t1 = c_i*alpha*(Time^gammaa)*exp(betaX)
  gr = as.vector(rep(0, length(par)))
  Gbetaa = (Xp_new*Delta)-(Xp_new*drop(t1))
  gr = (-1/n)*colSums(Gbetaa) + ( lambda_enet*(1-alpha_enet)*(betaT_p*(par_int^2)) )
  return(gr)
}

# weighted likelihood functions for the unpenalized parameters (b0,b_u)
lc1_weighted_elastic_unpenalized_min <- function(par, par_int, bT_p, Z_u, Z_p, pir, alpha_enet, lambda_enet) {
  n = nrow(Z_p)
  n_z = ncol(Z_p)
  b0 = par[1] 
  b_u = par[2:(ncol(Z_u)+1)]
  bp_nonzero = which(bT_p!=0)
  Zp_new = matrix(NA, nrow=n, ncol=n_z)
  for (i in 1:n_z) { Zp_new[,i] = Z_p[,i]*par_int[i]  }
  bZ = b0 + Z_u %*% b_u + Zp_new[,bp_nonzero, drop=FALSE] %*% bT_p[bp_nonzero]
  ll = (-1/n)*sum( pir*bZ - log(1+exp(expsafe(bZ))) ) + ( 0.5*lambda_enet*(1-alpha_enet)*sum((bT_p*par_int)^2) ) + lambda_enet*alpha_enet*sum(abs(bT_p))
  return(ll)
}

grad_weighted_lc1_elastic_unpenalized_min <- function(par, par_int, bT_p, Z_u , Z_p, pir, alpha_enet, lambda_enet) {
  n = nrow(Z_p)
  n_z = ncol(Z_p)
  b0 = par[1] 
  b_u = par[2:(ncol(Z_u)+1)]
  bp_nonzero = which(bT_p!=0)
  Zp_new = matrix(NA, nrow=n, ncol=n_z)
  for (i in 1:n_z) { Zp_new[,i] = Z_p[,i]*par_int[i]  }
  bZ = b0 + Z_u %*% b_u + Zp_new[,bp_nonzero, drop=FALSE] %*% bT_p[bp_nonzero]
  piZ = 1/( 1+exp(expsafe(-bZ)) )
  grad1 = (-1/n)*sum(pir-piZ) ## for b0 
  grad2 = (-1/n)*matrix(pir-piZ,1) %*% Zp_new ## for b_u
  return( c(grad1, grad2) )
}

# for the parameters (log_alpha, log_gamma, beta_u)
# Note: par = (log_alpha ,log_gamma, beta_u), par_int = betap_weights, betaT_p = beta_p tilda 
lc2_weighted_elastic_unpenalized_min <- function(par, par_int, betaT_p, X_u, X_p, c_i, alpha_enet, lambda_enet, Time, Delta) {
  n = nrow(X_u)
  n_x = ncol(X_p)
  log_alpha = par[1]
  log_gammaa = par[2]
  alpha = exp(log_alpha)
  gammaa = exp(log_gammaa)
  beta_u = par[3:(ncol(X_u)+2)]
  betap_nonzero = which(betaT_p!=0)
  
  Xp_new = matrix(NA, nrow=n, ncol=n_x)
  for (i in 1:n_x) { Xp_new[,i] = X_p[,i]*par_int[i]  }
  betaX = X_u %*% beta_u + Xp_new[,betap_nonzero, drop=FALSE] %*% betaT_p[betap_nonzero]
  
  ll1 = Delta *( log_alpha + log_gammaa + (gammaa-1)*log( pmax(Time,1e-15) ) + betaX )
  ll2 = c_i*alpha*(Time^gammaa)*exp(betaX)
  ll = (-1/n)* sum( ll1 - ll2 ) + ( 0.5*lambda_enet*(1-alpha_enet)*sum((betaT_p*par_int)^2) ) + lambda_enet*alpha_enet*sum(abs(betaT_p))  
  return(ll)
}

grad_lc2_weighted_elastic_unpenalized_min <- function(par, par_int, betaT_p, X_u, X_p, c_i, alpha_enet, lambda_enet, Time, Delta) {
  n = nrow(X_u)
  n_x = ncol(X_p)
  log_alpha = par[1]
  log_gammaa = par[2]
  alpha = exp(log_alpha)
  gammaa = exp(log_gammaa)
  beta_u = par[3:(ncol(X_u)+2)]
  betap_nonzero = which(betaT_p!=0)
  
  Xp_new = matrix(NA, nrow=n, ncol=n_x)
  for (i in 1:n_x) { Xp_new[,i] = X_p[,i]*par_int[i]  }
  betaX = X_u %*% beta_u + Xp_new[,betap_nonzero, drop=FALSE] %*% betaT_p[betap_nonzero]
  
  t1 = c_i*alpha*(Time^gammaa)*exp(betaX)
  grad1 = sum(Delta-t1)
  grad2 = sum( Delta+ (Delta*gammaa*log( pmax(Time,1e-15) )) - ( t1*gammaa*log( pmax(Time,1e-15) ) ) )
  grad3 = matrix(Delta-t1,1) %*% X_u
  return( (-1/n)*c(grad1, grad2, grad3) )
}

lc_3 <- function(par, Delta, a_i, b_i) {
  n = length(Delta)
  cnst1 = ( par*log(max(par,1e-15)) )- lgamma(max(par,1e-15))
  cnst2 = sum((Delta+par-1)*b_i)-sum(a_i*par) 
  return( (-1/n)*( (n*cnst1)+cnst2 ) )
}

grad_lc_3 <-function(par, Delta, a_i, b_i) {
  n = length(Delta)
  return( (-1/n)*( n*( log(max(par,1e-15))+1-digamma(par) )+ sum(b_i-a_i) ) ) }


# negative log-likelihood functions for the adaptive elastic-net case 
lc1_elastic_penalized_min  <- function( par, b0, b_u, Z_u, Z_p, pir, alpha, lambda ) {
  n = nrow(Z_p)
  n_z = ncol(Z_p)
  b_p = par
  bp_nonzero = which(b_p!=0)
  bZ = b0 + Z_u %*% b_u + Z_p[,bp_nonzero, drop=FALSE] %*% b_p[bp_nonzero]
  ll = (-1/n)* sum(pir*bZ - log(1+exp(bZ))) + ( 0.5*lambda*(1-alpha)*sum(b_p^2) )   # lambda*alpha*sum(abs(b_p)) is adding in the "lbfgs" function using "orthantwise_c = alpha_enet*lambda_enet"
  return(ll)
}

grad_lc1_elastic_penalized_min  <- function( par, b0, b_u, Z_u, Z_p, pir, alpha, lambda ) {
  n = nrow(Z_p)
  n_z = ncol(Z_p)
  b_p = par
  bp_nonzero = which(b_p!=0)
  bZ = b0 + Z_u %*% b_u + Z_p[,bp_nonzero, drop=FALSE] %*% b_p[bp_nonzero]
  piZ = 1/(1+exp(-bZ))
  grad = (-1/n)*matrix(pir-piZ,1) %*% Z_p + (lambda*(1-alpha)*b_p)
  return( grad) }


lc1_elastic_unpenalized_min <- function( par, b_p, Z_u, Z_p, pir, alpha, lambda ) {
  b0 = par[1] 
  b_u = par[2:(ncol(Z_u)+1)]
  n = nrow(Z_p)
  bp_nonzero = which(b_p!=0)
  bZ = b0 + Z_u %*% b_u + Z_p[,bp_nonzero, drop=FALSE] %*% b_p[bp_nonzero]
  ll = (-1/n)*sum( pir*bZ - log(1+exp(bZ)) ) + (0.5*lambda*(1-alpha)*sum(b_p^2)) + lambda*alpha*sum(abs(b_p))
  return(ll)
}

grad_lc1_elastic_unpenalized_min <- function( par, b_p, Z_u , Z_p, pir, alpha, lambda ) {
  b0 = par[1] 
  b_u = par[2:(ncol(Z_u)+1)]
  n = nrow(Z_p)
  bp_nonzero = which(b_p!=0)
  bZ = b0 + Z_u %*% b_u + Z_p[,bp_nonzero, drop=FALSE] %*% b_p[bp_nonzero]
  piZ = 1/(1+exp(-bZ))
  grad1 = (-1/n)*sum(pir-piZ) # for b0 
  grad2 = (-1/n)*matrix(pir-piZ,1) %*% Z_u # for b_u
  return( c(grad1, grad2) )
}


lc2_elastic_penalized_min <- function( par, alpha, gammaa, beta_u, X_u, X_p, c_i, alpha_enet, lambda_enet, Time, Delta ) {
  n = nrow(X_p)
  n_x = ncol(X_p)
  beta_p = par
  betap_nonzero = which(beta_p!=0)
  betaX = X_u %*% beta_u + X_p[,betap_nonzero, drop=FALSE] %*% beta_p[betap_nonzero]
  ll1 = Delta *( log(alpha) + log(gammaa) + (gammaa-1)*log( pmax(Time,1e-15) ) + betaX )
  ll2 = c_i*alpha*(Time^gammaa)*exp(betaX)
  ll = (-1/n)* sum( ll1 - ll2 ) + ( 0.5*lambda_enet*(1-alpha_enet)*sum(beta_p^2) )    # lambda*alpha*sum(abs(b_p)) is adding in the "lbfgs" function using "orthantwise_c"
  return(ll)
}

grad_lc2_elastic_penalized_min <- function( par, alpha, gammaa, beta_u, X_u, X_p, c_i, alpha_enet, lambda_enet, Time, Delta ) {
  n = nrow(X_p)
  n_x = ncol(X_p)
  beta_p = par
  betap_nonzero = which(beta_p!=0)
  betaX = X_u %*% beta_u + X_p[,betap_nonzero, drop=FALSE] %*% beta_p[betap_nonzero]
  t1 = c_i*alpha*(Time^gammaa)*exp(betaX)
  gr = as.vector(rep(0, length(par)))
  Gbetaa = (X_p*Delta)-(X_p*drop(t1))
  gr = (-1/n)*colSums(Gbetaa) + ( lambda_enet*(1-alpha_enet)*beta_p )
  return(gr)
}

# Note: par=(log_alpha, log_gamma, beta_u)
lc2_elastic_unpenalized_min <- function( par, beta_p, X_u, X_p, c_i, alpha_enet, lambda_enet, Time, Delta ) {
  n = nrow(X_u)
  log_alpha = par[1]
  log_gammaa = par[2]
  alpha = exp(log_alpha)
  gammaa = exp(log_gammaa)
  beta_u = par[3:(ncol(X_u)+2)]
  betap_nonzero = which(beta_p!=0)
  betaX = X_u %*% beta_u + X_p[,betap_nonzero, drop=FALSE] %*% beta_p[betap_nonzero]
  ll1 = Delta *( log_alpha + log_gammaa + (gammaa-1)*log( pmax(Time,1e-15) ) + betaX )
  ll2 = c_i*alpha*(Time^gammaa)*exp(betaX)
  ll = (-1/n)* sum( ll1 - ll2 ) + ( 0.5*lambda_enet*(1-alpha_enet)*sum(beta_p^2) ) + lambda_enet*alpha_enet*sum(abs(beta_p))   
  return(ll)
}

grad_lc2_elastic_unpenalized_min <- function(par, beta_p, X_u, X_p, c_i, alpha_enet, lambda_enet, Time, Delta) {
  n = nrow(X_u)
  log_alpha = par[1]
  log_gammaa = par[2]
  alpha = exp(log_alpha)
  gammaa = exp(log_gammaa)
  beta_u = par[3:(ncol(X_u)+2)]
  betap_nonzero = which(beta_p!=0)
  betaX = X_u %*% beta_u + X_p[,betap_nonzero, drop=FALSE] %*% beta_p[betap_nonzero]
  t1 = c_i*alpha*(Time^gammaa)*exp(betaX)
  grad1 = sum(Delta-t1)
  grad2 = sum( Delta+ (Delta*gammaa*log( pmax(Time,1e-15) )) - ( t1*gammaa*log( pmax(Time,1e-15) ) ) )
  grad3 = matrix(Delta-t1,1) %*% X_u
  return( (-1/n)*c(grad1, grad2, grad3) )
}

## The negative log-likelihood functions for the oracle case of EM algorithm

lc1_oracle_min <- function(par, Z_u, Z_p, pir) {
  b0 = par[1] 
  b_u = par[2:(ncol(Z_u)+1)]
  b_p = par[-c(1:(ncol(Z_u)+1))] #!b_p represents non-zero coef. of general b_p!
  n = nrow(Z_p)
  bZ = b0 + Z_u %*% b_u + Z_p %*% b_p
  ll = (-1/n)*sum( pir*bZ - log(1+exp(bZ)) ) 
  return(ll)
}

grad_lc1_oracle_min <- function(par, Z_u, Z_p, pir) {
  b0 = par[1] 
  b_u = par[2:(ncol(Z_u)+1)]
  b_p = par[-c(1:(ncol(Z_u)+1))]
  n = nrow(Z_p)
  bZ = b0 + Z_u %*% b_u + Z_p %*% b_p
  piZ = 1/(1+exp(-bZ))
  grad1 = (-1/n)*sum(pir-piZ) # for b0 
  grad2 = (-1/n)*matrix(pir-piZ,1) %*% cbind(Z_u,Z_p) # for b_u and b_p
  return( c(grad1, grad2) )
}

# Note: par=(log_alpha,log_gamma,beta_u, beta_p) and  
# "beta_p" represents non-zero coefficients in this function
lc2_oracle_min <- function(par, X_u, X_p, c_i, Time, Delta) {
  n = nrow(X_u)
  log_alpha = par[1]
  log_gammaa = par[2]
  alpha = exp(log_alpha)
  gammaa = exp(log_gammaa)
  beta_u = par[3:(ncol(X_u)+2)]
  beta_p = par[-c(1:(ncol(X_u)+2))]
  betaX = X_u %*% beta_u + X_p %*% beta_p
  ll1 = Delta *( log_alpha + log_gammaa + (gammaa-1)*log( pmax(Time,1e-15) ) + betaX )
  ll2 = c_i*alpha*(Time^gammaa)*exp(betaX)
  ll = (-1/n)* sum( ll1 - ll2 )   
  return(ll)
}

grad_lc2_oracle_min <- function(par, X_u, X_p, c_i, Time, Delta) {
  #!par=(log_alpha,log_gamma,beta_u, beta_p).  !beta_p represents non-zero coef.
  n = nrow(X_u)
  log_alpha = par[1]
  log_gammaa = par[2]
  alpha = exp(log_alpha)
  gammaa = exp(log_gammaa)
  beta_u = par[3:(ncol(X_u)+2)]
  beta_p = par[-c(1:(ncol(X_u)+2))]
  betaX = X_u %*% beta_u + X_p %*% beta_p
  t1 = c_i*alpha*(Time^gammaa)*exp(betaX)
  grad1 = sum(Delta-t1)
  grad2 = sum( Delta+ (Delta*gammaa*log( pmax(Time,1e-15) )) - ( t1*gammaa*log( pmax(Time,1e-15) ) ) )
  grad3 = matrix(Delta-t1,1) %*% X_u
  return( (-1/n)*c(grad1, grad2, grad3) )
}

# The following functions "C.stat" and "Cstat_beta" are created by using the similar functions in the Han Fu's github page: https://github.com/hanfu-bios/curemodels/blob/main/evaluation.R 
# which is related to study Fu et al. (2022).

# The function of the C-statistic for the cure rate
C.stat <- function(cure_cutoff, b0_hat, b_u_hat, b_p_hat, beta_u_hat, beta_p_hat, X_u, X_p, Z_u, Z_p, testing_time, testing_delta){
  C_csw_num = 0
  C_csw_denom = 0
  testing_n = length(testing_time)
  v = rep(0, testing_n)
  y = rep(999, testing_n)
  y[testing_time>cure_cutoff] = 0
  y[testing_time<=cure_cutoff & testing_delta==1] = 1
  v[y<2] = 1
  if(all(b_p_hat==0)) {
    p_hat = 1/(1+exp(-b0_hat-Z_u %*% b_u_hat))
  } else{
    p_hat = 1/(1+exp(-b0_hat- Z_u %*% b_u_hat - Z_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0])) 
  }
  temp = v*y + (1-v)*as.vector(p_hat)
  
  if(all(beta_p_hat==0)) {
    X_beta = X_u %*% beta_u_hat 
  } else{
    X_beta = X_u %*% beta_u_hat + X_p[,beta_p_hat!=0,drop=FALSE] %*% beta_p_hat[beta_p_hat!=0]
  }
  
  for(i in 1:testing_n)
    for(j in 1:testing_n){
      if (j==i | !testing_delta[i] | testing_time[i]>testing_time[j]) next
      I_ij = testing_time[i]<testing_time[j] | (testing_time[i]==testing_time[j] & !testing_delta[j])
      if (!I_ij) next
      if (X_beta[i]>X_beta[j]) C_csw_num = C_csw_num + temp[j]
      C_csw_denom = C_csw_denom + temp[j]
    }
  return(C_csw_num / C_csw_denom)
}


Cstat_beta <- function(beta_u_hat, beta_p_hat, X_u, X_p, time, delta){
  C_csw_num = 0
  C_csw_denom = 0
  n = length(time)
  if(all(beta_p_hat==0)) {
    X_beta = X_u %*% beta_u_hat 
  } else{
    X_beta = X_u %*% beta_u_hat + X_p[,beta_p_hat!=0,drop=FALSE] %*% beta_p_hat[beta_p_hat!=0] 
  }
  for(i in 1:n)
    for(j in 1:n){
      if (j==i | !delta[i] | time[i]>time[j]) next
      I_ij = time[i]<time[j] | (time[i]==time[j] & !delta[j])
      if (!I_ij) next
      if (X_beta[i]>X_beta[j]) C_csw_num = C_csw_num + 1
      C_csw_denom = C_csw_denom + 1
    }
  return(C_csw_num / C_csw_denom)
}

# Functions for the data generation of the Weibull Mixture Cure Frailty Model (MCFM)

# logit function
logit <- function(mX, vB) {
  w<- 1/(1+exp(-mX %*% vB))
  return(w)}

# the uncured rate pi(z)=logit(Z,b) and the cure rate is "1 - pi(z)"
cure_rate <- function(mZ, vB){
  return(1-logit (mZ, vB)) }


# "covariates_gen.R" functions is created by using the "data.gener.R" function of Han Fu github page https://github.com/hanfu-bios/curemodels/blob/main/data_generation.R
# Note: n=sample size, p=covariates size, p/nTrue=block size, rho=correlation, sd=standart sapma
covariates_gen<- function(n, p, nTrue, rho, sd){
  block_sz = round(p/nTrue) 
  corr_X_p = matrix(0, p, p)
  for (i in 1:nTrue) {
    corr_X_p[(block_sz*(i-1)+1):(block_sz*i),(block_sz*(i-1)+1):(block_sz*i)] = rho^abs(outer(1:block_sz, 1:block_sz, "-")) }
  Sigma_X_p = sd^2*corr_X_p
  X_p = mvnfast::rmvn(n, mu=rep(0,p), sigma = Sigma_X_p)
  return( list( Covariates=X_p, Covariance=Sigma_X_p ) )
}

# the cumulative distribution function of Weibull distribution as in Kizilaslan et al. (2024)
p_Weibull<- function(t, alpha, gammaa) {
  return(1- exp(-alpha*(t^gammaa) )  ) }

# the quantile function of Weibull distribution as in Kizilaslan et al. (2024)
q_Weibull<- function(p, alpha, gammaa) {
  return( ( (-1/alpha)*log(1-p) )^ (1/gammaa) ) }

# the random number generation from Weibull distribution as in Kizilaslan et al. (2024)
r_Weibull<- function(n, alpha, gammaa) {
  w<-c()
  for (i in 1:n) {
    u<- runif(1)
    w[i]<-q_Weibull(u, alpha, gammaa)  }
  return( w )   }

# random number generation from Weibull mixture cure model based on the covariates X, Z, b, beta coefficients and distribution parameters
r_Mixture_W<- function(alpha, gammaa, theta, b_coef, betaa_coef, mX, mZ) {
  n<- nrow(mX)
  w<- c()
  
  for (i in 1:n) {
    u<- runif(1) 
    pi_Z<- logit(mZ[i,],b_coef)
    
    if(u<pi_Z) {
       w1<- 1-(u/pi_Z)
       w2<- (w1^(-1/theta))-1
       w3<- (theta*w2)/(alpha*exp(mX[i,] %*% betaa_coef))
       w4<- w3^(1/gammaa)
       w4<- w3^(1/gammaa)
       w5<- ifelse(w4<1e-2,1e-2,w4)
       w[i]<- w5  }
    else{w[i]<-1e+10} }
  
  return(w)  } 

# survival data generation from the Weibull mixture cure model 
r_Mixture_W_cens <- function(alpha, gammaa,  theta, b_coef, betaa_coef, mX, mZ, censor.rate) {
  n<- nrow(mX)
  survival_time <- r_Mixture_W(alpha, gammaa, theta, b_coef, betaa_coef, mX, mZ) # from the Weibull mixture cure model 
  censor_time <- rexp(n, rate=censor.rate) 
  observed_time <- pmin(survival_time, censor_time)
  censor <- as.numeric(survival_time <= censor_time)
  return( list( observed_time = observed_time, censoring = censor, censor.rate = 1-(sum(censor)/n) ) ) }


## For getting average results from the simulation

# Function for getting the average results from the "EM.out" list variable in "EM_High_Cure_Adaptive_Enet_CV_fit" 

Results.for.criteria.penMCFM <- function(EM.out, true.parameters, M, b_p, beta_p, Zu, Zp, Zp_sigma, Xp, Xp_sigma, true.uncure.rate, criteria = c("AIC", "BIC", "Cstat") ) {
  
  EM.par <- EM.bias.par <- EM.mse.par <- EM.bias.par <- b_p.EM.hat <- b_p.bias<- beta_p.EM.hat <- beta_p.bias<- uncure.hat<- bias.uncure.hat<- mse.uncure.hat <- theta.est.res<- list()
  iteration.number <- MSE_bp <- ME_bp <- nonzeros.b_p <- Sensitivity.b_p <- Specificity.b_p<- FPR.b_p<- FDP.b_p <- C.stat <- c()
  MSE_betap <- ME_betap <- nonzeros.beta_p <- Sensitivity.beta_p <- Specificity.beta_p<- FPR.beta_p<- FDP.beta_p <- c()
  AIC.validation <- BIC.validation <- Cstat.validation <- c()
  bp.binary.est_matrix <-  matrix(NA, nrow = M, ncol=length(b_p)) 
  betap.binary.est_matrix <- matrix(NA, nrow = M, ncol=length(beta_p)) 
  
  for (k in 1:M) {
    if(criteria == "AIC")   {EM = EM.out[[k]][[2]]; C.stat[k] = EM.out[[k]][[10]][EM.out[[k]][[9]][1],3];  AIC.validation[k] = EM.out[[k]][[11]][1,1];  BIC.validation[k] = EM.out[[k]][[11]][2,1];  Cstat.validation[k] = EM.out[[k]][[11]][3,1] } #C.stat for each M based on the train.estimates and test.data in CV
    if(criteria == "BIC")   {EM = EM.out[[k]][[4]]; C.stat[k] = EM.out[[k]][[10]][EM.out[[k]][[9]][2],3];  AIC.validation[k] = EM.out[[k]][[11]][1,2];  BIC.validation[k] = EM.out[[k]][[11]][2,2];  Cstat.validation[k] = EM.out[[k]][[11]][3,2] }
    if(criteria == "Cstat") {EM = EM.out[[k]][[6]]; C.stat[k] = EM.out[[k]][[10]][EM.out[[k]][[9]][3],3];  AIC.validation[k] = EM.out[[k]][[11]][1,3];  BIC.validation[k] = EM.out[[k]][[11]][2,3];  Cstat.validation[k] = EM.out[[k]][[11]][3,3] }
    
    EM.par[[k]] <- c(EM$b0, EM$b_u, EM$b_p, EM$beta_u, EM$beta_p, EM$alpha, EM$gammaa, EM$theta )               # estimates of all the unknown parameters
    EM.bias.par[[k]] <- EM.par[[k]]-true.parameters
    EM.mse.par[[k]] <- EM.bias.par[[k]]^2
    iteration.number[k] <- EM$iteration
    theta.est.res[[k]] <- c(EM$theta, EM$theta-theta, (EM$theta-theta)^2)
    
    ## for b_p coefficients
    b_p.EM.hat[[k]] <- EM$b_p
    bp.binary.est_matrix[k,] <- sapply( b_p.EM.hat[[k]], function(x){ifelse(unlist(x) == 0, 0, 1) } )
    b_p.bias[[k]] <- b_p.EM.hat[[k]]-b_p
    MSE_bp[k]<- sum( b_p.bias[[k]]^2 )
    ME_bp[k]<- drop(t(b_p.bias[[k]]) %*% Zp_sigma[[k]] %*% b_p.bias[[k]])                                       # model error of b_p with true covariance matrix
    
    nonzeros.b_p[k]<- sum( EM$b_p!=0)
    Sensitivity.b_p[k]<-  sum(EM$b_p!=0 & b_p!=0 )/ sum(b_p!=0)
    Specificity.b_p[k]<- sum( EM$b_p==0 & b_p==0 ) / sum(b_p==0)
    FPR.b_p[k] <- sum(EM$b_p!=0 & b_p==0 )/ sum(b_p==0)
    FDP.b_p[k] <- ifelse( nonzeros.b_p[k]==0, 0, sum(EM$b_p!=0 & b_p==0 )/ nonzeros.b_p[k] )                    # false discovery proportion i.e. FP/(FP+TP) 
    
    ## for beta_p coefficients
    beta_p.EM.hat[[k]] <- EM$beta_p
    betap.binary.est_matrix[k,] <- sapply( beta_p.EM.hat[[k]], function(x){ifelse(unlist(x) == 0, 0, 1) } )
    beta_p.bias[[k]] <- beta_p.EM.hat[[k]]-beta_p
    MSE_betap[k]<- sum( beta_p.bias[[k]]^2 )
    ME_betap[k]<- drop(t(beta_p.bias[[k]]) %*% Xp_sigma[[k]] %*% beta_p.bias[[k]])                              # model error of beta_p with true covariance matrix
    
    nonzeros.beta_p[k]<- sum( EM$beta_p!=0)
    Sensitivity.beta_p[k]<- sum(EM$beta_p!=0 & beta_p!=0 )/ sum(beta_p!=0)
    Specificity.beta_p[k]<- sum( EM$beta_p==0 & beta_p==0 ) / sum(beta_p==0)
    FPR.beta_p[k] <- sum(EM$beta_p!=0 & beta_p==0 )/ sum(beta_p==0)
    FDP.beta_p[k] <- ifelse( nonzeros.beta_p[k]==0, 0, sum(EM$beta_p!=0 & beta_p==0 )/ nonzeros.beta_p[k] )     # false discovery proportion i.e. FP/(FP+TP) 
    
    Z<- as.matrix(cbind(Zu[[k]],Zp[[k]]))
    uncure.hat[[k]]<- logit(Z, c(EM$b0,EM$b_u,EM$b_p) )                                                         # for each individual in the sample based on their covariates values and their estimated regression coefficients
    bias.uncure.hat[[k]]<- uncure.hat[[k]]-true.uncure.rate[[k]]
    mse.uncure.hat[[k]]<- bias.uncure.hat[[k]]^2
    
  }
  
  Par.results <- colMeans(t( matrix(unlist(EM.par),nrow=length(EM.par[[1]]) )) )
  Bias.Par.results <- colMeans(t(matrix(unlist(EM.bias.par),nrow=length(EM.par[[1]]) )) )
  MSE.Par.results <- colMeans(t(matrix(unlist(EM.mse.par),nrow=length(EM.par[[1]]) )) )
  
  Uncure.results.Table <- c( Av_True_Uncure=mean(unlist(true.uncure.rate)), Est=mean(unlist(uncure.hat)), Bias=mean(unlist(bias.uncure.hat)), MSE=mean(unlist(mse.uncure.hat))  )
  Parameter.Results.Table <- rbind( True.parameters_b_beta_Weibull_Frailty = true.parameters, Par.results, Bias.Par.results, MSE.Par.results)
  Accuracy.Results.Table.bp <- c( Number_Non_zeros = mean(nonzeros.b_p), Sensitivity = mean(Sensitivity.b_p), Specificity = mean(Specificity.b_p), FPR = mean(FPR.b_p), FDP = mean(FDP.b_p), CStat = mean(C.stat)  )
  Accuracy.Results.Table.betap <- c( Number_Non_zeros = mean(nonzeros.beta_p), Sensitivity = mean(Sensitivity.beta_p), Specificity = mean(Specificity.beta_p), FPR = mean(FPR.beta_p), FDP= mean(FDP.beta_p), CStat = mean(C.stat)  )
  
  output = list()
  # for b_p coefficients
  output[[1]] = ME_bp                                                                                          # for model error ME comparison with oracle
  output[[2]] = MSE_bp                                                                                         # for MSE comparison with oracle
  output[[3]] = nonzeros.b_p
  output[[4]] = Sensitivity.b_p
  output[[5]] = Specificity.b_p
  output[[6]] = FPR.b_p
  output[[7]] = FDP.b_p
 
  # for beta_p coefficients
  output[[8]] = ME_betap                                                                                      # for model error ME comparison with oracle
  output[[9]] = MSE_betap                                                                                     # for MSE comparison with oracle
  output[[10]] = nonzeros.beta_p
  output[[11]] = Sensitivity.beta_p
  output[[12]] = Specificity.beta_p
  output[[13]] = FPR.beta_p
  output[[14]] = FDP.beta_p
  
  output[[15]] = C.stat
  output[[16]] = cbind(unlist(lapply(true.uncure.rate, mean)), unlist(lapply(uncure.hat, mean)), unlist(lapply(bias.uncure.hat, mean)), unlist(lapply(mse.uncure.hat, mean)) ) ##average wrt individuals
  output[[17]] = Uncure.results.Table
  output[[18]] = Accuracy.Results.Table.bp
  output[[19]] = Accuracy.Results.Table.betap
  output[[20]] = Parameter.Results.Table                                                                       # for all the parameters in the model
  output[[21]] = Parameter.Results.Table[,-c(which(true.parameters==0))]                                       # for only non-zero parameters in the model
  output[[22]] = cbind(AIC.valid = AIC.validation, BIC.valid = BIC.validation, Cstat.valid = Cstat.validation )
  output[[23]] = t(matrix(unlist(theta.est.res), ncol=M, nrow = 3))                                            # frailty distribution parameter theta estimate, and its bias and mse 
  output[[24]] = bp.binary.est_matrix                                        
  output[[25]] = betap.binary.est_matrix
  
  return(output)
  
}


uncure.estimate.result <- function( b0_estimate, bu_estimate, bp_estimate, Zu, Zp, true.uncure.rate){
  
  uncure.hat<- bias.uncure.hat<- mse.uncure.hat <- list()
  
  for (k in 1:length(Zp)) {
    
    Z<- as.matrix(cbind(Zu[[k]],Zp[[k]]))
    uncure.hat[[k]]<- logit(Z, c( b0_estimate[[k]], bu_estimate[[k]], bp_estimate[[k]]) )             # for each individual in the sample using the estimated regression coefficients and covariates values
    bias.uncure.hat[[k]]<- uncure.hat[[k]]-true.uncure.rate[[k]]
    mse.uncure.hat[[k]]<- bias.uncure.hat[[k]]^2
    
  }
  
  output = list()
  # average wrt each individual
  output = cbind( unlist(lapply(true.uncure.rate, mean)), unlist(lapply(uncure.hat, mean)), unlist(lapply(bias.uncure.hat, mean)), unlist(lapply(mse.uncure.hat, mean)) ) 
  
  return(output)
}


