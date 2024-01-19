#================================================================================================================
# EM algorithm and related functions for penMCFM
#
# author: Fatih Kızılaslan (fatih.kizilaslan@medisin.uio.no)
# date: 15-January-2024
#================================================================================================================

library(glmnet)
library(lbfgs)

# EM algorithm for the penMCFM method in the Kızılaslan, Swanson and Vitelli (2024).

# X_u, X_p, Z_u, Z_p should be scaled.
# alpha_enet and lambda_enet are the parameters in the elastic-net penalty 
# tol1, tol2, tol3: tolerance values for the complete penalized log-likelihood function partitions
# alpha, gammaa: Weibull distribution parameters
# theta: Frailty distribution (Gamma dist.) parameter

EM_frailty_High_Cure_Adaptive_Enet <- function(X_u, X_p, Z_u, Z_p, Time, Delta, b0, b_u, b_p, beta_u, beta_p, alpha, gammaa, theta, bp_weights, betap_weights, bT_p ,betaT_p, alpha_enet, lambda_enet, nIter, tol1, tol2, tol3=0.01 ) {
  
  N = length(Time)
  pir = rep(1, N)
  expNbZ = denominator_part = rep(0,N)
  llp1_0 = llp2_0 = llp3_0 = 0
  step = 1
  llp1 = llp2 =  llp3 =  numeric()
  conv1 = conv2 =  conv3 = FALSE
  
  repeat{
    
    ## E-step
    
    expNbZ[Delta==0] = exp( -( b0 + Z_u[Delta==0,,drop=FALSE] %*% b_u + Z_p[Delta==0,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]) )
    denominator_part[Delta==0] = (alpha*(Time[Delta==0]^gammaa)*exp( X_u[Delta==0,,drop=FALSE] %*% beta_u + X_p[Delta==0,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0] ) )/theta
    pir[Delta==0] = 1/( 1+expNbZ[Delta==0]*(1+denominator_part[Delta==0])^theta )
    
    denominator_ac = theta + (alpha*(Time^gammaa)*exp( X_u %*% beta_u + X_p %*% beta_p) )
    c_i = pir*( (Delta+theta)/denominator_ac)
    a_i = c_i + (1-pir)*( (Delta+theta)/theta )
    b_i = ( digamma(Delta+theta)-log(pmax(denominator_ac,1e-30)) )*pir + (digamma(Delta+theta)-log(max(theta, 1e-30)) ) * (1-pir)
    
    ## update penalized parameters
    
    if (!conv1){
      
      out_Lc1_penalized = lbfgs(lc1_weighted_elastic_penalized_min, grad_lc1_weighted_elastic_penalized_min, vars = bT_p, par_int = bp_weights,
                                b0 = b0, b_u = b_u, Z_u = Z_u, Z_p = Z_p, pir = pir, alpha = alpha_enet, lambda = lambda_enet,
                                invisible=1, orthantwise_c = alpha_enet*lambda_enet,
                                orthantwise_start = 0,
                                orthantwise_end = ncol(Z_p) )
      bT_p = out_Lc1_penalized$par
      b_p = out_Lc1_penalized$par * bp_weights  ## bp_adapted = bp_weighted*|bp_int| 
      
    }
    
    if (!conv2){
      
      out_Lc2_penalized = lbfgs(lc2_weighted_elastic_penalized_min, grad_lc2_weighted_elastic_penalized_min, vars = betaT_p, par_int = betap_weights,
                                alpha=alpha, gammaa=gammaa, beta_u=beta_u, X_u = X_u, X_p = X_p, c_i=c_i, alpha_enet = alpha_enet, lambda = lambda_enet,
                                Time=Time, Delta=Delta,
                                invisible=1, orthantwise_c = alpha_enet*lambda_enet,
                                orthantwise_start = 0,
                                orthantwise_end = ncol(X_p) )
      betaT_p = out_Lc2_penalized$par
      beta_p = out_Lc2_penalized$par * betap_weights
      
    }
    
    ## update unpenalized parameters
    
    if (!conv1){
      
      out_Lc1_unpenalized = lbfgs(lc1_weighted_elastic_unpenalized_min, grad_weighted_lc1_elastic_unpenalized_min, vars = c(b0,b_u), par_int = bp_weights,
                                  bT_p = bT_p, Z_u=Z_u, Z_p = Z_p, pir = pir, alpha_enet=alpha_enet, lambda_enet=lambda_enet,
                                  invisible = 1, orthantwise_c = 0 )
      b0 = out_Lc1_unpenalized$par[1]
      b_u = out_Lc1_unpenalized$par[2:(ncol(Z_u)+1)]
      llp1_1 = -out_Lc1_unpenalized$value 
      
    }
    
    llp1 = c(llp1, llp1_1)
    
    if (!conv2){
      
      out_Lc2_unpenalized = lbfgs(lc2_weighted_elastic_unpenalized_min, grad_lc2_weighted_elastic_unpenalized_min, vars = c(log(alpha),log(gammaa),beta_u), par_int = betap_weights,
                                  betaT_p = betaT_p, X_u=X_u, X_p = X_p, c_i = c_i, alpha_enet=alpha_enet, lambda_enet=lambda_enet,
                                  Time=Time, Delta=Delta,
                                  invisible = 1, orthantwise_c = 0 )
      alpha = exp( out_Lc2_unpenalized$par[1] )
      gammaa = exp( out_Lc2_unpenalized$par[2] )
      beta_u = out_Lc2_unpenalized$par[-(1:2)]
      llp2_1 = -out_Lc2_unpenalized$value
      
    }
    
    llp2 = c(llp2, llp2_1)
    
    if (!conv3){
      
      out_Lc3 = lbfgs(lc_3, grad_lc_3, vars= theta, Delta=Delta, a_i=a_i, b_i=b_i, invisible = 1, orthantwise_c = 0 )
      theta_new = out_Lc3$par
      llp3_1 = -out_Lc3$value
      
    }
    
    llp3 = c(llp3, llp3_1)
    theta_old = theta
    theta = theta_new
    
    
    if (!conv1 & abs(llp1_1- llp1_0)< tol1) conv1 <- TRUE
    if (!conv2 & abs(llp2_1- llp2_0)< tol2) conv2 <- TRUE
    if (!conv3 & (  (abs( (theta_new-theta_old)/ ((theta_old)+1) )< tol3) | (abs(llp3_1-llp3_0)< tol3) ) ) conv3 <- TRUE
    if (step > 1 & (conv1 & conv2 & conv3) | step >= nIter) {
      break
    }
    
    llp1_0 = llp1_1
    llp2_0 = llp2_1
    llp3_0 = llp3_1
    
    step = 1 + step
    
  }  
  
  ## output
  
  return( list(b_p = b_p, b_u = b_u, b0 = b0,
               beta_u = beta_u, beta_p= beta_p, alpha = alpha, gammaa = gammaa, theta = theta,
               logLik_Lc1 = llp1[(step-2):step], logLik_Lc2 = llp2[(step-2):step], logLik_Lc3 = llp3[(step-2):step], bp_weights = bp_weights, betap_weights = betap_weights,
               iteration = step,  nonzero_bp = sum(b_p!=0),  nonzero_betap = sum(beta_p!=0) ) )
}



# Initialization of the adaptive elastic-net for weights and lambda path

EM_frailty_High_Cure_Adaptive_Enet_Initial <- function(X_u, X_p, Z_u, Z_p, Time, Delta, alpha_enet, lambda_enet, multi_step_size) {
  N = length(Time)
  n_Zp = ncol(Z_p) 
  n_Xp = ncol(X_p)
  ## Initial step all b coef = 0
  b0 = 0; b_u = rep(0,ncol(Z_u)); b_p = rep(0,n_Zp)
  ## estimates of beta_u using Cox regression with glmnet package
  fit = glmnet(x = X_u, y = Surv(time = Time, event = Delta), family = "cox", lambda = 0) 
  beta_u = as.numeric( coef(fit) )
  beta_p = rep(0,n_Xp)
  
  ## Moment estimates of Weibull distribution parameters
  gammaa = uniroot(function(a) # MOM estimate of gammaa parameter
    log(gamma(1+2/a))-2*log(gamma(1+1/a))-log(var(Time)+(mean(Time))^2)+2*log(mean(Time)),
    c(0.01,10))$root
  alpha = ( gamma(1+1/gammaa)/mean(Time) )^(1/gammaa) # MOM estimate of alpha parameter
  
  theta = 1 ## assumption for frailty distribution parameter
  pir = rep(1, N)
  
  ## 1st step with the beginning assumptions
  expNbZ = denominator_part = rep(0,N)
  expNbZ[Delta==0] = exp( -(b0+Z_u[Delta==0,,drop=FALSE] %*% b_u + Z_p[Delta==0,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]) )
  denominator_part[Delta==0] = (alpha*(Time[Delta==0]^gammaa)*exp(X_u[Delta==0,,drop=FALSE] %*% beta_u + X_p[Delta==0,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0] ) )/theta
  pir[Delta==0] = 1/( 1+expNbZ[Delta==0]*(1+denominator_part[Delta==0])^theta )
  
  denominator_ac = theta + (alpha*(Time^gammaa)*exp(X_u %*% beta_u + X_p %*% beta_p ) )
  c_i = pir*( (Delta+theta)/denominator_ac)
  a_i = c_i + (1-pir)*( (Delta+theta)/theta )
  b_i = ( digamma(Delta+theta)-log(pmax(denominator_ac,1e-30)) )*pir + (digamma(Delta+theta)-log(max(theta, 1e-30)) ) * (1-pir)
  
  ## optimization of the unpenalized and penalized variables separately based on the first values of (b0,b_u,b_p) = beta_p = 0
  init_Lc1_penalized = lbfgs(lc1_elastic_penalized_min, grad_lc1_elastic_penalized_min, vars = b_p, 
                             b0 = b0, b_u = b_u, Z_u = Z_u, Z_p = Z_p, pir = pir, alpha = alpha_enet, lambda = lambda_enet,
                             invisible=1, orthantwise_c = alpha_enet*lambda_enet,
                             orthantwise_start = 0,
                             orthantwise_end = ncol(Z_p) )
  b_p = init_Lc1_penalized$par
  
  init_Lc2_penalized = lbfgs(lc2_elastic_penalized_min, grad_lc2_elastic_penalized_min, vars = beta_p,  
                             alpha=alpha, gammaa=gammaa, beta_u=beta_u, X_u = X_u, X_p = X_p, c_i=c_i, alpha_enet = alpha_enet, lambda = lambda_enet,
                             Time=Time, Delta=Delta,
                             invisible=1, orthantwise_c = alpha_enet*lambda_enet,
                             orthantwise_start = 0,
                             orthantwise_end = ncol(X_p) )
  beta_p = init_Lc2_penalized$par 
  
  
  init_Lc1_unpenalized = lbfgs(lc1_elastic_unpenalized_min, grad_lc1_elastic_unpenalized_min, vars = c(b0,b_u), 
                               b_p = b_p, Z_u=Z_u, Z_p = Z_p, pir = pir, alpha=alpha_enet, lambda=lambda_enet,
                               invisible = 1, orthantwise_c = 0 )
  b0 = init_Lc1_unpenalized$par[1]
  b_u = init_Lc1_unpenalized$par[2:(ncol(Z_u)+1)]
  
  init_Lc2_unpenalized = lbfgs(lc2_elastic_unpenalized_min, grad_lc2_elastic_unpenalized_min, vars = c(log(alpha),log(gammaa),beta_u),
                               beta_p= beta_p, X_u= X_u, X_p=X_p, c_i = c_i, alpha_enet=alpha_enet, lambda=lambda_enet, Time = Time, Delta = Delta,
                               invisible = 1, orthantwise_c = 0 )
  alpha = exp( init_Lc2_unpenalized$par[1] )
  gammaa = exp( init_Lc2_unpenalized$par[2] )
  beta_u = init_Lc2_unpenalized$par[-(1:2)]
  
  ## 2nd step with using the previous estimates
  expNbZ = denominator_part = rep(0,N)
  expNbZ[Delta==0] = exp( -(b0+Z_u[Delta==0,,drop=FALSE] %*% b_u + Z_p[Delta==0,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]) )
  denominator_part[Delta==0] = (alpha*(Time[Delta==0]^gammaa)*exp( X_u[Delta==0,,drop=FALSE] %*% beta_u + X_p[Delta==0,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0] ) )/theta
  pir[Delta==0] = 1/( 1+expNbZ[Delta==0]*(1+denominator_part[Delta==0])^theta )
  
  denominator_ac = theta + (alpha*(Time^gammaa)*exp(X_u %*% beta_u + X_p %*% beta_p) )
  c_i = pir*( (Delta+theta)/denominator_ac)
  a_i = c_i + (1-pir)*( (Delta+theta)/theta )
  b_i = ( digamma(Delta+theta)-log(pmax(denominator_ac,1e-30)) )*pir + (digamma(Delta+theta)-log(max(theta, 1e-30)) ) * (1-pir)
  
  ## optimization of un-penalized and penalized variables separately based on the previous estimates of (b0,b_u,b_p) and beta_p
  init_Lc1_penalized = lbfgs(lc1_elastic_penalized_min, grad_lc1_elastic_penalized_min, vars = b_p, 
                             b0 = b0, b_u = b_u, Z_u = Z_u, Z_p = Z_p, pir = pir, alpha = alpha_enet, lambda = lambda_enet,
                             invisible=1, orthantwise_c = alpha_enet*lambda_enet,
                             orthantwise_start = 0,
                             orthantwise_end = ncol(Z_p) )
  b_p = init_Lc1_penalized$par
  
  init_Lc2_penalized = lbfgs(lc2_elastic_penalized_min, grad_lc2_elastic_penalized_min, vars = beta_p,  
                             alpha=alpha, gammaa=gammaa, beta_u=beta_u, X_u = X_u, X_p = X_p, c_i=c_i, alpha_enet = alpha_enet, lambda = lambda_enet,
                             Time=Time, Delta=Delta,
                             invisible=1, orthantwise_c = alpha_enet*lambda_enet,
                             orthantwise_start = 0,
                             orthantwise_end = ncol(X_p) )
  beta_p = init_Lc2_penalized$par 
  
  
  init_Lc1_unpenalized = lbfgs(lc1_elastic_unpenalized_min, grad_lc1_elastic_unpenalized_min, vars = c(b0,b_u), 
                               b_p = b_p, Z_u=Z_u, Z_p = Z_p, pir = pir, alpha=alpha_enet, lambda=lambda_enet,
                               invisible = 1, orthantwise_c = 0 )
  b0 = init_Lc1_unpenalized$par[1]
  b_u = init_Lc1_unpenalized$par[2:(ncol(Z_u)+1)]
  
  init_Lc2_unpenalized = lbfgs(lc2_elastic_unpenalized_min, grad_lc2_elastic_unpenalized_min, vars = c(log(alpha),log(gammaa),beta_u),
                               beta_p= beta_p, X_u= X_u, X_p=X_p, c_i = c_i, alpha_enet=alpha_enet, lambda=lambda_enet, Time = Time, Delta = Delta,
                               invisible = 1, orthantwise_c = 0 )
  alpha = exp( init_Lc2_unpenalized$par[1] )
  gammaa = exp( init_Lc2_unpenalized$par[2] )
  beta_u = init_Lc2_unpenalized$par[-(1:2)]
  
  bT_p = b_p ; betaT_p = beta_p
  bp_weights = rep(1,length(b_p))
  betap_weights = rep(1,length(beta_p))
  
  m_step_size = 1
  
  repeat{
    
    expNbZ[Delta==0] = exp( -(b0+Z_u[Delta==0,,drop=FALSE] %*% b_u + Z_p[Delta==0,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]) )
    denominator_part[Delta==0] = (alpha*(Time[Delta==0]^gammaa)*exp( X_u[Delta==0,,drop=FALSE] %*% beta_u + X_p[Delta==0,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0] ) )/theta
    pir[Delta==0] = 1/( 1+expNbZ[Delta==0]*(1+denominator_part[Delta==0])^theta )
    
    denominator_ac = theta + (alpha*(Time^gammaa)*exp( X_u %*% beta_u + X_p %*% beta_p ) )
    c_i = pir*( (Delta+theta)/denominator_ac)
    a_i = c_i + (1-pir)*( (Delta+theta)/theta )
    b_i = ( digamma(Delta+theta)-log(pmax(denominator_ac,1e-30)) )*pir + (digamma(Delta+theta)-log(max(theta, 1e-30)) ) * (1-pir)
    
    ## optimization of the unpenalized and penalized variables for the weighted log-likelihood functions separately based on the last estimates of values of (b0,b_u,b_p) and beta_p
    out_Lc1_penalized = lbfgs(lc1_weighted_elastic_penalized_min, grad_lc1_weighted_elastic_penalized_min, vars = bT_p, par_int = bp_weights, 
                              b0 = b0, b_u = b_u, Z_u = Z_u, Z_p = Z_p, pir = pir, alpha = alpha_enet, lambda = lambda_enet,
                              invisible=1, orthantwise_c = alpha_enet*lambda_enet,
                              orthantwise_start = 0,
                              orthantwise_end = ncol(Z_p) )
    bT_p = out_Lc1_penalized$par # solving weighted equation wrt the b_tilda
    b_p = out_Lc1_penalized$par * bp_weights  # bp_adapted = bp_weighted*|bp_int| 
    bp_weights_new = abs(b_p) # new_weights od b_p for the next step
    
    out_Lc2_penalized = lbfgs(lc2_weighted_elastic_penalized_min, grad_lc2_weighted_elastic_penalized_min, vars = betaT_p, par_int = betap_weights, 
                              alpha=alpha, gammaa=gammaa, beta_u=beta_u, X_u = X_u, X_p = X_p, c_i=c_i, alpha_enet = alpha_enet, lambda = lambda_enet,
                              Time=Time, Delta=Delta,
                              invisible=1, orthantwise_c = alpha_enet*lambda_enet,
                              orthantwise_start = 0,
                              orthantwise_end = ncol(X_p) )
    betaT_p = out_Lc2_penalized$par
    beta_p = out_Lc2_penalized$par * betap_weights
    betap_weights_new = abs(beta_p) # new_weights of beta_p for the next step
    
    ## for the unpenalized parameters
    out_Lc1_unpenalized = lbfgs(lc1_weighted_elastic_unpenalized_min, grad_weighted_lc1_elastic_unpenalized_min, vars = c(b0,b_u), par_int = bp_weights,
                                bT_p = bT_p, Z_u = Z_u, Z_p = Z_p, pir = pir, alpha_enet = alpha_enet, lambda_enet = lambda_enet,
                                invisible = 1, orthantwise_c = 0 )
    b0 = out_Lc1_unpenalized$par[1]
    b_u = out_Lc1_unpenalized$par[2:(ncol(Z_u)+1)]
    
    out_Lc2_unpenalized = lbfgs(lc2_weighted_elastic_unpenalized_min, grad_lc2_weighted_elastic_unpenalized_min, vars = c(log(alpha),log(gammaa),beta_u), par_int = betap_weights,
                                betaT_p = betaT_p, X_u = X_u, X_p = X_p, c_i = c_i, alpha_enet = alpha_enet, lambda_enet = lambda_enet, 
                                Time = Time, Delta = Delta,
                                invisible = 1, orthantwise_c = 0 )
    alpha = exp( out_Lc2_unpenalized$par[1] )
    gammaa = exp( out_Lc2_unpenalized$par[2] )
    beta_u = out_Lc2_unpenalized$par[-(1:2)]
    
    bp_weights = bp_weights_new
    betap_weights = betap_weights_new
    
    m_step_size = m_step_size + 1
    if (m_step_size == multi_step_size){
      break
    }
  } 
  
  ## output 
  return( list(b_p = b_p, b_u = b_u, b0 = b0, beta_u = beta_u, beta_p = beta_p, alpha = alpha, gammaa = gammaa, theta = theta, bp_weights = bp_weights, betap_weights = betap_weights, bT_p = bT_p, betaT_p = betaT_p) ) 
  
}


# Cross-validation using the EM algorithm

EM_High_Cure_Adaptive_Enet_CV_fit <- function(X_u, X_p, Z_u, Z_p, Time, Delta, alpha_enet, nIter, n_lambda, tol, n_folds){
  
  lambda_list = lambda.path.glmnet.cox(Xu = scale.dummy.matrix(X_u), Xp=scale.dummy.matrix(X_p), Time = Time, Delta = Delta, nlambda=n_lambda, alpha_enet = alpha_enet)
  
  set.seed(data.seed)
  ## validation data: %20 of the data
  validation_i = sample(length(Time), size = length(Time)/5) 
  validation_Delta = Delta[validation_i]
  validation_Time = Time[validation_i]
  validation_X_u = scale.dummy.matrix(X_u[validation_i,])
  validation_X_p = scale.dummy.matrix(X_p[validation_i,])
  validation_Z_u = scale.dummy.matrix(Z_u[validation_i,]) 
  validation_Z_p = scale.dummy.matrix(Z_p[validation_i,]) 
  
  ## removing the validation data from the whole data 
  X_u_CV = X_u[-validation_i,]
  X_p_CV = X_p[-validation_i,]
  Z_u_CV = Z_u[-validation_i,]
  Z_p_CV = Z_p[-validation_i,]
  Time_CV = Time[-validation_i]
  Delta_CV = Delta[-validation_i]
  
  ## create folds for CV based on  X_u_CV, Z_u_CV,... , Delta_CV 
  set.seed(data.seed)  
  folds_i = sample(rep(1:n_folds, length.out = length(Time_CV)))
  
  train_out =  matrix(list(NA), nrow = n_lambda, ncol = n_folds )
  C_matrix = aic = bic = Obs_Loglik = iteration = non.zeros.bp =  non.zeros.betap = theta = lambda.update = matrix( NA, nrow = n_lambda, ncol = n_folds  )
  
  for (k in 1:n_folds ) {
    
    ## test data for kth fold
    test_i = which(folds_i == k)
    test_Delta = Delta_CV[test_i]
    test_Time = Time_CV[test_i]
    test_X_u = scale.dummy.matrix(X_u_CV[test_i,])
    test_X_p = scale.dummy.matrix(X_p_CV[test_i,]) 
    test_Z_u = scale.dummy.matrix(Z_u_CV[test_i,])
    test_Z_p = scale.dummy.matrix(Z_p_CV[test_i,]) 
    
    for (j in 1:n_lambda) {
      
      adaptive.initial = EM_frailty_High_Cure_Adaptive_Enet_Initial(X_u=scale.dummy.matrix(X_u_CV[-test_i,]), X_p=scale.dummy.matrix(X_p_CV[-test_i,]),  Z_u=scale.dummy.matrix(Z_u_CV[-test_i,]), Z_p=scale.dummy.matrix(Z_p_CV[-test_i,]),
                                                                    Time=Time_CV[-test_i], Delta=Delta_CV[-test_i], alpha_enet = alpha_enet , lambda_enet = lambda_list[j], 
                                                                    multi_step_size = 3) 
      
      train_out[[j,k]] = try ( EM_frailty_High_Cure_Adaptive_Enet(X_u=scale.dummy.matrix(X_u_CV[-test_i,]), X_p=scale.dummy.matrix(X_p_CV[-test_i,]),  Z_u=scale.dummy.matrix(Z_u_CV[-test_i,]), Z_p=scale.dummy.matrix(Z_p_CV[-test_i,]),
                                                                  Time=Time_CV[-test_i], Delta=Delta_CV[-test_i],
                                                                  b0 = adaptive.initial$b0 , b_u = adaptive.initial$b_u, b_p = adaptive.initial$b_p , beta_u = adaptive.initial$beta_u, beta_p = adaptive.initial$beta_p,
                                                                  alpha = adaptive.initial$alpha, gammaa = adaptive.initial$gammaa, theta = adaptive.initial$theta, 
                                                                  bp_weights = adaptive.initial$bp_weights, betap_weights = adaptive.initial$betap_weights,
                                                                  bT_p = adaptive.initial$bT_p, betaT_p = adaptive.initial$betaT_p, 
                                                                  alpha_enet = alpha_enet, lambda_enet = lambda_list[j],
                                                                  nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10 )  )
      lambda.update[j,k] = 0
      
      if (class( train_out[[j,k]] ) == "try-error") {
        adaptive.initial = EM_frailty_High_Cure_Adaptive_Enet_Initial(X_u=scale.dummy.matrix(X_u_CV[-test_i,]), X_p=scale.dummy.matrix(X_p_CV[-test_i,]),  Z_u=scale.dummy.matrix(Z_u_CV[-test_i,]), Z_p=scale(Z_p_CV[-test_i,]),
                                                                      Time=Time_CV[-test_i], Delta=Delta_CV[-test_i], alpha_enet = alpha_enet , lambda_enet = lambda_list[j]-2e-04,
                                                                      multi_step_size = 3) 
        train_out[[j,k]] = try ( EM_frailty_High_Cure_Adaptive_Enet(X_u=scale.dummy.matrix(X_u_CV[-test_i,]), X_p=scale.dummy.matrix(X_p_CV[-test_i,]),  Z_u=scale.dummy.matrix(Z_u_CV[-test_i,]), Z_p=scale(Z_p_CV[-test_i,]),
                                                                    Time=Time_CV[-test_i], Delta=Delta_CV[-test_i],
                                                                    b0 = adaptive.initial$b0 , b_u = adaptive.initial$b_u, b_p = adaptive.initial$b_p , beta_u = adaptive.initial$beta_u, beta_p = adaptive.initial$beta_p,
                                                                    alpha = adaptive.initial$alpha, gammaa = adaptive.initial$gammaa, theta = adaptive.initial$theta, 
                                                                    bp_weights = adaptive.initial$bp_weights, betap_weights = adaptive.initial$betap_weights,
                                                                    bT_p = adaptive.initial$bT_p, betaT_p = adaptive.initial$betaT_p, 
                                                                    alpha_enet = alpha_enet, lambda_enet = lambda_list[j] - 2e-04, 
                                                                    nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  )  )
        if ( class( train_out[[j,k]] ) != "try-error") lambda.update[j,k] = 1
      }
      
      if (class( train_out[[j,k]] ) == "try-error") {
        adaptive.initial = EM_frailty_High_Cure_Adaptive_Enet_Initial(X_u=scale.dummy.matrix(X_u_CV[-test_i,]), X_p=scale.dummy.matrix(X_p_CV[-test_i,]),  Z_u=scale.dummy.matrix(Z_u_CV[-test_i,]), Z_p=scale(Z_p_CV[-test_i,]),
                                                                      Time=Time_CV[-test_i], Delta=Delta_CV[-test_i], alpha_enet = alpha_enet , lambda_enet = lambda_list[j]+2e-04, 
                                                                      multi_step_size = 3) 
        train_out[[j,k]] = try ( EM_frailty_High_Cure_Adaptive_Enet(X_u=scale.dummy.matrix(X_u_CV[-test_i,]), X_p=scale.dummy.matrix(X_p_CV[-test_i,]),  Z_u=scale.dummy.matrix(Z_u_CV[-test_i,]), Z_p=scale(Z_p_CV[-test_i,]),
                                                                    Time=Time_CV[-test_i], Delta=Delta_CV[-test_i],
                                                                    b0 = adaptive.initial$b0 , b_u = adaptive.initial$b_u, b_p = adaptive.initial$b_p , beta_u = adaptive.initial$beta_u, beta_p = adaptive.initial$beta_p,
                                                                    alpha = adaptive.initial$alpha, gammaa = adaptive.initial$gammaa, theta = adaptive.initial$theta, 
                                                                    bp_weights = adaptive.initial$bp_weights, betap_weights = adaptive.initial$betap_weights,
                                                                    bT_p = adaptive.initial$bT_p, betaT_p = adaptive.initial$betaT_p, 
                                                                    alpha_enet = alpha_enet, lambda_enet = lambda_list[j] + 2e-04, 
                                                                    nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  )  )
        if ( class( train_out[[j,k]] ) != "try-error") lambda.update[j,k] = 1
      }
      
      if (class( train_out[[j,k]] ) == "try-error") {
        if ( class( train_out[[j,k]] ) != "try-error") lambda.update[j,k] = -1
        next
      }
      
      ## calculate the C-stat, AIC, BIC based on the test data using the estimates from train data
      C_matrix[j,k] = C.stat( cure_cutoff = 5, b0_hat=train_out[[j,k]]$b0, b_u_hat=train_out[[j,k]]$b_u, b_p_hat=train_out[[j,k]]$b_p, beta_u_hat=train_out[[j,k]]$beta_u, beta_p_hat=train_out[[j,k]]$beta_p,
                              X_u=test_X_u, X_p=test_X_p, Z_u=test_Z_u, Z_p=test_Z_p, testing_time=test_Time, testing_delta=test_Delta )
      
      Obs_Loglik[j,k] = -observed_LogLik_survival_data(par=train_out[[j,k]]$theta, b0=train_out[[j,k]]$b0, b_u=train_out[[j,k]]$b_u, b_p=train_out[[j,k]]$b_p, beta_u=train_out[[j,k]]$beta_u, beta_p=train_out[[j,k]]$beta_p, 
                                                       alpha=train_out[[j,k]]$alpha, gammaa=train_out[[j,k]]$gammaa, X_u=test_X_u, X_p=test_X_p, Z_u=test_Z_u, Z_p=test_Z_p,
                                                       Time=test_Time, Delta=test_Delta) #based on the test data using train estimates
      aic[j,k] = -2*Obs_Loglik[j,k] + 2*(train_out[[j,k]]$nonzero_bp + train_out[[j,k]]$nonzero_betap)
      bic[j,k] = -2*Obs_Loglik[j,k] + ( log(length(test_Time))* ( train_out[[j,k]]$nonzero_bp + train_out[[j,k]]$nonzero_betap ) )
      theta[j,k] = train_out[[j,k]]$theta
      iteration[j,k] = train_out[[j,k]]$iteration
      non.zeros.bp[j,k] = train_out[[j,k]]$nonzero_bp
      non.zeros.betap[j,k] = train_out[[j,k]]$nonzero_betap
      cat("completed step:", c(j,k),"\n")
      
    } 
    
    cat("Fold", k, "training finished\n")
    
  } 
  
  ## lambda selection based on C-stat criteria
  C_matrix.na.omit = na.omit(C_matrix)
  model_select_Cstat = which.max( rowMeans(C_matrix.na.omit, na.rm = T) ) 
  cat("Selected lambda using C_stat:", round( lambda_list[model_select_Cstat], 5),"\n" )  
  
  ## lambda selection based on AIC criteria
  aic.na.omit = na.omit(aic)
  model_select_aic = which.min( rowMeans(aic.na.omit, na.rm = T) ) 
  cat("Selected lambda using AIC:", round( lambda_list[model_select_aic], 5), "\n" )
  
  ## lambda selection based on BIC criteria
  bic.na.omit = na.omit(bic)
  model_select_bic = which.min( rowMeans(bic.na.omit, na.rm = T) ) 
  cat("Selected lambda using BIC:", round( lambda_list[model_select_bic], 5),"\n" )
  
  Av.logLik =  rowMeans( Obs_Loglik, na.rm = T) 
  
  ## Results based on the "train data + test data = all_data - validation_data" for the selected tuning parameter
  
  
  ## for selected "lambda_Cstat"
  
  adaptive.initial.Cstat = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                       Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                       alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_Cstat], multi_step_size = 3) 
  
  output.Cstat = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                          Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                          b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                          beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                          alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                          bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                          bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                          alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_Cstat],
                                                          nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )
  
  if (class( output.Cstat ) == "try-error") {
    
    adaptive.initial.Cstat = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                         Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                         alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_Cstat] - 2e-04, multi_step_size = 3) 
    
    output.Cstat = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                            Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                            b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                            beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                            alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                            bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                            bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                            alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_Cstat] - 2e-04,
                                                            nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )
    
  }
  
  if (class( output.Cstat ) == "try-error") {
    
    adaptive.initial.Cstat = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                         Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                         alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_Cstat] + 2e-04,  multi_step_size = 3) 
    
    output.Cstat = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                            Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                            b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                            beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                            alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                            bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                            bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                            alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_Cstat] + 2e-04, 
                                                            nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )
    
  }
  
  ## tuning alpha parameter based on the validation data
  
  Cstat_model_select_Cstat = C.stat( cure_cutoff = 5, b0_hat=output.Cstat$b0, b_u_hat=output.Cstat$b_u, b_p_hat=output.Cstat$b_p, beta_u_hat=output.Cstat$beta_u, beta_p_hat=output.Cstat$beta_p,
                                     X_u=scale.dummy.matrix(validation_X_u), X_p=scale.dummy.matrix(validation_X_p), Z_u=scale.dummy.matrix(validation_Z_u), Z_p=scale.dummy.matrix(validation_Z_p),
                                     testing_time=validation_Time, testing_delta=validation_Delta ) 
  
  Obs_Loglik_validation_model_select_Cstat = -observed_LogLik_survival_data(par=output.Cstat$theta, b0=output.Cstat$b0, b_u=output.Cstat$b_u, b_p=output.Cstat$b_p, beta_u=output.Cstat$beta_u, beta_p=output.Cstat$beta_p, 
                                                                            alpha=output.Cstat$alpha, gammaa=output.Cstat$gammaa, X_u=scale.dummy.matrix(validation_X_u), X_p=scale.dummy.matrix(validation_X_p), 
                                                                            Z_u=scale.dummy.matrix(validation_Z_u), Z_p=scale.dummy.matrix(validation_Z_p), Time=validation_Time, Delta=validation_Delta)
  AIC_model_select_Cstat = -2*Obs_Loglik_validation_model_select_Cstat + 2*(output.Cstat$nonzero_bp + output.Cstat$nonzero_betap) 
  BIC_model_select_Cstat = -2*Obs_Loglik_validation_model_select_Cstat + ( log(length(validation_Time))*(output.Cstat$nonzero_bp + output.Cstat$nonzero_betap) ) 
  
  cat( "output.Cstat completed\n")
  
  
  
  ## for selected "lambda_aic"
  
  
  adaptive.initial.aic = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                     Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                     alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_aic], multi_step_size = 3) 
  
  output.aic = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                        Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                        b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                        beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                        alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                        bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                        bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                        alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_aic], 
                                                        nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )
  
  if (class( output.aic ) == "try-error") {
    
    adaptive.initial.aic = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                       Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                       alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_aic] - 2e-04, multi_step_size = 3) 
    
    output.aic = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                          Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                          b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                          beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                          alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                          bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                          bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                          alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_aic] - 2e-04, 
                                                          nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )
    
  }
  
  if (class( output.aic ) == "try-error") {
    
    adaptive.initial.aic = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                       Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                       alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_aic] + 2e-04, multi_step_size = 5) 
    
    output.aic = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                          Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                          b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                          beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                          alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                          bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                          bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                          alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_aic] + 2e-04,
                                                          nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )    
    
  }
  
  ## tuning alpha parameter based on the validation data
  
  Cstat_model_select_aic = C.stat( cure_cutoff = 5, b0_hat=output.aic$b0, b_u_hat=output.aic$b_u, b_p_hat=output.aic$b_p, beta_u_hat=output.aic$beta_u, beta_p_hat=output.aic$beta_p,
                                   X_u=scale.dummy.matrix(validation_X_u), X_p=scale.dummy.matrix(validation_X_p), Z_u=scale.dummy.matrix(validation_Z_u), Z_p=scale(validation_Z_p), 
                                   testing_time=validation_Time, testing_delta=validation_Delta ) 
  
  Obs_Loglik_validation_model_select_aic = -observed_LogLik_survival_data(par=output.aic$theta, b0=output.aic$b0, b_u=output.aic$b_u, b_p=output.aic$b_p, beta_u=output.aic$beta_u, beta_p=output.aic$beta_p, 
                                                                          alpha=output.aic$alpha, gammaa=output.aic$gammaa, X_u=scale.dummy.matrix(validation_X_u), X_p=scale.dummy.matrix(validation_X_p),
                                                                          Z_u=scale.dummy.matrix(validation_Z_u),Z_p=scale.dummy.matrix(validation_Z_p), Time=validation_Time, Delta=validation_Delta)
  AIC_model_select_aic = -2*Obs_Loglik_validation_model_select_aic + 2*(output.aic$nonzero_bp + output.aic$nonzero_betap)  
  BIC_model_select_aic = -2*Obs_Loglik_validation_model_select_aic + ( log(length(validation_Time))*(output.aic$nonzero_bp+output.aic$nonzero_betap) )  
  
  cat( "output.aic completed\n")
  
  ## for selected "lambda_bic"
  
  
  adaptive.initial.bic = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                     Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                     alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_bic], multi_step_size = 3) 
  
  output.bic = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                        Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                        b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                        beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                        alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                        bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                        bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                        alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_bic],
                                                        nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )
  
  if (class( output.bic ) == "try-error") {
    
    adaptive.initial.bic = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                       Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                       alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_bic] - 2e-04, multi_step_size = 3) 
    
    output.bic = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                          Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                          b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                          beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                          alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                          bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                          bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                          alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_bic] - 2e-04, 
                                                          nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )
  }
  
  if (class( output.bic ) == "try-error") {
    
    adaptive.initial.bic = EM_frailty_High_Cure_Adaptive_Enet_Initial( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]),  Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                                       Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                                       alpha_enet = alpha_enet , lambda_enet = lambda_list[model_select_bic] + 2e-04, multi_step_size = 3) 
    
    output.bic = try( EM_frailty_High_Cure_Adaptive_Enet( X_u=scale.dummy.matrix(X_u[-validation_i,]),  X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), Z_p=scale.dummy.matrix(Z_p[-validation_i,]),
                                                          Time=Time[-validation_i], Delta=Delta[-validation_i],
                                                          b0 = adaptive.initial.Cstat$b0 , b_u = adaptive.initial.Cstat$b_u, b_p = adaptive.initial.Cstat$b_p , 
                                                          beta_u = adaptive.initial.Cstat$beta_u, beta_p = adaptive.initial.Cstat$beta_p,
                                                          alpha = adaptive.initial.Cstat$alpha, gammaa = adaptive.initial.Cstat$gammaa, theta = adaptive.initial.Cstat$theta, 
                                                          bp_weights = adaptive.initial.Cstat$bp_weights, betap_weights = adaptive.initial.Cstat$betap_weights,
                                                          bT_p = adaptive.initial.Cstat$bT_p, betaT_p = adaptive.initial.Cstat$betaT_p,
                                                          alpha_enet = alpha_enet, lambda_enet = lambda_list[model_select_bic] + 2e-04, 
                                                          nIter = nIter, tol1 = tol, tol2 = tol, tol3 = tol*10  ) )
  }
  
  
  ## tuning alpha parameter based on the validation data
  
  Cstat_model_select_bic = C.stat( cure_cutoff = 5, b0_hat=output.bic$b0, b_u_hat=output.bic$b_u, b_p_hat=output.bic$b_p, beta_u_hat=output.bic$beta_u, beta_p_hat=output.bic$beta_p,
                                   X_u=scale.dummy.matrix(validation_X_u), X_p=scale.dummy.matrix(validation_X_p), Z_u=scale.dummy.matrix(validation_Z_u), Z_p=scale.dummy.matrix(validation_Z_p), 
                                   testing_time=validation_Time, testing_delta=validation_Delta)
  
  Obs_Loglik_validation_model_select_bic = -observed_LogLik_survival_data(par=output.bic$theta, b0=output.bic$b0, b_u=output.bic$b_u, b_p=output.bic$b_p, beta_u=output.bic$beta_u, beta_p=output.bic$beta_p, 
                                                                          alpha=output.bic$alpha, gammaa=output.bic$gammaa, X_u=scale.dummy.matrix(validation_X_u), X_p=scale.dummy.matrix(validation_X_p),
                                                                          Z_u=scale.dummy.matrix(validation_Z_u), Z_p=scale.dummy.matrix(validation_Z_p), Time=validation_Time, Delta=validation_Delta)
  AIC_model_select_bic = -2*Obs_Loglik_validation_model_select_bic + 2*(output.bic$nonzero_bp + output.bic$nonzero_betap)  
  
  BIC_model_select_bic = -2*Obs_Loglik_validation_model_select_bic + ( log(length(validation_Time))*(output.bic$nonzero_bp + output.bic$nonzero_betap) )  
  
  cat( "output.bic completed\n")
  
  output = list()
  
  output[[1]] = cbind(lambda = lambda_list, aic, Average = rowMeans(aic,na.rm = TRUE) )
  output[[2]] = output.aic # based on the train data + test data sets 
  
  output[[3]] = cbind(lambda = lambda_list, bic, Average = rowMeans(bic,na.rm = TRUE) )
  output[[4]] = output.bic #based on the train data + test data sets 
  
  output[[5]] = cbind(lambda = lambda_list, C_matrix, Average = rowMeans(C_matrix,na.rm = TRUE) )
  output[[6]] = output.Cstat #based on the train data + test data sets 
  
  output[[7]] = lambda_list
  output[[8]] = lambda.update
  
  output[[9]] = c(model_select_aic, model_select_bic, model_select_Cstat)
  output[[10]] = cbind( AIC_Average = rowMeans(aic,na.rm = TRUE),  BIC_Average = rowMeans(bic,na.rm = TRUE), Cstat_Average = rowMeans(C_matrix,na.rm = TRUE), ObsLogLik_Average = Av.logLik ) #Average AIC, BIC, Cstat for each lambda based on test data from CV
  
  ## AIC, BIC and C-stat values for the selected lambda based on the validation data
  
  output[[11]] = rbind( AIC.validation = c( AIC_model_select_aic, AIC_model_select_bic,  AIC_model_select_Cstat ),
                        BIC.validation = c( BIC_model_select_aic, BIC_model_select_bic,  BIC_model_select_Cstat ),
                        C.stat.validation = c( Cstat_model_select_aic, Cstat_model_select_bic, Cstat_model_select_Cstat ) ) 
  colnames(output[[11]]) =  c("selected AIC model", "selected BIC model", "selected C-stat model") 
  
  output[[12]] = non.zeros.bp
  output[[13]] = non.zeros.betap
  output[[14]] = Obs_Loglik
  
  return(output)
  
}


# EM algorithm for the oracle case
# Oracle case: X_p and Z_p are the variables which have non-zero coefficients!!!!
# Note: X_u, X_p, Z_u, Z_p should be scaled.

EM_frailty_Cure_Oracle <- function( X_u, X_p, Z_u, Z_p, Time, Delta, nIter, tol1, tol2, tol3=0.01 ) {
  N = length(Time)
  n_Zu = ncol(Z_u); n_Zp = ncol(Z_p) 
  n_Xu = ncol(X_u); n_Xp = ncol(X_p)

  ## Initial step all b coef = 0
  b0 = 0; b_u = rep(0,n_Zu); b_p = rep(0,n_Zp)
  fit = glmnet(x = X_u, y = Surv(time = Time, event = Delta), family = "cox", lambda = 0) 
  beta_u = as.numeric( coef(fit) ); beta_p = rep(0,n_Xp)
  
  gammaa = uniroot(function(a) # MOM estimate of gammaa parameter 
    log(gamma(1+2/a))-2*log(gamma(1+1/a))-log(var(Time)+(mean(Time))^2)+2*log(mean(Time)),
    c(0.01,10))$root
  log_gammaa = log(max(gammaa, 1e-15))
  alpha = ( gamma(1+1/gammaa)/mean(Time) )^(1/gammaa) # MOM estimate of alpha parameter 
  log_alpha = log(max(alpha, 1e-15))
  
  theta = 1 # starting value of the frailty distribution parameter
  
  pir = rep(1, N)
  expNbZ = denominator_part = rep(0,N)
  expNbZ[Delta==0] =  exp( -(b0+Z_u[Delta==0,,drop=FALSE] %*% b_u + Z_p[Delta==0,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]) )
  denominator_part[Delta==0] = (alpha*(Time[Delta==0]^gammaa)*exp( X_u[Delta==0,,drop=FALSE] %*% beta_u + X_p[Delta==0,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]) )/theta
  pir[Delta==0] = 1/( 1+expNbZ[Delta==0]*(1+denominator_part[Delta==0])^theta )
  
  denominator_ac = theta + (alpha*(Time^gammaa)*exp(X_u %*% beta_u + X_p %*% beta_p ) )
  c_i = pir*( (Delta+theta)/denominator_ac)
  a_i = c_i + (1-pir)*( (Delta+theta)/theta )
  b_i = ( digamma(Delta+theta)-log(pmax(denominator_ac,1e-100)) )*pir + (digamma(Delta+theta)-log(max(theta, 1e-100)) ) * (1-pir)
  
  init_fit_lc1  = lbfgs(lc1_oracle_min, grad_lc1_oracle_min, vars = c(b0,b_u,b_p),
                        Z_u = Z_u, Z_p = Z_p, pir = pir,
                        invisible = 1, orthantwise_c = 0 )
  b0 = init_fit_lc1$par[1]
  b_u = init_fit_lc1$par[2:(n_Zu+1)]
  b_p = init_fit_lc1$par[-c(1:(n_Zu+1))]
  llp1_0 = -init_fit_lc1$value
  
  init_fit_Lc2 = optim(fn=lc2_oracle_min, gr=grad_lc2_oracle_min, par = c(log(alpha),log(gammaa),beta_u,beta_p),
                       X_u = X_u, X_p = X_p, c_i = c_i, Time = Time, Delta = Delta) 
  
  if( is.na(init_fit_Lc2$value) == TRUE ) {
    init_fit_Lc2 = lbfgs(lc2_oracle_min, grad_lc2_oracle_min, vars = c(log(alpha),log(gammaa),beta_u,beta_p),
                         X_u = X_u, X_p = X_p, c_i = c_i, Time = Time, Delta = Delta,
                         invisible = 1, orthantwise_c = 0 )
  }
  
  alpha = exp( init_fit_Lc2$par[1] )
  gammaa = exp( init_fit_Lc2$par[2] )
  beta_u = init_fit_Lc2$par[3:(n_Xu+2)]
  beta_p = init_fit_Lc2$par[-c(1:(n_Xu+2))]
  llp2_0 = -init_fit_Lc2$value
  
  init_fit_Lc3 = lbfgs(lc_3, grad_lc_3, vars= theta, Delta=Delta, a_i=a_i, b_i=b_i, invisible = 1, orthantwise_c = 0 )
  theta = init_fit_Lc3$par
  llp3_0 = ifelse(-init_fit_Lc3$value !="NaN", -init_fit_Lc3$value, 0)

  step = 1
  llp1 = llp2 = llp3 = numeric()
  conv1 = conv2 = conv3 = FALSE
  
  repeat{
    
    ## E-step
    expNbZ[Delta==0] =  exp( -(b0+Z_u[Delta==0,,drop=FALSE] %*% b_u + Z_p[Delta==0,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]) )
    denominator_part[Delta==0] = (alpha*(Time[Delta==0]^gammaa)*exp(X_u[Delta==0,,drop=FALSE] %*% beta_u + X_p[Delta==0,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]) )/theta
    pir[Delta==0] = 1/( 1+expNbZ[Delta==0]*(1+denominator_part[Delta==0])^theta )
    
    denominator_ac = theta + (alpha*(Time^gammaa)*exp(X_u %*% beta_u + X_p %*% beta_p ) )
    c_i = pir*( (Delta+theta)/denominator_ac)
    a_i = c_i + (1-pir)*( (Delta+theta)/theta )
    b_i = ( digamma(Delta+theta)-log(pmax(denominator_ac,1e-100)) )*pir + (digamma(Delta+theta)-log(max(theta, 1e-100)) ) * (1-pir)
    
    ## update the penalized parameters
    if (!conv1){
      out_Lc1 = lbfgs(lc1_oracle_min, grad_lc1_oracle_min, vars = c(b0,b_u,b_p),
                      Z_u=Z_u, Z_p = Z_p, pir = pir,
                      invisible = 1, orthantwise_c = 0 )
      b0 = out_Lc1$par[1]
      b_u = out_Lc1$par[2:(n_Zu+1)]
      b_p = out_Lc1$par[-c(1:(n_Zu+1))]
      llp1_1 = -out_Lc1$value
    }
    llp1 = c(llp1, llp1_1)
    
    if (!conv2){
      out_Lc2 = optim(fn=lc2_oracle_min, gr=grad_lc2_oracle_min, par = c(log(alpha),log(gammaa),beta_u,beta_p),
                      X_u = X_u, X_p = X_p, c_i = c_i, Time = Time, Delta = Delta) 
      
      if( is.na(out_Lc2$value)==TRUE ) {
        out_Lc2 = lbfgs(lc2_oracle_min, grad_lc2_oracle_min, vars = c(log(alpha),log(gammaa),beta_u,beta_p),
                        X_u = X_u, X_p = X_p, c_i = c_i, Time = Time, Delta = Delta,
                        invisible = 1, orthantwise_c = 0 )
      }
      alpha = exp( out_Lc2$par[1] )
      gammaa = exp( out_Lc2$par[2] )
      beta_u = init_fit_Lc2$par[3:(n_Xu+2)]
      beta_p = init_fit_Lc2$par[-c(1:(n_Xu+2))]
      llp2_1 = -out_Lc2$value
    }
    llp2 = c(llp2, llp2_1)
    
    if (!conv3){
      out_Lc3 = lbfgs(lc_3, grad_lc_3, vars= theta, Delta=Delta, a_i=a_i, b_i=b_i, invisible = 1, orthantwise_c = 0 )
      theta_new = out_Lc3$par
      llp3_1 = -out_Lc3$value
    }
    llp3 = c(llp3, llp3_1)
    theta_old = theta
    theta = theta_new

    if (!conv1 & abs(llp1_1- llp1_0)< tol1) conv1 <- TRUE
    if (!conv2 & abs(llp2_1- llp2_0)< tol2) conv2 <- TRUE
    if (!conv3 & (  (abs( (theta_new-theta_old)/ ((theta_old)+1) )< tol3) | (abs(llp3_1-llp3_0)< tol3) ) ) conv3 <- TRUE
    if (step > 1 & (conv1 & conv2 & conv3) | step >= nIter) {
      break
    }
    llp1_0 = llp1_1
    llp2_0 = llp2_1
    llp3_0 = llp3_1
    step = 1 + step
  }

  ## output 
  return( list(b_p = drop(b_p), b_u = b_u, b0 = b0,
               beta_u = beta_u, beta_p = beta_p, alpha = alpha, gammaa = gammaa, theta = theta, 
               logLik_Lc1 = llp1[(step-2):step], logLik_Lc2 = llp2[(step-2):step], logLik_Lc3 = llp3[(step-2):step], 
               iteration = step ) )  
}



