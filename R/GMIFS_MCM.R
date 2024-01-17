#================================================================================================================
# Functions for the GMIFS method for Weibull MCM from Fu et al. (2022) 
# All the codes were adapted by using the codes 
# from Han Fu's GitHub page: https://github.com/hanfu-bios/curemodels/blob/main/gmifs.R

# author: Fatih Kızılaslan (fatih.kizilaslan@medisin.uio.no)
# date: 17-January-2024
#================================================================================================================

library(optimg)

weibull.cure <- function(X_u, X_p, W_u, W_p, time, delta, epsilon, tol, nIter){  
  
  # X_u: N by nonp, non-penalized covariate matrix associated with incidence
  # X_p: N by J, penalized covariate matrix associated with incidence
  # W_u: N by nonp, non-penalized covariate matrix associated with latency
  # W_p:  N by M, penalized covariate matrix associated with incidence
  # time: vector of length N, observed survival time
  # delta: vector of length N, censoring status (not censored = 0)
  # epsilon: incremental size
  # tol: difference between log-likelihood
  
  N = length(time)
  J = ncol(X_p)                       # number of penalized incidence covariates
  M = ncol(W_p)                       # number of penalized latency covariates
  X_p = cbind(X_p, -X_p) # N by 2J
  W_p = cbind(W_p, -W_p) # N by 2M
  
  ## initialization
  step = 1
  b_p = rep(0, 2*J)                   # penalized incidence 
  beta_p = rep(0, 2*M)                # penalized latency  
  alpha = uniroot(function(a)         # Moment estimates of Weibull distribution parameters
    log(gamma(1+2/a))-2*log(gamma(1+1/a))-log(var(time)+(mean(time))^2)+2*log(mean(time)),
    c(0.01,10))$root
  log_alpha = log(max(alpha, 1e-15))
  lambda = gamma(1+1/alpha)/mean(time)
  log_lambda = log(max(lambda, 1e-15))
  
  b_u = rep(0, ncol(X_u))
  beta_u = rep(0, ncol(W_u))
  itct = 0
  init = optimg::optimg(par = c(log_alpha, log_lambda, b_u, itct, beta_u), fn = negloglik, 
                        b_p = b_p, beta_p = beta_p, X_u = X_u, X_p = X_p, W_u = W_u, W_p = W_p, time = time, delta = delta)
  
  log_alpha = init$par[1]
  alpha = exp(log_alpha)
  log_lambda = init$par[2]
  lambda = exp(log_lambda)
  b_u = init$par[3:(ncol(X_u)+2)]                             # unpenalized incidence 
  itct = init$par[ncol(X_u)+3]
  beta_u = init$par[(ncol(X_u)+4):(ncol(X_u)+ncol(W_u)+3)]    # unpenalized latency 
  LL0 = -init$value
  
  b_p_path <- beta_p_path <- alpha_path <- lambda_path <- b_u_path <- beta_u_path <- itct_path <- NULL
  logLikelihood <- numeric()
  
  repeat{
    
    ## update penalized parameters
    upd = update(alpha, lambda, b_p, beta_p, b_u, itct, beta_u, X_u , X_p, W_u, W_p, time, delta, epsilon)
    b_p = upd$b_p
    beta_p = upd$beta_p
    b_p_path = rbind(b_p_path, b_p)
    beta_p_path = rbind(beta_p_path, beta_p)
    
    ## update other parameters
    out = optim(par = c(log_alpha, log_lambda, b_u, itct, beta_u), fn = negloglik, gr = gradient,
                b_p = b_p, beta_p = beta_p, X_u = X_u, X_p = X_p, W_u = W_u, W_p = W_p,
                time = time, delta = delta, method="BFGS")
    
    log_alpha = out$par[1]
    alpha = exp(log_alpha)
    log_lambda = out$par[2]
    lambda = exp(log_lambda)
    b_u = out$par[3:(ncol(X_u)+2)]                            # unpenalized incidence 
    itct = out$par[ncol(X_u)+3]
    beta_u = out$par[(ncol(X_u)+4):(ncol(X_u)+ncol(W_u)+3)]  # unpenalized latency  
    
    alpha_path = c(alpha_path, alpha)
    lambda_path = c(lambda_path, lambda)
    b_u_path = rbind(b_u_path, b_u)
    itct_path = c(itct_path, itct)
    beta_u_path = rbind(beta_u_path, beta_u)
    
    LL1 = -out$value
    logLikelihood = c(logLikelihood, LL1)

    if ( step > 1 &  abs( (LL1 - LL0)/ (LL0+1) ) < tol | step >= nIter  ) { 
      break
    }
    LL0 <- LL1
    step <- 1 + step
  }

  cat("step=", step, ", |LL1-LL0|=", abs(LL1 - LL0), ", |LL1-LL0|/(LL0+1)", abs( (LL1 - LL0)/(LL0+1) ), "\n")
  
  ## output 
  b_p_path = b_p_path[,1:J] - b_p_path[,(J+1):(2*J)]
  beta_p_path = beta_p_path[,1:M] - beta_p_path[,(M+1):(2*M)]
  output <- list(b_p_path = b_p_path, beta_p_path = beta_p_path, alpha_path = alpha_path,
                 lambda_path = lambda_path, b_u_path = b_u_path, itct_path = itct_path,
                 beta_u_path = beta_u_path, logLikelihood = logLikelihood)
  return(output)
}



negloglik <- function(theta, b_p, beta_p, X_u , X_p, W_u, W_p, time, delta){
  N = length(time)
  alpha = exp(theta[1])
  lambda = exp(theta[2])
  b_u = theta[3:(ncol(X_u)+2)]                             # unpenalized incidence 
  itct = theta[ncol(X_u)+3]
  beta_u = theta[(ncol(X_u)+4):(ncol(X_u)+ncol(W_u)+3)]    # unpenalized latency
  b_nonzero = which(b_p!=0)
  beta_nonzero = which(beta_p!=0)
  logC_b = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  logC_beta = W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero]
  temp1 = 1/(1+exp(-logC_b)) 
  logC_lbeta = -(lambda * time)^alpha * exp(logC_beta)
  ll1 = delta * (log(temp1) + theta[2] + theta[1] + 
                   (alpha-1)*log(pmax(lambda * time,rep(1e-100,N))) + 
                   logC_beta + logC_lbeta)
  ll2 = (1-delta) * log(pmax(1 - temp1 * (1-exp(logC_lbeta)),rep(1e-100,N)))
  return(-sum(ll1+ll2))
}

gradient <- function(theta, b_p, beta_p, X_u , X_p, W_u, W_p, time, delta){
  N = length(time)
  alpha = exp(theta[1])
  lambda = exp(theta[2])
  b_u = theta[3:(ncol(X_u)+2)]                                # unpenalized incidence 
  itct = theta[ncol(X_u)+3]
  beta_u = theta[(ncol(X_u)+4):(ncol(X_u)+ncol(W_u)+3)]       # unpenalized latency
  b_nonzero = which(b_p!=0)
  beta_nonzero = which(beta_p!=0)
  C_b = exp(itct + X_u %*% b_u + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
  C_beta = exp(W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
  C_lbeta = exp(-(lambda * time)^alpha * C_beta)
  temp1 = log (pmax(lambda*time,rep(1e-100,N)))
  temp2 = (lambda * time)^alpha * C_beta
  temp3 = 1 / (1 + C_b*C_lbeta)
  temp4 = (1-delta) * C_b * C_lbeta * temp2 * temp3
  temp5 = delta * (1-temp2)
  temp6 = 1/(1+C_b)
  temp7 = temp6*(delta - (1-delta) * C_b * (1-C_lbeta) * temp3)
  grad1 = sum(delta * (1/alpha + temp1 - temp2 * temp1) - temp4 * temp1) * alpha
  grad2 = alpha * sum(temp5 - temp4)
  grad3 = matrix(temp7,1) %*% X_u
  grad4 = sum(temp7)                                           # intercept
  grad5 = matrix(temp5 - temp4, 1) %*% W_u
  return(-c(grad1, grad2, grad3, grad4, grad5))
}

update <- function(alpha, lambda, b_p, beta_p, b_u, itct, beta_u, X_u , X_p, W_u, W_p, time, delta, epsilon){
  b_nonzero = which(b_p!=0)
  beta_nonzero = which(beta_p!=0)
  C_b = exp(itct + X_u %*% b_u + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
  C_beta = exp(W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
  C_lbeta = exp(-(lambda * time)^alpha * C_beta)
  temp2 = (lambda * time)^alpha * C_beta
  temp3 = 1 / (1 + C_b*C_lbeta)
  temp4 = (1-delta) * C_b * C_lbeta * temp2 * temp3
  temp5 = delta * (1-temp2)
  temp6 = 1/(1+C_b)
  ## update b_p
  grad_b = matrix(temp6*(delta - (1-delta) * C_b * (1-C_lbeta) * temp3),1) %*% X_p # length J
  j_b = which.max(grad_b)
  b_p[j_b] = b_p[j_b] + epsilon
  ## update beta_p
  grad_beta = matrix(temp5 - temp4,1) %*% W_p                  # length M
  j_beta = which.max(grad_beta)
  beta_p[j_beta] = beta_p[j_beta] + epsilon
  
  return(list(b_p = b_p, beta_p = beta_p))
}


## Applying the GMIFS method to MCM by splitting the data into train, test, and validation datasets
gmifs_Vfit <- function(X_u, X_p, W_u, W_p, time, delta, cure_cutoff, nIter, tol, n_folds, epsilon){
  
  set.seed(data.seed)
  ## validation data: %20 of the data
  validation_i = sample(length(time), size = length(time)/5)        
  validation_delta = delta[validation_i]
  validation_time = time[validation_i]
  validation_X_u = scale.dummy.matrix(X_u[validation_i,])
  validation_X_p = scale.dummy.matrix(X_p[validation_i,])
  validation_W_u = scale.dummy.matrix(W_u[validation_i,]) 
  validation_W_p = scale.dummy.matrix(W_p[validation_i,])           
  
  ## removing the validation data from the whole data 
  X_u_CV = X_u[-validation_i,]
  X_p_CV = X_p[-validation_i,]
  W_u_CV = W_u[-validation_i,]
  W_p_CV = W_p[-validation_i,]
  time_CV = time[-validation_i]
  delta_CV = delta[-validation_i]
  
  ## create folds for CV based on  X_u_CV, Z_u_CV,... , Delta_CV 
  set.seed(data.seed)  
  folds_i = sample(rep(1:n_folds, length.out = length(Time_CV)))
  
  train_out =  matrix(list(NA), nrow = n_lambda, ncol = n_folds )
  C_matrix = aic = bic = Obs_Loglik = iteration = non.zeros.bp =  non.zeros.betap = theta = lambda.update = matrix( NA, nrow = n_lambda, ncol = n_folds  )
  Cstat = matrix(NA, nIter, n_folds)
  
  for (k in 1:n_folds) {
    
    test_i = which(folds_i == k)
    test_delta = delta_CV[test_i]
    test_time = time_CV[test_i]
    test_X_u = scale.dummy.matrix(X_u_CV[test_i,])
    test_X_p = scale.dummy.matrix(X_p_CV[test_i,])
    test_W_u = scale.dummy.matrix(W_u_CV[test_i,])
    test_W_p = scale.dummy.matrix(W_p_CV[test_i,])
    
    train_out = weibull.cure(X_u=scale.dummy.matrix(X_u_CV[-test_i,]), X_p=scale.dummy.matrix(X_p_CV[-test_i,]),
                             W_u=scale.dummy.matrix(W_u_CV[-test_i,]), W_p=scale.dummy.matrix(W_p_CV[-test_i,]),
                             time=time_CV[-test_i], delta=delta_CV[-test_i],
                             epsilon = epsilon, nIter=nIter, tol=tol)
    
    Cstat[1:length(train_out$itct_path),k] = 
      sapply(1:length(train_out$itct_path), 
             function(x) 
               C.stat(cure_cutoff=cure_cutoff, 
                      b0_hat=train_out$itct_path[x],  b_u_hat=train_out$b_u_path[x,], b_p_hat= train_out$b_p_path[x,],
                      beta_u_hat=train_out$beta_u_path[x,], beta_p_hat=train_out$beta_p_path[x,],
                      X_u=test_W_u, X_p=test_W_p, Z_u=test_X_u, Z_p=test_X_p, testing_time=test_time, testing_delta=test_delta)  )
    cat("Fold", k, "training finished\n")
  }

  c_stat = rowMeans(Cstat, na.rm = T)
  model_select = which.max(c_stat)
  print(model_select)
  
  ## Results based on the "train+test sets = all_data-validation_data" using the selected model, namely model_select = "number of iteration"
  output = try( weibull.cure(X_u=scale.dummy.matrix(X_u[-validation_i,]), X_p=scale.dummy.matrix(X_p[-validation_i,]),
                             W_u=scale.dummy.matrix(W_u[-validation_i,]), W_p=scale.dummy.matrix(W_p[-validation_i,]),
                             time=time[-validation_i], delta=delta[-validation_i], epsilon = epsilon, nIter=model_select, tol=tol) )
  
  # if we have error, try again using tol = tol*6 
  if (class( output ) == "try-error") {
    output = try ( weibull.cure(X_u=scale.dummy.matrix(X_u[-validation_i,]), X_p=scale.dummy.matrix(X_p[-validation_i,]),
                                W_u=scale.dummy.matrix(W_u[-validation_i,]), W_p=scale.dummy.matrix(W_p[-validation_i,]),
                                time=time[-validation_i], delta=delta[-validation_i], epsilon = epsilon, nIter=model_select, tol=tol*6) ) }
  
  # if we have error, try again using tol = tol*9 
  if (class( output ) == "try-error") {
    output = try ( weibull.cure(X_u=scale.dummy.matrix(X_u[-validation_i,]), X_p=scale.dummy.matrix(X_p[-validation_i,]),
                                W_u=scale.dummy.matrix(W_u[-validation_i,]), W_p=scale.dummy.matrix(W_p[-validation_i,]),
                                time=time[-validation_i], delta=delta[-validation_i], epsilon = epsilon, nIter=model_select, tol=tol*9) ) }
  
  # if we have error, try again using tol = tol*90 
  if (class( output ) == "try-error") {
    output = try ( weibull.cure(X_u=scale.dummy.matrix(X_u[-validation_i,]), X_p=scale.dummy.matrix(X_p[-validation_i,]),
                                W_u=scale.dummy.matrix(W_u[-validation_i,]), W_p=scale.dummy.matrix(W_p[-validation_i,]),
                                time=time[-validation_i], delta=delta[-validation_i], epsilon = epsilon, nIter=model_select, tol=tol*90) ) }
  
  ## Results based on only the validation_data
  validation_out = try( weibull.cure(X_u=validation_X_u, X_p=validation_X_p, W_u=validation_W_u, W_p=validation_W_p, time = validation_time, delta = validation_delta, epsilon = epsilon, nIter=model_select, tol=tol) )
 
   #if we have error, try again using tol = tol*6
  if (class( validation_out ) == "try-error") {
    validation_out = try( weibull.cure(X_u=validation_X_u, X_p=validation_X_p, W_u=validation_W_u, W_p=validation_W_p, time = validation_time, delta = validation_delta, epsilon = epsilon, nIter=model_select, tol=tol*6) )
  }
  
  #if we have error, try again using tol = tol*9
  if (class( validation_out ) == "try-error") {
    validation_out = try( weibull.cure(X_u=validation_X_u, X_p=validation_X_p, W_u=validation_W_u, W_p=validation_W_p, time = validation_time, delta = validation_delta, epsilon = epsilon, nIter=model_select, tol=tol*9) )
  }
  
  #if we have error, try again using tol = tol*90
  if (class( validation_out ) == "try-error") {
    validation_out = try( weibull.cure(X_u=validation_X_u, X_p=validation_X_p, W_u=validation_W_u, W_p=validation_W_p, time = validation_time, delta = validation_delta, epsilon = epsilon, nIter=model_select, tol=tol*90) )
  }
  
  nIter_validation = length(validation_out$alpha_path)                 # Note: It can be less than "model_select" when the convergency criteria satisfied.
  if(nIter_validation == 1){
    Cstat.validation = C.stat(cure_cutoff=cure_cutoff, b0_hat=as.numeric(validation_out$itct_path),  b_u_hat=as.numeric(validation_out$b_u_path), b_p_hat=as.numeric(validation_out$b_p_path), beta_u_hat=as.numeric(validation_out$beta_u_path),
                              beta_p_hat=as.numeric(validation_out$beta_p_path), X_u=validation_W_u, X_p=validation_W_p, Z_u=validation_X_u, Z_p=validation_X_p, testing_time = validation_time, testing_delta = validation_delta)
  }
  if(nIter_validation != 1){
    Cstat.validation = C.stat(cure_cutoff=cure_cutoff, b0_hat=validation_out$itct_path[nIter_validation ],  b_u_hat=validation_out$b_u_path[nIter_validation, ], b_p_hat=validation_out$b_p_path[nIter_validation,], beta_u_hat=validation_out$beta_u_path[nIter_validation, ],
                              beta_p_hat=validation_out$beta_p_path[nIter_validation, ], X_u=validation_W_u, X_p=validation_W_p, Z_u=validation_X_u, Z_p=validation_X_p, testing_time = validation_time, testing_delta = validation_delta)
  }
  
  list(out = output, Selected_model = model_select, Average_Cstat = c_stat,  C.stat.validation = Cstat.validation )
}



