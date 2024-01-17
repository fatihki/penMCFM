#================================================================================================================
# Functions for the GMIFS and penCox.1se methods for penMCFM
#
# author: Fatih Kızılaslan (fatih.kizilaslan@medisin.uio.no)
# date: 17-January-2024
#================================================================================================================

library(glmnet)

# par=(log alpha, log gammaa, log theta, beta_u, b0, b_u)
observed_negloglik <- function( par, b_p, beta_p, X_u, X_p, Z_u, Z_p, Time, Delta){
  n = length(Time)
  alpha = exp(par[1])
  gammaa = exp(par[2])
  theta = exp(par[3])
  beta_u = par[4:(ncol(X_u)+3)]
  b0 = par[ncol(X_u)+4]
  b_u = par[(ncol(X_u)+5):(ncol(X_u)+ncol(Z_u)+4)]
  
  bp_nonzero = which(b_p!=0)
  betap_nonzero = which(beta_p!=0)
  
  betaX = X_u %*% beta_u + X_p[,betap_nonzero, drop=FALSE] %*% beta_p[betap_nonzero]
  bZ = b0 + Z_u %*% b_u + Z_p[,bp_nonzero, drop=FALSE] %*% b_p[bp_nonzero]
  
  A = 1+ ( (alpha*(Time^gammaa)*exp(betaX))/theta )
  B0 =  (( A^theta) + exp(bZ) )/( 1 + exp(bZ) )
  B = log( pmax(B0,rep(1e-100,n)) )
  
  ll1 = Delta*( bZ - log(1+exp(bZ)) + par[1] + par[2] + log(pmax(Time^(gammaa-1),rep(1e-100,n))) + betaX - (theta+1)*log( pmax(A,rep(1e-100,n)) ) )
  ll2 = (1-Delta)* log( pmax(B,rep(1e-100,n)) )
  
  return( -sum(ll1+ll2) )
}


gradient_observed_negloglik <- function( par, b_p, beta_p, X_u, X_p, Z_u, Z_p, Time, Delta){
  n = length(Time)
  alpha = exp(par[1])
  gammaa = exp(par[2])
  theta = exp(par[3])
  beta_u = par[4:(ncol(X_u)+3)]
  b0 = par[ncol(X_u)+4]
  b_u = par[(ncol(X_u)+5):(ncol(X_u)+ncol(Z_u)+4)]
  
  bp_nonzero = which(b_p!=0)
  betap_nonzero = which(beta_p!=0)
  
  betaX = X_u %*% beta_u + X_p[,betap_nonzero, drop=FALSE] %*% beta_p[betap_nonzero]
  bZ = b0 + Z_u %*% b_u + Z_p[,bp_nonzero, drop=FALSE] %*% b_p[bp_nonzero]
  
  A = 1+ ( (alpha*(Time^gammaa)*exp(betaX))/theta )
  B0 =  (( A^theta) + exp(bZ) )/( 1 + exp(bZ) )
  B = log( pmax(B0,rep(1e-100,n)) )
  
  c1 = 1-(1/A)
  c2 = (A^theta) + exp(bZ)
  c3 = log( pmax(Time^gammaa,rep(1e-100,n)) )
  c4 = log(pmax(A,rep(1e-100,n)))
  c5 = (1+exp(bZ))*c2
  
  gr1 = sum( Delta * ( 1- (theta+1) * c1 ) + (1-Delta) * ( ( theta * (A^(theta-1)) * (A-1) * (1/c2) ) - (theta * c1) ) )
  gr2 = sum( Delta * ( 1 + c3 - ((theta+1) * c1 * c3) ) + (1-Delta) * theta * (A-1) * c3 * ( ( A^(theta-1) /c2 ) - (1/(A^2)) ) )
  gr3 = sum( Delta * (-c4 + (theta+1) * c1 ) + (1-Delta) * theta * ( ( (A^theta) * (c4-c1) * (1/c2) ) - c4 + c1 ) )
  gr4_1 = Delta * (1-(theta+1) * c1) + (1-Delta) * theta * (A-1) * ( (A^(theta-1)/c2) - (1/A) )
  gr4 = matrix(gr4_1,1) %*% X_u # for beta_p
  gr5_1 = Delta * (1/(1+exp(bZ))) + (1-Delta) * ( (1-A^theta) * exp(bZ) *(1/c5) )
  gr5 = sum(gr5_1) # for b0
  gr6 = matrix(gr5_1,1) %*% Z_u # for b_u
  
  return( -c(gr1, gr2, gr3, gr4, gr5, gr6) )
}


update.penMCFM <- function(alpha, gammaa, theta, b_p, beta_p, b0, b_u, beta_u, X_u , X_p, Z_u, Z_p, Time, Delta, epsilon){
  bp_nonzero = which(b_p!=0)
  betap_nonzero = which(beta_p!=0)
  
  betaX = X_u %*% beta_u + X_p[,betap_nonzero, drop=FALSE] %*% beta_p[betap_nonzero]
  bZ = b0 + Z_u %*% b_u + Z_p[,bp_nonzero, drop=FALSE] %*% b_p[bp_nonzero]
  
  A = 1+ ( (alpha*(Time^gammaa)*exp(betaX))/theta )
  c1 = 1-(1/A)
  c2 = (A^theta) + exp(bZ)
  c5 = (1+exp(bZ))*c2
  
  gr4_1 = Delta * (1-(theta+1) * c1) + (1-Delta) * theta * (A-1) * ( (A^(theta-1)/c2) - (1/A) )
  gr5_1 = Delta * ( 1/(1+exp(bZ)) ) + (1-Delta) * ( (1-(A^theta)) * exp(bZ) *(1/c5) )
  ## update beta_p
  grad_betap = matrix(gr4_1,1) %*% X_p
  j_betap = which.max(grad_betap)
  beta_p[j_betap] = beta_p[j_betap] + epsilon
  ## update b_p
  grad_bp = matrix(gr5_1,1) %*% Z_p
  j_bp = which.max(grad_bp)
  b_p[j_bp] = b_p[j_bp] + epsilon
  
  return(list(b_p = b_p, beta_p = beta_p))
}



weibull.cure.penMCFM <- function(X_u, X_p, Z_u, Z_p, Time, Delta, epsilon, tol, nIter){  
  # epsilon: incremental size
  # tol: difference between log-likelihood
  
  n = length(Time)
  P1 = ncol(Z_p)                                             # number of penalized incidence covariates (i.e. cure part)
  P2 = ncol(X_p)                                             # number of penalized latency covariates (i.e. survival part)
  Z_p = cbind(Z_p, -Z_p)                                     # n x 2P1
  X_p = cbind(X_p, -X_p)                                     # n x 2P2
  
  ## initialization
  step = 1
  b_p = rep(0, 2*P1)                                         # penalized incidence 
  beta_p = rep(0, 2*P2)                                      # penalized latency 
  ## Moment estimates of Weibull distribution parameters
  gammaa = uniroot(function(a)                               # MOM estimate of gammaa parameter
    log(gamma(1+2/a))-2*log(gamma(1+1/a))-log(var(Time)+(mean(Time))^2)+2*log(mean(Time)),
    c(0.01,10))$root
  alpha = ( gamma(1+1/gammaa)/mean(Time) )^(1/gammaa)        # MOM estimate of alpha parameter
  log_alpha = log(max(alpha, 1e-15))
  log_gammaa = log(max(gammaa, 1e-15))
  theta = 1
  log_theta = log(theta)
  
  b_u = rep(0, ncol(Z_u))
  beta_u = rep(0, ncol(X_u))
  b0 = 0
  init = optim(par = c(log_alpha, log_gammaa, log_theta, beta_u, b0, b_u), fn = observed_negloglik, gr = gradient_observed_negloglik,
               b_p = b_p, beta_p = beta_p, X_u = X_u, X_p = X_p, Z_u = Z_u, Z_p = Z_p,
               Time = Time, Delta = Delta, method ="BFGS")
  
  log_alpha = init$par[1]
  alpha = exp(log_alpha)
  log_gammaa = init$par[2]
  gammaa = exp(log_gammaa)
  log_theta = init$par[3]
  theta = exp(log_theta)
  beta_u = init$par[4:(ncol(X_u)+3)]
  b0 = init$par[ncol(X_u)+4]
  b_u = init$par[(ncol(X_u)+5):(ncol(X_u)+ncol(Z_u)+4)]
  LL0 = -init$value
  
  b_p_path <- beta_p_path <- alpha_path <- gammaa_path <- theta_path <- b_u_path <- beta_u_path <- b0_path <- NULL
  logLikelihood <- numeric()
  
  ## loop 
  
  repeat{
    ## update penalized parameters
    upd = update.penMCFM(alpha, gammaa, theta, b_p, beta_p, b0, b_u, beta_u, X_u , X_p, Z_u, Z_p, Time, Delta, epsilon)
    b_p = upd$b_p
    beta_p = upd$beta_p
    b_p_path = rbind(b_p_path, b_p)
    beta_p_path = rbind(beta_p_path, beta_p)
    
    ## update other parameters
    out = optim(par = c(log_alpha, log_gammaa, log_theta, beta_u, b0, b_u), fn = observed_negloglik, gr = gradient_observed_negloglik,
                b_p = b_p, beta_p = beta_p, X_u = X_u, X_p = X_p, Z_u = Z_u, Z_p = Z_p,
                Time = Time, Delta = Delta, method ="BFGS")
    
    log_alpha = out$par[1]
    alpha = exp(log_alpha)
    log_gammaa = out$par[2]
    gammaa = exp(log_gammaa)
    log_theta = out$par[3]
    theta = exp(log_theta)
    beta_u = out$par[4:(ncol(X_u)+3)]
    b0 = out$par[ncol(X_u)+4]
    b_u = out$par[(ncol(X_u)+5):(ncol(X_u)+ncol(Z_u)+4)]
    
    alpha_path = c(alpha_path, alpha)
    gammaa_path = c(gammaa_path, gammaa)
    theta_path = c(theta_path, theta)
    b0_path = c(b0_path, b0)
    b_u_path = rbind(b_u_path, b_u)
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
  b_p_path = b_p_path[,1:P1] - b_p_path[,(P1+1):(2*P1)]
  beta_p_path = beta_p_path[,1:P2] - beta_p_path[,(P2+1):(2*P2)]
  output <- list(b_p_path = b_p_path, beta_p_path = beta_p_path, alpha_path = alpha_path,
                 gammaa_path = gammaa_path, theta_path = theta_path,  b0_path = b0_path,
                 b_u_path = b_u_path, beta_u_path = beta_u_path, logLikelihood = logLikelihood)
  return( output )
}


gmifs_penMCFM_fit <- function(X_u, X_p, Z_u, Z_p, Time, Delta, cure_cutoff, nIter, tol, n_folds, epsilon){
  #Cross-validation
  set.seed(data.seed)
  folds_i = sample(rep(1:n_folds, length.out = length(Time)))
  Cstat = matrix(NA, nIter, n_folds)
  for (k in 1:n_folds) {
    test_i = which(folds_i == k)
    test_delta = Delta[test_i]
    test_time = Time[test_i]
    test_X_u = scale.dummy.matrix(X_u[test_i,])
    test_X_p = scale.dummy.matrix(X_p[test_i,])
    test_Z_u = scale.dummy.matrix(Z_u[test_i,])
    test_Z_p = scale.dummy.matrix(Z_p[test_i,])
    
    train_out = weibull.cure.penMCFM(X_u=scale.dummy.matrix(X_u[-test_i,]), X_p=scale.dummy.matrix(X_p[-test_i,]),
                                     Z_u=scale.dummy.matrix(Z_u[-test_i,]), Z_p=scale.dummy.matrix(Z_p[-test_i,]),
                                     Time=Time[-test_i], Delta=Delta[-test_i],
                                     epsilon = epsilon, nIter=nIter, tol=tol)
    
    Cstat[1:length(train_out$b0_path),k] = sapply(1:length(train_out$b0_path), function(x) C.stat(cure_cutoff=cure_cutoff, 
                                                                                                  b0_hat=train_out$b0_path[x],  b_u_hat=train_out$b_u_path[x,], b_p_hat= train_out$b_p_path[x,],
                                                                                                  beta_u_hat=train_out$beta_u_path[x,], beta_p_hat=train_out$beta_p_path[x,],
                                                                                                  X_u=test_X_u, X_p=test_X_p, Z_u=test_Z_u, Z_p=test_Z_p, testing_time=test_time, testing_delta=test_delta)  )
    
    cat("Fold", k, "training finished\n")
  }
  c_stat = rowMeans(Cstat, na.rm = T)
  model_select = which.max(c_stat)
  print(model_select)
  output = weibull.cure.penMCFM(X_u=scale.dummy.matrix(X_u), X_p=scale.dummy.matrix(X_p), Z_u=scale.dummy.matrix(Z_u), 
                                Z_p=scale.dummy.matrix(Z_p), Time = Time, Delta = Delta, epsilon = epsilon, nIter=model_select, tol=tol)
  
  return( list(out=output, Selected_model=model_select, Average_Cstat=c_stat) )
}

## Applying the GMIFS method to penMCFM by splitting the data into train, test, and validation datasets
gmifs_penMCFM_Vfit<- function(X_u, X_p, Z_u, Z_p, Time, Delta, cure_cutoff, nIter, tol, n_folds, epsilon){
  
  set.seed(data.seed)
  validation_i = sample(length(Time), size = length(Time)/5) # validation data is %20 of the data
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
  
  Cstat = matrix(NA, nIter, n_folds)
  
  for (k in 1:n_folds) {
    
    test_i = which(folds_i == k)
    test_delta = Delta[test_i]
    test_time = Time[test_i]
    test_X_u = scale.dummy.matrix(X_u_CV[test_i,])
    test_X_p = scale.dummy.matrix(X_p_CV[test_i,])
    test_Z_u = scale.dummy.matrix(Z_u_CV[test_i,])
    test_Z_p = scale.dummy.matrix(Z_p_CV[test_i,])
    
    train_out = weibull.cure.penMCFM(X_u=scale.dummy.matrix(X_u_CV[-test_i,]), X_p=scale.dummy.matrix(X_p_CV[-test_i,]),
                                     Z_u=scale.dummy.matrix(Z_u_CV[-test_i,]), Z_p=scale.dummy.matrix(Z_p_CV[-test_i,]),
                                     Time=Time_CV[-test_i], Delta=Delta_CV[-test_i],
                                     epsilon = epsilon, nIter=nIter, tol=tol)
    
    Cstat[1:length(train_out$b0_path),k] = sapply(1:length(train_out$b0_path), function(x) C.stat(cure_cutoff=cure_cutoff, 
                                                                                                  b0_hat=train_out$b0_path[x],  b_u_hat=train_out$b_u_path[x,], b_p_hat= train_out$b_p_path[x,],
                                                                                                  beta_u_hat=train_out$beta_u_path[x,], beta_p_hat=train_out$beta_p_path[x,],
                                                                                                  X_u=test_X_u, X_p=test_X_p, Z_u=test_Z_u, Z_p=test_Z_p, testing_time=test_time, testing_delta=test_delta)  )
    
    cat("Fold", k, "training finished\n")
  }
  
  c_stat = rowMeans(Cstat, na.rm = T)
  model_select = which.max(c_stat)
  print(model_select)
  
  ## Results based on the "train+test sets = all_data-validation_data" using the selected model, namely model_select = "number of iteration"
  output = weibull.cure.penMCFM(X_u=scale.dummy.matrix(X_u[-validation_i,]), X_p=scale.dummy.matrix(X_p[-validation_i,]), Z_u=scale.dummy.matrix(Z_u[-validation_i,]), 
                                Z_p=scale.dummy.matrix(Z_p[-validation_i,]), Time = Time[-validation_i], Delta = Delta[-validation_i], epsilon = epsilon, nIter=model_select, tol=tol)
  
  
  ## Results based on only the validation_data and its output
  validation_out = weibull.cure.penMCFM(X_u=validation_X_u, X_p=validation_X_p, Z_u=validation_Z_u, Z_p=validation_Z_p, Time = validation_Time, Delta = validation_Delta, epsilon = epsilon, nIter=model_select, tol=tol)
  nIter_validation = length(validation_out$alpha_path)          # Note: It can be less than "model_select" when the convergency criteria satisfied.
  if(nIter_validation == 1){
    Cstat.validation = C.stat(cure_cutoff=cure_cutoff, b0_hat=as.numeric(validation_out$b0_path),  b_u_hat=as.numeric(validation_out$b_u_path), b_p_hat=as.numeric(validation_out$b_p_path), beta_u_hat=as.numeric(validation_out$beta_u_path),
                              beta_p_hat=as.numeric(validation_out$beta_p_path), X_u=validation_X_u, X_p=validation_X_p, Z_u=validation_Z_u, Z_p=validation_Z_p, testing_time = validation_Time, testing_delta = validation_Delta)
  }
  if(nIter_validation != 1){
    Cstat.validation = C.stat(cure_cutoff=cure_cutoff, b0_hat=validation_out$b0_path[nIter_validation],  b_u_hat=validation_out$b_u_path[nIter_validation, ], b_p_hat= validation_out$b_p_path[nIter_validation, ], beta_u_hat=validation_out$beta_u_path[nIter_validation, ], 
                              beta_p_hat=validation_out$beta_p_path[nIter_validation, ], X_u=validation_X_u, X_p=validation_X_p, Z_u=validation_Z_u, Z_p=validation_Z_p, testing_time = validation_Time, testing_delta = validation_Delta) 
  }
  
  return( list(out = output, Selected_model = model_select, Average_Cstat = c_stat, C.stat.validation = Cstat.validation ) )
}


## Cox penalized regression with glmnet package for the latency part of the model
cox_glmnet_fit <- function(X_u, X_p, Time, Delta, alpha_enet, n_folds){
  
  set.seed(data.seed)
  ## validation data: %20 of the data
  validation_i = sample(length(Time), size = length(Time)/5)
  validation_Delta = Delta[validation_i]
  validation_Time = Time[validation_i]
  validation_X_u = scale.dummy.matrix(X_u[validation_i,])
  validation_X_p = scale.dummy.matrix(X_p[validation_i,])
  
  ## removing the validation data from the whole data 
  X_u_CV = X_u[-validation_i,]
  X_p_CV = X_p[-validation_i,]
  Time_CV = Time[-validation_i]
  Delta_CV = Delta[-validation_i]
  
  ## create folds for CV based on  X_u_CV, Z_u_CV,... , Delta_CV
  set.seed(data.seed)  
  folds_i = sample(rep(1:n_folds, length.out = length(Time_CV)))
  
  cvfit = cv.glmnet(x = cbind(scale.dummy.matrix(X_u_CV), scale.dummy.matrix(X_p_CV)), 
                    y = Surv(time = Time_CV, event = Delta_CV),
                    family = "cox",
                    penalty.factor = c( rep(0, ncol(X_u)), rep(1, ncol(X_p)) ),  
                    type.measure = "C",
                    nfolds = n.fold,
                    foldid = folds_i,
                    alpha = alpha_enet)
  
  coef.cvfit.lambda.min = as.numeric(coef(cvfit, s = "lambda.min")) # Note: "lambda.min" was used in Fu et al. (2022) study 
  beta_u.lambda.min = coef.cvfit.lambda.min[ 1:ncol(X_u) ]
  beta_p.lambda.min = coef.cvfit.lambda.min[ -(1:ncol(X_u)) ]
  
  coef.cvfit.lambda.1se = as.numeric(coef(cvfit, s = "lambda.1se")) 
  beta_u.lambda.1se = coef.cvfit.lambda.1se[1:ncol(X_u)]
  beta_p.lambda.1se = coef.cvfit.lambda.1se[-(1:ncol(X_u))]

  ## C-statistics value based on the validation set for tuning alpha parameter  using estimates from CV
  Cstat_model_lambda.min = Cstat_beta(beta_u_hat = beta_u.lambda.min ,  beta_p_hat = beta_p.lambda.min, X_u=scale.dummy.matrix(validation_X_u), X_p=scale.dummy.matrix(validation_X_p), time = validation_Time, delta = validation_Delta)
  Cstat_model_lambda.1se = Cstat_beta(beta_u_hat = beta_u.lambda.1se ,  beta_p_hat = beta_p.lambda.1se, X_u=scale.dummy.matrix(validation_X_u), X_p=scale.dummy.matrix(validation_X_p), time = validation_Time, delta = validation_Delta)
  
  return( list(beta_u.min = beta_u.lambda.min, beta_p.min = beta_p.lambda.min, Cstat.validat.min = Cstat_model_lambda.min, 
               beta_u.1se = beta_u.lambda.1se, beta_p.1se = beta_p.lambda.1se, Cstat.validat.1se = Cstat_model_lambda.1se ) )
}


