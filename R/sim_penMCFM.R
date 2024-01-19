#================================================================================================================
# Main program using EM algorithm for penMCFM, and ... methods in Kizilaslan et al. (2024)
#
# author: Fatih Kızılaslan (fatih.kizilaslan@medisin.uio.no)
# date: 15-January-2024
#================================================================================================================

rm(list=ls())
start_time <- Sys.time()

library(fastDummies)
library(glmnet)
library(mvtnorm)
library(survival) 
library(optimg)

source("functions.R")
source("EM_penMCFM.R")
source("GMIFS_MCM.R")
source("GMIFS_penCox_penMCFM.R")

Xu <- Xp <- Xp_sigma <- Z <- Zu <- Zp <- Zp_sigma <- t <- delta <- censoring.rate <- true.uncure.rate <- list()
EM.out <- gmifs.penMCFM.out <-  gmifs.MCM.out <-  penCox.out <- list()
em.oracle <- bp_oracle_est <- betap_oracle_est <- uncure.hat.oracle <- list()

alpha <- 1.25; gammaa <- 2.5                      # Weibull distribution parameters
theta <- 0.5                                      # frailty distribution parameter 
n <- 500                                          # sample size
M <- 100                                          # number of replication size 

Xu.cont.cov.size<- 10                             # unpenalized continuous covariates size for X
Xp.cont.cov.size<- 1000                           # penalized continuous covariates size for X
non.zero.Xp<- 20

#Zu.cont.cov.size<- ?                             # unpenalized continuous covariates size for Z_u
Zp.cont.cov.size<- 1000                           # penalized continuous covariates size for Z
non.zero.Zp<- 20                                  # total number of non-zero coefficients

rho.Xu <- rho.Zu <- 0                             # correlation value for the unpenalized covariates
rho.Xp <- 0                                       # correlation value for the penalized covariates, it is chosen among the values (0, 0.2, 0.5)
rho.Zp <- 0 

# location of the nonzero coefficients 
nonzero.location.Zp <- nonzero.location.Xp <- floor( seq(from = 1, to = Zp.cont.cov.size, length.out = non.zero.Zp) ) 
nonzero.location.Zp1 <- nonzero.location.Xp1 <- nonzero.location.Zp[ 1:(non.zero.Zp/2) ] 

b0 <- -2                                           # Intercept
b_u <- c(-1,1)                                     # unpenalized coefficients related to Z_u
b_p <- rep(0, Zp.cont.cov.size)                    # penalized coefficients related to Z_p
b_p[nonzero.location.Zp] <- 1

set.seed(170723)
beta_u<- runif(Xu.cont.cov.size, min=-3, max=3)    # unpenalized coefficients related to X_p
beta_p<- rep(0,Xp.cont.cov.size)
beta_p[nonzero.location.Xp]<- 1 

SNR <- 1                                           # SNR is a signal values for the nonzero coefficients, it is chosen among the values (0.5, 1, 1.5, 2, 2.5)
alpha.enet <- 1                                    # alpha parameter in the elastic-net penalty, it is chosen among the values (0.1, 0.5, 0.9, 1)
n.fold <- 4                                        # fold size for CV
size_lambda <- 50                                  # size of the sequence for the lambda parameter in the elastic-net penalty
tolerance <- 1e-04                                 # tolerance limit for the difference in cost function values between consecutive two iterations
data.seed <- 10723                                 # seed number    

b_p <- b_p*SNR;               bp_true <- b_p
beta_p <- beta_p*SNR;         betap_true <- beta_p


# main for loop for the EM algorithms 
for (k in 1:M) {
  
  ## setting seed for the algorithms
  s <- data.seed + (k*25)
  set.seed(s)

  ## Covariates and coefficients generation
  Xu[[k]]<- covariates_gen(n, p=Xu.cont.cov.size, nTrue=Xu.cont.cov.size, rho=rho.Xu, sd=1)$Covariates
  # set.seed(s)                                     # We use the same covariates for X_p and Z_p
  # X_pen_cov <- covariates_gen(n, p=Xp.cont.cov.size, nTrue=10, rho=rho.Xp, sd=1)
  # Xp[[k]] <- X_pen_cov $Covariates                # We will use the same as Zp.
  # Xp_sigma[[k]] <- Z_pen_cov$Covariance           # Covariance matrix for the penalized variables
  
  ## setting categorical covariate for Z_u which has 3 levels
  categoric.variables <- sample( LETTERS[1:3], n, replace=TRUE, prob=c(0.4,0.35,0.25 )) 
  dummy.variables.Z <- dummy_columns(categoric.variables, remove_first_dummy = TRUE, remove_selected_columns=TRUE)
  Zu[[k]] <- cbind( rep(1,n), dummy.variables.Z )   # for intercept + dummy variables 
  
  set.seed(s)
  Z_pen_cov <- covariates_gen(n, p=Zp.cont.cov.size, nTrue=10, rho=rho.Zp, sd=1)
  Zp[[k]] <- Z_pen_cov$Covariates
  Zp_sigma[[k]] <- Z_pen_cov$Covariance             # Covariance matrix for the penalized variables
  
  Z <- as.matrix( cbind(Zu[[k]], Zp[[k]]) )
  Xp[[k]] <- Zp[[k]] ; Xp_sigma[[k]] <- Zp_sigma[[k]]
  X <- as.matrix( cbind(Xu[[k]], Xp[[k]]) )  
  
  all.b.coef <- c(b0,b_u,b_p)
  all.beta.coef <- c(beta_u,beta_p)
  
  true.parameters <- c(all.b.coef, all.beta.coef, alpha, gammaa, theta)
  true.uncure.rate[[k]] <- logit(Z, all.b.coef)      # for each individual in the sample
  
  ## survival data generation
  set.seed(s)
  data.Exp <- r_Mixture_W_cens(alpha, gammaa, theta, all.b.coef, all.beta.coef, mX=X, mZ=Z, censor.rate=0.5)
  t[[k]] <- data.Exp$observed_time 
  delta[[k]] <- data.Exp$censoring 
  censoring.rate[[k]] <- data.Exp$censor.rate 
  
  ## EM algorithm for penMCFM
  EM.out[[k]] <- EM_High_Cure_Adaptive_Enet_CV_fit( X_u=Xu[[k]], X_p=Xp[[k]], Z_u=as.matrix(Zu[[k]][,-1]), Z_p=Zp[[k]], Time=t[[k]], Delta=delta[[k]],
                                                  alpha_enet=alpha.enet, nIter=500, n_lambda=size_lambda, tol=tolerance, n_folds=n.fold )
  
  ## GMIFS method for penMCFM
  gmifs.penMCFM.out[[k]] <- gmifs_penMCFM_Vfit( X_u=as.matrix(Xu[[k]]), X_p=as.matrix(Xp[[k]]),  Z_u=as.matrix(Zu[[k]][,-1]), Z_p=Zp[[k]], Time=t[[k]], Delta=delta[[k]],
                                              cure_cutoff=5, nIter=10000/2, tol = 1e-05, n_folds=n.fold, epsilon = 0.001 )
  
  ## GMIFS method for Weibull MCM from Fu et al. (2022) study using their code in GitHub page: https://github.com/hanfu-bios/curemodels/blob/main/gmifs.R
  gmifs.MCM.out[[k]] <- gmifs_Vfit( X_u=as.matrix(Zu[[k]][,-1]), X_p=Zp[[k]], W_u=Xu[[k]], W_p=Xp[[k]], time=t[[k]], delta=delta[[k]], 
                                  cure_cutoff=5, nIter=1000, tol = 1e-05, n_folds=n.fold, epsilon = 0.001 )
  
  ## penCox method (both using lambda.min and lambda.1se) for MCFM
  penCox.out[[k]] <- cox_glmnet_fit( X_u=as.matrix(Xu[[k]]), X_p=as.matrix(Xp[[k]]), Time=t[[k]], Delta=delta[[k]], alpha_enet=1, n_folds=n.fold ) 
  
  
  ## EM algorithm for the oracle case
  Xp_oracle <- Xp[[k]][,nonzero.location.Xp]
  Zp_oracle <- Zp[[k]][,nonzero.location.Zp]  
  em.oracle[[k]] <- EM_frailty_Cure_Oracle( X_u=scale(Xu[[k]]),  X_p=scale(Xp_oracle), Z_u=as.matrix(Zu[[k]][,-1]), Z_p=scale(Zp_oracle), Time=t[[k]],
                                          Delta = delta[[k]], nIter = 500,  tol1 = tolerance, tol2 = tolerance, tol3 = tolerance*10 )
  
  bp_oracle <- rep(0,Zp.cont.cov.size)
  bp_oracle[nonzero.location.Zp] <- em.oracle[[k]]$b_p                                         # all coefficient with zeros
  bp_oracle_est[[k]] <- bp_oracle
  uncure.hat.oracle[[k]] <- logit(Z, c(em.oracle[[k]]$b0, em.oracle[[k]]$b_u, bp_oracle_est[[k]]) )
  
}

end_time <- Sys.time()
show(end_time - start_time)


## Average results of the used methods 

## Average results for the oracle estimates based on the obtained "em.oracle" list variable
bp_oracle_est<- betap_oracle_est<- bias.uncure.hat.oracle<- mse.uncure.hat.oracle<- list()
MSE_bp_oracle<- ME_bp_oracle<- MSE_betap_oracle<- ME_betap_oracle<- c()

for (k in 1:M) {
  
  bp_oracle <- rep(0,Zp.cont.cov.size)
  bp_oracle[nonzero.location.Zp] <- em.oracle[[k]]$b_p                                         # all coefficient with zeros
  bp_oracle_est[[k]] <- bp_oracle
  bp_oracle_bias <- (bp_oracle_est[[k]]-b_p)
  MSE_bp_oracle[k] <- sum(bp_oracle_bias^2)
  ME_bp_oracle[k] <- drop(t(bp_oracle_bias) %*% Zp_sigma[[k]] %*% bp_oracle_bias)              # model error of b_p with true covariance matrix
  
  betap_oracle <- rep(0,Xp.cont.cov.size)
  betap_oracle[nonzero.location.Xp] <- em.oracle[[k]]$beta_p                                   # all coefficient with zeros
  betap_oracle_est[[k]] <- betap_oracle
  betap_oracle_bias <- (betap_oracle_est[[k]]-beta_p)
  MSE_betap_oracle[k] <- sum(betap_oracle_bias^2)
  ME_betap_oracle[k] <- drop(t(betap_oracle_bias) %*% Xp_sigma[[k]] %*% betap_oracle_bias)     # model error of beta_p with true covariance matrix
  
  bias.uncure.hat.oracle[[k]]<- uncure.hat.oracle[[k]]-true.uncure.rate[[k]]
  mse.uncure.hat.oracle[[k]]<- bias.uncure.hat.oracle[[k]]^2
}

Oracle.Uncure.results<- c( Av_True_Uncure = mean(unlist(true.uncure.rate)), Est = mean(unlist(uncure.hat.oracle)),
                           Bias = mean(unlist(bias.uncure.hat.oracle)), MSE = mean(unlist(mse.uncure.hat.oracle))  )


## Average results wrt C-statistic, AIC and BIC criteria based on the obtained "EM.out" list variable (using "Results.for.criteria.penMCFM.R")
Cstat.out = Results.for.criteria.penMCFM(EM.out, true.parameters, M, b_p, beta_p, Zu, Zp, Zp_sigma, Xp, Xp_sigma=Zp_sigma, true.uncure.rate, criteria="Cstat")
AIC.out = Results.for.criteria.penMCFM(EM.out, true.parameters, M, b_p, beta_p, Zu, Zp, Zp_sigma, Xp, Xp_sigma=Zp_sigma, true.uncure.rate, criteria="AIC") 
BIC.out = Results.for.criteria.penMCFM(EM.out, true.parameters, M, b_p, beta_p, Zu, Zp, Zp_sigma, Xp, Xp_sigma=Zp_sigma, true.uncure.rate, criteria="BIC") 

## The C-statistics values for the selected model of three methods (wrt C-statistic, AIC and BIC criteria) based on the validation dataset 
Cstat.validation.selected.aic <- Cstat.validation.selected.bic <-  Cstat.validation.selected.Cstat <- c()

for (i in 1:M) {
  
  Cstat.validation.selected.aic[i] = EM.out[[i]][[11]][3,1]
  Cstat.validation.selected.bic[i] = EM.out[[i]][[11]][3,2]
  Cstat.validation.selected.Cstat[i] = EM.out[[i]][[11]][3,3]
  
}

## for b_p coefficients
RME_bp_Cstat = Cstat.out[[1]]/ME_bp_oracle 
ERR_bp_Cstat = Cstat.out[[2]]/MSE_bp_oracle 

RME_bp_AIC =  AIC.out[[1]]/ME_bp_oracle                                                       # relative model error wrt b_p for our.AIC.model/oracle.model
ERR_bp_AIC =  AIC.out[[2]]/MSE_bp_oracle                                                      # estimation error wrt b_p for our.AIC.model/oracle.model

RME_bp_BIC =  BIC.out[[1]]/ME_bp_oracle 
ERR_bp_BIC =  BIC.out[[2]]/MSE_bp_oracle 

## for beta_p coefficients
RME_betap_Cstat = Cstat.out[[8]]/ME_betap_oracle 
ERR_betap_Cstat = Cstat.out[[9]]/MSE_betap_oracle 

RME_betap_AIC =  AIC.out[[8]]/ME_betap_oracle                                                 # relative model error wrt beta_p for our.AIC.model/oracle.model
ERR_betap_AIC =  AIC.out[[9]]/MSE_betap_oracle                                                # estimation error wrt beta_p for our.AIC.model/oracle.model

RME_betap_BIC =  BIC.out[[8]]/ME_betap_oracle 
ERR_betap_BIC =  BIC.out[[9]]/MSE_betap_oracle 


## Average Results for GMIFS and penCox methods
b_u.penMCFM.hat <- beta_u.penMCFM.hat <- b_p.penMCFM.hat <- b_p.bias.penMCFM<- beta_p.penMCFM.hat <- beta_p.bias.penMCFM<- list()
b0.penMCFM.hat <- iteration.penMCFM <- MSE_bp.penMCFM<- ME_bp.penMCFM<- nonzeros.b_p.penMCFM <- Sensitivity.b_p.penMCFM <- Specificity.b_p.penMCFM<- FPR.b_p.penMCFM<- FDP.b_p.penMCFM <- c()
ME_betap.penMCFM <- MSE_betap.penMCFM <- nonzeros.beta_p.penMCFM <- Sensitivity.beta_p.penMCFM <- Specificity.beta_p.penMCFM<- FPR.beta_p.penMCFM<- FDP.beta_p.penMCFM <- Cstat.penMCFM <-c()

b_u.MCM.hat <- beta_u.MCM.hat <- b_p.MCM.hat <- b_p.bias.MCM <- beta_p.MCM.hat <- beta_p.bias.MCM <- list()
b0.MCM.hat <- iteration.MCM <- MSE_bp.MCM<- ME_bp.MCM<- nonzeros.b_p.MCM <- Sensitivity.b_p.MCM <- Specificity.b_p.MCM<- FPR.b_p.MCM<- FDP.b_p.MCM <- c()
ME_betap.MCM <- MSE_betap.MCM <- nonzeros.beta_p.MCM <- Sensitivity.beta_p.MCM <- Specificity.beta_p.MCM<- FPR.beta_p.MCM<- FDP.beta_p.MCM <-   Cstat.MCM <-c()

beta_u.penCox.min.hat <- beta_p.penCox.min.hat<- beta_p.bias.penCox.min<- beta_u.penCox.1se.hat <- beta_p.penCox.1se.hat<- beta_p.bias.penCox.1se<- list()
nonzeros.beta_p.penCox.min <- Sensitivity.beta_p.penCox.min <- Specificity.beta_p.penCox.min <- FPR.beta_p.penCox.min <- FDP.beta_p.penCox.min <- Cstat.penCox.min <- c()
nonzeros.beta_p.penCox.1se <- Sensitivity.beta_p.penCox.1se <- Specificity.beta_p.penCox.1se <- FPR.beta_p.penCox.1se <- FDP.beta_p.penCox.1se <- Cstat.penCox.1se <- c()
MSE_betap.penCox.min<- ME_betap.penCox.min<- MSE_betap.penCox.1se<- ME_betap.penCox.1se<- c()

Cstat.penMCFM.validation.set <- Cstat.MCM.validation.set <- Cstat.penCox.min.validation.set<- Cstat.penCox.1se.validation.set<-c()

for (k in 1:M) {
  
  ## GMIFS method results for penMCFM
  iteration.penMCFM[k] <- length(gmifs.penMCFM.out[[k]]$out$alpha_path)
  if(iteration.penMCFM[k]==1){
    b_p.penMCFM.hat[[k]] <- gmifs.penMCFM.out[[k]]$out$b_p_path
    beta_p.penMCFM.hat[[k]] <- gmifs.penMCFM.out[[k]]$out$beta_p_path
  }else{
    b_p.penMCFM.hat[[k]] <- gmifs.penMCFM.out[[k]]$out$b_p_path[iteration.penMCFM[k],]
    beta_p.penMCFM.hat[[k]] <- gmifs.penMCFM.out[[k]]$out$beta_p_path[iteration.penMCFM[k],]
  }
  
  b0.penMCFM.hat[k] <- gmifs.penMCFM.out[[k]]$out$b0_path[iteration.penMCFM[k]]
  b_u.penMCFM.hat[[k]] <- gmifs.penMCFM.out[[k]]$out$b_u_path[iteration.penMCFM[k],]
  beta_u.penMCFM.hat[[k]] <- gmifs.penMCFM.out[[k]]$out$beta_u_path[iteration.penMCFM[k],]
  
  # b_p coefficients
  b_p.bias.penMCFM[[k]] <- b_p.penMCFM.hat[[k]]-b_p
  MSE_bp.penMCFM[k]<- sum( b_p.bias.penMCFM[[k]]^2 )
  ME_bp.penMCFM[k]<- drop(t(b_p.bias.penMCFM[[k]]) %*% Zp_sigma[[k]] %*% b_p.bias.penMCFM[[k]])                                    # model error of b_p with true covariance matrix
  
  nonzeros.b_p.penMCFM[k] <- sum(b_p.penMCFM.hat[[k]]!=0)
  Sensitivity.b_p.penMCFM[k] <- sum( b_p.penMCFM.hat[[k]]!=0 & b_p!=0 )/ sum(b_p!=0)
  Specificity.b_p.penMCFM[k]<- sum( b_p.penMCFM.hat[[k]]==0 & b_p==0 ) / sum(b_p==0)
  FPR.b_p.penMCFM[k] <- sum(b_p.penMCFM.hat[[k]]!=0 & b_p==0 )/ sum(b_p==0)
  FDP.b_p.penMCFM[k] <- ifelse( nonzeros.b_p.penMCFM[k]==0, 0, sum(b_p.penMCFM.hat[[k]]!=0 & b_p==0 )/ nonzeros.b_p.penMCFM[k] )   # false discovery proportion i.e. FP/(FP+TP) 
  
  # beta_p coefficients
  beta_p.bias.penMCFM[[k]] <- beta_p.penMCFM.hat[[k]]-beta_p
  MSE_betap.penMCFM[k]<- sum( beta_p.bias.penMCFM[[k]]^2 )
  ME_betap.penMCFM[k]<- drop(t(beta_p.bias.penMCFM[[k]]) %*% Xp_sigma[[k]] %*% beta_p.bias.penMCFM[[k]])                           # model error of beta_p with true covariance matrix
  
  nonzeros.beta_p.penMCFM[k] <- sum(beta_p.penMCFM.hat[[k]]!=0) 
  Sensitivity.beta_p.penMCFM[k]<- sum(beta_p.penMCFM.hat[[k]]!=0 & beta_p!=0 )/ sum(beta_p!=0)
  Specificity.beta_p.penMCFM[k]<- sum( beta_p.penMCFM.hat[[k]]==0 & beta_p==0 ) / sum(beta_p==0)
  FPR.beta_p.penMCFM[k] <- sum(beta_p.penMCFM.hat[[k]]!=0 & beta_p==0 )/ sum(beta_p==0)
  FDP.beta_p.penMCFM[k] <- ifelse( nonzeros.beta_p.penMCFM[k]==0, 0, sum(beta_p.penMCFM.hat[[k]]!=0 & beta_p==0 )/ nonzeros.beta_p.penMCFM[k] )  
  
  Cstat.penMCFM[k] <- C.stat(cure_cutoff=5, b0_hat= b0.penMCFM.hat[k],  b_u_hat= b_u.penMCFM.hat[[k]], b_p_hat=  b_p.penMCFM.hat[[k]],
                             beta_u_hat= beta_u.penMCFM.hat[[k]], beta_p_hat= beta_p.penMCFM.hat[[k]],
                             X_u=as.matrix(Xu[[k]]), X_p=as.matrix(Xp[[k]]),  Z_u=as.matrix(Zu[[k]][,-1]), Z_p=Zp[[k]], testing_time=t[[k]], testing_delta=delta[[k]] ) 
  Cstat.penMCFM.validation.set[k] <- gmifs.penMCFM.out[[k]]$C.stat.validation
  
  ## GMIFS method results for MCM
  iteration.MCM[k] <- length(gmifs.MCM.out[[k]]$out$alpha_path)
  if(iteration.MCM[k] ==1 ){
    b_p.MCM.hat[[k]] <- gmifs.MCM.out[[k]]$out$b_p_path
    beta_p.MCM.hat[[k]] = gmifs.MCM.out[[k]]$out$beta_p_path
  }else{
    b_p.MCM.hat[[k]] <- gmifs.MCM.out[[k]]$out$b_p_path[iteration.MCM[k],]
    beta_p.MCM.hat[[k]] = gmifs.MCM.out[[k]]$out$beta_p_path[iteration.MCM[k],]
  }
  b0.MCM.hat[k] = gmifs.MCM.out[[k]]$out$itct_path[iteration.MCM[k]]
  b_u.MCM.hat[[k]] = gmifs.MCM.out[[k]]$out$b_u_path[iteration.MCM[k],]
  beta_u.MCM.hat[[k]] = gmifs.MCM.out[[k]]$out$beta_u_path[iteration.MCM[k],]
  
  # b_p coefficients
  b_p.bias.MCM[[k]] <- b_p.MCM.hat[[k]]-b_p
  MSE_bp.MCM[k]<- sum( b_p.bias.MCM[[k]]^2 )
  ME_bp.MCM[k]<- drop(t(b_p.bias.MCM[[k]]) %*% Zp_sigma[[k]] %*% b_p.bias.MCM[[k]])                                                  # model error of b_p with true covariance matrix
  
  nonzeros.b_p.MCM[k] = sum(b_p.MCM.hat[[k]]!=0)
  Sensitivity.b_p.MCM[k] = sum( b_p.MCM.hat[[k]]!=0 & b_p!=0 )/ sum(b_p!=0)
  Specificity.b_p.MCM[k]<- sum( b_p.MCM.hat[[k]]==0 & b_p==0 ) / sum(b_p==0)
  FPR.b_p.MCM[k] <- sum(b_p.MCM.hat[[k]]!=0 & b_p==0 )/ sum(b_p==0)
  FDP.b_p.MCM[k] <- ifelse( nonzeros.b_p.MCM[k]==0, 0, sum(b_p.MCM.hat[[k]]!=0 & b_p==0 )/ nonzeros.b_p.MCM[k] )  
  
  # beta_p coefficients
  beta_p.bias.MCM[[k]] <- beta_p.MCM.hat[[k]]-beta_p
  MSE_betap.MCM[k]<- sum( beta_p.bias.MCM[[k]]^2 )
  ME_betap.MCM[k]<- drop(t(beta_p.bias.MCM[[k]]) %*% Xp_sigma[[k]] %*% beta_p.bias.MCM[[k]])                                         # model error of beta_p with true covariance matrix
  
  nonzeros.beta_p.MCM[k] <- sum(beta_p.MCM.hat[[k]]!=0) 
  Sensitivity.beta_p.MCM[k]<- sum(beta_p.MCM.hat[[k]]!=0 & beta_p!=0 )/ sum(beta_p!=0)
  Specificity.beta_p.MCM[k]<- sum( beta_p.MCM.hat[[k]]==0 & beta_p==0 ) / sum(beta_p==0)
  FPR.beta_p.MCM[k] <- sum(beta_p.MCM.hat[[k]]!=0 & beta_p==0 )/ sum(beta_p==0)
  FDP.beta_p.MCM[k] <- ifelse( nonzeros.beta_p.MCM[k]==0, 0, sum(beta_p.MCM.hat[[k]]!=0 & beta_p==0 )/ nonzeros.beta_p.MCM[k] )  
  
  Cstat.MCM[k] <- C.stat(cure_cutoff=5, b0_hat= b0.MCM.hat[k],  b_u_hat= b_u.MCM.hat[[k]], b_p_hat = b_p.MCM.hat[[k]],
                         beta_u_hat= beta_u.MCM.hat[[k]], beta_p_hat= beta_p.MCM.hat[[k]],
                         X_u=as.matrix(Xu[[k]]), X_p=as.matrix(Xp[[k]]),  Z_u=as.matrix(Zu[[k]][,-1]), Z_p=Zp[[k]], testing_time=t[[k]], testing_delta=delta[[k]] ) 
  Cstat.MCM.validation.set[k] <- gmifs.MCM.out[[k]]$C.stat.validation
  
  ## Penalized Cox regression results which consider only latency part of the model!
  # Results based on lambda.min
  beta_u.penCox.min.hat[[k]] <- penCox.out[[k]]$beta_u.min 
  beta_p.penCox.min.hat[[k]] <- penCox.out[[k]]$beta_p.min
  beta_p.bias.penCox.min[[k]] <- beta_p.penCox.min.hat[[k]]-beta_p
  MSE_betap.penCox.min[k]<- sum( beta_p.bias.penCox.min[[k]]^2 )
  ME_betap.penCox.min[k]<- drop(t(beta_p.bias.penCox.min[[k]]) %*% Xp_sigma[[k]] %*% beta_p.bias.penCox.min[[k]])                     # model error beta_p with true covariance matrix
  
  nonzeros.beta_p.penCox.min[k] <- sum(beta_p.penCox.min.hat[[k]]!=0) 
  Sensitivity.beta_p.penCox.min[k]<- sum(beta_p.penCox.min.hat[[k]]!=0 & beta_p!=0 )/ sum(beta_p!=0)
  Specificity.beta_p.penCox.min[k]<- sum( beta_p.penCox.min.hat[[k]]==0 & beta_p==0 ) / sum(beta_p==0)
  FPR.beta_p.penCox.min[k] <- sum(beta_p.penCox.min.hat[[k]]!=0 & beta_p==0 )/ sum(beta_p==0)
  FDP.beta_p.penCox.min[k] <- ifelse( nonzeros.beta_p.penCox.min[k]==0, 0, sum(beta_p.penCox.min.hat[[k]]!=0 & beta_p==0 )/ nonzeros.beta_p.penCox.min[k] )  
  Cstat.penCox.min[k] <- Cstat_beta(beta_u_hat=beta_u.penCox.min.hat[[k]], beta_p_hat= beta_p.penCox.min.hat[[k]],  X_u=as.matrix(Xu[[k]]), X_p=as.matrix(Xp[[k]]),  time=t[[k]], delta=delta[[k]])
  Cstat.penCox.min.validation.set[k] <- penCox.out[[k]]$Cstat.validat.min
  
  # Results based on lambda.1se
  beta_u.penCox.1se.hat[[k]] <- penCox.out[[k]]$beta_u.1se 
  beta_p.penCox.1se.hat[[k]] <- penCox.out[[k]]$beta_p.1se 
  beta_p.bias.penCox.1se[[k]] <- beta_p.penCox.1se.hat[[k]]-beta_p
  MSE_betap.penCox.1se[k]<- sum( beta_p.bias.penCox.1se[[k]]^2 )
  ME_betap.penCox.1se[k]<- drop(t(beta_p.bias.penCox.1se[[k]]) %*% Xp_sigma[[k]] %*% beta_p.bias.penCox.1se[[k]])                      # model error beta_p with true covariance matrix
  
  nonzeros.beta_p.penCox.1se[k] <- sum(beta_p.penCox.1se.hat[[k]]!=0) 
  Sensitivity.beta_p.penCox.1se[k]<- sum(beta_p.penCox.1se.hat[[k]]!=0 & beta_p!=0 )/ sum(beta_p!=0)
  Specificity.beta_p.penCox.1se[k]<- sum( beta_p.penCox.1se.hat[[k]]==0 & beta_p==0 ) / sum(beta_p==0)
  FPR.beta_p.penCox.1se[k] <- sum(beta_p.penCox.1se.hat[[k]]!=0 & beta_p==0 )/ sum(beta_p==0)
  FDP.beta_p.penCox.1se[k] <- ifelse( nonzeros.beta_p.penCox.1se[k]==0, 0, sum(beta_p.penCox.1se.hat[[k]]!=0 & beta_p==0 )/ nonzeros.beta_p.penCox.1se[k] )  
  Cstat.penCox.1se[k] <- Cstat_beta(beta_u_hat=beta_u.penCox.1se.hat[[k]], beta_p_hat= beta_p.penCox.1se.hat[[k]],  X_u=as.matrix(Xu[[k]]), X_p=as.matrix(Xp[[k]]),  time=t[[k]], delta=delta[[k]])
  Cstat.penCox.1se.validation.set[k] <- penCox.out[[k]]$Cstat.validat.1se
  
}


# b_p
RME_bp_penMCFM = ME_bp.penMCFM/ME_bp_oracle                                                                                             # relative model error of b_p for penMCFM.model/oracle.model
ERR_bp_penMCFM = MSE_bp.penMCFM/MSE_bp_oracle                                                                                           # estimation error of b_p for penMCFM.model/oracle.model

RME_bp_MCM = ME_bp.MCM/ME_bp_oracle 
ERR_bp_MCM = MSE_bp.MCM/MSE_bp_oracle 

# beta_p
RME_betap_penMCFM = ME_betap.penMCFM/ME_betap_oracle                                                                                    # relative model error of beta_p for penMCFM.model/oracle.model
ERR_betap_penMCFM = MSE_betap.penMCFM/MSE_betap_oracle                                                                                  # estimation error of beta_p for penMCFM.model/oracle.model

RME_betap_MCM = ME_betap.MCM/ME_betap_oracle 
ERR_betap_MCM = MSE_betap.MCM/MSE_betap_oracle 

RME_betap_penCox.min = ME_betap.penCox.min/ME_betap_oracle 
ERR_betap_penCox.min = MSE_betap.penCox.min/MSE_betap_oracle 

RME_betap_penCox.1se = ME_betap.penCox.1se/ME_betap_oracle 
ERR_betap_penCox.1se = MSE_betap.penCox.1se/MSE_betap_oracle 

# uncure rate estimate results
uncure_penMCFM <- uncure.estimate.result( b0_estimate=b0.penMCFM.hat, bu_estimate=b_u.penMCFM.hat, bp_estimate=b_p.penMCFM.hat, Zu, Zp, true.uncure.rate)
uncure_MCM <- uncure.estimate.result( b0_estimate=b0.MCM.hat, bu_estimate=b_u.MCM.hat, bp_estimate=b_p.MCM.hat, Zu, Zp, true.uncure.rate)

# Combining all the methods
Uncure.results <- rbind(Cstat = Cstat.out[[17]], AIC = AIC.out[[17]], BIC = BIC.out[[17]], Oracle = Oracle.Uncure.results, penMCFM = colMeans(uncure_penMCFM) , MCM = colMeans(uncure_MCM) )
Accuracy.results.bp <- rbind( Cstat = Cstat.out[[18]], AIC = AIC.out[[18]], BIC = BIC.out[[18]],  penMCFM = bp.penMCFM , MCM = bp.MCM  )
Accuracy.results.betap <- rbind( Cstat = Cstat.out[[19]], AIC = AIC.out[[19]], BIC = BIC.out[[19]], penMCFM = betap.penMCFM , MCM = betap.MCM, penCox.min = betap.penCox.min, penCox.1se = betap.penCox.1se  )

RME.ERR.results.bp <- rbind( Cstat = c(mean(RME_bp_Cstat), mean(ERR_bp_Cstat)), AIC = c(mean(RME_bp_AIC), mean(ERR_bp_AIC)), BIC = c(mean(RME_bp_BIC), mean(ERR_bp_BIC)), 
                             penMCFM = c(mean(RME_bp_penMCFM), mean(ERR_bp_penMCFM)),  MCM = c(mean(RME_bp_MCM), mean(ERR_bp_MCM))   )
colnames(RME.ERR.results.bp) = c("RME","ERR")

RME.ERR.results.betap <- rbind( Cstat = c(mean(RME_betap_Cstat), mean(ERR_betap_Cstat)), AIC = c(mean(RME_betap_AIC), mean(ERR_betap_AIC)), BIC = c(mean(RME_betap_BIC), mean(ERR_betap_BIC)),
                                penMCFM = c(mean(RME_betap_penMCFM), mean(ERR_betap_penMCFM)),  MCM = c(mean(RME_betap_MCM), mean(ERR_betap_MCM)),
                                penCox.min = c(mean(RME_betap_penCox.min), mean(ERR_betap_penCox.min)),   penCox.1se = c(mean(RME_betap_penCox.1se), mean(ERR_betap_penCox.1se)) )
colnames(RME.ERR.results.betap) = c("RME","ERR")

Cstat.Average.Validation.set <- cbind(AIC= mean(Cstat.validation.selected.aic), BIC=mean(Cstat.validation.selected.bic), 
                                      Cstat= mean(Cstat.validation.selected.Cstat), 
                                      penMCFM=mean(Cstat.penMCFM.validation.set), penMCM=mean(Cstat.MCM.validation.set),  
                                      penCox1se=mean(Cstat.penCox.1se.validation.set), penCoxmin=mean(Cstat.penCox.min.validation.set) )
rownames(Cstat.Average.Validation.set ) <- c("Cstat.validation")


print( round(Uncure.results,5) )
print( round(Accuracy.results.bp,5) )
print( round(RME.ERR.results.bp,5) )
print( round(Accuracy.results.betap,5) )
print( round(RME.ERR.results.betap,5) )
print(round(Cstat.Average.Validation.set, 5))

cat("\n average.censoring.rate:", round(mean(unlist(censoring.rate)),5), "___alpha.enet:", alpha.enet, "___SNR:", SNR,   "___sample size:", n, 
    "\n rhoZp:", rho.Zp, "___penalized Zp cov size:", Zp.cont.cov.size, "___non-zero Zp penalized cov:",  non.zero.Zp,
    "\n rhoXp:", rho.Xp, "___penalized Xp cov size:", Xp.cont.cov.size, "___non-zero Xp penalized cov:",  non.zero.Xp)


# To save results as data frames for b_p and beta_p separately, intended for drawing related figures
# Note 1: We save AIC, BIC, Cstat results as row data, so the length of row is "3*n" 
# Note 2: "..._bp_penMCFM" represents GMIFS method for penMCFM, and "...bp_MCM" represents GMIFS method for MCM

# for b_p coefficients
df_bp = data.frame( cbind(RME_bp = c(RME_bp_AIC, RME_bp_BIC, RME_bp_Cstat, RME_bp_penMCFM, RME_bp_MCM ), 
                          ERR_bp = c(ERR_bp_AIC, ERR_bp_BIC, ERR_bp_Cstat, ERR_bp_penMCFM, ERR_bp_MCM), 
                          nonzeros_bp = c(AIC.out[[3]], BIC.out[[3]], Cstat.out[[3]], nonzeros.b_p.penMCFM, nonzeros.b_p.MCM ),
                          Sensitivity_bp = c(AIC.out[[4]], BIC.out[[4]], Cstat.out[[4]], Sensitivity.b_p.penMCFM, Sensitivity.b_p.MCM ),
                          Specificity_bp = c(AIC.out[[5]], BIC.out[[5]], Cstat.out[[5]], Specificity.b_p.penMCFM, Specificity.b_p.MCM ),
                          FPR_bp = c(AIC.out[[6]], BIC.out[[6]], Cstat.out[[6]], FPR.b_p.penMCFM, FPR.b_p.MCM ),
                          FDP_bp = c(AIC.out[[7]], BIC.out[[7]], Cstat.out[[7]], FDP.b_p.penMCFM, FDP.b_p.MCM ),
                          Cstat = c(AIC.out[[15]], BIC.out[[15]], Cstat.out[[15]], Cstat.penMCFM, Cstat.MCM ),
                          uncure.hat = c( AIC.out[[16]][,2], BIC.out[[16]][,2], Cstat.out[[16]][,2], uncure_penMCFM[,2],  uncure_MCM[,2] ),
                          uncure.bias = c( AIC.out[[16]][,3], BIC.out[[16]][,3], Cstat.out[[16]][,3], uncure_penMCFM[,3], uncure_MCM[,3] ),
                          uncure.mse = c( AIC.out[[16]][,4], BIC.out[[16]][,4], Cstat.out[[16]][,4], uncure_penMCFM[,4], uncure_MCM[,4] ),
                          SNR = rep(SNR, length(c(RME_bp_AIC, RME_bp_BIC, RME_bp_Cstat, RME_bp_penMCFM, RME_bp_MCM)) ),
                          alpha.Enet = rep( alpha.enet, length(c(RME_bp_AIC, RME_bp_BIC, RME_bp_Cstat, RME_bp_penMCFM, RME_bp_MCM)) ),
                          CorrXpZp = rep( rho.Zp ,length(c(RME_bp_AIC, RME_bp_BIC, RME_bp_Cstat, RME_bp_penMCFM, RME_bp_MCM)) )
) )



# Defining the "AIC, BIC, Cstat, penMCFM, MCM" methods as factor variable "A,B,C,D, E", respectively
df_bp$method = as.factor( c(rep("A",length(RME_bp_AIC)), rep("B",length(RME_bp_BIC)), rep("C",length(RME_bp_Cstat)),  rep("D",length(RME_bp_penMCFM)), rep("E",length(RME_bp_MCM)) ) )

# for beta_p coefficients
df_betap = data.frame( cbind(RME_betap = c(RME_betap_AIC, RME_betap_BIC, RME_betap_Cstat, RME_betap_penMCFM, RME_betap_MCM, RME_betap_penCox.1se, RME_betap_penCox.min ), 
                             ERR_betap = c(ERR_betap_AIC, ERR_betap_BIC, ERR_betap_Cstat, ERR_betap_penMCFM, ERR_betap_MCM, ERR_betap_penCox.1se, ERR_betap_penCox.min), 
                             nonzeros_betap = c(AIC.out[[10]], BIC.out[[10]], Cstat.out[[10]], nonzeros.beta_p.penMCFM, nonzeros.beta_p.MCM, nonzeros.beta_p.penCox.1se, nonzeros.beta_p.penCox.min),
                             Sensitivity_betap = c(AIC.out[[11]], BIC.out[[11]], Cstat.out[[11]], Sensitivity.beta_p.penMCFM, Sensitivity.beta_p.MCM, Sensitivity.beta_p.penCox.1se, Sensitivity.beta_p.penCox.min),
                             Specificity_betap = c(AIC.out[[12]], BIC.out[[12]], Cstat.out[[12]], Specificity.beta_p.penMCFM, Specificity.beta_p.MCM, Specificity.beta_p.penCox.1se, Specificity.beta_p.penCox.min),
                             FPR_betap = c(AIC.out[[13]], BIC.out[[13]], Cstat.out[[13]], FPR.beta_p.penMCFM, FPR.beta_p.MCM, FPR.beta_p.penCox.1se, FPR.beta_p.penCox.min),
                             FDP_betap = c(AIC.out[[14]], BIC.out[[14]], Cstat.out[[14]], FDP.beta_p.penMCFM, FDP.beta_p.MCM, FDP.beta_p.penCox.1se, FDP.beta_p.penCox.min),
                             CstatValidation = c(Cstat.validation.selected.aic, Cstat.validation.selected.bic, Cstat.validation.selected.Cstat, Cstat.penMCFM.validation.set, Cstat.MCM.validation.set ,  Cstat.penCox.1se.validation.set, Cstat.penCox.min.validation.set), 
                             SNR = rep( SNR, length(c(RME_betap_AIC, RME_betap_BIC, RME_betap_Cstat, RME_betap_penMCFM, RME_betap_MCM ,RME_betap_penCox.1se, RME_betap_penCox.min)) ),
                             alpha.Enet = rep( alpha.enet, length( c(RME_betap_AIC, RME_betap_BIC, RME_betap_Cstat, RME_betap_penMCFM, RME_betap_MCM, RME_betap_penCox.1se, RME_betap_penCox.min) ) ),
                             CorrXpZp = rep( rho.Zp, length(c(RME_betap_AIC, RME_betap_BIC, RME_betap_Cstat, RME_betap_penMCFM, RME_betap_MCM , RME_betap_penCox.1se, RME_betap_penCox.min)) )
) )

# Defining the "AIC, BIC, Cstat, penMCFM, MCM, penCox.1se and penCox.min" methods as factor variable "A, B, C, D, E, F, G", respectively
df_betap$method = as.factor( c(rep("A",length(RME_betap_AIC)), rep("B",length(RME_betap_BIC)), rep("C",length(RME_betap_Cstat)), rep("D",length(RME_betap_penMCFM)), rep("E",length(RME_betap_MCM)), rep("F",length(RME_betap_penCox.1se)), rep("G",length(RME_betap_penCox.min))    ) )  

# saving results as data frame wrt alpha.enet, rho.Zp and SNR values for b_p
save( df_bp, file = paste0("df_bp","_alpha",alpha.enet,"_rhoZp",rho.Zp,"_SNR", SNR, ".RData"))

# saving results wrt alpha.enet, rho.Zp and SNR values for beta_p
save( df_betap, file = paste0("df_betap","_alpha",alpha.enet,"_rhoZp",rho.Zp,"_SNR", SNR, ".RData" ) )



