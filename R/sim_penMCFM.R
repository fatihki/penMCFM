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

Xu <- Xp <- Xp_sigma <- Z <- Zu <- Zp <- Zp_sigma <- t <- delta <- censoring.rate <- list()
true.uncure.rate <- EM.out <- em.oracle <- bp_oracle_est <- betap_oracle_est <- uncure.hat.oracle <-list()

alpha <- 1.25; gammaa <- 2.5 ### Weibull distribution parameters
theta <- 0.5                  ## frailty distribution parameter 
n <- 500                    ## sample size

Xu.cont.cov.size<- 10     ## unpenalized covariates size for X
Xp.cont.cov.size<- 1000   ## penalized covariates size for X
non.zero.Xp<- 20

Zu.cont.cov.size<- 3     ## unpenalized covariates size for Z  ????for the categorical variable
Zp.cont.cov.size<- 1000  ## penalized covariates size for Z
non.zero.Zp<- 20         ## total number of non-zero coefficients

rho.Xu <- rho.Zu <- 0     ## correlation value for the covariates, it is chosen among the values (0, 0.2, 0.5).

# location of the nonzero coefficients 
nonzero.location.Zp <- nonzero.location.Xp <- floor( seq(from = 1, to = Zp.cont.cov.size, length.out = non.zero.Zp) ) 
nonzero.location.Zp1 <- nonzero.location.Xp1 <- nonzero.location.Zp[ 1:(non.zero.Zp/2) ] 

b0 <- -2.                                           ## Intercept
b_u <- c(-1,1)                                      ## unpenalized coefficients related to Z_u
b_p <- rep(0, Zp.cont.cov.size)                     ## penalized coefficients related to Z_p
b_p[nonzero.location.Zp] <- 1

set.seed(170723)
beta_u<- runif(Xu.cont.cov.size, min=-3, max=3)   ## unpenalized coefficients related to X_p
beta_p<- rep(0,Xp.cont.cov.size)
beta_p[nonzero.location.Xp]<- 1 

SNR <- 1                                          # SNR is a signal values for the nonzero coefficients, it is chosen among the values (0.5, 1, 1.5, 2, 2.5).
alpha.enet <- 1                                   # alpha parameter in the elastic-net penalty, it is chosen among the values (0.1, 0.5, 0.9, 1).
n.fold <- 4                                       # fold size for CV
size_lambda <- 50                                 # size of the sequence for the lambda parameter in the elastic-net penalty
tolerance <- 1e-04                                # tolerance limit for the consecutive two iterations
data.seed <- 10723                                # seed number    

b_p <- b_p*SNR;               bp_true <- b_p
beta_p <- beta_p*SNR;         betap_true <- beta_p



for (k in M) {
  
  # setting seed for the algorithms
  s <- data.seed + (k*25)
  set.seed(s)

  # Covariate and coefficients generation
  Xu[[k]]<- covariates_gen(n, p=Xu.cont.cov.size, nTrue=Xu.cont.cov.size, rho=rho.Xu, sd=1)$Covariates
  # set.seed(s)
  # X_pen_cov <- covariates_gen(n, p=Xp.cont.cov.size, nTrue=10, rho=rho.Xp, sd=1)
  # Xp[[k]] <- X_pen_cov $Covariates # We use the same as Zp.
  # Xp_sigma[[k]] <- Z_pen_cov$Covariance ## Covariance matrix for the penalized variables
  
  categoric.variables <- sample( LETTERS[1:3], n, replace=TRUE, prob=c(0.4,0.35,0.25)) 
  dummy.variables.Z <- dummy_columns(categoric.variables, remove_first_dummy = TRUE, remove_selected_columns=TRUE)
  Zu[[k]] <- cbind(rep(1,n),dummy.variables.Z) ## for intercept+dummy variables 
  
  set.seed(s)
  Z_pen_cov <- covariates_gen(n, p=Zp.cont.cov.size, nTrue=10, rho=rho.Zp, sd=1)
  Zp[[k]] <- Z_pen_cov$Covariates
  Zp_sigma[[k]] <- Z_pen_cov$Covariance ## Covariance matrix for the penalized variables
  
  Z <- as.matrix( cbind(Zu[[k]], Zp[[k]]) )
  Xp[[k]] <- Zp[[k]]
  X <- as.matrix( cbind(Xu[[k]], Xp[[k]]) )  
  
  all.b.coef <- c(b0,b_u,b_p)
  all.beta.coef <- c(beta_u,beta_p)
  
  true.parameters <- c(all.b.coef, all.beta.coef, alpha, gammaa, theta)
  true.uncure.rate[[k]] <- logit(Z, all.b.coef)   # for each individual in the sample
  
  # survival data generation
  set.seed(s)
  data.Exp <- r_Mixture_W_cens(alpha, gammaa, theta, all.b.coef, all.beta.coef, mX=X, mZ=Z, censor.rate=0.5)
  t[[k]] <- data.Exp$observed_time 
  delta[[k]] <- data.Exp$censoring 
  censoring.rate[[k]] <- data.Exp$censor.rate 
  
  # EM algorithm for penMCFM
  EM.out[[k]] <- EM_High_Cure_Adaptive_Enet_CV_fit(X_u=Xu[[k]], X_p=Xp[[k]], Z_u=as.matrix(Zu[[k]][,-1]), Z_p=Zp[[k]], Time=t[[k]], Delta=delta[[k]],
                                                  alpha_enet=alpha.enet, nIter=500, n_lambda = size_lambda, grid_size = grid_size,  tol=tolerance, n_folds=n.fold )
  
  ## EM algorithm for oracle case
  Xp_oracle <- Xp[[k]][,nonzero.location.Xp]
  Zp_oracle <- Zp[[k]][,nonzero.location.Zp]  
  em.oracle[[k]] <- EM_frailty_Cure_Oracle(X_u=scale(Xu[[k]]),  X_p=scale(Xp_oracle), Z_u=as.matrix(Zu[[k]][,-1]), Z_p=scale(Zp_oracle), Time=t[[k]],
                                          Delta = delta[[k]], nIter = 500,  tol1 = tolerance, tol2 = tolerance, tol3 = tolerance*10)
  bp_oracle <- rep(0, Zp.cont.cov.size)
  bp_oracle[nonzero.location.Zp] <- em.oracle[[k]]$b_p  # nonzero coefficients from the EM algorithm
  bp_oracle_est[[k]] <- bp_oracle
  
  betap_oracle <- rep(0,Xp.cont.cov.size)
  betap_oracle[nonzero.location.Xp] <- em.oracle[[k]]$beta_p  # nonzero coefficients from the EM algorithm
  betap_oracle_est[[k]] <- betap_oracle
  
  uncure.hat.oracle[[k]] <- logit(Z, c(em.oracle[[k]]$b0, em.oracle[[k]]$b_u, bp_oracle_est[[k]]) )
}

end_time <- Sys.time()
show(end_time - start_time)
