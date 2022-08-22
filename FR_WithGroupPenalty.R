rm(list=ls())
source("./Util_FR.R")
library(expm)
library(fda)
library(vars)
library(caret)
library(foreach)
library(doParallel)
library(Rcpp)
library(gss)
library(fda.usc)
library(refund)
sourceCpp("./Rcpp_FR.cpp")
# registerDoParallel(24)
set.seed(7)

# Basic simulation setting
current_A_setting <- c('Random')
current_A_setting <- c('Exp')

simulation_settings <- list(c(totalT_train=50, n=5, R=5, p=3, spectral_norm=1),
                            c(totalT_train=100, n=20, R=20, p=3, spectral_norm=1),
                            c(totalT_train=200, n=40, R=50, p=3, spectral_norm=1))

current_xepsilon_decay_setting <- c(T)
noise_range <- 0.2 # functional noise, Z is U(-noise_range, noise_range)
sampling_noise <- 0 # i.i.d. sampling noise (no sampling noise)
no_experiments <- 500

# Simulation starts
for(current_simulation_setting in simulation_settings){
  totalT_train <- current_simulation_setting[1] # number of curves
  totalT_test <- round(totalT_train*0.5)
  totalT <- totalT_train + totalT_test
  n <- m <- current_simulation_setting[2] # number of sampling points per curve
  R <- current_simulation_setting[3] # rank of transition operator A and predictor X and functional coefficient beta
  p <- current_simulation_setting[4] # number of scalar predictors
  spectral_norm <- current_simulation_setting[5] # the spectral norm of the transition operator A and functional coefficient beta
  
  # Simulation
  simulation_result <- foreach(experiment_index=1:no_experiments, .packages=c('fda','fda.usc','refund','caret'))%dopar%{
    print(paste('Experiment', experiment_index))
    set.seed(experiment_index)
    
    #### Simulate sampling points and transition operator A
    sample_pointsX <- seq(0,1,length.out=m) # sampling points for covariate X() at s1,s2,...,sm 
    sample_pointsY <- seq(0,1,length.out=n) # sampling points for response Y() r1,r2,...,rn
    # sample_points <- sort(runif(n,0,1))
    
    #### Simulate A(s,t) and B(t)
    if(current_A_setting=='Random'){
      Lambda <- randomLambda(R, spectral_norm=spectral_norm)
      beta_coef_nonzero <- runif(R,-1,1) # coefficient function for scalar predictor z
      beta_coef_nonzero <- spectral_norm*beta_coef_nonzero/sqrt(sum(beta_coef_nonzero)^2)
      beta_coef <- cbind(beta_coef_nonzero, matrix(0, nrow=R, ncol=p-1)) # only the first beta(t) is not zero (i.e. all other scalar covariates are insignificant)
    }
    if(current_A_setting=='Exp'){
      Lambda <- sqrt(3)*spectral_norm*exp_ts_vector(R) #*sqrt(3)
      beta_coef <- c(sqrt(3)*spectral_norm, rep(0, p-1))
    }
    if(current_A_setting=='SinCos'){
      Lambda <- sqrt(3)*spectral_norm*cossin_ts_vector(R) #*sqrt(3)
      beta_coef <- c(sqrt(3)*spectral_norm, rep(0, p-1))
    }
    
    #### Simulate Functional regression
    X_coef <- Epsilon_coef <- scalar_Z <- c() # X(s) and Y(s) coefficient
    for(t in 1:totalT){
      # data generation
      if(current_xepsilon_decay_setting){ # variance decay with i=1,2,...,R
        epsilon_coef <- runif(R, min=-noise_range, max=noise_range)/(1:R) # epsilon(s)
        x_coef <- runif(R, min=-1, max=1)/(1:R)
      }else{ # constant variance across i=1,2,...,R
        epsilon_coef <- runif(R, min=-noise_range, max=noise_range)/sqrt(R) # epsilon(s)
        x_coef <- runif(R, min=-1, max=1)/sqrt(R) #*sqrt(3)
      }
      scalar_z <- runif(p, min=-1, max=1)/sqrt(3) # scalar predictor z
      # combine
      scalar_Z <- cbind(scalar_Z, scalar_z)
      X_coef <-  cbind(X_coef, c(x_coef))
      Epsilon_coef <- cbind(Epsilon_coef, c(epsilon_coef))
    }
    Y_coef_best <- Lambda%*%X_coef+beta_coef%*%scalar_Z # best estimation of Y_coef given X_coef

    #### Sample Covariate X() at (s1,s2,...,sm), Sample Response Y() at (r1,r2,...,rn)
    X_sample <- Y_sample <- c()
    for(t in 1:totalT){
      # covariate X
      tmpX <- functionValue(sample_pointsX, X_coef[,t], R=R)
      tmpX <- tmpX+rnorm(m, sd=sampling_noise) # add iid sampling noise to each observation point (default sampling_noise=0)
      X_sample <- cbind(X_sample, tmpX)
      # response Y
      tmpY <- functionValue_all(s=sample_pointsY, basis_coef=Y_coef_best[,t], R=R, type=current_A_setting) + # signal part
        functionValue(sample_pointsY, Epsilon_coef[,t], R=R)
      tmpY <- tmpY+rnorm(n, sd=sampling_noise) # add iid sampling noise to each observation point (default sampling_noise=0)
      Y_sample <- cbind(Y_sample, tmpY)
    }
    Y_test <- Y_sample[,(totalT_train+1):totalT] 
    X_test <- X_sample[,(totalT_train+1):totalT] # later used for prediction
    Z_test <- scalar_Z[,(totalT_train+1):totalT]
    
    # Oracle prediction one step ahead
    Y_pred_oracle <- c() # signal part
    Epsilon_noise <- c() # noise part
    for(t in 1:totalT){
      tmp <- functionValue_all(s=sample_pointsY, basis_coef=Y_coef_best[,t], R=R, type=current_A_setting)
      Y_pred_oracle <- cbind(Y_pred_oracle, tmp)
      tmp <- functionValue(sample_pointsY, Epsilon_coef[,t], R=R)
      Epsilon_noise <- cbind(Epsilon_noise, tmp)
    }
    Y_oracle_prediction <- Y_pred_oracle[,(totalT_train+1):(totalT)] # later used for prediction
    constzero_prederr <- c(mean(apply(Y_oracle_prediction-0, 2, mse)),mean(apply(Y_oracle_prediction-0, 2, mae)))
    noise_part <- c(mean(apply(Epsilon_noise-0, 2, mse)),mean(apply(Epsilon_noise-0, 2, mae)))
    
    #### RKHS estimation
    K1 <- outer(sample_pointsY, sample_pointsY, kernel) # kernel matrix K(ri,rj) for Y
    K1_sqr <- sqrtm(K1) # K1^(1/2)
    K1_inv <- solve(K1)
    K2 <- outer(sample_pointsX, sample_pointsX, kernel) # kernel matrix K(si,sj) for X
    K2_sqr <- sqrtm(K2) # K2^(1/2)
    
    Y <- Y_sample[,1:totalT_train] # observed response matrix (training)
    X <- X_sample[,1:totalT_train] # observed covariate matrix (training)
    Z <- scalar_Z[,1:totalT_train] # observed scalar covariate matrix (training)
    
    #### 5-fold cross-validation to select lambda
    cv_folds <- 5
    cv_indices <- createFolds(1:totalT_train, k=cv_folds)
    converge_tol <- 1e-10
    # make the RKHS penalty for A and Beta the same
    if(current_A_setting=='Random'){
      lambda_penalty_RKHS_range <- 10^seq(-18,-5,1)
    }else{
      lambda_penalty_RKHS_range <- 10^seq(-8,-1,1)
    }
    lambda_penalty_groupL2_range <- 10^seq(-6,2,1)
    lambda_penalty_RKHS <- rep(lambda_penalty_RKHS_range, each=length(lambda_penalty_groupL2_range)) # penalty for the RKHS norm of A and H
    lambda_penalty_groupL2 <- c() # penalty for the group L2 norm of H
    for(times in 1:length(lambda_penalty_RKHS_range)){
      lambda_penalty_groupL2 <- c(lambda_penalty_groupL2, lambda_penalty_groupL2_range)
      lambda_penalty_groupL2_range <- rev(lambda_penalty_groupL2_range)
    }
    lambda_penalty0 <- cbind(lambda_penalty_RKHS, lambda_penalty_RKHS, lambda_penalty_groupL2)
    lambda_penalty <- as.list(data.frame(t(lambda_penalty0)))
    # allow the RKHS penalty for A and Beta to be different
    # lambda_penalty0 <- arrange_lambda_penalty(lambda_penalty_RKHS_range, lambda_penalty_RKHS_range, lambda_penalty_groupL2_range)
    # lambda_penalty <- as.list(data.frame(t(lambda_penalty0)))
    
    rkhs_time1 <- proc.time()
    cv_error <- CV_coordinate_descent(Y, X, Z, K1, K1_sqr, K1_inv, K2, K2_sqr, lambda_penalty, cv_indices, converge_tol)
    lambda_pen_cv_selected <- lambda_penalty[[which.min(apply(cv_error,2,mean))]]

    rkhs_est <- coordinate_descent(Y, X, Z, K1, K1_sqr, K1_inv, K2, K2_sqr, lambda_pen_cv_selected[1], lambda_pen_cv_selected[2], lambda_pen_cv_selected[3], converge_tol)
    A_est <- rkhs_est$A_est
    H_est <- rkhs_est$H_est
    Y_forecast_rkhs <- A_est%*%X_test/m + H_est%*%Z_test # prediction
    rkhs_prederr <- c(mean(apply(Y_oracle_prediction-Y_forecast_rkhs, 2, mse)), mean(apply(Y_oracle_prediction-Y_forecast_rkhs, 2, mae)))
    rkhs_time2 <- proc.time()
    
    #### PFFR estimation
    pffr_time1 <- proc.time()
    if(n==5){
      nbasisXY_pffr <- 10 # for n=5, PFFR gives error due to lack of data if nbasis=20
    }else{
      nbasisXY_pffr <- 20 # PFFR may not work for dimension higher than 25. (super slow and may collapse)
    }
    X_trans <- t(X)
    Y_trans <- t(Y)
    Z1 <- Z[1,]; Z2 <- Z[2,]; Z3 <- Z[3,]
    pffr_result <- pffr(Y_trans~Z1+Z2+Z3+ff(X_trans, xind=sample_pointsX, splinepars=list(bs="ps",m=list(c(2,1),c(2,1)), k=c(nbasisXY_pffr,nbasisXY_pffr))),
                        yind=sample_pointsY, bs.yindex=list(bs="ps", k=nbasisXY_pffr,m=c(2,1)), bs.int=list(bs="ps",k=nbasisXY_pffr,m=c(2,1)))
    Y_forecast_pffr <- t(predict(pffr_result, newdata=list(X_trans=t(X_test), Z1=Z_test[1,], Z2=Z_test[2,], Z3=Z_test[3,]), type='response'))
    pffr_prederr <- c(mean(apply(Y_oracle_prediction-Y_forecast_pffr, 2, mse)), mean(apply(Y_oracle_prediction-Y_forecast_pffr, 2, mae)))
    pffr_time2 <- proc.time()
  }
}

