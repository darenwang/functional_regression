rm(list=ls())
source("./Util_FR.R")

####need to install packages like fda and refund
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
set.seed(7)

# Basic simulation setting
current_A_setting <- c('Random')
current_A_setting <- c('Exp')

simulation_settings <- list(c(totalT_train=50, n=5, R=5), # main text
                            c(totalT_train=100, n=20, R=20),
                            c(totalT_train=200, n=40, R=50))

current_xepsilon_decay_setting <- c(T)
spectral_norm_settings <- c(0.5,1,2)
noise_range <- 0.2 # functional noise, Z is U(-noise_range, noise_range)
sampling_noise <- 0 # i.i.d. sampling noise 
no_experiments <- 1

# Simulation starts
for(current_simulation_setting in simulation_settings){ # basic simulation setting
  totalT_train <- current_simulation_setting[1] # number of curves
  totalT_test <- round(totalT_train*0.5)
  totalT <- totalT_train + totalT_test
  n <- m <- current_simulation_setting[2] # number of sampling points per curve
  R <- current_simulation_setting[3] # rank of transition operator A
  
  for(spectral_norm in spectral_norm_settings){
    #### Simulation starts
    simulation_result <- foreach(experiment_index=1:no_experiments, .packages=c('fda','fda.usc','refund','caret'))%dopar%{ # .errorhandling='remove', 
      print(experiment_index)
      set.seed(experiment_index)
      
      #### Simulate sampling points and transition operator A
      sample_pointsX <- seq(0,1,length.out=m) # sampling points for covariate X() at s1,s2,...,sm 
      sample_pointsY <- seq(0,1,length.out=n) # sampling points for response Y() r1,r2,...,rn
      
      ##### Simulate data 
      if(current_A_setting=='Diag'){
        Lambda <- diag(rep(spectral_norm, R))
      }
      if(current_A_setting=='Random'){
        Lambda <- randomLambda(R, spectral_norm=spectral_norm)
      }
      if(current_A_setting=='Exp'){
        Lambda <- sqrt(3)*spectral_norm*exp_ts_vector(R) #*sqrt(4)
      }
      if(current_A_setting=='SinCos'){
        Lambda <- sqrt(3)*spectral_norm*cossin_ts_vector(R) #*sqrt(4)
      }
      
      # Simulate Functional regression exactly
      X_coef <- Epsilon_coef <- c() # X(s) and Y(s) coefficient
      for(t in 1:totalT){
        if(current_xepsilon_decay_setting){ # variance decay with i=1,2,...,R
          epsilon_coef <- runif(R, min=-noise_range, max=noise_range)/(1:R) # epsilon(s)
          x_coef <- runif(R, min=-1, max=1)/(1:R)
        }else{ # constant variance across i=1,2,...,R
          epsilon_coef <- runif(R, min=-noise_range, max=noise_range)/sqrt(R) # epsilon(s)
          x_coef <- runif(R, min=-1, max=1)/sqrt(R) #*sqrt(3)
        }
        X_coef <-  cbind(X_coef, c(x_coef))  
        Epsilon_coef <- cbind(Epsilon_coef, c(epsilon_coef))
      }
      Y_coef_best <- Lambda%*%X_coef  # best estimation of Y_coef given X_coef (Y_ceof in terms of exp(-t))
      # Sample Observed Covariate X() at (s1,s2,...,sm), Sample Observed Response Y() at (r1,r2,...,rn)
      X_sample <- Y_sample <- c()
      for(t in 1:totalT){
        # covariate X
        tmpX <- functionValue(sample_pointsX, X_coef[,t], R=R)
        tmpX <- tmpX+rnorm(m, sd=sampling_noise) # add iid sampling noise to each observation point (default sampling_noise=0)
        X_sample <- cbind(X_sample, tmpX)
        # response Y
        tmpY <- functionValue_all(s=sample_pointsY, basis_coef=Y_coef_best[,t], R=R, type=current_A_setting) + # signal part
          functionValue(sample_pointsY, Epsilon_coef[,t], R=R) # functional noise part
        tmpY <- tmpY+rnorm(n, sd=sampling_noise) # add iid sampling noise to each observation point (default sampling_noise=0)
        Y_sample <- cbind(Y_sample, tmpY)
      }
      Y_test <- Y_sample[,(totalT_train+1):totalT]
      X_test <- X_sample[,(totalT_train+1):totalT] # later used for prediction
      # Oracle prediction one step ahead
      Y_pred_oracle <- c()
      for(t in 1:totalT){
        tmp <- functionValue_all(s=sample_pointsY, basis_coef=Y_coef_best[,t], R=R, type=current_A_setting)
        Y_pred_oracle <- cbind(Y_pred_oracle, tmp)
      }
      Y_oracle_prediction <- Y_pred_oracle[,(totalT_train+1):(totalT)] # later used for prediction
      constzero_prederr <- c(mean(apply(Y_oracle_prediction-0, 2, mse)),mean(apply(Y_oracle_prediction-0, 2, mae)))
      # print(constzero_prederr)

      ##### RKHS preparation
      K1 <- outer(sample_pointsY, sample_pointsY, kernel) # kernel matrix K(ri,rj) for Y
      K1_sqr <- sqrtm(K1) # K1^(1/2)
      K2 <- outer(sample_pointsX, sample_pointsX, kernel) # kernel matrix K(si,sj) for X
      K2_sqr <- sqrtm(K2) # K2^(1/2)
      
      Y <- Y_sample[,1:totalT_train] # observed response matrix (training)
      X <- X_sample[,1:totalT_train] # observed covariate matrix (training)
      # 5-fold cross-validation to select lambda
      lambda_penalty <- 10^c(seq(-15,2,0.5))
      cv_folds <- 5
      cv_indices <- createFolds(1:totalT_train, k=cv_folds)
      cv_data <- list()
      for(cv_fold_index in 1:cv_folds){
        cv_index_test <- sort(cv_indices[[cv_fold_index]])
        cv_index_train <- sort(unlist(cv_indices)[!unlist(cv_indices)%in%cv_index_test])
        cv_data[[cv_fold_index]] <- list(Y_train=Y[,cv_index_train], X_train=X[,cv_index_train],
                                         Y_test=Y[,cv_index_test], X_test=X[,cv_index_test])
      }
      # RKHS estimation using SVD Rcpp
      rkhs_time1 <- proc.time()
      cv_error <- CV_estEfunc2_Rcpp(K1, K2, K1_sqr, K2_sqr, Y, X, lambda_penalty, cv_data) # 5-fold cross-validation to select lambda
      lambda_pen_cv_selected <- lambda_penalty[which.min(apply(cv_error,2,mean))]
      A_est <- estE_func2_Rcpp(K1, K2, K1_sqr, K2_sqr, Y, X, lambda_pen_cv_selected) # estimation based on the best lambda
      Y_forecast_rkhs <- A_est%*%X_test/m # prediction
      rkhs_prederr <- c(mean(apply(Y_oracle_prediction-Y_forecast_rkhs, 2, mse)), mean(apply(Y_oracle_prediction-Y_forecast_rkhs, 2, mae)))
      rkhs_time2 <- proc.time()
  
      ##### FDA (fda.usc)
      fda_time1 <- proc.time()
      domain <- c(0,1)
      nbasis_fda <- 20
      nbasisXY_fda <- nbasis_fda # for smooth X and Y 
      basisY <- basisX <- create.bspline.basis(rangeval=domain, nbasis=nbasisXY_fda, norder=4) # cubic splines
      functionalY <- Data2fd(argvals=sample_pointsY, y=Y, basisobj=basisY)
      functionalX <- Data2fd(argvals=sample_pointsX, y=X, basisobj=basisX)
      nbasisBeta_fda <- nbasis_fda # for beta function
      basiss <- create.bspline.basis(rangeval=domain, nbasis=nbasisBeta_fda)
      basist <- create.bspline.basis(rangeval=domain, nbasis=nbasisBeta_fda)
      lambda_penalty_fda <- 10^c(seq(-15,2,0.5))
      cv_error_fda <- CV_fda(Y, X, lambda_penalty_fda, cv_indices, nbasisXY_fda, nbasisBeta_fda, domain, sample_pointsY, sample_pointsX) # 5-fold cross-validation to select lambda
      lambda_pen_cv_selected_fda <- lambda_penalty_fda[which.min(apply(cv_error_fda,2,mean))]
      # penalization fda estimation
      Lfdobjt <- Lfdobjs <- vec2Lfd(c(0,0), domain) # roughness penalty
      fda_pen <- fregre.basis.fr(x=functionalX, y=functionalY, basis.s=basiss, basis.t=basist,
                                 lambda.s=lambda_pen_cv_selected_fda, lambda.t=lambda_pen_cv_selected_fda, Lfdobj.s=Lfdobjs, Lfdobj.t=Lfdobjt)
      functionalX_test <- Data2fd(argvals=sample_pointsX, y=X_test, basisobj=basisX)
      functionalY_forecast <- predict(fda_pen, new.fdataobj=functionalX_test)
      Y_forecast_fda <- eval.fd(sample_pointsY, functionalY_forecast)
      fda_prederr <- c(mean(apply(Y_oracle_prediction-Y_forecast_fda, 2, mse)), mean(apply(Y_oracle_prediction-Y_forecast_fda, 2, mae)))
      fda_time2 <- proc.time()
      
      ##### PFFR refund
      pffr_time1 <- proc.time()
      if(n==5){
        nbasisXY_pffr <- 15 # for n=5, PFFR gives error due to lack of data if nbasis=20
      }else{
        nbasisXY_pffr <- 20 # PFFR may not work for dimension higher than 25. (super slow and may collapse)
      }
      X_trans <- t(X)
      Y_trans <- t(Y)
      pffr_result <- pffr(Y_trans~ff(X_trans, xind=sample_pointsX, splinepars=list(bs="ps",m=list(c(2,1),c(2,1)), k=c(nbasisXY_pffr,nbasisXY_pffr))),
                          yind=sample_pointsY, bs.yindex=list(bs="ps",k=nbasisXY_pffr,m=c(2,1)), bs.int=list(bs="ps",k=nbasisXY_pffr,m=c(2,1)))
      Y_forecast_pffr <- t(predict(pffr_result, newdata=list(X_trans=t(X_test)), type='response'))
      pffr_prederr <- c(mean(apply(Y_oracle_prediction-Y_forecast_pffr, 2, mse)), mean(apply(Y_oracle_prediction-Y_forecast_pffr, 2, mae)))
      pffr_time2 <- proc.time()
    }
  }
}

