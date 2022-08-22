mse <- function(x){
  # mean sqaured error
  mean(x^2)
}

rmse <- function(x){
  # root mean sqaured error
  sqrt(mean(x^2))
}

mae <- function(x){
  # mean absolute error
  mean(abs(x))
}

rank_calculation <- function(x, cut_off=1e-1){
  # calculate the rank of the final estimated W
  W <- x$par
  if(is.null(W)){
    W <- x$par1
  }
  sum(svd(W)$d>cut_off)
}

rank_calculation_FAR2 <- function(x, cut_off=1e-1){
  # calculate the rank of the final estimated W
  W1 <- x$par1
  W2 <- x$par2
  return(c(sum(svd(W1)$d>cut_off), sum(svd(W2)$d>cut_off)))
}

truncate_A <- function(A, d){
  # approximate matrix A using first d svd
  tmp <- svd(A)
  A_approx <- tmp$u[,1:d]%*%diag(tmp$d[1:d])%*%t(tmp$v[,1:d])
  return(A_approx)
}

randomLambda <- function(R, spectral_norm=0.9){
  # function to simulate a R matrix for the R-rank transition operator A
  Lambda <- matrix(rnorm(R^2), R, R)
  Lambda_svd <- svd(Lambda)
  Lambda <- Lambda_svd$u%*%diag(Lambda_svd$d/max(Lambda_svd$d)*spectral_norm)%*%t(Lambda_svd$v)
  return(Lambda)
}

randomLambda_Aue <- function(R, sigma_l, spectral_norm){
  # function to simulate a R matrix for the R-rank transition operator A based on Aue (2014)
  sd <- sigma_l%*%t(sigma_l)
  Lambda <- matrix(rnorm(R^2, c(sd)), R, R)
  Lambda_svd <- svd(Lambda)
  Lambda <- Lambda_svd$u%*%diag(Lambda_svd$d/max(Lambda_svd$d)*spectral_norm)%*%t(Lambda_svd$v)
  return(Lambda)
}

cosine_fun <- function(s,k=1){
  # cosine function for RKHS basis
  if(k==1){
    return(cos(0*s))
  }else{
    return(sqrt(2)*cos((k-1)*pi*s))
  }
}

cosine_fun_expt <- function(s,k=1){
  # cosine * exp(-t)
  if(k==1){
    return(cos(0*s)*exp(-s))
  }else{
    return(sqrt(2)*cos((k-1)*pi*s)*exp(-s))
  }
}

cosine_fun_cost <- function(s,k=1){
  # cosine * exp(-t)
  if(k==1){
    return(cos(0*s)*exp(-s))
  }else{
    return(sqrt(2)*cos((k-1)*pi*s)*cos(pi/3*s))
  }
}

exp_ts_vector <- function(R){ # A(s,t)=exp(-(s+t))=exp(-s)*exp(-t)
  integral_seq <- c()
  for(k in 1:R){
    integral_seq <- c(integral_seq, integrate(cosine_fun_expt, lower=0, upper=1, k=k, rel.tol=1e-10)$value)
  }
  return(integral_seq)
}

cossin_ts_vector <- function(R){ # A(s,t)=sin(pi/5*s)*cos(pi/3*t)
  integral_seq <- c()
  for(k in 1:R){
    integral_seq <- c(integral_seq, integrate(cosine_fun_cost, lower=0, upper=1, k=k, rel.tol=1e-10)$value)
  }
  return(integral_seq)
}

A_transition <- function(x1,x2,Lambda){
  # function to calculate the true transition operator at A(x1,x2)
  R <- dim(Lambda)[1]
  u_x1 <- u_x2 <- c()
  for(r in 1:R){
    u_x1 <- c(u_x1, cosine_fun(x1,k=r))
    u_x2 <- c(u_x2, cosine_fun(x2,k=r))
  }
  value <- c(t(u_x1)%*%Lambda%*%t(t(u_x2)))
  return(value)
}

A_transition_vec <- Vectorize(A_transition, vectorize.args=c('x1','x2'))

beta_transition <- function(x,beta_coef){
  # function to calculate the true beta(t) function at beta(x)
  R <- dim(beta_coef)[1]
  u_x <- c()
  for(r in 1:R){
    u_x <- c(u_x, cosine_fun(x,k=r))
  }
  value <- c(t(beta_coef)%*%u_x)
}

functionValue <- function(s, x_coef, R=6){
  # function to recover the true function at time s
  Xt_value <- rep(0, length(s))
  for(k in 1:R){
    tmp <- cosine_fun(s, k=k)*c(x_coef[k])
    Xt_value <- Xt_value + tmp
  }
  return(Xt_value)
}

functionValue_all <- function(s, basis_coef, R=NULL, type='Random'){
  # function to recover the true function at time s
  if(type%in%c('Diag','Random')){
    f_value <- rep(0, length(s))
    for(k in 1:R){
      tmp <- cosine_fun(s, k=k)*c(basis_coef[k]) # basis_coef is a vector
      f_value <- f_value + tmp
    }
  }
  if(type=='Exp'){
    f_value <- exp(-s)*c(basis_coef) # basis_coef is a scalar
  }
  if(type=='SinCos'){
    f_value <- sin(2*pi*s)*c(basis_coef) # basis_coef is a scalar
  }
  return(f_value)
}

kernel <- function(x1,x2){
  # function for reproducing kernel
  # return value k(x1,x2) = 1+k1(x1)k1(x2)+k2(x1)*k2(x2)-k4(abs(x1-x2)) [see Gu(2013) in Chapter 2.3.3]
  val <- 1+(x1-0.5)*(x2-0.5) #1+k1(x1)k1(x2)
  val <- val+((x1-0.5)^2-1/12)/2*((x2-0.5)^2-1/12)/2 #k2(x1)*k2(x2)
  wk <- abs(x1-x2)
  val <- val-((wk-0.5)^4-(wk-0.5)^2/2+7/240)/24
  return(val)
}

tr <- function(W){
  # trace function
  sum(diag(W))
}

estE_func1 <- function(K1, K2, K1_sqr, K2_sqr, Y, X, lambda){
  m <- dim(X)[1]
  n <- dim(Y)[1]
  S1 <- 1/m*(t(X)%*%K2_sqr) %x% K1_sqr # Kronecker product
  Z <- t(S1)%*%c(Y)
  STS_svd <- svd(t(S1)%*%S1)
  STS_svd$d[STS_svd$d<0] <- 0
  est_e <- STS_svd$v%*%diag(1/(STS_svd$d+lambda))%*%t(STS_svd$v)%*%Z # solve(t(S1)%*%S1+lambda_pen_cv_selected*diag(1,m*n))%*%Z
  E_est <- matrix(est_e, nrow=n, ncol=m)
  A_est <- K1_sqr%*%E_est%*%K2_sqr
  return(A_est)
}

CV_estEfunc1 <- function(K1, K2, K1_sqr, K2_sqr, Y, X, lambda_penalty, cv_indices){
  m <- dim(X)[1]
  cv_error <- c()
  totalT_train <- dim(Y)[2]
  cv_folds <- length(cv_indices)
  for(cv_fold_index in 1:cv_folds){
    print(cv_fold_index)
    cv_index_test <- sort(cv_indices[[cv_fold_index]])
    cv_index_train <- sort(unlist(cv_indices)[!unlist(cv_indices)%in%cv_index_test])
    # CV Training
    Y_tmp <- Y[,cv_index_train]
    X_tmp <- X[,cv_index_train]
    S1_tmp <- 1/m*(t(X_tmp)%*%K2_sqr) %x% K1_sqr # Kronecker product
    Z_tmp <- t(S1_tmp)%*%c(Y_tmp)
    STS_svd_tmp <- svd(t(S1_tmp)%*%S1_tmp)
    STS_svd_tmp$d[STS_svd_tmp$d<0] <- 0
    # CV Testing prediction
    Y_true_tmp <- Y[,cv_index_test]
    
    # 5-fold CV
    tmp_cv_error <- c()
    for(lambda_pen in lambda_penalty) {
      # print(lambda_pen)
      # Estimation
      est_e_tmp <- STS_svd_tmp$v%*%diag(1/(STS_svd_tmp$d+lambda_pen))%*%t(STS_svd_tmp$v)%*%Z_tmp # solve(t(S1_tmp)%*%S1_tmp+lambda_pen*diag(1,m*n))%*%Z_tmp
      E_est_tmp <- matrix(est_e_tmp, nrow=n, ncol=m)
      A_est_tmp <- K1_sqr%*%E_est_tmp%*%K2_sqr
      
      # Prediction
      Y_forecast_rkhs_tmp <- A_est_tmp%*%X[,cv_index_test]/m
      tmp_cv_error <- c(tmp_cv_error, mean(apply(Y_true_tmp-Y_forecast_rkhs_tmp, 2, mse)))
    }
    cv_error <- rbind(cv_error, tmp_cv_error)
  }
  return(cv_error)
}

estE_func2 <- function(K1, K2, K1_sqr, K2_sqr, Y, X, lambda){
  m <- dim(X)[1]
  n <- dim(Y)[1]
  SS1 <- 1/m^2*K2_sqr%*%X%*%t(X)%*%K2_sqr
  SS1_svd <- svd(SS1)
  K1_svd <- svd(K1)
  SS1_svd$d[SS1_svd$d<0] <- 0
  K1_svd$d[K1_svd$d<0] <- 0
  KYX <- K1_sqr%*%Y%*%t(X)%*%K2_sqr
  E_est2 <- matrix(0,n,m)
  for(i in 1:m){
    for(j in 1:n){
      E_est2 <- E_est2+(1/m)*(1/(lambda+SS1_svd$d[i]*K1_svd$d[j]))*(K1_svd$u[,j]%*%t(K1_svd$u[,j]))%*%KYX%*%(SS1_svd$u[,i]%*%t(SS1_svd$u[,i]))
    }
  }
  A_est2 <- K1_sqr%*%E_est2%*%K2_sqr
  return(A_est2)
}

CV_estEfunc2 <- function(K1, K2, K1_sqr, K2_sqr, Y, X, lambda_penalty, cv_indices){
  m <- dim(X)[1]
  n <- dim(Y)[1]
  cv_error <- c()
  cv_folds <- length(cv_indices)
  K1_svd <- svd(K1)
  K1_svd$d[K1_svd$d<0] <- 0
  for(cv_fold_index in 1:cv_folds){
    print(cv_fold_index)
    cv_index_test <- sort(cv_indices[[cv_fold_index]])
    cv_index_train <- sort(unlist(cv_indices)[!unlist(cv_indices)%in%cv_index_test])
    # CV Training
    Y_tmp <- Y[,cv_index_train]
    X_tmp <- X[,cv_index_train]
    
    SS1_tmp <- 1/m^2*K2_sqr%*%X_tmp%*%t(X_tmp)%*%K2_sqr
    SS1_svd_tmp <- svd(SS1_tmp)
    SS1_svd_tmp$d[SS1_svd_tmp$d<0] <- 0
    KYX_tmp <- K1_sqr%*%Y_tmp%*%t(X_tmp)%*%K2_sqr
    
    # CV Testing prediction
    Y_true_tmp <- Y[,cv_index_test]
    
    # 5-fold CV
    tmp_cv_error <- c()
    for(lambda_pen in lambda_penalty) {
      # print(lambda_pen)
      # Estimation
      E_est2_tmp <- matrix(0,n,m)
      for(i in 1:m){
        for(j in 1:n){
          E_est2_tmp <- E_est2_tmp+(1/m)*(1/(lambda_pen+SS1_svd_tmp$d[i]*K1_svd$d[j]))*(K1_svd$u[,j]%*%t(K1_svd$u[,j]))%*%KYX_tmp%*%(SS1_svd_tmp$u[,i]%*%t(SS1_svd_tmp$u[,i]))
        }
      }
      A_est2_tmp <- K1_sqr%*%E_est2_tmp%*%K2_sqr
      
      # Prediction
      Y_forecast_rkhs_tmp2 <- A_est2_tmp%*%X[,cv_index_test]/m
      tmp_cv_error <- c(tmp_cv_error, mean(apply(Y_true_tmp-Y_forecast_rkhs_tmp2, 2, mse)))
    }
    cv_error <- rbind(cv_error, tmp_cv_error)
  }
  return(cv_error)
}

GCV_estEfunc2 <- function(K1, K2, K1_sqr, K2_sqr, Y, X, lambda_penalty, alpha=1.4){
  # This function follows the formulation in SunDu2018JASA
  m <- dim(X)[1]
  n <- dim(Y)[1]
  totalT <- dim(Y)[2]
  K1_svd <- svd(K1)
  K1_svd$d[K1_svd$d<0] <- 0

  SS1 <- 1/m^2*K2_sqr%*%X%*%t(X)%*%K2_sqr # SS1 = S4%*%t(S4)
  SS1_svd <- svd(SS1)
  SS1_svd$d[SS1_svd$d<0] <- 0
  S4 <- (1/m)*K2_sqr%*%X
  KYX <- K1_sqr%*%Y%*%t(S4)
  
  GCV_numerator <- GCV_denominator <- c()
  for(lambda_pen in lambda_penalty) {
    # print(lambda_pen)
    estY <- matrix(0,n,totalT)
    traceA <- 0
    for(i in 1:m){
      for(j in 1:n){
        estY <- estY+(1/(lambda_pen+SS1_svd$d[i]*K1_svd$d[j]))*(K1_sqr%*%K1_svd$u[,j]%*%t(K1_svd$u[,j]))%*%KYX%*%(SS1_svd$u[,i]%*%t(SS1_svd$u[,i])%*%S4)
        traceA <- traceA + (1/(lambda_pen+SS1_svd$d[i]*K1_svd$d[j]))*tr(SS1_svd$u[,i]%*%t(SS1_svd$u[,i])%*%SS1)*tr(K1_svd$u[,j]%*%t(K1_svd$u[,j])%*%K1)
      }
    }
    # GCV
    GCV_numerator <- c(GCV_numerator, sum(Y-estY)^2/n/totalT)
    GCV_denominator <- c(GCV_denominator, ((n*totalT-alpha*traceA)/n/totalT)^2)
  }
  GCV <- GCV_numerator/GCV_denominator
  return(GCV)
}

CV_fda <- function(Y, X, lambda_penalty, cv_indices, nbasisXY, nbasisBeta, domain, sample_pointsY, sample_pointsX){
  # smoothing parameters
  basisY <- basisX <- create.bspline.basis(rangeval=domain, nbasis=nbasisXY, norder=4) # cubic splines
  basiss <- create.bspline.basis(rangeval=domain, nbasis=nbasisBeta)
  basist <- create.bspline.basis(rangeval=domain, nbasis=nbasisBeta)
  Lfdobjt <- Lfdobjs <- vec2Lfd(c(0,0), domain) # roughness penalty
  
  # CV
  cv_error <- c()
  cv_folds <- length(cv_indices)
  for(cv_fold_index in 1:cv_folds){
    cat(cv_fold_index)
    cv_index_test <- sort(cv_indices[[cv_fold_index]])
    cv_index_train <- sort(unlist(cv_indices)[!unlist(cv_indices)%in%cv_index_test])
    # CV Training
    functionalY_tmp <- Data2fd(argvals=sample_pointsY, y=Y[,cv_index_train], basisobj=basisY)
    functionalX_tmp <- Data2fd(argvals=sample_pointsX, y=X[,cv_index_train], basisobj=basisX)
  
    # CV Testing prediction
    Y_test <- Y[,cv_index_test]
    
    # 5-fold CV
    tmp_cv_error <- c()
    for(lambda_pen in lambda_penalty) {
      # Penalization fda estimation
      fda_pen <- fregre.basis.fr(x=functionalX_tmp, y=functionalY_tmp, basis.s=basiss, basis.t=basist,
                                 lambda.s=lambda_pen, lambda.t=lambda_pen, Lfdobj.s=Lfdobjs, Lfdobj.t=Lfdobjt)
      # Prediction
      functionalX_test <- Data2fd(argvals=sample_pointsX, y=X[,cv_index_test], basisobj=basisX)
      functionalY_forecast <- predict(fda_pen, new.fdataobj=functionalX_test)
      Y_forecast_fda <- eval.fd(sample_pointsY, functionalY_forecast)

      tmp_cv_error <- c(tmp_cv_error, mean(apply(Y_test-Y_forecast_fda, 2, mse)))
    }
    cv_error <- rbind(cv_error, tmp_cv_error)
  }
  return(cv_error)
}

### Functions used for coordinate descents
one_dim_h <- function(hk, h, k, n, Qkl, Zl, K1_invk, lambda2, lambda3){
  h[k] <- hk
  value <- sum((Qkl-hk*Zl)^2)+lambda2*(2*hk*sum(K1_invk*h)-K1_invk[k]*hk^2)+lambda3/sqrt(n)*sqrt(sum(h^2))
  return(value)
}

beta_h <- function(h, n, Ql, Zl, K1_inv, lambda2, lambda3){
  value <- sum((Ql-h%*%t(Zl))^2)+lambda2*t(h)%*%K1_inv%*%h+lambda3/sqrt(n)*sqrt(sum(h^2))
  return(value)
}

beta_H <- function(H, n, Q, Z, K1_inv, lambda2, lambda3){
  value <- sum((Q-H%*%Z)^2)+lambda2*tr(t(H)%*%K1_inv%*%H)+lambda3/sqrt(n)*sum(sqrt(diag(t(H)%*%H)))
  return(value)
}

beta_EH <- function(E, H, n, m, Y, X, Z, K1, K1_sqr, K1_inv, K2_sqr, lambda1, lambda2, lambda3){
  value <- sum((Y-1/m*K1_sqr%*%E%*%K2_sqr%*%X-H%*%Z)^2)+lambda1*tr(t(E)%*%E)+lambda2*tr(t(H)%*%K1_inv%*%H)+lambda3/sqrt(n)*sum(sqrt(diag(t(H)%*%H)))
  return(value)
}

coordinate_descent <- function(Y, X, Z, K1, K1_sqr, K1_inv, K2, K2_sqr, lambda1, lambda2, lambda3, converge_tol=1e-8, max_iter=1e5, E=NULL, H=NULL){
  # function to conduct coordinate descent and optimize E and H (and return A_est and H_est)
  # H=K1%*%B, A=K1_sqr%*%E%*%K2_sqr
  n <- dim(Y)[1]
  m <- dim(X)[1]
  p <- dim(Z)[1]
  if(is.null(E)){
    E <- matrix(0,n,m) # initial values of E
  }
  if(is.null(H)){
    H <- matrix(0,n,p) # initial values of H
  }
  ZtZ <- Z%*%t(Z)
  YtZ <- Y%*%t(Z)
  XtZ <- X%*%t(Z)
  
  # optimization starts
  optim_err0 <- beta_EH(E, H, n, m, Y, X, Z, K1, K1_sqr, K1_inv, K2_sqr, lambda1, lambda2, lambda3)
  optim_improve0 <- 1; optim_err0_seq <- c()
  while(optim_improve0>converge_tol){ # first coordinate descent to solve (E, H)
    optim_err0_seq <- c(optim_err0_seq, optim_err0)
    # fix H, optimize E
    Q <- Y-H%*%Z
    E <- estEE_func2_Rcpp(K1, K2, K1_sqr, K2_sqr, Y=Q, X, lambda=lambda1)
    A <- K1_sqr%*%E%*%K2_sqr
    
    # fix E, optimize H
    Q <- Y-1/m*A%*%X
    QtZ <- YtZ-1/m*A%*%XtZ
    optim_err1 <- beta_H(H, n, Q, Z, K1_inv, lambda2, lambda3)
    optim_improve1 <- 1; optim_err1_seq <- c()
    while(optim_improve1>converge_tol){ # First coordinate descent (beta_1, beta_2, ..., beta_p)
      optim_err1_seq <- c(optim_err1_seq, optim_err1)
      for(scalar_index in 1:p){
        zero_check_vec <- QtZ[,scalar_index]-H[,-scalar_index]%*%ZtZ[-scalar_index,scalar_index]
        zero_check <- 2*sqrt(n)/lambda3*sqrt(sum(zero_check_vec^2))
        if(zero_check<1){ # Check group penalty
          H[,scalar_index] <- 0
          next
        }else{
          h <- H[,scalar_index]
          Ql <- Q-H[,-scalar_index]%*%Z[-scalar_index,]
          optim_err2 <- beta_h(h, n, Ql, Zl=Z[scalar_index,], K1_inv, lambda2, lambda3)
          optim_improve2 <- 1; optim_err2_seq <- c()
          while(optim_improve2>converge_tol){ # Second coordinate descent (within beta_l)
            optim_err2_seq <- c(optim_err2_seq, optim_err2)
            for(k in 1:n){
              hk <- optimize(one_dim_h, interval=c(-1e4,1e4), h, k, n, Qkl=Ql[k,], Zl=Z[scalar_index,], K1_invk=K1_inv[k,], lambda2, lambda3)
              h[k] <- hk$minimum
            }
            optim_err2_new <- beta_h(h, n, Ql, Zl=Z[scalar_index,], K1_inv, lambda2, lambda3)
            # optim_improve2 <- 1 - optim_err2_new/optim_err2 # relative error
            optim_improve2 <- optim_err2-optim_err2_new # absolute error
            optim_err2 <- optim_err2_new
          }
          # print(optim_err2_seq)
          H[,scalar_index] <- h
        }
      }
      optim_err1_new <- beta_H(H, n, Q, Z, K1_inv, lambda2, lambda3)
      # optim_improve1 <- 1 - optim_err1_new/optim_err1 # relative error
      optim_improve1 <- optim_err1-optim_err1_new # absolute error
      optim_err1 <- optim_err1_new
    }
    optim_err0_new <- beta_EH(E, H, n, m, Y, X, Z, K1, K1_sqr, K1_inv, K2_sqr, lambda1, lambda2, lambda3)
    # optim_improve0 <- 1 - optim_err0_new/optim_err0 # relative error
    optim_improve0 <- optim_err0-optim_err0_new # absolute error
    optim_err0 <- optim_err0_new
  }
  opt_result <- list(A_est=A, H_est=H, E_est=E, optim_err0_seq=optim_err0_seq)
}

CV_coordinate_descent <- function(Y, X, Z, K1, K1_sqr, K1_inv, K2, K2_sqr, lambda_penalty, cv_indices, converge_tol=1e-8, max_iter=1e5){
  # function to conduct CV  to select lambda_penalty=c(lambda1,lambda2,lambda3)
  cv_error <- c()
  cv_folds <- length(cv_indices)
  m <- dim(X)[1]
  print('CV starts')
  for(cv_fold_index in 1:cv_folds){
    print(paste('CV fold', cv_fold_index))
    cv_index_test <- sort(cv_indices[[cv_fold_index]])
    cv_index_train <- sort(unlist(cv_indices)[!unlist(cv_indices)%in%cv_index_test])
    # CV Training
    Y_tmp <- Y[,cv_index_train]
    X_tmp <- X[,cv_index_train]
    Z_tmp <- Z[,cv_index_train]
  
    # CV Testing prediction
    tmp_cv_error <- c()
    for(lambda_pen_index in 1:length(lambda_penalty)){
      lambda_pen <- lambda_penalty[[lambda_pen_index]]
      if(lambda_pen_index==1){
        rkhs_result_tmp <- coordinate_descent(Y_tmp, X_tmp, Z_tmp, K1, K1_sqr, K1_inv, K2, K2_sqr, lambda_pen[1], lambda_pen[2], lambda_pen[3], converge_tol)
      }else{
        # warm start from previous rkhs_result_tmp
        E_init <- rkhs_result_tmp$E_est
        H_init <- rkhs_result_tmp$H_est
        rkhs_result_tmp <- coordinate_descent(Y_tmp, X_tmp, Z_tmp, K1, K1_sqr, K1_inv, K2, K2_sqr, lambda_pen[1], lambda_pen[2], lambda_pen[3], converge_tol,
                                              E=E_init, H=H_init)
      }
      # Estimation
      A_est_tmp <- rkhs_result_tmp$A_est
      H_est_tmp <- rkhs_result_tmp$H_est
      
      # Prediction
      Y_forecast_rkhs_tmp <- A_est_tmp%*%X[,cv_index_test]/m + H_est_tmp%*%Z[,cv_index_test]
      tmp_cv_error <- c(tmp_cv_error, mean(apply(Y[,cv_index_test]-Y_forecast_rkhs_tmp, 2, mse)))
    }
    cv_error <- rbind(cv_error, tmp_cv_error)
  }
  return(cv_error)
}

arrange_lambda_penalty <- function(lambda_penalty_RKHS_rangeA, lambda_penalty_RKHS_rangeB, lambda_penalty_groupL2_range){
  # Function for arranging the lambda penalty combinations for better warm start
  kA <- length(lambda_penalty_RKHS_rangeA)
  kB <- length(lambda_penalty_RKHS_rangeB)
  kL2 <- length(lambda_penalty_groupL2_range)
  # RKHS penalty for A kernel
  lambda_penalty_RKHS_A <- lambda_penalty_RKHS_B <- lambda_penalty_RKHS_groupL2 <- c()
  lambda_penalty_RKHS_A <- rep(lambda_penalty_RKHS_rangeA, each=kB*kL2)
  
  # RKHS penalty for Beta
  lambda_penalty_RKHS_rangeB <- rep(lambda_penalty_RKHS_rangeB, each=kL2)
  for(times in 1:kA){
    lambda_penalty_RKHS_B <- c(lambda_penalty_RKHS_B, lambda_penalty_RKHS_rangeB)
    lambda_penalty_RKHS_rangeB <- rev(lambda_penalty_RKHS_rangeB)
  }
  
  # RKHS penalty for L2 group
  for(times in 1:(kA*kB)){
    lambda_penalty_RKHS_groupL2 <- c(lambda_penalty_RKHS_groupL2, lambda_penalty_groupL2_range)
    lambda_penalty_groupL2_range <- rev(lambda_penalty_groupL2_range)
  }
  lambda_penalty <- cbind(lambda_penalty_RKHS_A, lambda_penalty_RKHS_B, lambda_penalty_RKHS_groupL2)
  return(lambda_penalty)      
}

coordinate_descent_Beta <- function(Y, Z, K1_inv, lambda2, lambda3, converge_tol=1e-8, max_iter=1e5, H=NULL){
  # function to conduct coordinate descent and optimize E and H (and return A_est and H_est)
  # H=K1%*%B, A=K1_sqr%*%E%*%K2_sqr
  n <- dim(Y)[1]
  p <- dim(Z)[1]
  if(is.null(H)){
    H <- matrix(0,n,p) # initial values of H
  }
  ZtZ <- Z%*%t(Z)
  YtZ <- Y%*%t(Z)

  # optimization starts
  optim_err1 <- beta_H(H, n, Y, Z, K1_inv, lambda2, lambda3)
  optim_improve1 <- 1; optim_err1_seq <- c()
  while(optim_improve1>converge_tol){ # First coordinate descent (beta_1, beta_2, ..., beta_p)
    optim_err1_seq <- c(optim_err1_seq, optim_err1)
    for(scalar_index in 1:p){
      zero_check_vec <- YtZ[,scalar_index]-H[,-scalar_index]%*%ZtZ[-scalar_index,scalar_index]
      zero_check <- 2*sqrt(n)/lambda3*sqrt(sum(zero_check_vec^2))
      if(zero_check<1){ # Check group penalty
        H[,scalar_index] <- 0
        next
      }else{
        h <- H[,scalar_index]
        Ql <- Y-H[,-scalar_index]%*%Z[-scalar_index,]
        optim_err2 <- beta_h(h, n, Ql, Zl=Z[scalar_index,], K1_inv, lambda2, lambda3)
        optim_improve2 <- 1; optim_err2_seq <- c()
        while(optim_improve2>converge_tol){ # Second coordinate descent (within beta_l)
          optim_err2_seq <- c(optim_err2_seq, optim_err2)
          for(k in 1:n){
            hk <- optimize(one_dim_h, interval=c(-1e4,1e4), h, k, n, Qkl=Ql[k,], Zl=Z[scalar_index,], K1_invk=K1_inv[k,], lambda2, lambda3)
            h[k] <- hk$minimum
          }
          optim_err2_new <- beta_h(h, n, Ql, Zl=Z[scalar_index,], K1_inv, lambda2, lambda3)
          # optim_improve2 <- 1 - optim_err2_new/optim_err2 # relative error
          optim_improve2 <- optim_err2-optim_err2_new # absolute error
          optim_err2 <- optim_err2_new
        }
        # print(optim_err2_seq)
        H[,scalar_index] <- h
      }
    }
    optim_err1_new <- beta_H(H, n, Y, Z, K1_inv, lambda2, lambda3)
    # optim_improve1 <- 1 - optim_err1_new/optim_err1 # relative error
    optim_improve1 <- optim_err1-optim_err1_new # absolute error
    optim_err1 <- optim_err1_new
  }
  opt_result <- list(H_est=H, optim_err1_seq=optim_err1_seq)
}