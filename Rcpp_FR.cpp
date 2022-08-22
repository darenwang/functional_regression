// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;

// [[Rcpp::export]]
List svd_Rcpp(arma::mat X, bool dc){
  arma::mat u, v;
  arma::vec d;
  if(dc){
    arma::svd(u, d, v, X, "dc"); // divide and conquer svd (only works for n < 40)
  }else{
    arma::svd(u, d, v, X, "std"); // standard svd
  }
  List result = List::create(Named("u")=u, Named("v")=v, Named("d")=d);
  return result;
}

// [[Rcpp::export]]
arma::mat estE_func1_Rcpp(arma::mat K1, arma::mat K2, arma::mat K1_sqr, arma::mat K2_sqr, arma::mat Y, arma::mat X, double lambda){
  int m = X.n_rows;
  int n = Y.n_rows;
  // Rprintf("%f%f", m, n);
  
  double m_float = X.n_rows;
  arma::mat S1 = arma::kron(1/m_float*(X.t()*K2_sqr), K1_sqr); // Kronecker product
  arma::mat Z = S1.t()*Y.as_col();
  List STS_svd = svd_Rcpp(S1.t()*S1, false);
  // NumericVector d;
  // d = STS_svd["d"];
  // d = pmax(0, d) + lambda;
  arma::vec d = STS_svd["d"];
  arma::vec dd = 1/(lambda+d);
  arma::mat v = STS_svd["v"];
  
  arma::vec est_e = v*arma::diagmat(dd)*v.t()*Z; //solve(t(S1)%*%S1+lambda_pen_cv_selected*diag(1,m*n))%*%Z
  arma::mat E_est = reshape(est_e, n, m);
  arma::mat A_est = K1_sqr*E_est*K2_sqr;
  return A_est;
}

// [[Rcpp::export]]
arma::mat estE_func2_Rcpp(arma::mat K1, arma::mat K2, arma::mat K1_sqr, arma::mat K2_sqr, arma::mat Y, arma::mat X, double lambda){
  int m = X.n_rows;
  int n = Y.n_rows;
  // Rprintf("%f%f", m, n);
  
  double m_float = X.n_rows;
  arma::mat SS1 = 1/m_float/m_float*K2_sqr*X*X.t()*K2_sqr;
  List SS1_svd = svd_Rcpp(SS1, false);
  List K1_svd = svd_Rcpp(K1, false);
  arma::mat KYX = K1_sqr*Y*X.t()*K2_sqr;
  arma::vec SS1_svd_d = SS1_svd["d"];
  arma::vec K1_svd_d = K1_svd["d"];
  arma::mat SS1_svd_u = SS1_svd["u"];
  arma::mat K1_svd_u = K1_svd["u"];
  
  arma::mat E_est(n, m, arma::fill::zeros);
  for(int iter1 = 0; iter1 < m; ++iter1){
    for(int iter2 = 0; iter2 < n; ++iter2){
      E_est = E_est+(1/m_float)*(1/(lambda+SS1_svd_d(iter1)*K1_svd_d(iter2))*(K1_svd_u.col(iter2)*K1_svd_u.col(iter2).t())*KYX*(SS1_svd_u.col(iter1)*SS1_svd_u.col(iter1).t()));
    }
  }
  arma::mat A_est = K1_sqr*E_est*K2_sqr;
  return A_est;
}

// [[Rcpp::export]]
arma::mat CV_estEfunc2_Rcpp(arma::mat K1, arma::mat K2, arma::mat K1_sqr, arma::mat K2_sqr, arma::mat Y, arma::mat X, arma::vec lambda_penalty, List cv_data){
  int m = X.n_rows;
  int n = Y.n_rows;
  double m_float = X.n_rows;
  int number_lambda = lambda_penalty.size();
  int cv_folds = cv_data.length();
  arma::mat cv_error(cv_folds, number_lambda);
  List K1_svd = svd_Rcpp(K1, false);
  arma::vec K1_svd_d = K1_svd["d"];
  arma::mat K1_svd_u = K1_svd["u"];
  arma::vec cv_index_test, cv_index_train;
  arma::mat SS1_tmp, KYX_tmp;
  for(int cv_fold_index = 0; cv_fold_index < cv_folds; ++cv_fold_index){
    Rprintf("%u", cv_fold_index);
    List cv_data_tmp = cv_data[cv_fold_index];
    
    // CV Training
    arma::mat Y_tmp = cv_data_tmp["Y_train"];
    arma::mat X_tmp = cv_data_tmp["X_train"];
    
    SS1_tmp = 1/m_float/m_float*K2_sqr*X_tmp*X_tmp.t()*K2_sqr;
    List SS1_svd_tmp = svd_Rcpp(SS1_tmp, false);
    arma::vec SS1_svd_d_tmp = SS1_svd_tmp["d"];
    arma::mat SS1_svd_u_tmp = SS1_svd_tmp["u"];
    KYX_tmp = K1_sqr*Y_tmp*X_tmp.t()*K2_sqr;
    // CV Testing prediction
    arma::mat Y_test_tmp = cv_data_tmp["Y_test"];
    arma::mat X_test_tmp = cv_data_tmp["X_test"];
      
    // 5-fold CV
    arma::vec tmp_cv_error(number_lambda);
    double lambda_pen;
    for(int iter_lambda = 0; iter_lambda < number_lambda; ++iter_lambda){
      lambda_pen = lambda_penalty(iter_lambda);
      // Estimation
      arma::mat E_est_tmp(n, m, arma::fill::zeros);
      for(int iter1 = 0; iter1 < m; ++iter1){
        for(int iter2 = 0; iter2 < n; ++iter2){
          E_est_tmp = E_est_tmp+(1/m_float)*(1/(lambda_pen+SS1_svd_d_tmp(iter1)*K1_svd_d(iter2))*(K1_svd_u.col(iter2)*K1_svd_u.col(iter2).t())*KYX_tmp*(SS1_svd_u_tmp.col(iter1)*SS1_svd_u_tmp.col(iter1).t()));
        }
      }
      arma::mat A_est_tmp = K1_sqr*E_est_tmp*K2_sqr;
      
      // Prediction
      arma::mat Y_forecast_rkhs_tmp = A_est_tmp*X_test_tmp/m_float;
      tmp_cv_error(iter_lambda) = arma::mean(arma::mean(arma::square(Y_test_tmp-Y_forecast_rkhs_tmp), 0)); // mse
    }
    cv_error.row(cv_fold_index) = tmp_cv_error.t();
  }
  return cv_error;
}

// [[Rcpp::export]]
arma::mat estEE_func2_Rcpp(arma::mat K1, arma::mat K2, arma::mat K1_sqr, arma::mat K2_sqr, arma::mat Y, arma::mat X, double lambda){
  int m = X.n_rows;
  int n = Y.n_rows;
  // Rprintf("%f%f", m, n);
  
  double m_float = X.n_rows;
  arma::mat SS1 = 1/m_float/m_float*K2_sqr*X*X.t()*K2_sqr;
  List SS1_svd = svd_Rcpp(SS1, false);
  List K1_svd = svd_Rcpp(K1, false);
  arma::mat KYX = K1_sqr*Y*X.t()*K2_sqr;
  arma::vec SS1_svd_d = SS1_svd["d"];
  arma::vec K1_svd_d = K1_svd["d"];
  arma::mat SS1_svd_u = SS1_svd["u"];
  arma::mat K1_svd_u = K1_svd["u"];
  
  arma::mat E_est(n, m, arma::fill::zeros);
  for(int iter1 = 0; iter1 < m; ++iter1){
    for(int iter2 = 0; iter2 < n; ++iter2){
      E_est = E_est+(1/m_float)*(1/(lambda+SS1_svd_d(iter1)*K1_svd_d(iter2))*(K1_svd_u.col(iter2)*K1_svd_u.col(iter2).t())*KYX*(SS1_svd_u.col(iter1)*SS1_svd_u.col(iter1).t()));
    }
  }
  return E_est;
}