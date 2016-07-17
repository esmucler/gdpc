#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
List getMatrixBeta(const arma::mat & Z, const arma::vec & f, const int & k, const int & sel) {
  // This function finds the optimal beta and alpha and the mse corresponding to Z, f and k. 
  int m = Z.n_rows;
  int N = Z.n_cols;
  arma::mat beta = mat(m, k + 2);
  arma::mat Fmat = mat(k + 1, N);
  arma::mat FF = mat(k + 2, k + 2);
  arma::mat invFF = mat(k + 2, k + 2);
  for ( int i = 0; i < N; i++){
    Fmat.col(i) = f.subvec(i, i + k);
  }
  Fmat.insert_rows(k + 1, ones(1, N));
  FF = Fmat * Fmat.t();
  double condition_FF = rcond(FF);
  if(condition_FF > 1e-10){
    invFF = inv_sympd(FF);
  } else{
    invFF = pinv(FF);
  }
  arma::mat Proj = Fmat.t() * invFF * Fmat;
  beta = Z * Fmat.t() * invFF.t();
  arma::mat res = Z - beta * Fmat;
  double mse = (accu(pow(Z, 2)) - trace(Z * Proj * Z.t() ))/ N;
  double crit = 0;
  if (sel == 1){
    crit = N * log(mse) + m * (k + 2) * 2;  
  } else {
    crit = N * log(mse) + m * (k + 2) * log(N);  
  }
  
  List ret;
  ret["mse"] = mse / m;
  ret["beta"] = beta;
  ret["res"] = res;
  ret["crit"] = crit;
  return(ret);
}

arma::mat getMatrixC(const arma::subview_row<double> & rowZ, const double & alpha, const int & k) {
  //This function constructs the matrix C correspoding to rowZ, alpha and k 
  int N = rowZ.n_elem;
  arma::mat C = zeros(N + k, k + 1);
  arma::vec liml = vec(2);
  liml[0] = 0;
  arma::vec limu = vec(2);
  limu[0] = k;
  for ( int t = 1; t <= N + k ; t++){
    liml[1] = t - N;
    limu[1] = t - 1;
    for (int q = 1; q <= k + 1; q++){
      if( (q >= max(liml) + 1) & (q<= min(limu) + 1)){
        C(t - 1, q - 1) = rowZ(t - q) - alpha;
      }
    }
  }
  return (C);
}

arma::mat getMatrixD(const arma::subview_row<double> & rowbeta, const int & N, const int & k) {
  //This function constructs the matrix D correspoding to rowbeta, N and k 
  arma::mat beta_mat = zeros(N + k, N + k);
  arma::vec liml = vec(2);
  liml[1] = 1;
  arma::vec limu = vec(2);
  limu[1] = N;
  for ( int t = 1; t <= N + k ; t++){
    liml[0] = t - k;
    limu[0] = t;
    for (int r = max(liml); r <= min(limu); r++){
      for (int q = r; q <= r + k; q++){
        beta_mat(t - 1, q - 1) = beta_mat(t - 1, q - 1) + rowbeta(q - r) * rowbeta(t - r);
      }
    }
  }
  return (beta_mat);
}

// [[Rcpp::export]]
arma::vec getF(const arma::mat & Z, const arma::mat & beta, const int & k) {
  //This functions finds the optimal f corresponding to Z, beta, alpha and k.
  int m = Z.n_rows;
  int N = Z.n_cols;
  arma::vec f = zeros(N + k);
  arma::mat D = zeros(N + k, N + k);
  for (int j=0 ; j < m ; j++){
    D = D + getMatrixD(beta(j, span(0, k)), N, k);
    f = f + getMatrixC(Z.row(j), beta(j, k + 1), k) * (beta(j, span(0, k))).t();
  }
  
  double condition_D = rcond(D);
  if (condition_D > 1e-10) {
    f = solve(D,f);
  } else {
    f = pinv(D) * f;
  }
  return(f);
}

// [[Rcpp::export]]
arma::mat getFitted(arma::vec & f_fin, const arma::vec & f_ini, const arma::mat & beta, const arma::vec & alpha, const int & k) {
  // This function finds the fitted values associated with f and beta and alpha
  int N = f_fin.n_elem;
  if (k > 0){
    f_fin.insert_rows(0, f_ini);
  }
  arma::mat Fmat = mat(k + 1, N);
  for ( int i = 0; i < N; i++){
    Fmat.col(i) = f_fin.subvec(i, i + k);
  }
  Fmat.insert_rows(k + 1, ones(1, N));
  arma::mat betalpha = fliplr(beta);
  betalpha.insert_cols(k + 1, alpha);
  arma::mat fit = Fmat.t() * betalpha.t();
  return(fit);
}

// [[Rcpp::export]]
arma::vec getFini(const arma::mat & Z, const int & k){
  // Get initial estimator: ordinary principal component with k leads
  int N = Z.n_cols;
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::vec f_ini = vec(N + k);
  arma::mat Z_trans = Z.t();
  arma::rowvec mean_Zt = mean(Z_trans);
  Z_trans.each_row() -= mean_Zt; 
  svd_econ(U, s, V, Z_trans, "right");
  f_ini.rows(0, N - 1) = Z_trans * V.col(0);
  if (k != 0) {
    f_ini.rows(N, N + k - 1) = zeros(k, 1) + f_ini(N - 1);
  }
  f_ini = (f_ini - mean(f_ini)) / stddev(f_ini);
  return(f_ini);
}