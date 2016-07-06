#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
List betaf(arma::mat & Z, arma::rowvec & f, int & k, int & sel) {
  // This function finds the beta, alpha and mse corresponding to Z, f and k. 
  int m = Z.n_rows;
  int N = Z.n_cols;
  arma::mat beta = mat(m,k+1);
  arma::vec alpha = vec(m);
  arma::mat Fmat = mat(N, k+1);
  arma::mat FF = mat(k+2,k+2);
  arma::mat invFF = mat(k+2,k+2);
  for ( int i=0; i<N; i++){
    Fmat.row(i) = f.subvec(i,i+k);
  }
  Fmat.insert_cols(k+1,ones(N,1));
  FF = Fmat.t()*Fmat;
  double condition_FF = cond( FF );
  if(condition_FF<1e10){
    invFF = inv_sympd(FF);
  } else{
    invFF = pinv(FF);
  }
  arma::mat Proj = Fmat*invFF*Fmat.t();
  beta = Z * Fmat * invFF.t();
  arma::mat res = Z - beta * Fmat.t();
  alpha = beta.col(k + 1);
  beta.shed_col(k+1);
  double mse = (accu(pow(Z,2)) - trace(Z * Proj * Z.t() ))/ N;
  double crit = 0;
  if (sel == 1){
    crit = N * log( mse ) + m * (k+2) * 2;  
  } else {
    crit = N * log( mse ) + m * (k+2) * log(N);  
  }
  
  List ret;
  ret["mse"] = mse;
  ret["alpha"] = alpha;
  ret["beta"] = beta;
  ret["res"] = res;
  ret["crit"] = crit;
  return(ret);
}

// [[Rcpp::export]]
arma::mat matrix_C(arma::rowvec & betav, double & alfa, int & k) {
  //This function constructs the matrix C correspoding to betav, alpha and k 
  int N = betav.n_elem;
  arma::mat C = zeros(N + k, k + 1);
  for ( int t = 1; t <= N+k ; t++){
    arma::vec liml = vec(2);
    liml[0] = 0;
    liml[1] = t-N;
    arma::vec limu = vec(2);
    limu[0] = k;
    limu[1] = t-1;
    for (int q = 1; q <= k+1; q++){
      if( (q>= max(liml)+1) & (q<=min(limu)+1) ){
        C(t-1,q-1)=betav(t-q)-alfa;
      }
    }
  }
  return (C);
}

// [[Rcpp::export]]
arma::mat matrix_D(arma::rowvec & betav, int & N, int & k) {
  //This function constructs the matrix C correspoding to betav, N and k 
  arma::mat beta_mat = zeros(N + k, N + k);
  for ( int t = 1; t <= N+k ; t++){
    arma::vec liml = vec(2);
    liml[0] = t-k;
    liml[1] = 1;
    arma::vec limu = vec(2);
    limu[0] = t;
    limu[1] = N;
    for (int r=max(liml); r<=min(limu); r++){
      for (int q = r; q <= r+k; q++){
        beta_mat(t-1,q-1) = beta_mat(t-1,q-1)+betav(q-r)*betav(t-r);
      }
    }
  }
  return (beta_mat);
}

// [[Rcpp::export]]
arma::vec matrix_ff(arma::mat & Z, arma::mat & beta, arma::vec & alpha, int & k) {
  int m = Z.n_rows;
  int N = Z.n_cols;
  arma::vec f = zeros(N+k);
  arma::mat D = zeros(N+k, N+k);
  arma::rowvec betaj = zeros<rowvec>(k+1);
  arma::rowvec zetaj = zeros<rowvec>(N);
  for (int j=0; j<m; j++){
    //There has to be a better way to do this...
    betaj = beta.row(j);
    zetaj = Z.row(j);
    D = D + matrix_D(betaj, N, k);
    f = f + matrix_C(zetaj, alpha(j), k) * betaj.t();
  }
  
  double condition_D = cond( D );
  if(condition_D<1e10){
    f = solve(D,f);
  } else{
    f = pinv(D) * f;
  }
  return(f);
}