#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
List betaf(arma::mat & Z, arma::vec & f, int & k, int & sel) {
  // This function finds the optimal beta and alpha and the mse (only in N) corresponding to Z, f and k. 
  int m = Z.n_rows;
  int N = Z.n_cols;
  arma::mat beta = mat(m,k+1);
  arma::vec alpha = vec(m);
  arma::mat Fmat = mat(k+1, N);
  arma::mat FF = mat(k+2,k+2);
  arma::mat invFF = mat(k+2,k+2);
  for ( int i=0; i<N; i++){
    Fmat.col(i) = f.subvec(i,i+k);
  }
  Fmat.insert_rows(k+1,ones(1,N));
  FF = Fmat*Fmat.t();
  double condition_FF = cond( FF );
  if(condition_FF<1e10){
    invFF = inv_sympd(FF);
  } else{
    invFF = pinv(FF);
  }
  arma::mat Proj = Fmat.t()*invFF*Fmat;
  beta = Z * Fmat.t() * invFF.t();
  arma::mat res = Z - beta * Fmat;
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

arma::mat matrix_C(arma::rowvec & betav, double & alfa, int & k) {
  //This function constructs the matrix C correspoding to betav, alpha and k 
  int N = betav.n_elem;
  arma::mat C = zeros(N + k, k + 1);
  arma::vec liml = vec(2);
  liml[0] = 0;
  arma::vec limu = vec(2);
  limu[0] = k;
  for ( int t = 1; t <= N+k ; t++){
    liml[1] = t-N;
    limu[1] = t-1;
    for (int q = 1; q <= k+1; q++){
      if( (q>= max(liml)+1) & (q<=min(limu)+1) ){
        C(t-1,q-1)=betav(t-q)-alfa;
      }
    }
  }
  return (C);
}

arma::mat matrix_D(arma::rowvec & betav, int & N, int & k) {
  //This function constructs the matrix D correspoding to betav, N and k 
  arma::mat beta_mat = zeros(N + k, N + k);
  arma::vec liml = vec(2);
  liml[1] = 1;
  arma::vec limu = vec(2);
  limu[1] = N;
  for ( int t = 1; t <= N+k ; t++){
    liml[0] = t-k;
    limu[0] = t;
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
  //This functions finds the optimal f corresponding to Z, beta, alpha and k.
  int m = Z.n_rows;
  int N = Z.n_cols;
  arma::vec f = zeros(N + k);
  arma::mat D = zeros(N + k, N + k);
  arma::rowvec betaj = zeros<rowvec>(k+1);
  arma::rowvec zetaj = zeros<rowvec>(N);
  for (int j=0 ; j < m ; j++){
    //There has to be a better way to do this...
    betaj = beta.row(j);
    zetaj = Z.row(j);
    D = D + matrix_D(betaj, N, k);
    f = f + matrix_C(zetaj, alpha(j), k) * betaj.t();
  }
  
  double condition_D = cond(D);
  if (condition_D<1e10) {
    f = solve(D,f);
  } else {
    f = pinv(D) * f;
  }
  return(f);
}

// [[Rcpp::export]]
arma::mat fits(arma::vec & f_fin, arma::vec & f_ini,arma::mat beta, arma::vec alpha, int & k) {
  // This function finds the fitted values associated with f and beta and alpha
  int N = f_fin.n_elem;
  f_fin.insert_rows(0, f_ini);
  arma::mat Fmat = mat(k+1, N);
  for ( int i=0; i<N; i++){
    Fmat.col(i) = f_fin.subvec(i,i+k);
  }
  Fmat.insert_rows(k+1,ones(1,N));
  arma::mat betalpha = fliplr(beta);
  betalpha.insert_cols(k+1, alpha);
  arma::mat fit = Fmat.t() * betalpha.t();
  return(fit);
}