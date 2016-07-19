// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// getFitted
arma::mat getFitted(arma::vec& f_fin, const arma::vec& f_ini, const arma::mat& beta, const arma::vec& alpha, const int& k);
RcppExport SEXP gdpc_getFitted(SEXP f_finSEXP, SEXP f_iniSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec& >::type f_fin(f_finSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type f_ini(f_iniSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    __result = Rcpp::wrap(getFitted(f_fin, f_ini, beta, alpha, k));
    return __result;
END_RCPP
}
// getFini
arma::vec getFini(const arma::mat& Z, const int& k);
RcppExport SEXP gdpc_getFini(SEXP ZSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    __result = Rcpp::wrap(getFini(Z, k));
    return __result;
END_RCPP
}
// gdpc_priv
List gdpc_priv(const arma::mat& Z, const int& k, const double& tol, const int& niter_max, const int& sel);
RcppExport SEXP gdpc_gdpc_priv(SEXP ZSEXP, SEXP kSEXP, SEXP tolSEXP, SEXP niter_maxSEXP, SEXP selSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int& >::type niter_max(niter_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type sel(selSEXP);
    __result = Rcpp::wrap(gdpc_priv(Z, k, tol, niter_max, sel));
    return __result;
END_RCPP
}
