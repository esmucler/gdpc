# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

getFitted <- function(f_fin, f_ini, beta, alpha, k) {
    .Call(`_gdpc_getFitted`, f_fin, f_ini, beta, alpha, k)
}

getFini <- function(Z, k) {
    .Call(`_gdpc_getFini`, Z, k)
}

gdpc_priv <- function(Z, k, f_ini, passf_ini, tol, niter_max, sel) {
    .Call(`_gdpc_gdpc_priv`, Z, k, f_ini, passf_ini, tol, niter_max, sel)
}

