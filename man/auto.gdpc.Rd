\name{auto.gdpc}
\alias{auto.gdpc}
\title{Automatic Fitting of Generalized Dynamic Principal Components}
\description{
Computes Generalized Dynamic Principal Components. The number of components can be supplied by the user or chosen automatically so that a given proportion of variance is explained. The number of lags is chosen automatically using an AIC or a BIC type criterion.
}


\usage{auto.gdpc(Z, crit = 'AIC', normalize = TRUE, auto_comp = TRUE, expl_var = 0.9,
          num_comp = 5, tol = 1e-4, k_max = 10, niter_max = 500, ncores = 1)}
\arguments{
  \item{Z}{Data matrix. Each column is a different time series.}
  \item{crit}{A string specifying the criterion to be used. Either 'AIC' or 'BIC'. Default is 'AIC'.}
  \item{normalize}{Logical. If TRUE the data is standardized to zero mean and unit variance. Default is TRUE.}
  \item{auto_comp}{Logical. If TRUE compute components until the proportion of explained variance is equal to expl_var, otherwise use num_comp components. Default is TRUE.}
  \item{expl_var}{A number between 0 and 1. Desired proportion of explained variance (only used if auto_comp==TRUE). Default is 0.9.}
  \item{num_comp}{Integer. Number of components to be computed (only used if auto_comp==FALSE). Default is 5.}
  \item{tol}{Relative precision. Default is 1e-4.}
  \item{k_max}{Integer. Maximum possible number of lags. Default is 10.}
  \item{niter_max}{Integer. Maximum number of iterations. Default is 500.}
  \item{ncores}{Integer. Number of cores to be used for parallel computations. Default is 1.}
}

\details{

Suppose the data matrix consists of \eqn{m} series of length \eqn{T}. If crit = 'AIC' the number of lags is chosen among \eqn{0,\dots, k_{max}} as the value \eqn{k_{opt}} that minimizes the following AIC type criterion
\deqn{ AIC = T \log(trace(\Sigma)) + 2 m (k+2), }
where \eqn{\Sigma = (res^{T} res) / T} and \eqn{res} is the matrix of residuals of the fit using \eqn{k} lags. 

If crit = 'BIC' the number of lags is chosen among \eqn{0,\dots, k_{max}} as the value \eqn{k_{opt}} that minimizes the following BIC type criterion
\deqn{ BIC = T \log(trace(\Sigma)) + m (k+2) \log(T), }
where \eqn{\Sigma = (res^{T} res) / T} and \eqn{res} is the matrix of residuals of the fit using \eqn{k} lags. 


}

\value{
  An object of class \code{gdpcs}, that is, a list of length equal to the number of computed components. The i-th entry of this list is an object of class \code{gdpc}, that is, a list with entries
  \item{f}{Coordinates of the i-th Principal Component corresponding to the periods \eqn{1,\dots,T}.}
  \item{initial_f}{Coordinates of the i-th Principal Component corresponding to the periods \eqn{-k+1,\dots,0}. Only for the case \eqn{k>0}.}
  \item{beta}{Matrix of loadings corresponding to f. Column number \eqn{k} is the vector of \eqn{k-1} lag loadings.}
  \item{alpha}{Vector of intercepts corresponding to f.}
  \item{mse}{Mean squared error of the residuals of the fit using the first i components.}
  \item{crit}{The AIC or BIC of the fitted model, according to what the user specified.}
  \item{k_opt}{Number of lags chosen.}
  \item{expart}{Proportion of the variance explained by the first i components.}
  \item{res}{Matrix of residuals of the fit using the first i components.}
    \item{fitted}{Matrix of fitted values of the fit using the first i components.}
      \item{call}{The matched call.}
}
\author{
  Daniel Peña, Ezequiel Smucler, Victor Yohai
}
\references{
  Peña and Yohai (2016) "Generalized Dynamic Principal Components". Upcoming on Journal of the American Statistical Association.
}

\seealso{
\code{\link{gdpc}}, \code{\link{plot.gdpc}}, \code{\link{plot.gdpcs}}
}
\examples{
N <- 200 #length of series
m <- 200 #number of series
set.seed(1234)
f <- rnorm(N + 1)
x <- matrix(0, N, m)
u <- matrix(rnorm(N * m), N, m)
for (i in 1:m) {
    x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:N] + 10 * cos(2 * pi * (i/m)) * f[2:(N + 1)] + u[, i]
}
#Choose number of lags using the AIC type criterion.
#Data is not normalized.
#k_max=3 to keep computation time low
autofit <- auto.gdpc(x, normalize = FALSE, k_max = 3)
autofit
resid <- residuals(autofit[[1]]) #Get residuals of the fit
fit_val <- fitted(autofit[[1]]) #Get fitted values
plot.gdpcs(autofit, which_comp = 1)
}

\keyword{ts}