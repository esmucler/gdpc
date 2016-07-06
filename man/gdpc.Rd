\name{gdpc}
\alias{gdpc}
\title{Generalized Dynamic Principal Components}
\description{
Function to compute a single Dynamic Principal Component with a given number of lags.
}


\usage{gdpc(Z, k, f_ini = NULL, tol = 1e-4, niter_max = 500, crit = 'AIC')}
\arguments{
  \item{Z}{Data matrix. Each column is a different time series.}
  \item{k}{Integer. Number of lags to use.}
  \item{f_ini}{(Optional) Numeric vector. Starting point for the iterations. If no argument is passed, the standard first Principal Component with k lags is used.}
  \item{tol}{Relative precision. Default is 1e-4.}
  \item{niter_max}{Integer. Maximum number of iterations. Default is 500.}
  \item{crit}{A string specifying the criterion to be used to evaluate the fitted model. Either 'AIC' or 'BIC'. Default is 'AIC'.}
}

\details{

See \code{\link{auto.gdpc}} for the definition of the AIC or BIC type criterion that is part of the output of this function.

}

\value{
  A list with entries:
 \item{f}{Coordinates of the first Principal Component corresponding to the periods \eqn{1,\dots,T}.}
  \item{initial_f}{Coordinates of the first Principal Component corresponding to the periods \eqn{-k+1,\dots,0}. Only for the case \eqn{k>0}.}
 \item{beta}{Matrix beta corresponding to f.}
 \item{alpha}{Vector alpha corresponding to f.}
 \item{mse}{Mean squared error.}
 \item{crit}{The criterion specified by the user.}
 \item{res}{Matrix of residuals.}
 }
\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai
}
\references{
Peña and Yohai (2016) "Generalized Dynamic Principal Components". Upcoming on Journal of the American
Statistical Association.
}


\seealso{
\code{\link{auto.gdpc}}
}

\examples{
N <- 200 #length of series
m <- 500 #number of series
set.seed(1234)
f <- rnorm(N + 1)
x <- matrix(0, N, m)
u <- matrix(rnorm(N * m), N, m)
for (i in 1:m) {
    x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:N] + 10 * cos(2 * pi * (i/m)) * f[2:(N + 1)] + u[, i]
}
fit <- gdpc(x, k = 1) #find first DPC with one lag
par(mfrow = c(1, 2)) #plot loadings
plot(fit$beta[, 1], xlab = "variable", ylab = "loadings", type = "l")
title("Zero lag loadings") 
plot(fit$beta[, 2], xlab = "variable", ylab = "loadings", type = "l")
title("One lag loadings")
}

\keyword{ts}
