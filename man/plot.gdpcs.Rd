\name{plot.gdpcs}
\alias{plot.gdpcs}
\title{Plot Generalized Dynamic Principal Components}
\description{
Plots a \code{gdpcs} object.
}


\usage{
\method{plot}{gdpcs}(x, which_comp = 1, \dots)
}
\arguments{
  \item{x}{An object of class \code{gdpcs}, usually the result of \link{auto.gdpc}.}
  \item{which_comp}{Numeric vector indicating which components to plot. Default is 1.}
  \item{\dots}{Additional arguments to be passed to the plotting functions.}
}


\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai
}
\references{
Peña and Yohai (2016) "Generalized Dynamic Principal Components". Upcoming on Journal of the American
Statistical Association.
}


\seealso{
\code{\link{gdpc}}, \code{\link{auto.gdpc}}, \code{\link{plot.gdpc}}
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
#k_max=2 to keep computation time low
autofit <- auto.gdpc(x, normalize = FALSE, k_max = 2, auto_comp = FALSE, num_comp = 2) 
plot(autofit, which_comp = c(1,2), xlab = '', ylab = '')
}

\keyword{ts}