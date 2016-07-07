\name{plot.gdpc}
\alias{plot.gdpc}
\alias{plot}
\title{Plot Dynamic Principal Components}
\description{
Plots a \code{gdpc} object.
}


\usage{
\method{plot}{gdpc}(x, which = 'Loadings', which_load = 0, \dots)
}
\arguments{
  \item{x}{An object of class \code{gdpc}, the result of \link{gdpc} or one of the entries of the result of \link{auto.gdpc}.}
  \item{which}{String. Indicates what to plot, either 'Component' or 'Loadings'. Default is 'Loadings'.}
  \item{which_load}{Lag number indicating which loadings should be plotted. Only used if which = 'Loadings'. Default is 0.}
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
\code{\link{gdpc}}, \code{\link{auto.gdpc}}
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
plot(autofit[[1]], which = 'Loadings', which_load = 0)
}

\keyword{ts}
\keyword{ts}