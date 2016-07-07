\name{ipi91}
\alias{ipi91}
\docType{data}
\title{
Industrial Production Index (IPI) of France, Germany, Italy, United Kingdom, USA and Japan
}
\description{
Six series corresponding to the Industrial Production Index (IPI) of France, Germany, Italy, United Kingdom, USA and Japan. Monthly data from January 1991 to December 2012.
}
\usage{data(ipi91)}
\format{
  A data frame with 264 observations on the following 6 variables.
 \describe{
    \item{\code{France}}{IPI of France.}
    \item{\code{Germany}}{IPI of Germany.}
    \item{\code{Italy}}{IPI of Italy.}
    \item{\code{United Kingdom}}{IPI of United Kingdom.}
    \item{\code{USA}}{IPI of USA.}
    \item{\code{Japan}}{IPI of Japan.}
  }
}

\references{
Peña and Yohai (2016) "Generalized Dynamic Principal Components". Upcoming on Journal of the American
Statistical Association.
}

\examples{
data(ipi91)

\dontrun{
#Compute first GDPC with twelve lags; this may take a bit.
gdpc_ipi <- gdpc(ipi91, 12, niter_max = 1000)
plot(gdpc_ipi, which='Component') #plot the component
}
}
\keyword{datasets}