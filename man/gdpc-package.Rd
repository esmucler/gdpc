\name{gdpc-package}
\alias{gdpc-package}
\docType{package}

\title{Generalized Dynamic Principal Components}

\description{
Functions to compute the Generalized Dynamic Principal Components proposed in Peña and Yohai (2016).
}
\details{
\tabular{ll}{
Package: \tab gdpc\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2016-7-1\cr
License: \tab GPL (>= 2)\cr
Imports: \tab Rcpp, MASS, parallel, doParallel, foreach\cr
LinkingTo: \tab Rcpp, RcppArmadillo\cr
NeedsCompilation: \tab yes\cr
}

Index:
\preformatted{
gdpc    Function to compute a single Dynamic Principal Component with
        a given number of lags.
}
\preformatted{
auto.gdpc   Function to compute Dynamic Principal Components, number of 
            components and lags can be chosen automatically.
}
\preformatted{
ipi91   Data frame. Six series corresponding to the Industrial Production
        Index (IPI) of France, Germany, Italy, United Kingdom, USA and 
        Japan. Monthly data from January 1991 to December 2012.
}
}
\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai

Maintainer: Ezequiel Smucler <ezequiels.90@gmail.com>
}

\references{
Peña and Yohai (2016) "Generalized Dynamic Principal Components". Upcoming on Journal of the American
Statistical Association.
}


\examples{
data(ipi91)
\dontrun{
  #Compute GDPC, number of components and number of lags is chosen automatically.
  #This might take a bit.
  ipi_autogdpc <- auto.gdpc(ipi91)
  }
}
\keyword{ts}
\keyword{package}