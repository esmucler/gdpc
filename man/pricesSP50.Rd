\name{pricesSP50}
\alias{pricesSP50}
\docType{data}
\title{
Stock Prices of the First 50 Components of S&P500
}
\description{
Fifty series corresponding to the stock prices of the first 50 components of the Standard&Poor's 500 index. Five hundred daily observations starting 1/1/2010.
}
\usage{data(pricesSP50)}
\format{
  A matrix time series with 500 observations on the stock prices of the first 50 components of the Standard&Poor's 500 index.
}



\examples{
data(pricesSP50)
\dontrun{
#This may take a bit.
gdpc_SP50 <- auto.gdpc(pricesSP50, normalize = 2)
}
}
\keyword{datasets}