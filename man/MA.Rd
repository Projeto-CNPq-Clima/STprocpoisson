% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MA.R
\docType{data}
\name{MA}
\alias{MA}
\title{Border Map of Maranhão State, Brazil}
\format{
An \code{sf} object with the following details:
\describe{
\item{geometry}{A polygon representing the boundary of Maranhão.}
\item{CRS}{The coordinate reference system is EPSG:4326 (WGS84).}
}
}
\source{
Official cartographic data from Brazilian government agencies.
}
\usage{
data(MA)
}
\description{
The \code{MA} dataset provides spatial data representing the border of Maranhão state, Brazil.
It is useful for spatial analyses and modeling in combination with other data from the
\code{STprocpoisson} package. This dataset serves as a geographical base map.
}
\examples{
library(STprocpoisson)
data(MA)
plot(MA, main = "Maranhão State Border")
}
\keyword{datasets}
