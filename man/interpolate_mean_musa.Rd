% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolate_mean_musa.R
\name{interpolate_mean_musa}
\alias{interpolate_mean_musa}
\title{Interpolate Accumulated Mean Values at an Unobserved Location}
\usage{
interpolate_mean_musa(results, data, sites, Xw, Sites1, gridt, Xwr)
}
\arguments{
\item{results}{A list containing the output from \code{STModelMusaOkumotoMCMC}, including:
\itemize{
\item \code{MW}: Samples for parameter W.
\item \code{Malpha}: Samples for parameter alpha.
\item \code{Mb}: Samples for parameter phi.
\item \code{Mv}: Samples for parameter sigma^2.
\item \code{MPsi}: Samples for parameter Psi.
}}

\item{data}{A matrix of occurrence times for the event of interest. Each column corresponds to a monitoring station.}

\item{sites}{A matrix of geographic coordinates where the process was observed.}

\item{Xw}{Covariates associated with the scale parameter at the location.}

\item{Sites1}{A vector representing the geographic coordinates of the location where interpolation is to be performed.}

\item{gridt}{A vector of times at which the accumulated mean is to be estimated.}

\item{Xwr}{Covariates associated with the scale parameter at the location where interpolation is to be performed.}
}
\value{
A matrix (\code{MatMean}) containing the interpolated accumulated mean values at the times specified in \code{gridt}.
}
\description{
This function interpolates the accumulated mean values for a location where no data was observed,
using outputs from the STModelMusaOkumotoMCMC model. The interpolation is performed over a vector of times.
}
