% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_mean_surfaceGOELSA.R
\name{compute_mean_surfaceGOELSA}
\alias{compute_mean_surfaceGOELSA}
\title{Compute Interpolated Mean Surface for Spatiotemporal Model}
\usage{
compute_mean_surfaceGOELSA(
  resultsSA,
  sites,
  X,
  Z,
  M,
  DNO,
  CovXNO,
  CovZNO,
  CovMNO,
  tau
)
}
\arguments{
\item{resultsSA}{A list containing the output from \code{STModelGoelMCMCSA}, including:
\itemize{
\item \code{MW}: Samples for parameter W.
\item \code{MWT}: Acceptance indicators for parameter W.
\item \code{MMj}: Samples for parameter M.
\item \code{MMT}: Acceptance indicators for parameter M.
\item \code{Mvw}: Samples for parameter sigma^2_w.
\item \code{Mvm}: Samples for parameter sigma^2_m.
\item \code{MBeta}: Samples for parameter Beta.
\item \code{Mbw}: Samples for parameter phi_w.
\item \code{MbwT}: Acceptance indicators for phi_w.
\item \code{Mbm}: Samples for parameter phi_m.
\item \code{MbmT}: Acceptance indicators for phi_m.
\item \code{MPsi}: Samples for parameter Psi.
\item \code{Mdelta}: Samples of parameter delta obtained during the MCMC procedure (iteration - burnin).
\item \code{Mtheta}: Samples of parameter theta obtained during the MCMC procedure (iteration - burnin).
\item \code{Mf}: Samples of parameter f obtained during the MCMC procedure (iteration - burnin).
}}

\item{sites}{A matrix with geographic coordinates of the monitoring stations.}

\item{X}{Covariates for the W parameter of the Goel intensity.}

\item{Z}{Covariates for the eta parameter of the Goel intensity.}

\item{M}{Covariates for the gama parameter of the Goel intensity.}

\item{DNO}{A grid of points where interpolation is to be performed.}

\item{CovXNO}{Covariates for the W parameter at the grid points.}

\item{CovZNO}{Covariates for the eta parameter at the grid points.}

\item{CovMNO}{Covariates for the gama parameter at the grid points.}

\item{tau}{A vector of temporal points for which the mean surface is computed.}
}
\value{
A list: (\code{Surface}) containing the interpolated mean values at the grid points
for each temporal point in \code{tau}. The first column contains the mean values at the initial
time step, and subsequent columns contain the differences between consecutive time steps.
}
\description{
This function interpolates the mean values on a grid of points (\code{DNO}) for a spatiotemporal
nonhomogeneous Poisson model with Goel intensity. The function utilizes the MCMC outputs
from \code{STModelGoelMCMCSA} and applies a Gaussian process-based interpolation.
}
