% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STModelGoelMCMC.R
\name{STModelGoelMCMC}
\alias{STModelGoelMCMC}
\title{STModelGoelMCMC: A Bayesian Space-Time Model Using MCMC}
\usage{
STModelGoelMCMC(
  data,
  sites,
  X = cbind(as.matrix(rep(1, ncol(data))), as.matrix(sites)),
  Z = X,
  M = Z,
  prior = list(c3 = (-2 * log(0.05)/max(dist(sites))) * 0.1, d3 = 0.1, BB1 = diag(100,
    ncol(Z)), AA1 = as.matrix(rep(0, ncol(Z))), BB2 = diag(100, ncol(Z)), AA2 =
    as.matrix(rep(0, ncol(Z))), aa1 = 2.01, bb1 = 1.005, V = diag(100, ncol(X)), MM1 =
    as.matrix(rep(0, ncol(X)))),
  iteration,
  burnin
)
}
\arguments{
\item{data}{A matrix of failure times, where each column represents a station.}

\item{sites}{A matrix of geographic coordinates for the stations (e.g., longitude and latitude).}

\item{X}{A matrix of covariates associated with the parameter W. Default is a column of ones and the coordinates from \code{sites}.}

\item{Z}{A matrix of covariates associated with the parameter beta. Default is the same as \code{X}.}

\item{M}{A matrix of covariates associated with the parameter alpha. Default is the same as \code{Z}.}

\item{prior}{A list specifying the hyperparameters for the prior distributions:
\itemize{
\item \code{c3}, \code{d3}: Shape and rate parameters for the prior of phi_w.
\item \code{BB1}, \code{AA1}: Mean vector and covariance matrix for the prior of gamma.
\item \code{BB2}, \code{AA2}: Mean vector and covariance matrix for the prior of eta. Defaults are equal to \code{BB1} and \code{AA1}.
\item \code{aa1}, \code{bb1}: Shape and scale parameters for the inverse-Gamma prior of sigma^2_w.
\item \code{V}, \code{MM1}: Covariance matrix and mean vector for the prior of Psi.
}}

\item{iteration}{The total number of MCMC iterations.}

\item{burnin}{The number of burn-in iterations to discard.}
}
\value{
A list containing the following components:
\itemize{
\item \code{Mgama}: MCMC chain for the parameter gamma.
\item \code{MgamaT}: A vector of zeros and ones indicating acceptance (1) or rejection (0) of gamma proposals in the Metropolis-Hastings algorithm.
\item \code{Meta}: MCMC chain for the parameter eta.
\item \code{MetaT}: A vector of zeros and ones indicating acceptance (1) or rejection (0) of eta proposals in the Metropolis-Hastings algorithm.
\item \code{Mv}: MCMC chain for the parameter sigma^2.
\item \code{Mb}: MCMC chain for the parameter phi.
\item \code{MbT}: A vector of zeros and ones indicating acceptance (1) or rejection (0) of phi proposals in the Metropolis-Hastings algorithm.
\item \code{MW}: MCMC chain for the parameter W.
\item \code{MWT}: A vector of zeros and ones indicating acceptance (1) or rejection (0) of W proposals in the Metropolis-Hastings algorithm.
\item \code{MPsi}: MCMC chain for the parameter Psi.
\item \code{MPsiT}: A vector of zeros and ones indicating acceptance (1) or rejection (0) of Psi proposals in the Metropolis-Hastings algorithm.
}
}
\description{
This function implements a Bayesian space-time model using MCMC for failure time data across geographical locations.
It estimates parameters associated with different components of the model, including covariates, spatial dependencies, and prior distributions.
}
