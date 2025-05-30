% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STModelWeibullMCMC.R
\name{STModelWeibullMCMC}
\alias{STModelWeibullMCMC}
\title{Spatiotemporal Nonhomogeneous Poisson Model with Weibull Intensity (MCMC)}
\usage{
STModelWeibullMCMC(
  data,
  sites,
  X = cbind(as.matrix(rep(1, ncol(data))), as.matrix(sites)),
  Z = X,
  prior = list(aa1 = 0.001, bb1 = 0.001, aa2 = 2.01, bb2 = 1.005, c1 = (-2 *
    log(0.05)/max(dist(sites))) * 0.1, d1 = 0.1, c2 = (-2 * log(0.05)/max(dist(sites))) *
    0.1, d2 = 0.1, A1 = as.matrix(rep(0, ncol(X))), B1 = diag(100, ncol(X)), A =
    as.matrix(rep(0, ncol(Z))), B = diag(100, ncol(Z))),
  iteration,
  burnin
)
}
\arguments{
\item{data}{A matrix representing the occurrence times of the event of interest at each monitoring station. Each column corresponds
to the occurrence times of a specific station. Dimensions: mxn, where \code{m} is the maximum number of occurrences
across the stations, and \code{n} is the number of monitoring stations.}

\item{sites}{A matrix with the geographic coordinates of the monitoring stations. Dimensions: nx2.}

\item{X}{A list of covariates for the scale parameter of the Weibull intensity function. Defaults to a matrix combining
a column of ones (intercept) and the coordinates in \code{sites}.}

\item{Z}{A list of covariates for the shape parameter of the Weibull intensity function. Defaults to the same
covariates as \code{X}.}

\item{prior}{A list of hyperparameters for the prior distributions of the model parameters:
\describe{
\item{sigma^2_w}{Inverse Gamma distribution parameters \code{aa1} and \code{bb1}. Default: \code{aa1=0.001}, \code{bb1=0.001}.}
\item{sigma^2_m}{Inverse Gamma distribution parameters \code{aa2} and \code{bb2}. Default: \code{aa2=0.001}, \code{bb2=0.001}.}
\item{phi_w}{Gamma distribution parameters \code{c1} and \code{d1}. Default: \code{c1=(-2*log(0.05)/max(dist(sites)))*0.1}, \code{d1=0.1}.}
\item{phi_m}{Gamma distribution parameters \code{c2} and \code{d2}. Default: \code{c2=c1}, \code{d2=d1}.}
\item{Psi}{Multivariate normal distribution parameters \code{A1} (mean vector) and \code{B1} (covariance matrix). Default: \code{A1=0}, \code{B1=diag(100, ncol(X))}.}
\item{Beta}{Multivariate normal distribution parameters \code{A} (mean vector) and \code{B} (covariance matrix). Default: \code{A=0}, \code{B=diag(100, ncol(Z))}.}
}}

\item{iteration}{Number of MCMC iterations.}

\item{burnin}{Number of burn-in iterations to discard before the chains converge.}
}
\value{
A list containing:
\describe{
\item{MMj}{Samples of parameter M obtained during the MCMC procedure (\code{iteration - burnin}).}
\item{MMT}{A binary vector indicating whether each proposed value of M was accepted (1) or rejected (0). Used to compute the acceptance rate for M.}
\item{MW}{Samples of parameter W obtained during the MCMC procedure (\code{iteration - burnin}).}
\item{MWT}{A binary vector indicating whether each proposed value of W was accepted (1) or rejected (0). Used to compute the acceptance rate for W.}
\item{MPsi}{Samples of parameter Psi obtained during the MCMC procedure (\code{iteration - burnin}).}
\item{MBeta}{Samples of parameter Beta obtained during the MCMC procedure (\code{iteration - burnin}).}
\item{Mvw}{Samples of parameter sigma^2_w obtained during the MCMC procedure (\code{iteration - burnin}).}
\item{Mbw}{Samples of parameter phi_w obtained during the MCMC procedure (\code{iteration - burnin}).}
\item{MbwT}{A binary vector indicating whether each proposed value of phi_w was accepted (1) or rejected (0). Used to compute the acceptance rate for phi_w.}
\item{Mvm}{Samples of parameter sigma^2_m obtained during the MCMC procedure (\code{iteration - burnin}).}
\item{Mbm}{Samples of parameter phi_m obtained during the MCMC procedure (\code{iteration - burnin}).}
\item{MbmT}{A binary vector indicating whether each proposed value of phi_m was accepted (1) or rejected (0). Used to compute the acceptance rate for phi_m.}
}
}
\description{
Performs a Markov Chain Monte Carlo (MCMC) procedure to estimate the parameters of a spatiotemporal
nonhomogeneous Poisson model. This model is designed for analyzing extreme rainfall,
as proposed by Fidel Ernesto Castro Morales & Daniele Torres Rodrigues.
}
