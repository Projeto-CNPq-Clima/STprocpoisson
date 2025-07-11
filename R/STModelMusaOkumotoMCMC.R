#' Spatiotemporal Nonhomogeneous Poisson Model with Musa Okumoto Intensity (MCMC)
#'
#' Performs a Markov Chain Monte Carlo (MCMC) procedure to estimate the parameters of a spatiotemporal
#' nonhomogeneous Poisson model. This model is designed for analyzing extreme rainfall,
#' as proposed by Fidel Ernesto Castro Morales & Daniele Torres Rodrigues.
#'
#' @param data A matrix representing the occurrence times of the event of interest at each monitoring station. Each column corresponds
#' to the occurrence times of a specific station. Dimensions: mxn, where `m` is the maximum number of occurrences
#' across the stations, and `n` is the number of monitoring stations.
#'
#' @param sites A matrix with the geographic coordinates of the monitoring stations. Dimensions: nx2.
#'
#' @param X A list of covariates for the scale parameter of the Musa Okumoto intensity function. Defaults to a matrix combining
#' a column of ones (intercept) and the coordinates in `sites`.
#'
#'
#' @param prior A list of hyperparameters for the prior distributions of the model parameters:
#'   \describe{
#'     \item{Psi}{Multivariate normal distribution parameters `M` (mean vector) and `V` (covariance matrix). Default: `A1=0`, `B1=diag(100, ncol(X))`.}
#'     \item{alpha}{Gamma distribution parameters `c2` and `d2`. Default:`c2=1e-05`,`d2=0.001`.}
#'     \item{sigma^2_w}{Inverse Gamma distribution parameters `c4` and `d4`. Default: `c4=2.01`, `d4=1.005`.}
#'     \item{phi_w}{Gamma distribution parameters `c3` and `d3`. Default: `c3=(-2*log(0.05)/max(dist(sites)))*0.1`, `d3=0.1`.}
#'
#'   }
#'
#' @param iteration Number of MCMC iterations.
#'
#' @param burnin Number of burn-in iterations to discard before the chains converge.
#'
#' @return A list containing:
#' \describe{
#'   \item{MW}{Samples of parameter W obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{MWT}{A binary vector indicating whether each proposed value of W was accepted (1) or rejected (0). Used to compute the acceptance rate for W.}
#'   \item{Malpha}{Samples of parameter alpha obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{MalphaT}{A binary vector indicating whether each proposed value of alpha was accepted (1) or rejected (0). Used to compute the acceptance rate for alpha.}
#'   \item{Mb}{Samples of parameter phi_w obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{MbT}{A binary vector indicating whether each proposed value of phi_w was accepted (1) or rejected (0). Used to compute the acceptance rate for phi_w.}
#'   \item{Mv}{Samples of parameter sigma^2_w obtained during the MCMC procedure (`iteration - burnin`).}'
#'   \item{MPsi}{Samples of parameter Psi obtained during the MCMC procedure (`iteration - burnin`).}
#'     }
#'
#' @export

STModelMusaOkumotoMCMC <- function(data, sites, X = cbind(as.matrix(rep(1, ncol(data))), (1 / 100) * sites),
                                   prior = list(
                                     V = diag(100, ncol(X)),
                                     M = as.matrix(rep(0, ncol(X))),
                                     c2 = 1e-05,
                                     d2 = 0.001,
                                     c4 = 2.01,
                                     d4 = 1.005,
                                     c3 = (-2 * log(0.05) / max(dist(sites))) * 0.1,
                                     d3 = 0.1
                                   ), iteration, burnin) {
  b <- 1
  v <- 1
  alpha <- 1
  W <- as.matrix(rep(0, ncol(data)))
  data <- as.matrix(data)
  sites <- as.matrix(sites)

  # Hiperparametros



  c2 <- prior$c2
  d2 <- prior$d2

  c4 <- prior$c4
  d4 <- prior$d4


  c3 <- prior$c3
  d3 <- prior$d3

  SU2 <- 1000
  SU3 <- 28.01968
  SU5 <- 0.001118574

  Psi <- as.matrix(rep(0, ncol(X)))

  V <- prior$V
  M <- prior$M

  n <- ncol(data)
  m <- nrow(data)
  tempdados <- is.na(data)
  nj <- m - apply(tempdados, 2, sum)
  # p=ncol(X)



  Tt <- array(NA, dim = c(1, n))
  for (y in 1:n) {
    Tt[1, y] <- data[nj[y], y]
  }

  # pp=ncol(X)
  Tt <- t(Tt)

  Malpha <- NULL
  MalphaT <- NULL

  MW <- NULL
  MWT <- NULL

  MPsi <- NULL

  Mv <- NULL
  Mb <- NULL
  MbT <- NULL


  #############################
  ## Programa principal
  ############################ 3
  for (j in 1:iteration) {
    if (j <= burnin) {
      SIGMA <- gSigma(b, v, sites)
      DELTA <- gCorr(b, sites)

      temp <- amostrarWMUSA(W, sites, X, Psi, b, v, nj, SU5, Tt, alpha)
      W <- as.matrix(temp[[1]])
      MWT <- c(MWT, temp[[2]])

      temp <- amostraralphaMUSA(alpha, W, Tt, c2, d2, data, SU2)
      alpha <- temp[[1]]
      MalphaT <- c(MalphaT, temp[[2]])

      temp <- amostrarb(W, v, b, sites, c3, d3, X, Psi, SU3)
      b <- temp[[1]]
      MbT <- c(MbT, temp[[2]])


      aa1 <- 0.5 * ncol(data) + c4
      bb1 <- 0.5 * t(W - X %*% Psi) %*% solve(DELTA) %*% (W - X %*% Psi) + d4
      v <- 1 / rgamma(1, aa1, bb1)


      AA <- solve(solve(V) + t(X) %*% solve(SIGMA) %*% X) %*% (solve(V) %*% M + t(X) %*% solve(SIGMA) %*% W)
      BB <- solve(solve(V) + t(X) %*% solve(SIGMA) %*% X)

      Psi <- as.matrix(MASS::mvrnorm(1, AA, BB))


      SU2 <- sintonizarMUSA(burnin, 0.30, SU2, MalphaT, j)
      SU3 <- sintonizarMUSA(burnin, 0.44, SU3, MbT, j)
      SU5 <- sintonizarNMUSA(burnin, 0.25, SU5, MWT, j)
      print(j)
    } else {
      SIGMA <- gSigma(b, v, sites)
      DELTA <- gCorr(b, sites)

      temp <- amostrarWMUSA(W, sites, X, Psi, b, v, nj, SU5, Tt, alpha)
      W <- as.matrix(temp[[1]])
      MW <- rbind(MW, t(W))
      MWT <- c(MWT, temp[[2]])

      temp <- amostraralphaMUSA(alpha, W, Tt, c2, d2, data, SU2)
      alpha <- temp[[1]]
      Malpha <- c(Malpha, alpha)
      MalphaT <- c(MalphaT, temp[[2]])

      temp <- amostrarb(W, v, b, sites, c3, d3, X, Psi, SU3)
      b <- temp[[1]]
      Mb <- c(Mb, b)
      MbT <- c(MbT, temp[[2]])

      aa1 <- 0.5 * ncol(data) + c4
      bb1 <- 0.5 * t(W - X %*% Psi) %*% solve(DELTA) %*% (W - X %*% Psi) + d4
      v <- 1 / rgamma(1, aa1, bb1)
      Mv <- c(Mv, v)


      AA <- solve(solve(V) + t(X) %*% solve(SIGMA) %*% X) %*% (solve(V) %*% M + t(X) %*% solve(SIGMA) %*% W)
      BB <- solve(solve(V) + t(X) %*% solve(SIGMA) %*% X)

      Psi <- MASS::mvrnorm(1, AA, BB)
      MPsi <- rbind(MPsi, t(Psi))




      print(j)
    }
  }
  resul <- list(MW, MWT, Malpha, MalphaT, Mb, MbT, Mv, MPsi)
  names(resul) <- c("MW", "MWT", "Malpha", "MalphaT", "Mb", "MbT", "Mv", "MPsi")
  return(resul)
}
