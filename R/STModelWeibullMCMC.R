#' Spatiotemporal Nonhomogeneous Poisson Model with Weibull Intensity (MCMC)
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
#' @param X A list of covariates for the scale parameter of the Weibull intensity function. Defaults to a matrix combining
#' a column of ones (intercept) and the coordinates in `sites`.
#'
#' @param Z A list of covariates for the shape parameter of the Weibull intensity function. Defaults to the same
#' covariates as `X`.
#'
#' @param prior A list of hyperparameters for the prior distributions of the model parameters:
#'   \describe{
#'     \item{sigma^2_w}{Inverse Gamma distribution parameters `aa1` and `bb1`. Default: `aa1=0.001`, `bb1=0.001`.}
#'     \item{sigma^2_m}{Inverse Gamma distribution parameters `aa2` and `bb2`. Default: `aa2=0.001`, `bb2=0.001`.}
#'     \item{phi_w}{Gamma distribution parameters `c1` and `d1`. Default: `c1=(-2*log(0.05)/max(dist(sites)))*0.1`, `d1=0.1`.}
#'     \item{phi_m}{Gamma distribution parameters `c2` and `d2`. Default: `c2=c1`, `d2=d1`.}
#'     \item{Psi}{Multivariate normal distribution parameters `A1` (mean vector) and `B1` (covariance matrix). Default: `A1=0`, `B1=diag(100, ncol(X))`.}
#'     \item{Beta}{Multivariate normal distribution parameters `A` (mean vector) and `B` (covariance matrix). Default: `A=0`, `B=diag(100, ncol(Z))`.}
#'   }
#'
#' @param iteration Number of MCMC iterations.
#'
#' @param burnin Number of burn-in iterations to discard before the chains converge.
#'
#' @return A list containing:
#' \describe{
#'   \item{MMj}{Samples of parameter M obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{MMT}{A binary vector indicating whether each proposed value of M was accepted (1) or rejected (0). Used to compute the acceptance rate for M.}
#'   \item{MW}{Samples of parameter W obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{MWT}{A binary vector indicating whether each proposed value of W was accepted (1) or rejected (0). Used to compute the acceptance rate for W.}
#'   \item{MPsi}{Samples of parameter Psi obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{MBeta}{Samples of parameter Beta obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{Mvw}{Samples of parameter sigma^2_w obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{Mbw}{Samples of parameter phi_w obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{MbwT}{A binary vector indicating whether each proposed value of phi_w was accepted (1) or rejected (0). Used to compute the acceptance rate for phi_w.}
#'   \item{Mvm}{Samples of parameter sigma^2_m obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{Mbm}{Samples of parameter phi_m obtained during the MCMC procedure (`iteration - burnin`).}
#'   \item{MbmT}{A binary vector indicating whether each proposed value of phi_m was accepted (1) or rejected (0). Used to compute the acceptance rate for phi_m.}
#' }
#'
#' @export

STModelWeibullMCMC <- function(data, sites, X = cbind(as.matrix(rep(1, ncol(data))), as.matrix(sites)), Z = X,
                               prior = list(
                                 aa1 = 0.001,
                                 bb1 = 0.001,
                                 aa2 = 2.01,
                                 bb2 = 1.005,
                                 c1 = (-2 * log(0.05) / max(dist(sites))) * 0.1,
                                 d1 = 0.1,
                                 c2 = (-2 * log(0.05) / max(dist(sites))) * 0.1,
                                 d2 = 0.1,
                                 A1 = as.matrix(rep(0, ncol(X))),
                                 B1 = diag(100, ncol(X)),
                                 A = as.matrix(rep(0, ncol(Z))),
                                 B = diag(100, ncol(Z))), iteration, burnin) {
  # Valores iniciais
  bw <- 1
  vw <- 1
  bm <- 1
  vm <- 1

  M <- as.matrix(rep(log(0.89), ncol(data)))
  W <- as.matrix(rep(0, ncol(data)))

  # hiperparametros

  # hiperparametros
  aa1 <- prior$aa1
  bb1 <- prior$bb1
  aa2 <- prior$aa2
  bb2 <- prior$bb2
  c1 <- prior$c1
  d1 <- prior$d1
  c2 <- prior$c2
  d2 <- prior$d2
  A1 <- prior$A1
  B1 <- prior$B1
  A <- prior$A
  B <- prior$B


  SU1 <- 0.000001
  SU2 <- 0.000001
  SU3 <- 100
  SU4 <- 100
  Psi <- as.matrix(rep(0, ncol(X)))
  Beta <- as.matrix(rep(0, ncol(Z)))
  n <- ncol(data)
  m <- nrow(data)
  tempdados <- is.na(data)
  nj <- m - apply(tempdados, 2, sum)

  Tt <- array(NA, dim = c(1, n))
  for (y in 1:n) {
    Tt[1, y] <- data[nj[y], y]
  }


  MMj <- NULL
  MMT <- NULL

  MW <- NULL
  MWT <- NULL

  MPsi <- NULL
  MBeta <- NULL

  Mvw <- NULL
  Mbw <- NULL
  MbwT <- NULL

  Mvm <- NULL
  Mbm <- NULL
  MbmT <- NULL


  for (j in 1:iteration) {
    if (j <= burnin) {
      temp <- amostrarW(W, M, sites, X, Psi, bw, vw, nj, Tt, SU1)
      W <- as.matrix(temp[[1]])
      MWT <- c(MWT, temp[[2]])

      temp <- amostrarM(W, M, sites, Z, Beta, bm, vm, data, nj, Tt, SU2)
      M <- as.matrix(temp[[1]])
      MMT <- c(MMT, temp[[2]])

      RR <- gCorr(bw, sites)
      aaW <- (n / 2) + aa1
      bbW <- 0.5 * t(W - X %*% Psi) %*% solve(RR) %*% (W - X %*% Psi) + bb1

      vw <- 1 / rgamma(1, shape = aaW, rate = bbW)

      RR1 <- gCorr(bm, sites)
      aaM <- (n / 2) + aa2
      bbM <- 0.5 * t(M - Z %*% Beta) %*% solve(RR1) %*% (M - Z %*% Beta) + bb2

      vm <- 1 / rgamma(1, shape = aaM, rate = bbM)

      varPsi <- solve(solve(B) + t(X) %*% solve(gSigma(bw, vw, sites)) %*% X)
      medPsi <- (t(A) %*% solve(B) + t(W) %*% solve(gSigma(bw, vw, sites)) %*% X) %*% varPsi
      Psi <- as.matrix(MASS::mvrnorm(1, medPsi, varPsi))


      varBeta <- solve(solve(B) + t(Z) %*% solve(gSigma(bm, vm, sites)) %*% Z)
      medBeta <- (t(A) %*% solve(B) + t(M) %*% solve(gSigma(bm, vm, sites)) %*% Z) %*% varBeta
      Beta <- as.matrix(MASS::mvrnorm(1, medBeta, varBeta))

      temp <- amostrarb(W, vw, bw, sites, c1, d1, X, Psi, SU3)
      bw <- temp[[1]]
      MbwT <- c(MbwT, temp[[2]])

      temp <- amostrarb(M, vm, bm, sites, c2, d2, Z, Beta, SU4)
      bm <- temp[[1]]
      MbmT <- c(MbmT, temp[[2]])


      if ((j %% 50) == 0) {
        SU1 <- sintonizarN(burnin, 0.15, SU1, MWT, j)
        SU2 <- sintonizarN(burnin, 0.15, SU2, MMT, j)
        SU3 <- sintonizar(burnin, 0.44, SU3, MbwT, j)
        SU4 <- sintonizar(burnin, 0.44, SU4, MbmT, j)
      } else {

      }





    } else {
      temp <- amostrarW(W, M, sites, X, Psi, bw, vw, nj, Tt, SU1)
      W <- as.matrix(temp[[1]])
      MW <- rbind(MW, t(W))
      MWT <- c(MWT, temp[[2]])

      temp <- amostrarM(W, M, sites, Z, Beta, bm, vm, data, nj, Tt, SU2)
      M <- as.matrix(temp[[1]])
      MMj <- rbind(MMj, t(M))
      MMT <- c(MMT, temp[[2]])

      RR <- gCorr(bw, sites)
      aaW <- (n / 2) + aa1
      bbW <- 0.5 * t(W - X %*% Psi) %*% solve(RR) %*% (W - X %*% Psi) + bb1

      vw <- 1 / rgamma(1, shape = aaW, rate = bbW)
      Mvw <- c(Mvw, vw)

      RR1 <- gCorr(bm, sites)
      aaM <- (n / 2) + aa2
      bbM <- 0.5 * t(M - Z %*% Beta) %*% solve(RR1) %*% (M - Z %*% Beta) + bb2

      vm <- 1 / rgamma(1, shape = aaM, rate = bbM)
      Mvm <- c(Mvm, vm)

      varPsi <- solve(solve(B) + t(X) %*% solve(gSigma(bw, vw, sites)) %*% X)
      medPsi <- (t(A) %*% solve(B) + t(W) %*% solve(gSigma(bw, vw, sites)) %*% X) %*% varPsi
      Psi <- as.matrix(MASS::mvrnorm(1, medPsi, varPsi))
      MPsi <- rbind(MPsi, t(Psi))

      varBeta <- solve(solve(B) + t(Z) %*% solve(gSigma(bm, vm, sites)) %*% Z)
      medBeta <- (t(A) %*% solve(B) + t(M) %*% solve(gSigma(bm, vm, sites)) %*% Z) %*% varBeta
      Beta <- as.matrix(MASS::mvrnorm(1, medBeta, varBeta))
      MBeta <- rbind(MBeta, t(Beta))

      temp <- amostrarb(W, vw, bw, sites, c1, d1, X, Psi, SU3)
      bw <- temp[[1]]
      Mbw <- c(Mbw, bw)
      MbwT <- c(MbwT, temp[[2]])

      temp <- amostrarb(M, vm, bm, sites, c2, d2, Z, Beta, SU4)
      bm <- temp[[1]]
      Mbm <- c(Mbm, bm)
      MbmT <- c(MbmT, temp[[2]])



      resul <- list(MMj, MMT, MW,MWT,MPsi,MBeta,Mvw,Mbw,MbwT, Mvm,Mbm, MbmT)
      names(resul) <- c("MMj", "MMT", "MW","MWT","MPsi","MBeta","Mvw","Mbw","MbwT","Mvm","Mbm","MbmT")
    }
  }
  return(resul)
}
