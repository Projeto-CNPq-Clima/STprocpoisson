#' Compute Interpolated Mean Surface for Spatiotemporal Model
#'
#' This function interpolates the mean values on a grid of points (`DNO`) for a spatiotemporal
#' nonhomogeneous Poisson model with Goel intensity. The function utilizes the MCMC outputs
#' from `STModelGoelMCMCSA` and applies a Gaussian process-based interpolation.
#'
#' @param resultsSA A list containing the output from `STModelGoelMCMCSA`, including:
#'   - `MW`: Samples for parameter W.
#'   - `MWT`: Acceptance indicators for parameter W.
#'   - `MMj`: Samples for parameter M.
#'   - `MMT`: Acceptance indicators for parameter M.
#'   - `Mvw`: Samples for parameter sigma^2_w.
#'   - `Mvm`: Samples for parameter sigma^2_m.
#'   - `MBeta`: Samples for parameter Beta.
#'   - `Mbw`: Samples for parameter phi_w.
#'   - `MbwT`: Acceptance indicators for phi_w.
#'   - `Mbm`: Samples for parameter phi_m.
#'   - `MbmT`: Acceptance indicators for phi_m.
#'   - `MPsi`: Samples for parameter Psi.
#'   - `Mdelta`: Samples of parameter delta obtained during the MCMC procedure (iteration - burnin).
#'   - `Mtheta`: Samples of parameter theta obtained during the MCMC procedure (iteration - burnin).
#'   - `Mf`: Samples of parameter f obtained during the MCMC procedure (iteration - burnin).
#'
#' @param sites A matrix with geographic coordinates of the monitoring stations.
#' @param X Covariates for the W parameter of the Goel intensity.
#' @param Z Covariates for the eta parameter of the Goel intensity.
#' @param M Covariates for the gama parameter of the Goel intensity.
#' @param DNO A grid of points where interpolation is to be performed.
#' @param CovXNO Covariates for the W parameter at the grid points.
#' @param CovZNO Covariates for the eta parameter at the grid points.
#' @param CovMNO Covariates for the gama parameter at the grid points.
#' @param tau A vector of temporal points for which the mean surface is computed.
#'
#' @return A list: (`Surface`) containing the interpolated mean values at the grid points
#' for each temporal point in `tau`. The first column contains the mean values at the initial
#' time step, and subsequent columns contain the differences between consecutive time steps.
#' @export

compute_mean_surfaceGOELSA <- function(resultsSA, sites, X, Z, M, DNO, CovXNO, CovZNO, CovMNO, tau) {

  X<-as.matrix(X)
  Z<-as.matrix(Z)
  M<-as.matrix(M)

  jj <- nrow(DNO)
  res <- rbind(DNO, as.matrix(sites))
  tt <- nrow(res)

  Xr <- rbind(CovXNO, X)
  Zr <- rbind(CovXNO, Z)
  Mr <- rbind(CovZNO, M)

  MW <- resultsSA$MW
  Mvw <- resultsSA$Mv
  Mbw <- resultsSA$Mb
  Mgama <- resultsSA$Mgama
  Meta <- resultsSA$Meta
  MPsi<- resultsSA$MPsi
  Mdelta<-resultsSA$Mdelta
  Mtheta<-resultsSA$Mtheta
  Mf<-resultsSA$Mf

  Eta<-exp(CovZNO%*%t(Mgama))

  Gama<-exp(CovMNO%*%t(Meta))

  MWNO <- array(NA, dim = c(nrow(MPsi), jj))
  MfmNO <- NULL


  for (h in 1:nrow(MW)) {
    WM <- as.matrix(MW[h, ])
    sigma <- gSigma(Mbw[h], Mvw[h], res)
    Psih <- as.matrix(MPsi[h, ])

    Mpro <- Xr %*% as.matrix(Psih)

    A1 <- as.matrix(Mpro[(jj + 1):tt, ])
    A2 <- as.matrix(Mpro[1:jj, ])

    SSA1 <- sigma[(jj + 1):tt, (jj + 1):tt]
    SSA2 <- sigma[1:jj, 1:jj]
    SSA12 <- sigma[1:jj, (jj + 1):tt]

    A2est <- A2 + SSA12 %*% solve(SSA1) %*% (as.matrix(WM) - A1)
    SSA2est <- SSA2 - SSA12 %*% solve(SSA1) %*% t(SSA12)
    WNO <- as.matrix(MASS::mvrnorm(1, A2est, SSA2est))
    MWNO[h, ] <- t(WNO)


  }


  Meanmf <- NULL
  for (i in 1:length(tau)) {
    for (h in 1:nrow(MW)) {
      MfmNO <- rbind(MfmNO, t(as.matrix(mfGOELSA(MWNO[h, ], Gama[h],Eta[h], tau[i],Mdelta[h],Mf[h],Mtheta[h]))))
    }

    Meanmf <- cbind(Meanmf, as.matrix(apply(MfmNO, 2, mean)))
    MfmNO <- NULL
  }


  Surface <- array(NA, dim = c(nrow(DNO), length(tau)))
  Surface[, 1] <- Meanmf[, 1]


  for (kk in 2:length(tau)) {
    minterest <- Meanmf[, kk] - Meanmf[, (kk - 1)]
    Surface[, kk] <- minterest
  }

  output <- list(Surface,MWNO)
  names(output) <- c("Surface","MWNO")

  return(output)

}
