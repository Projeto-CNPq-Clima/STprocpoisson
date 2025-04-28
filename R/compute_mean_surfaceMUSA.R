#' Compute Interpolated Mean Surface for Spatiotemporal Model
#'
#' This function interpolates the mean values on a grid of points (`DNO`) for a spatiotemporal
#' nonhomogeneous Poisson model with Weibull intensity. The function utilizes the MCMC outputs
#' from `STModelMusaOkMCMC` and applies a Gaussian process-based interpolation.
#'
#' @param results A list containing the output from `STModelMusaOkumotoMCMC`, including:
#'   - `MW`: Samples for parameter W.
#'   - `MWT`: Acceptance indicators for parameter W.
#'   - `Malpha`: Samples for parameter alpha.
#'   - `MalphaT`: Acceptance indicators for parameter alpha.
#'   - `Mb`: Samples for parameter phi.
#'   - `MbT`: Acceptance indicators for phi.
#'   - `Mv`: Samples for parameter sigma^2.
#'   - `MPsi`: Samples for parameter Psi.
#' @param sites A matrix with geographic coordinates of the monitoring stations.
#' @param X Covariates for the scale parameter of the Weibull intensity.
#' @param DNO A grid of points where interpolation is to be performed.
#' @param CovXNO Covariates for the scale parameter at the grid points.
#' @param tau A vector of temporal points for which the mean surface is computed.
#'
#' @return A list: (`Surface`) containing the interpolated mean values at the grid points
#' for each temporal point in `tau`. The first column contains the mean values at the initial
#' time step, and subsequent columns contain the differences between consecutive time steps.
#' @export

compute_mean_surfaceMUSA <- function(results, sites, X, DNO, CovXNO, tau) {

  X<-as.matrix(X)

  jj <- nrow(DNO)
  res <- rbind(DNO, as.matrix(sites))
  tt <- nrow(res)

  Xr <- rbind(CovXNO, X)

  MW <- results$MW
  Malpha <- results$Malpha
  Mb <- results$Mb
  Mv <- results$Mv
  MPsi <- results$MPsi

  MWNO <- array(NA, dim = c(nrow(MPsi), jj))
  MfmNO <- NULL


  for (h in 1:nrow(MW)) {
    WM <- as.matrix(MW[h, ])
    sigma <- gSigma(Mb[h], Mv[h], res)
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
      MfmNO <- rbind(MfmNO, t(as.matrix(mfMUSA(MWNO[h, ], Malpha[h], tau[i]))))
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
