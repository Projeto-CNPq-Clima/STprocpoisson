#' Interpolate Accumulated Mean Values at an Unobserved Location
#'
#' This function interpolates the accumulated mean values for a location where no data was observed,
#' using outputs from the STModelWeibullMCMC model. The interpolation is performed over a vector of times.
#'
#' @param resultsSA A list containing the output from `STModelWeibullMCMCSA`, including:
#'   - `MW`: Samples for parameter W.
#'   - `MMj`: Samples for parameter M.
#'   - `Mvw`: Samples for variance sigma^2_w.
#'   - `Mvm`: Samples for variance sigma^2_m.
#'   - `Mbw`: Samples for spatial decay phi_w.
#'   - `Mbm`: Samples for spatial decay phi_m.
#'   - `MPsi`: Samples for regression coefficients Psi.
#'   - `MBeta`: Samples for regression coefficients Beta.
#' @param data A matrix of occurrence times for the event of interest. Each column corresponds to a monitoring station.
#' @param sites A matrix of geographic coordinates where the process was observed.
#' @param Sites1 A vector representing the geographic coordinates of the location where interpolation is to be performed.
#' @param gridt A vector of times at which the accumulated mean is to be estimated.
#' @param Xmr Covariates associated with the shape parameter at the location.
#' @param Xwr Covariates associated with the scale parameter at the location.
#' @param Xw aaaa
#' @param Xm aaaaaa
#'
#' @return A matrix (`MatMean`) containing the interpolated accumulated mean values at the times specified in `gridt`.
#'
#' @export
interpolate_meanSA <- function(resultsSA, data, sites,Xm,Xw, Sites1, gridt, Xmr, Xwr) {
  data<-as.matrix(data)
  sites<-as.matrix(sites)
  Sites1<-as.matrix(Sites1)
  Xm<-as.matrix(Xm)
  Xw<-as.matrix(Xw)
  Xmr<-as.matrix(Xmr)
  Xwr<-as.matrix(Xwr)

  Stotal <- rbind(sites, t(as.matrix(Sites1)))
  MatMean <- NULL
  n <- nrow(sites)

  for (i in 1:length(resultsSA$MW[, 1])) {
    SIGMAWtotal <- gSigma(resultsSA$Mbw[i], resultsSA$Mvw[i], Stotal)
    SIGMAWA12 <- t(as.matrix(SIGMAWtotal[1:n, (n + 1)]))
    SIGMAWA1 <- SIGMAWtotal[1:n, 1:n]
    A2estw <- t(as.matrix(Xwr)) %*% as.matrix(resultsSA$MPsi[i, ]) +
      SIGMAWA12 %*% solve(SIGMAWA1) %*% (as.matrix(resultsSA$MW[i, ]) - Xw %*% as.matrix(resultsSA$MPsi[i, ]))
    SIGMAA2estw <- SIGMAWtotal[(n + 1), (n + 1)] - SIGMAWA12 %*% solve(SIGMAWA1) %*% t(SIGMAWA12)
    WNO <- rnorm(1, A2estw, sd = sqrt(SIGMAA2estw))

    Gama <- exp(WNO)

    SIGMAMtotal <- gSigma(resultsSA$Mbm[i], resultsSA$Mvm[i], Stotal)
    SIGMAMA12 <- t(as.matrix(SIGMAMtotal[1:n, (n + 1)]))
    SIGMAMA1 <- SIGMAMtotal[1:n, 1:n]
    A2estm <- t(as.matrix(Xmr)) %*% as.matrix(resultsSA$MBeta[i, ]) +
      SIGMAMA12 %*% solve(SIGMAMA1) %*% (as.matrix(resultsSA$MMj[i, ]) - Xm %*% as.matrix(resultsSA$MBeta[i, ]))
    SIGMAA2estm <- SIGMAMtotal[(n + 1), (n + 1)] - SIGMAMA12 %*% solve(SIGMAMA1) %*% t(SIGMAMA12)
    MNO <- rnorm(1, A2estm, sd = sqrt(SIGMAA2estm))

    Eta <- exp(MNO)
    Delta <- resultsSA$Mdelta[i]
    Ff <- resultsSA$Mf[i]
    Theta <- resultsSA$Mtheta[i]

    vect <- NULL
    for (j in 1:length(gridt)) {

      #delta,eta,gama,f,ttheta,t
      Mean <- mfSA(Delta,Eta,Gama,Ff,Theta,gridt[j])
      vect <- c(vect, Mean)
    }

    MatMean <- rbind(MatMean, t(as.matrix(vect)))
  }

  return(MatMean)
}
