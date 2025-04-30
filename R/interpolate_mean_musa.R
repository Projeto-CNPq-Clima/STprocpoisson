#' Interpolate Accumulated Mean Values at an Unobserved Location
#'
#' This function interpolates the accumulated mean values for a location where no data was observed,
#' using outputs from the STModelMusaOkumotoMCMC model. The interpolation is performed over a vector of times.
#'
#' @param results A list containing the output from `STModelMusaOkumotoMCMC`, including:
#'   - `MW`: Samples for parameter W.
#'   - `Malpha`: Samples for parameter alpha.
#'   - `Mb`: Samples for parameter phi.
#'   - `Mv`: Samples for parameter sigma^2.
#'   - `MPsi`: Samples for parameter Psi.
#' @param data A matrix of occurrence times for the event of interest. Each column corresponds to a monitoring station.
#' @param sites A matrix of geographic coordinates where the process was observed.
#' @param Xw Covariates associated with the scale parameter at the location.
#' @param Sites1 A vector representing the geographic coordinates of the location where interpolation is to be performed.
#' @param gridt A vector of times at which the accumulated mean is to be estimated.
#' @param Xwr Covariates associated with the scale parameter at the location where interpolation is to be performed.
#'
#'
#' @return A matrix (`MatMean`) containing the interpolated accumulated mean values at the times specified in `gridt`.
#'
#' @export
interpolate_mean_musa<- function(results, data, sites,Xw, Sites1, gridt, Xwr) {

  data<-as.matrix(data)
  sites<-as.matrix(sites)
  Sites1<-as.matrix(Sites1)
  Xw<-as.matrix(Xw)
  Xwr<-as.matrix(Xwr)

  Stotal <- rbind(sites, t(as.matrix(Sites1)))
  MatMean <- NULL
  n <- nrow(sites)

  for (i in 1:length(results$MW[, 1])) {
    SIGMAWtotal <- gSigma(results$Mb[i], results$Mv[i], Stotal)
    SIGMAWA12 <- t(as.matrix(SIGMAWtotal[1:n, (n + 1)]))
    SIGMAWA1 <- SIGMAWtotal[1:n, 1:n]
    A2estw <- t(as.matrix(Xwr)) %*% as.matrix(results$MPsi[i, ]) +
      SIGMAWA12 %*% solve(SIGMAWA1) %*% (as.matrix(results$MW[i, ]) - Xw %*% as.matrix(results$MPsi[i, ]))
    SIGMAA2estw <- SIGMAWtotal[(n + 1), (n + 1)] - SIGMAWA12 %*% solve(SIGMAWA1) %*% t(SIGMAWA12)
    WNO <- rnorm(1, A2estw, sd = sqrt(SIGMAA2estw))

    Gama <- WNO

    Alpha<-results$Malpha[i]

    vect <- NULL
    for (j in 1:length(gridt)) {
      Mean <- mfMUSA(Gama,Alpha,gridt[j])
      vect <- c(vect, Mean)
    }

    MatMean <- rbind(MatMean, t(as.matrix(vect)))
  }

  return(MatMean)
}
