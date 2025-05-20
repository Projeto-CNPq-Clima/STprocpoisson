#' Fit and Plot Nonhomogeneous Poisson Weibull Model
#'
#' This function fits a nonhomogeneous Poisson Weibull model to data from a specified monitoring station,
#' generates a plot of the observed and estimated cumulative mean functions, and returns relevant metrics.
#'
#' @param data A matrix of occurrence times for the event of interest. Each column corresponds to a monitoring station.
#' @param results A list containing the output from `STModelWeibullMCMC`, which includes:
#'   - `MW`: Samples for parameter W.
#'   - `Mbeta`: Samples for parameter Beta.
#'   - `Malpha`: Samples for parameter Alpha.
#' @param l An integer specifying the index of the monitoring station to analyze.
#'
#' @return A list containing:
#' \describe{
#'   \item{retemp}{The occurrence times of the events of interest.}
#'   \item{medy1}{The estimated cumulative mean function at the event times.}
#'   \item{mtnp}{The observed cumulative mean function at the event times.}
#' }
#'
#' @export

fit_weibull_model <- function(data, results, l) {

  Mean <- NULL
  MatMean <- NULL
  for (i in 1:length(results$MW[, l])) {
    Gama <- results$MW[i, l]
    Eta <- results$MMj[i, l]

    gridt <- data[1:(length(data[, l]) - sum(is.na(data[, l]))), l]

    Mean <- mfWEIBULL(Gama, Eta, gridt)

    MatMean <- rbind(MatMean, t(as.matrix(Mean)))
  }

 retemp <- mtnp(gridt)[,1]

  # Summary statistics for the estimated function
  medy1 <- apply(MatMean, 2, mean)
  per25 <- apply(MatMean, 2, quantile, probs = 0.025)
  per975 <- apply(MatMean, 2, quantile, probs = 0.975)
  infy <- min(per25)
  supy <- max(per975)
  infx <- min(retemp)
  supx <- max(retemp)
  xx <- c(gridt, rev(gridt))
  yy <- c(per25, rev(per975))

  # Plot the results
  plot(xx, yy,
    type = "n", xlab = "Time", ylab = "Cumulative Mean Function",
    xlim = c(infx, supx), ylim = c(infy, supy), cex.lab = 1.5
  )
  polygon(xx, yy, col = "gray", border = NA)
  par(new = TRUE)
  plot(mtnp(data[, l]),
    type = "l", xlim = c(infx, supx), ylim = c(infy, supy),
    xlab = "", ylab = "", main = paste("Station", l), col = "blue", lwd = 2
  )
  par(new = TRUE)
  plot(retemp, medy1,
    type = "l", xlim = c(infx, supx), ylim = c(infy, supy),
    xlab = "", ylab = "", lwd = 2, col = "red"
  )

  # Return the results
  return(list(
    retemp = retemp,
    medy1 = medy1,
    mtnp = mtnp(data[, l])[, 2]
  ))
}
