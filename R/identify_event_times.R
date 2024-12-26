#' Identify Event Times Based on Threshold
#'
#' This function identifies the times at which the values in each column of a matrix exceed a user-defined threshold.
#'
#' @param series A matrix of dimensions Txn, where each column represents the time series of a random variable from a monitoring station.
#' @param threshold A numeric value representing the user-defined threshold to identify times of interest.
#' @return It returns a matrix where each column corresponds to a location, and each row contains the times of occurrence for that location. The matrix dimensions are m x n, where:
#' - `m`: Maximum number of occurrences across all locations.
#' - `n`: Number of locations.
#'
#' @examples
#' library(STprocpoisson)
#' data(series)
#' Times<-identify_event_times(series,20)
#'
#'
#' @export

identify_event_times <- function(series, threshold) {
  limia <- threshold
  tamanho <- NULL
  for (i in 1:ncol(series)) {
    temp <- length(transforma(series[, i], limia))
    tamanho <- c(tamanho, temp)
  }
  MATRIZ <- array(NA, dim = c(max(tamanho), ncol(series)))
  for (i in 1:ncol(series)) {
    temp <- as.matrix(transforma(series[, i], limia))
    MATRIZ[1:tamanho[i], i] <- temp
  }
  MATRIZ
  return(MATRIZ)
}
