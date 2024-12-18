#' Identify Event Times Based on Threshold
#'
#' This function identifies the times at which the values in each column of a matrix exceed a user-defined threshold.
#'
#' @param series A matrix of dimensions Txn, where each column represents the time series of a random variable from a monitoring station.
#' @param threshold A numeric value representing the user-defined threshold to identify times of interest.
#' @return A matrix where each column contains the times (indices) at which the values of the respective column in `series` exceeded the threshold. Rows are padded with NA where no event occurred.
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
  data <- MATRIZ
  return(data)
}
