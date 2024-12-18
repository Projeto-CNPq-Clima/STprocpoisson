# This function converts the time series matrix into a matrix with the times of occurrence of the events of interest.
# The arguments of this function are:
# mat: an array containing in each column the time series corresponding to each location.
# lim: The threshold of interest.

#' Title
#'
#' @param mat aaaa
#' @param lim aaaa
#'
#' @return aaa
transforma <- function(mat, lim) {
  n <- length(mat)

  res <- NULL

  for (i in 1:n) {
    if ((mat[i] >= lim)) {
      temp <- i
      res <- c(res, temp)
    }
  }

  res
}
