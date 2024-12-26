#' Geographic Coordinates of Rain Gauge Stations in Maranhão and Piauí, Brazil
#'
#' This dataset contains the geographic coordinates of 20 rain gauge stations located in the states of Maranhão and Piauí, Brazil. These stations were used to monitor daily precipitation data over a 10-year period, from January 1, 2013, to December 31, 2022.
#'
#' @format A data frame with 20 rows and 2 variables:
#' \describe{
#'   \item{latitude}{The latitude of the station, in decimal degrees.}
#'   \item{longitude}{The longitude of the station, in decimal degrees.}
#' }
#'
#' @details
#' The rain gauge stations are part of a network managed by the National Institute of Meteorology (INMET) and the National Water and Basic Sanitation Agency (ANA). These stations are distributed across Maranhão and Piauí, northeastern Brazil, covering a region known for its diverse climatic and geographic conditions.
#'
#' The coordinates in this dataset can be used for spatial analysis or visualization of the study area. They complement the `series` dataset, which contains daily precipitation data recorded at these locations.
#'
#' @source National Institute of Meteorology (INMET) and National Water and Basic Sanitation Agency (ANA).
#'
#' @examples
#' # Load the dataset
#' data(sites)
#'
#' # View the structure of the data
#' str(sites)
#'
#' # Plot the geographic locations of the rain gauge stations
#' plot(sites$longitude, sites$latitude, pch = 16, col = "blue",
#'      xlab = "Longitude", ylab = "Latitude",
#'      main = "Rain Gauge Stations in Maranhão and Piauí")
"sites"
