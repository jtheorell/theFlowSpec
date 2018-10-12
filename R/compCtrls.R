#' Compensation controls
#'
#' This is a file containing information from 14 compensation controls,
#' converted to a data frame with the function readFlowSetDf:
#' 12 single-stained bead populations, 1 unstained bead and 1 PBMC control,
#' experimentally treated the same way as the "fullPanel" file. 2000 events
#' from each file are included. The cells have been pre-gated based on SSC/FSC
#' characteristics, and the PBMC file has had its singlets removed with FSC-A/FSC-H,
#' to decrease inter-sample variability. For all, all 18 channels where left open when acquired on
#' the BD Fortessa. Date of  acquisition: 29.07.2015.
#'
#'
#' @docType data
#'
#' @usage data(compCtrls)
#'
#' @format An object of class \code{"data.frame"};
#'
#' @keywords datasets
"compCtrls"
