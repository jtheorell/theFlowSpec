#' Plotting all variables against a single variable
#'
#'
#' This function is useful both when setting appropriate gates and when the adjustments of the compensation are done
#' @importFrom grDevices col2rgb colorRampPalette densCols dev.off png
#' @importFrom graphics mtext par plot title
#' @param rawData This can be either be a data frame, a FlowFrame or a FlowSet.
#' @param specMat Here, the matrix generated by the secMatCalc function is added, that is used for the
#' @return The unmixed data. It will be returned in the format it was imported as.
#'
#' @export specUnmix
specUnmix <- function(rawData, specMat){

}