#' Correct defects in spectral unmixing by compensation
#'
#'
#' This function aims at solving the common problem of imperfect unmixing and compensation in flow cytometry using traditional fluorescence compensation.
#' @param unmixData Unmixed data.
#' @param correcMat A correction matrix. If the corrections are unknown, as with new data, this can be generated with \code{\link{creEmptyCorrMat}}.
#' @param transform If transformation should be performed. Defaults to TRUE. Might be good to turn off once the correction phase of the analysis is over, and the transofrmaiton should be optimized.
#' @param coFactors Used if transform==TRUE. See \code{\link{transformData}}.
#' @seealso \code{\link{creEmptyCorrMat}}, \code{\link{specUnmix}}, \code{\link{transformData}}
#' @examples
#' #Load some data
#' data(fullPanel)
#'
#' #Load a spectral matrix
#' data(specMat)
#'
#' #Select the rows that should be compensated
#' fluoExprs <- fullPanel[,grepl("[VBRY]", colnames(fullPanel))]
#'
#' #Unmix this dataset
#' dataUnmixed <- specUnmix(fluoExprs, specMat)
#'
#' #Create an empty correction matrix
#' correcMat <- creEmptyCorrMat(dataUnmixed)
#'
#' #Now correct the data with this.
#' correcData <- correctUnmix(dataUnmixed, correcMat)
#'
#' #In this first run, no correction will be included, but in later iterations, the correcMat can be changed to include the corrections needed.
#'
#' #For example, Qdot705 is "overcompensated" to AF700 in the spectral unmixing process. For this reason, a correction is included:
#' correcMat["Qdot705","AF700"] <- 0.1
#'
#' #And now, a new correction is made
#' correcData <- correctUnmix(dataUnmixed, correcMat)
#'
#' #And so on, and so forth.
#'
#' @export correctUnmix
correctUnmix <- function(unmixData, correcMat, transform=TRUE, coFactors=rep(1500, ncol(unmixData))){


  corrUnmixed <- specUnmix(unmixData, correcMat)

  print("Now, transformation is started")
  if(transform==TRUE){
    transformedData <- transformData(corrUnmixed, coFactors=coFactors)
    return(transformedData)
  } else {
    return(corrUnmixed)
  }

}
