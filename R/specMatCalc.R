#' Calculating the matrix used for spectral unmixing
#'
#'
#' This algoritm takes the (pre-gated) single-stained controls and negative controls, including an autofluorescence control and estimates the unmixing for all fluorescent variables.
#' @param compControls A dataframe with two components (and many columns):
#' \describe{
#'            \item{Name}{A column that associates each observation with a certain sample}
#'            \item{Channels}{All the channels detecting fluosrescent information in the file.}
#'            \item{Negative}{The file number in the compControls list/FlowSet negative for the fluorochrome/marker in question. For autofluorescense subtraction, the value here is set to "NA".}
#'          }
#' @param compScheme A dataframe with three columns:
#' \describe{
#'            \item{Fluorochrome}{The name of the fluorochrome/marker present in the Name column of the compControls file}
#'            \item{Positive}{The file number in the compControls list/FlowSet positive for the fluorochrome/marker in question.}
#'            \item{Negative}{The file number in the compControls list/FlowSet negative for the fluorochrome/marker in question. For autofluorescense subtraction, the value here is set to "NA".}
#'          }
#' @param Ids A vector that connects each row in the dataframe with a certain single-stained sample.
#' @return A data frame with each row representing a fluorochrome and each column representing a detector
#' @examples
#' #Load suitable compensation controls. Use the readFlowSetDf function.
#' #NB! The compensation controls need to be cleaned up, so that doublets do not
#' #cause any harm.
#' data(compCtrls)
#'
#' #Now construct the compScheme file
#' allControls <- unique(compCtrls$id)
#' fluorochromes <- allControls[-grep("stained", allControls)]
#'
#' #Now pair each fluorochrome with a positive and a negative sample name.
#' #NB! The negative sample needs to have exacely the same autofluorescence as the
#' #stained sample (not ture for the autofluorescence control, that instead needs to have
#' #the same autofluorescence as the actual samples that will be compensated).
#'
#' Ids <- c(fluorochromes, "Autofluorescence")
#' PosSamples <- c(fluorochromes, "PBMC_unstained")
#' NegSamples <- c(rep("Unstained", times=length(fluorochromes)), "0")
#' compScheme <- cbind(Ids, PosSamples, NegSamples)
#'
#' #Now, select only the fluorescence channels, marked by one of the lasers
#' compControls <- compCtrls[,grepl("[VBRY]", colnames(compCtrls))]
#' ids <- compCtrls$id
#'
#'
#' @export specMatCalc
specMatCalc <- function(compControls, compScheme, ids){


  specMat <- list()
  for(i in 1:nrow(compScheme)){

    posSample <- compControls[ids==compScheme[2,i],]
    negSample <- compControls[ids==compScheme[3,i],]


  }


  if(center=="peak"){
    #The peak of the data is defined
    if(length(x)<500){
      nBreaks <- 10
    } else {
      nBreaks <- length(x)/50
    }
    histdata <- hist(responseVector, breaks=nBreaks, plot=FALSE)
    zeroPosition <- histdata$mids[match(max(histdata$counts), histdata$counts)]

    #And the position for this this peak is subtracted from all points
    responseVector <- responseVector-zeroPosition
  }
}
