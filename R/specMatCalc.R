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
#' NegSamples <- c(rep("Unstained", times=length(fluorochromes)), 0)
#' compScheme <- cbind(Ids, PosSamples, NegSamples)
#'
#' #Now, select only the fluorescence channels, marked by one of the lasers
#' compControls <- compCtrls[,grepl("[VBRY]", colnames(compCtrls))]
#' ids <- compCtrls$id
#'
#' #And run the function
#' spectralMatrix <- specMatCalc(compControls, compScheme, ids)
#'
#' @export specMatCalc
specMatCalc <- function(compControls, compScheme, ids){

  #Now, to save some computational time, the negative samples, that are reused, are first calculated
  negSamples <- unique(compScheme[,3])

  negPeaksList=list()
  for(i in 1:length(negSamples)){
    if(negSamples[i]=="0"){
      negPeaksList[[i]] <- rep(0, times=ncol(compControls))
    } else {
      negSample <-  compControls[ids==negSamples[i],]
      negPeaksList[[i]] <- specMatCalcCoCoFunction(negSample)
    }

  }
  specMat <- list()
  for(i in 1:nrow(compScheme)){


    posSample <- compControls[ids==compScheme[i,2],]
    posPeaks <- specMatCalcCoCoFunction(posSample)
    negPeaks <- negPeaksList[[which(negSamples==compScheme[i,3])]]

    peaks <- posPeaks-negPeaks
    peaksFraction <- peaks/max(peaks)
    specMat[[i]] <-peaksFraction

  }

  specDf <- do.call("rbind", specMat)

  row.names(specDf) <- compScheme[,1]
  return(specDf)
}
