#' Setup the file used for batch normalisation
#'
#' This function is used to generate the batch normalization values
#' that ideally are based on an identical control that is included
#' in each experiment.
#' @param ctrlFrame A control flowFrame that all the batches will be normalized
#' to.
#' @param normNames The variables that should be normalized. Default is all that
#' do not fall into one of the two excl categories below.
#' @param exclColNmStr A vector of strings that are not accepted in any
#' normNames.
#' @param volThresh Defining how small the smaller of the two peaks can be to be
#' considered a true peak. It is a fraction of the volume of the larger peak.
#' Default is 0.05, i.e. if the volume of the second peak is 5 percent or larger
#' than the volume of the first peak, it is considered a peak.
#' @param lowVertThresh Defining how low the turnpoint needs to be between the
#' peaks for these to be considered well-separated. 25 percent of the lower
#' peak is the default.
#' @param transCoFacs This vector of named values define the values for the
#' transformation during the normalization. This is only applied internally, so
#' transformation needs to be performed afterwards, preferrably with individual
#' values for each channel. In the "default" case, the function defines the file
#' as a CyTOF file, and applies the transformation value 8, if >5 percent of the
#' values are 0. Otherwise, the value 256 is applied. NB. The entries need to be
#' named in the same way as the normNames to secure that the right factor is
#' added to each variable.
#' @param nonLevel1Var If any of the markers are known to be unreliable,
#' these will be excluded from being used as filters for other markers.
#' @importFrom BiocGenerics ncol
#' @export batchNormCtrlSetup
batchNormCtrlSetup <- function(ctrlFrame, normNames = colnames(ctrlFrame),
                               exclColNmStr = c("ime"), volThresh = 0.05,
                               lowVertThresh = 0.25, transCoFacs = "default",
                               nonLevel1Var){

  # First, the nrow of the file is evaluated, so that it
  #contains more than 5000 events.
  if (BiocGenerics::nrow(ctrlFrame) < 5000) {
    stop("The control file contains less than 5000
              events, which is the (completely arbitrary) lower border for
              useful normalization")
  }

  # Here, zero trimming and downsampling is performed
  transCoFacs <- fixTransNamesAndCoFacs(focusFrame = ctrlFrame,
                                        transNames = normNames,
                                        transCoFacs = transCoFacs,
                                        exclColNmStr = exclColNmStr)
  ctrlFrameTrim <- oneDZeroTrim(focusFrame = ctrlFrame, trimFrac = 0.1,
                                nRowOut = 100000)

  # Here, the first of two functions identifying the peaks comes:
  newNormVals <- lapply(names(transCoFacs), function(i)
    theFlowSpec:::bNormCSetupCo1(ctrlVar = exprs(ctrlFrameTrim[, i])[,1],
                   transCoFac = transCoFacs[i],
                   lowVertThresh = lowVertThresh,
                   volThresh = volThresh*2))

  names(newNormVals) <- names(transCoFacs)

  if(length(which(lapply(newNormVals, `[[`, 1)==TRUE))==0){
    stop("This dataset does not contain any markers with a well-separated
         bimodal distribution, which complicates normalisation. Try increasing
         the low vertex threshold, which will allow for variables with two
         good peaks, that are a poorly separated, to be used for further
         analysis. Or consider other strategies whatsoever.")
  }

  #Now, the variables that passed this first test are ranked according to their
  #separability. However, if there are concerns about specific markers,
  #they can also be excluded here, for further refinement in the second
  #round

  WellDefined1 <- unlist(lapply(newNormVals, `[[`, "WellDefined"))

  newNormValsWDef1 <- newNormVals[WellDefined1]

  if(missing(nonLevel1Var)==FALSE){
    newNormValsWDef1 <- newNormValsWDef1[-which(names(newNormValsWDef1) %in%
                                                  nonLevel1Var)]
  }
  ctrlVarsUnRanked <- names(newNormValsWDef1)
  LowVertexes <- unlist(lapply(newNormValsWDef1, `[[`, "LowVertex"))
  PeakDistDivided <- 1/unlist(lapply(newNormValsWDef1, `[[`, "PeakDist"))

  ctrlVarsRanked <- ctrlVarsUnRanked[order(LowVertexes, PeakDistDivided)]

  #And now, the variables that did not pass the first test for well-definedness
  #will now go thourhg a second round, where dividing the events into a group
  #positive for each of the factors below

  ctrlFrameTrans <- arcTrans(flowObj = ctrlFrame, transCoFacs = transCoFacs,
                             transNames = normNames)

  newNormVals <- lapply(names(transCoFacs), function(i)
    theFlowSpec:::bNormCSetupCo2(newNormVal = newNormVals[[i]],
                                 transCoFac = transCoFacs[i],
                                 ctrlVarName = i,
                                 ctrlVarsRankedFrameTrans =
                                   ctrlFrameTrans[, c(i, ctrlVarsRanked)],
                                 lowVertThresh = lowVertThresh,
                                 volThresh = volThresh,
                                 nonLevel1Var = nonLevel1Var))

  names(newNormVals) <- names(transCoFacs)

  return(newNormVals)

}

######################
#####################
bNormCSetupCo1 <- function(ctrlVar, transCoFac, lowVertThresh, volThresh) {

  ctrlVarTrans <- asinh(ctrlVar/transCoFac)
  # And here, it is made clear if it is meaningful to talk about two peaks, or
  # if it is indeed only one.
  ctrlVarPeaks <- peakIdenti(ctrlVarTrans, nPeaks = 3, volThresh = volThresh,
                             returnStats = TRUE)

  # Here, the variables are classified, according to the complexity of the peak
  # result. A variable passes if it contains two well-separated peaks. If
  #one or three peaks are detected, or if the two identified peaks are not
  #readily separated, the variable will undergo a second attempt to further
  #stratify the data, to increase the separation of the populations.

  ctrlVarPeaksUnTrans <- sinh(ctrlVarPeaks$PeakPos)*transCoFac

  if(length(ctrlVarPeaksUnTrans) == 2){

    localNormVals <- list("WellDefined" = NA,
                          "LowPeakCtrl" = ctrlVarPeaksUnTrans[1],
                          "HighPeakCtrl" = ctrlVarPeaksUnTrans[2],
                          "DensVol" = ctrlVarPeaks$DensVol,
                          "LowVertex" = ctrlVarPeaks$LowVertex,
                          "PeakDist" = abs(ctrlVarPeaksUnTrans[1]-
                                             ctrlVarPeaksUnTrans[2]),
                          "FilterVar" = NA,
                          "OptFilter" = NA)

    localNormVals$WellDefined <-
      ifelse(ctrlVarPeaks$LowVertex/min(ctrlVarPeaks$Height) <=
               lowVertThresh, TRUE, FALSE)
  } else if(length(ctrlVarPeaksUnTrans) == 1){
    localNormVals <- list("WellDefined" = FALSE,
           "LowPeakCtrl" = ctrlVarPeaksUnTrans[1],
           "HighPeakCtrl" = NA,
           "FilterVar" = NA)
  } else {
    localNormVals <- list("WellDefined" = FALSE,
                          "LowPeakCtrl" = ctrlVarPeaksUnTrans[1],
                          "HighPeakCtrl" =
                            ctrlVarPeaksUnTrans[c(2,3)][
                              which.max(ctrlVarPeaks$DensVol[c(2,3)])],
                          "DensVol" = c(ctrlVarPeaks$DensVol[1],
                                        ctrlVarPeaks$DensVol[c(2,3)][
                                          which.max(ctrlVarPeaks$DensVol[c(2,3)])]),
                          "FilterVar" = NA)
    }

  return(localNormVals)

}

bNormCSetupCo2 <- function(newNormVal, ctrlVar, ctrlVarName, transCoFac,
                           ctrlVarsRankedFrameTrans, lowVertThresh, volThresh,
                           nonLevel1Var){

  if(newNormVal$WellDefined && ctrlVarName %in% nonLevel1Var == FALSE){
    return(newNormVal)
  } else {
    #Here, a loop is applied, where the ideal variables are applied in their
    #ranked order, to see if the gating of these markers can enhance the
    #separation.
    wellSep <- FALSE
    n <- 2

    while(wellSep == FALSE && n <
          BiocGenerics:::ncol(ctrlVarsRankedFrameTrans)){

      ctrlWhileTrimPos <- oneDZeroTrim(focusFrame =
                                         ctrlVarsRankedFrameTrans[,n],
                                       trimFrac = 0.1, nRowOut = 100000,
                                       returnFlowFrame = FALSE,
                                       returnPositions = TRUE)

      ctrlWhileTrim <- ctrlVarsRankedFrameTrans[ctrlWhileTrimPos,]

      focusFilt1And2 <- madFilter(flowObj = ctrlWhileTrim,
                                  gateVar = n, nGates = 2,
                                  returnSepFilter = TRUE)

      #Here, the two gated subsets are investigated further,
      ctrlVarPeaksLocal <- list()
      for(i in seq_len(ncol(focusFilt1And2))){
         ctrlVarFilt <- ctrlWhileTrim[which(focusFilt1And2[,i] == 1),1]
         ctrlVarTrim <- oneDZeroTrim(focusFrame = ctrlVarFilt, trimFrac = 0.1,
                                       nRowOut = 100000)
         ctrlVarPeaksLocal[[i]] <- peakIdenti(exprs(ctrlVarTrim)[,1],
                                              nPeaks = 3,
                                              volThresh = volThresh,
                                              returnStats = TRUE)
      }

      optimalFilter <- which.min(unlist(lapply(ctrlVarPeaksLocal, "[[",
                                               "LowVertex")))[1]
      ctrlVarPeaksOpt <- ctrlVarPeaksLocal[[optimalFilter]]

      if(length(ctrlVarPeaksOpt$PeakPos) == 2 &&
         ctrlVarPeaksOpt$LowVertex/min(ctrlVarPeaksOpt$Height) <=
         lowVertThresh){

        optFilter <- c("Low","High")[optimalFilter]

        ctrlVarPeaksOptUnTrans <- sinh(ctrlVarPeaksOpt$PeakPos)*transCoFac

        localNormVals <- list("WellDefined" = TRUE,
             "LowPeakCtrl" = ctrlVarPeaksOptUnTrans[1],
             "HighPeakCtrl" = ctrlVarPeaksOptUnTrans[2],
             "DensVol" = ctrlVarPeaksOpt$DensVol,
             "LowVertex" = ctrlVarPeaksOpt$LowVertex,
             "PeakDist" = abs(ctrlVarPeaksOptUnTrans[1]-
                                ctrlVarPeaksOptUnTrans[2]),
             "FilterVar" = colnames(ctrlVarsRankedFrameTrans[,n]),
             "OptFilter" = optFilter)

        wellSep <- TRUE
      }


      n <- n + 1

    }

  }

  print(paste0("Done with ", ctrlVarName))

  if(wellSep){
    return(localNormVals)
  } else {
    return(newNormVal)
  }

}
