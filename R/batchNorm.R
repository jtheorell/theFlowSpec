#' Batch normalisation
#'
#'
#' This function is intended to be used, when files from the same individual has
#' been acquired on multiple time points, and the files from the different
#' time points should be compared. The function works most optimally if all the
#' data came from a sample that was processed and frozen the same day, i.e. a
#' true technical control, but it also works reasonably well with data from the
#' same donor from different time points, as it is robust to outlier populations
#' changing between the dates. It is not recommended to apply this function to
#' files with less than 5000 events.
#' @param normSet A flowSet that should be normalized. If the intCtrl should
#' be normalized and used downstream, it needs to be included here.
#' @param intCtrl A control flowFrame that was acquired the same day as the
#' normSet.
#' @param batchNormCtrlFile A file generated with the batchNormCtrlSetup
#' function, using an universal external control.
#' @param normNames The variables that should be normalized. Default is all that
#' do not fall into one of the two excl categories below.
#' @param exclColNmStr A vector of strings that are not accepted in any
#' normNames.
#' @param exclCategorical Logical: should normNames with less than 10
#' unique values be excluded? Defaults to TRUE.
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
#'
#' @return The normSet normalized. If the intCtrl was not among the normSet
#' files, it will be added as the final object to the list.
#'
#' @importFrom BiocGenerics nrow
#' @export batchNorm
#' @export batchNormCoFunc2
batchNorm <- function(normSet, intCtrlFrame, batchNormCtrlFile,
                      normNames = colnames(intCtrlFrame),
                      exclColNmStr = c("ime"), volThresh = 0.05,
                      lowVertThresh = 0.25,
                      transCoFacs = "default") {

    # First, the nrow of the control file is evaluated, so that it
    #contains more than 5000 events.
    if (BiocGenerics::nrow(intCtrlFrame) < 5000){
        stop("The internal control contains less than 5000
              events, which is considered the lower border for
              meaningful normalization")
    }

    # Here, zero trimming and downsampling is performed
    transCoFacs <- fixTransNamesAndCoFacs(focusFrame = intCtrlFrame,
                                          transNames = normNames,
                                          transCoFacs = transCoFacs,
                                          exclColNmStr = exclColNmStr)

    intCtrlFrameTrans <- arcTrans(flowObj = intCtrlFrame,
                                  transCoFacs = transCoFacs,
                                  transNames = normNames)

    # Now, the procedure that identified the most clearly separated subsets with
    #the external control file is repeated for the internal control file
    intNormVals <- lapply(names(transCoFacs), function(i){
            if(is.na(batchNormCtrlFile[[i]]$FilterVar)){
                focCtrlFrameTrans <- intCtrlFrameTrans[,i]
            } else {
                focCtrlFrameTrans <-
                    intCtrlFrameTrans[,c(i,batchNormCtrlFile[[i]]$FilterVar)]
            }
        batchNormCoFunc1(batchNormVals = batchNormCtrlFile[[i]],
                         focCtrlFrameTrans = focCtrlFrameTrans,
                         transCoFac = transCoFacs[i],
                         lowVertThresh = lowVertThresh,
                         volThresh = volThresh,
                         varName = names(transCoFacs[i]))
    })


    names(intNormVals) <- names(transCoFacs)

    # Here depending on if the data is from CyTOF or flow cytometry the analysis
    #takes different paths. In the case of CyTOF, the zero values are always the
    #lower peak and they should not be moved.
    exprsInt <- exprs(intCtrlFrame)
    dimExprsInt <- dim(exprsInt)
    if(sum(exprsInt == 0, na.rm = TRUE) >0.05*dimExprsInt[1]*dimExprsInt[2] &&
       sum(exprsInt < 0, na.rm = TRUE) == 0){
        message("The data is supposed to be CyTOF-data, and normalized
                accordingly")
      intNormVals <- lapply(intNormVals, function(x) {
            if(length(x) == 2){
                x[c(1,2)] <- c(0,0)
            } else {
                x[c(1,3)] <- c(0,0)
            }
            return(x)
        })
    } else {
        message("The data is supposed to be flow cytometry data, is normalized
                accordingly.")
        #Here, no changes are made, as the peaks are thought to represent the
        #true peaks.
    }

    # If two peaks are present in both files, the data is first normalized to
    #its internal peaks, and then re-organized with the information from the
    #other peaks. If, on the other hand there is only one peak to normalize to,
    #only centering is performed with this peak in mind (this will have no
    #effect at all on CyTOF data).

    normSetNorm <- fsApply(normSet, batchNormCoFunc2,
                                newNormVals = intNormVals,
                                normNames = names(transCoFacs))
    return(normSetNorm)
}

batchNormCoFunc1 <- function(batchNormVals, focCtrlFrameTrans, transCoFac,
                             lowVertThresh, volThresh, varName){

    #This function is all in all very similar to the batch norm control setup
    #co-functions, so look there for annotation.
    if(BiocGenerics:::ncol(focCtrlFrameTrans) == 2){
        focCtrlTrimPos <- oneDZeroTrim(focCtrlFrameTrans[,2], trimFrac = 0.1,
                                       nRowOut = 100000,
                                       returnFlowFrame = FALSE,
                                       returnPositions = TRUE)

        focCtrlTrim <- focCtrlFrameTrans[focCtrlTrimPos,]

        focusFilt1And2 <- madFilter(flowObj = focCtrlTrim,
                                    gateVar = 2, nGates = 2,
                                    returnSepFilter = TRUE)
        #Here, it is made sure that the gating indeed returned two filtered
        #populations.
        if(ncol(focusFilt1And2) == 2){
          lowOrHighFilt <- ifelse(batchNormVals$OptFilter == "Low", 1, 2)
          focCtrlFilt <-
            focCtrlTrim[which(focusFilt1And2[,lowOrHighFilt] == 1),1]
          focCtrlTrim <-
            oneDZeroTrim(focusFrame = focCtrlFilt, trimFrac = 0.1,
                         nRowOut = 100000)
          focCtrlPeaks <- peakIdenti(exprs(focCtrlTrim)[,1], nPeaks = 2,
                                     volThresh = volThresh,
                                     returnStats = TRUE)
        } else {
          focCtrlFrameTrim <- oneDZeroTrim(focCtrlFrameTrans[,1],
                                           trimFrac = 0.1,
                                           nRowOut = 100000)
          focCtrlPeaks <- peakIdenti(exprs(focCtrlFrameTrim)[,1],
                                     nPeaks = 1,
                                     volThresh = volThresh,
                                     returnStats = TRUE)
        }
    } else {
      focCtrlFrameTrim <- oneDZeroTrim(focCtrlFrameTrans, trimFrac = 0.1,
                                     nRowOut = 100000)
      nPeaks <- ifelse(is.na(batchNormVals$HighPeakCtrl), 1, 2)

      focCtrlPeaks <- peakIdenti(exprs(focCtrlFrameTrim)[,1],
                               nPeaks = nPeaks,
                               volThresh = volThresh,
                               returnStats = TRUE)
  }

    intPeaksUnTrans <- sinh(focCtrlPeaks$PeakPos)*transCoFac
    # Now, two different methods are used, depending on if two peaks were found
    # in both files or not.
    if (is.na(batchNormVals$HighPeakCtrl) == FALSE &&
              length(intPeaksUnTrans) == 2) {

        print(paste0("Done with ", varName))

        return(c("IntPeakLow" = intPeaksUnTrans[1],
                   "IntPeakHigh" = intPeaksUnTrans[2],
                   "ExtPeakLow" = batchNormVals$LowPeakCtrl,
                   "ExtPeakHigh" = batchNormVals$HighPeakCtrl))

    } else {
        # In this case, we are writing the function in a slightly convoluted
        #way, to allow for the situation where there is two values for one of
        #the control and the sample. If so, the dominant peak are selected,
        #based on the density volume of the peaks as defined in the peakIdenti
        #function.
        if (is.na(batchNormVals$HighPeakCtrl)) {
            if(length(focCtrlPeaks$PeakPos) == 2){
                intPeak <- focCtrlPeaks$PeakPos[which.max(focCtrlPeaks$DensVol)]
                extPeak <- batchNormVals$LowPeakCtrl
            } else {
                intPeak <- focCtrlPeaks$PeakPos
                extPeak <- batchNormVals$LowPeakCtrl
            }

        } else{
            intPeak <- focCtrlPeaks$PeakPos
            extPeak <- unlist(batchNormVals[2:3])[
                which.max(batchNormVals$DensVol)]
        }

        intPeakUnTrans <- sinh(intPeak)*transCoFac

        print(paste0("Done with ", varName))

        return(c("IntPeak" = intPeakUnTrans,
                 "ExtPeak" = extPeak))

    }

}

batchNormCoFunc2 <- function(focusNormFrame, newNormVals, normNames) {
    #First, the expression part of the FCS file is extracted
    focusNormExprs <- exprs(focusNormFrame)
    #This function has an internal loop, so that each FCSfile is handled
    #separately, reducing calculation time and memory usage.
    for(i in seq_along(newNormVals)){
        focusName <- normNames[i]
        focusVals <- newNormVals[[i]]
        focusData <- focusNormExprs[,focusName]

        # First, it is evaluated wheter the data should be centered only or
        # also normalized

        if (length(focusVals) == 2) {
            # This means that centering should be performed
            centerFile1 <- focusData - focusVals[1]
            normalizedFile <- centerFile1 + focusVals[2]
        } else {
            # First, the data is normalized to its internal peaks
            normalizedFile1 <- (focusData - focusVals[1]) /
                (focusVals[2] - focusVals[1])

            # Then, the normalization procedure is reversed with the peak values
            # from the control file
            normalizedFile <- normalizedFile1 * (focusVals[4] -
                                                     focusVals[3]) +
                focusVals[3]
        }

        #Here, the data is put back into the exprs frame where it came from.
        #To reduce memory usage, the same frame is reused.
        focusNormExprs[,focusName] <- normalizedFile
    }

    exprs(focusNormFrame) <- focusNormExprs

    return(focusNormFrame)

}

