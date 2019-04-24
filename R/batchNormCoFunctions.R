#' @export batchNormCoFunc1
#' @export batchNormCoFunc2
batchNormCoFunc1 <- function(intVar, extVar, transCoFac, volThresh) {

    #Now, the normalization files are transformed, using the transformation
    #values above and the base hyperbolic arc sinus function
    
    intCtrlTrans <- asinh(intVar/transCoFac)
    extCtrlTrans <- asinh(extVar/transCoFac)
    
    # And here, it is made clear if it is meaningful to talk about two peaks, or
    # if it is indeed only one.
    peakResultInt <- peakIdenti(intCtrlTrans, volThresh = volThresh)

    peakResultExt <- peakIdenti(extCtrlTrans, volThresh = volThresh)

    # Now, two different methods are used, depending on if two peaks were found
    # in both files or not
    if (length(peakResultInt) == 2 && length(peakResultExt) == 2) {

        intPeakUnTrans <- sinh(peakResultInt)*transCoFac
        extPeakUnTrans <- sinh(peakResultExt)*transCoFac
        
        data.frame("LowPeakFile" = intPeakUnTrans[1],
                   "HighPeakFile" = intPeakUnTrans[2],
                   "LowPeakCtrl" = extPeakUnTrans[1],
                   "HighPeakCtrl" = extPeakUnTrans[2])

    } else {
        # In this case, we are writing the function in a slightly convoluted 
        #way, to allow for the situation where there is two values for one of 
        #the control and the sample.
        if (length(peakResultInt) == 2) {
            intPeak <-
                peakResultInt[which.min(abs(peakResultInt - peakResultExt))]
            extPeak <- peakResultExt
        } else if (length(peakResultExt) == 2) {
            extPeak <-
                peakResultExt[which.min(abs(peakResultExt - peakResultInt))]
            intPeak <- peakResultInt
        } else {
            extPeak <- peakResultExt
            intPeak <- peakResultInt
        }

        intPeakUnTrans <- sinh(intPeak)*transCoFac
        extPeakUnTrans <- sinh(extPeak)*transCoFac
        
        data.frame("IntPeak" = intPeakUnTrans, "ExtPeak" = extPeakUnTrans)

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
            centerFile1 <- focusData - focusVals[, 1]
            normalizedFile <- centerFile1 + focusVals[, 2]
        } else {
            # First, the data is normalized to its internal peaks
            normalizedFile1 <- (focusData - focusVals[, 1]) /
                (focusVals[, 2] - focusVals[, 1])
            
            # Then, the normalization procedure is reversed with the peak values
            # from the control file
            normalizedFile <- normalizedFile1 * (focusVals[, 4] -
                                                     focusVals[, 3]) +
                focusVals[, 3]
        }
        
        #Here, the data is put back into the exprs frame where it came from. 
        #To reduce memory usage, the same frame is reused.
        focusNormExprs[,focusName] <- normalizedFile
    }
    
    exprs(focusNormFrame) <- focusNormExprs
    
    return(focusNormFrame)

}
