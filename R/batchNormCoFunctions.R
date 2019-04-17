batchNormCoFunc1 <- function(normVar, controlVar, volRatio=volRatio){
    
    #And here, it is made clear if it is meaningful to talk about two peaks, or
    #if it is indeed only one.
    peakResultNorm <- twoPeakIdentification(normVar, volRatio = volRatio)
    
    peakResultCtrl <- twoPeakIdentification(controlVar, volRatio = volRatio)
    
       
    #Now, two different methods are used, depending on if two peaks were found
    #in both files or not
    if(peakResultNorm[[1]] == 2 && peakResultCtrl[[1]] == 2){
        
        data.frame("LowPeakFile" = peakResultNorm[[2]][1], 
                   "HighPeakFile" = peakResultNorm[[2]][2], 
                   "LowPeakCtrl" = peakResultCtrl[[2]][1], 
                   "HighPeakCtrl" = peakResultCtrl[[2]][2])
        
    } else {
        #In this case, we are writing the function in a slightly convoluted way,
        #to allow for the situation where there is two values for one of the
        #control and the sample. 
        if(peakResultNorm[[1]] == 2){
            normPeak <- 
                peakResultNorm[[2]][which.min(abs(peakResultNorm[[2]]-
                                                      peakResultCtrl[[2]]))]
            ctrlPeak <- peakResultCtrl[[2]]
        } else if(peakResultCtrl[[1]] == 2){
            ctrlPeak <- 
                peakResultCtrl[[2]][which.min(abs(peakResultCtrl[[2]]-
                                                      peakResultNorm[[2]]))]
            normPeak <- peakResultNorm[[2]]
        } else {
            ctrlPeak <- peakResultCtrl[[2]]
            normPeak <- peakResultNorm[[2]]
        }
        
        data.frame("PeakFile" = normPeak, "PeakCtrl" = ctrlPeak)
        
    }
    
}

batchNormCoFunc2 <- function(newNormValsFocus, focusColNormFile){
    #First, it is evaluated wheter the data should be centered only or 
    #also normalized
    if(length(newNormValsFocus) == 2){
        #This means that centering should be performed
        centerFile1 <- focusColNormFile - newNormValsFocus[,1]
        normalizedFile <- centerFile1 + newNormValsFocus[,2]
    } else {
        #First, the data is normalized to its internal peaks
        normalizedFile1 <- (focusColNormFile - newNormValsFocus[,1]) /
            (newNormValsFocus[,2] - newNormValsFocus[,1])

        #Then, the normalization procedure is reversed with the peak values from
        #the control file
        normalizedFile <- normalizedFile1*(newNormValsFocus[,4] - 
                                               newNormValsFocus[,3]) + 
            newNormValsFocus[,3]
    }
    
    return(normalizedFile)

}

