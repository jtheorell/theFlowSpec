#' Peak identification for higher-level functions. 
#' 
#' This function is primarily thought to be used internally to define peaks
#' in data. 
#' One function is borrowed from package vulcan, namely the vulcan::densityauc, 
#' which is neat, but the package is large and significantly increases the 
#' installation time, and the importing is thus discarded. 
#' @param markerData The data that the peaks should be identified for
#' @param volThresh The cutoff ratio of the volume for each secondary peak, under
#' which it is not considered to be a peak
#' @param distThresh The cutoff under which two peaks are considered one, as
#' they are too close to each other. This value between 0 and 1 corresponds to a
#' fraction from the 10th to the 90th percentile of the data range that the 
#' peaks must be separated by to count. Defaults to 0.1 or 10 percent of the 
#' distance. 
#' @param adjust The value deciding the accuracy of the density calculation. The
#' higher the value, the lower the sensitivity for small aberrations in the 
#' density. 
#' @param nPeaks The number of peaks that should be exported. If n+1 
#' fulfilling the volRatio criterion are found, the peaks most separated in 
#' space are chosen. 
#' @importFrom zoo rollmean
#' @export peakIdenti
peakIdenti <- function(markerData, volThresh = 0.05, distThresh = 0.1, 
                       adjust = 2, nPeaks = 2) {
    
    Da = density(markerData, adjust = adjust, na.rm = TRUE)
    DeltaY = diff(Da$y)
    Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    
    #Here, the top peak decides if the peaks are odd or even numbers in Turns.
    maxPeakPos <- which.max(Da$y[Turns])
    if(maxPeakPos %% 2 == 0){
        allPeakMaxima <- Turns[seq_along(Turns) %% 2 == 0]

    } else {
        allPeakMaxima <- Turns[seq_along(Turns) %% 2 != 0]
    }
    #Either, there is only one peak, and in that case, no further analyses are
    #needed. If not, a set of functions will find which peaks to choose
    if(length(allPeakMaxima) == 1){
        return(Da$x[allPeakMaxima])
    }
    #Now, the volume for each peak is calculated.
    #Each peak is defined as the region from the preceding to the succeding 
    #deflection point
    densVols <- vector()
    
    for(i in seq_along(allPeakMaxima)){
        focusTurnPos <- which(Turns == allPeakMaxima[i])
        if(focusTurnPos == 1){
            peakRegion <- c(1, Turns[2])
        } else if(focusTurnPos == length(Turns)){
            peakRegion <- c(Turns[focusTurnPos-1], length(Da$x))
        } else {
            peakRegion <- c(Turns[focusTurnPos-1], Turns[focusTurnPos+1])
        }
        peakRegVals <- c(Da$x[peakRegion[1]], Da$x[peakRegion[2]])
        
        #Now calculate the volume of the region in question
        xt <- diff(Da$x[Da$x > peakRegVals[1] & Da$x < peakRegVals[2]])
        yt <- rollmean(Da$y[Da$x > peakRegVals[1] & Da$x < peakRegVals[2]], 
                       2)
        densVols[i] <- sum(xt * yt)
    }
    
    
    #Here, the volumes are sorted by size
    allPeakMaximaSorted <- allPeakMaxima[order(densVols, decreasing = TRUE)]
    #Now, all peaks with a volume smaller than volThresh are excluded
    
    peakMaximaReal <- allPeakMaximaSorted[which(densVols > volThresh)]
    
    #Now, if this returns only one value, this value is returned here
    if(length(peakMaximaReal) == 1){
        return(Da$x[peakMaximaReal])
    }
    #Here, it is investigaded if the peaks are separated by more than 
    #a set fraction of the total width of the data. This is done 
    #iteratively, first considering the largest and the second largest 
    #peak, then considering if the third peak is separated from both the
    #first and the second peak, and so on. 
    
    threshDist <- length(Da$x)*distThresh
    peakVector <- peakMaximaReal[1]
    j <- 2
    for(i in seq(2, length(peakMaximaReal))){
        if(all(vapply(seq(1,i-1), function(x) 
            abs(peakMaximaReal[x] - peakMaximaReal[i]) > threshDist, 
            TRUE))){
            peakVector[j] <- peakMaximaReal[i]
            j <- j+1
        }
    }
    
    if(nPeaks < length(peakVector)){
        if(nPeaks == 1){
            return(Da$x[peakVector[1]])
        }
        
        mostSepVals <- vapply(peakVector, function(x) 
            sum(abs(x - peakVector)), 1)
        peakVectorReduced <- 
            sort(peakVector[
                order(mostSepVals, 
                      decreasing = TRUE)][seq_len(nPeaks)])
        peakVectorSorted <- sort(peakVectorReduced)
            
    } else {
        peakVectorSorted <- sort(peakVector)
    }
    
    return(Da$x[peakVectorSorted])
}
