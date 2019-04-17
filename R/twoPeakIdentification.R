#' @export twoPeakIdentification
twoPeakIdentification <- function(markerData, volRatio = volRatio) {
    Da = density(markerData)
    DeltaY = diff(Da$y)
    Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    # Identify the two highest turnpoints
    topTwo <- Turns[order(Da$y[Turns], decreasing = TRUE)][1:2]

    # And here, it is made clear if it is meaningful to talk about two peaks, or
    # if it is indeed only one.
    # First, it is determined if there are two peaks whatsoever, and if so, by
    # how much these are separated. Currently, they need to be separated by at
    # least 10% of the total range from the first to th 99th percentile to count
    # as two peaks.
    robRange <- quantile(markerData, c(0.1, 0.9))
    tenth <- abs(unname((robRange[2] - robRange[1]) * 0.1))
    if (is.na(topTwo[2]) == FALSE && abs(Da$x[topTwo][1] -
        Da$x[topTwo][2]) > tenth) {
        # Now, if two peaks are present, identify the lowest turnpoint between
        # the topTwo
        betweenTurns <- Turns[Turns < max(topTwo) & Turns > min(topTwo)]
        minTurn <- betweenTurns[which.min(Da$y[betweenTurns])]
        if (topTwo[1] < topTwo[2]) {
            startHigh <- Da$x[1]
            stopHigh <- Da$x[minTurn]
            startLow <- Da$x[minTurn]
            stopLow <- Da$x[length(Da$x)]
            topTurns <- c(topTwo[1], topTwo[2])
        } else {
            startLow <- Da$x[1]
            stopLow <- Da$x[minTurn]
            startHigh <- Da$x[minTurn]
            stopHigh <- Da$x[length(Da$x)]
            topTurns <- c(topTwo[2], topTwo[1])
        }

        volHigh <- length(markerData[which(markerData > startHigh & markerData
        < stopHigh)])
        volLow <- length(markerData[which(markerData > startLow & markerData
        < stopLow)])

        # And here, it is determined if the second peak counts, which it only
        # does if it has a volume that is at least a set percentage of the
        # largest one. Default being 5%.
        if (volHigh > 0 && volLow > 0 && volLow / volHigh > volRatio) {
            peakResult <- list("nPeaks" = 2, "peakData" = Da$x[topTwo])

        } else {
            peakResult <- list("nPeaks" = 1, "peakData" =  Da$x[topTwo[1]])
        }
    } else {
        peakResult <- list("nPeaks" = 1, "peakData" = Da$x[topTwo[1]])
    }

}
